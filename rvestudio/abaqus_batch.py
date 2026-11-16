from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path
from time import perf_counter
from typing import Callable, Iterable
import json
import math

import numpy as np
import pandas as pd

from .abaqus import AbaqusSettings, prepare_abaqus_run, run_abaqus, _safe_name
from .abaqus_homogenization import (
    PeriodicHomogenizationSettings,
    prepare_periodic_homogenization,
    run_periodic_homogenization,
)
from .database import _read_table, _write_table, build_case_from_row
from .homogenization import homogenize_rve

BACKEND_ANALYTICAL = "Analytical"
BACKEND_STRUCTURAL = "Abaqus - Structural FE"
BACKEND_CH21 = "Abaqus - Periodic CH21"


def _active_rows(df: pd.DataFrame, case_limit: int | None = None, case_ids: Iterable[str] | None = None) -> pd.DataFrame:
    active = df.copy()
    if "Active" in active.columns:
        col = active["Active"]
        if col.dtype == bool:
            active = active[col]
        else:
            truth = col.astype(str).str.strip().str.lower().isin({"1", "true", "yes", "y", "on"})
            active = active[truth]
    if case_ids:
        wanted = {str(x).strip() for x in case_ids if str(x).strip()}
        active = active[active["CaseID"].astype(str).isin(wanted)]
    if case_limit:
        active = active.head(int(case_limit))
    return active


def _term_name(i: int, j: int) -> str:
    return f"C{i + 1}{j + 1}"


def _flatten_ch21_result(case_id: str, result: dict, analytical_ch: np.ndarray | None, elapsed: float) -> dict:
    C = np.asarray(result.get("CH", result.get("CH_energy")), float).reshape(6, 6)
    row: dict = {
        "CaseID": str(case_id), "Status": "OK", "Backend": BACKEND_CH21, "Error": "",
        "RunDir": str(result.get("run_dir", "")), "Runtime_s": float(elapsed),
        "PositiveDefinite": bool(result.get("positive_definite", False)),
        "StiffnessConditionNumber": float(result.get("stiffness_condition_number", math.nan)),
        "MacroReactionSymmetryError": float(result.get("macro_reaction_symmetry_error", math.nan)),
        "EnergyReactionRelativeError": float(result.get("energy_vs_reaction_relative_error", math.nan)),
        "StrainAmplitude": float(result.get("strain_amplitude", math.nan)),
    }
    for i in range(6):
        for j in range(i, 6):
            key = _term_name(i, j)
            fe = float(C[i, j]); row[key] = fe
            if analytical_ch is not None:
                aa = float(analytical_ch[i, j]); row[f"Analytical_{key}"] = aa
                row[f"DifferencePct_{key}"] = 100.0 * (fe - aa) / fe if abs(fe) > 1e-14 else math.nan
    for key, value in result.get("engineering_constants", {}).items():
        row[key] = float(value)
    # Plane-stress shell properties comparable to the Chapter-3/FE Q tables.
    for key, value in result.get("Q_thesis_terms", {}).items():
        row[f"Periodic_{key}"] = float(value)
    for key, value in result.get("thesis_FE_prescribed_Q_reference", {}).items():
        row[f"ThesisFE_{key}"] = float(value)
    for key, value in result.get("periodic_Q_vs_thesis_FE_difference_pct", {}).items():
        row[f"PeriodicVsThesisFE_DifferencePct_{key}"] = float(value)
    row["RVEArea_mm2"] = float(result.get("rve_area", math.nan))
    row["RVEVolume_mm3"] = float(result.get("rve_volume", math.nan))
    domain = result.get("periodic_domain", {}) or {}
    audit = result.get("pbc_audit", {}) or {}
    row["PBCDomainMode"] = str(domain.get("domain_mode", ""))
    row["PrimitiveCellCount"] = int(domain.get("primitive_cell_count", 1) or 1)
    row["SupercellWidth_mm"] = float(domain.get("rectangle_width", math.nan))
    row["SupercellHeight_mm"] = float(domain.get("rectangle_height", math.nan))
    row["PBCAuditPass"] = bool(audit.get("passed", False))
    row["PBCUnpairedBoundaryNodes"] = int(audit.get("unpaired_boundary_nodes", 0) or 0)
    row["PBCRigidBodyRank"] = int(audit.get("rigid_body_constraint_rank", 0) or 0)
    return row


def _flatten_structural_result(case_id: str, result: dict, elapsed: float) -> dict:
    rf = list(result.get("top_reaction_force", [math.nan, math.nan, math.nan]))
    uu = list(result.get("max_abs_top_displacement", [math.nan, math.nan, math.nan]))
    return {
        "CaseID": str(case_id), "Status": "OK", "Backend": BACKEND_STRUCTURAL, "Error": "",
        "RunDir": str(result.get("run_dir", "")), "Runtime_s": float(elapsed),
        "LoadCase": str(result.get("load_case", "")),
        "AppliedDisplacement_mm": float(result.get("applied_displacement", math.nan)),
        "ReactionX_N": float(rf[0]) if len(rf) > 0 else math.nan,
        "ReactionY_N": float(rf[1]) if len(rf) > 1 else math.nan,
        "ReactionZ_N": float(rf[2]) if len(rf) > 2 else math.nan,
        "MaxU1_mm": float(uu[0]) if len(uu) > 0 else math.nan,
        "MaxU2_mm": float(uu[1]) if len(uu) > 1 else math.nan,
        "MaxU3_mm": float(uu[2]) if len(uu) > 2 else math.nan,
        "StructuralStiffness_N_per_mm": float(result.get("structural_stiffness_N_per_mm", math.nan)),
        "MaxMises": float(result.get("max_mises", math.nan)),
    }


def _analytical_ch(model, mats, lams) -> np.ndarray:
    return np.asarray(homogenize_rve(model, mats, lams).CH, float)


def _run_periodic_case(
    row: dict,
    *,
    batch_dir: Path,
    abaqus_command: str,
    settings: PeriodicHomogenizationSettings,
    resume: bool,
    compare_analytical: bool,
    case_progress: Callable[[str, int, int, str], None] | None = None,
) -> dict:
    cid = str(row.get("CaseID", "Case")); started = perf_counter()
    try:
        model, mats, lams = build_case_from_row(row)
        model.name = cid
        analytical = _analytical_ch(model, mats, lams) if compare_analytical else None
        files = prepare_periodic_homogenization(
            batch_dir, model, mats, lams, settings,
            run_name=f"Case_{cid}", resume=resume,
        )
        result = run_periodic_homogenization(
            files, abaqus_command, settings,
            (lambda i, n, name: case_progress(cid, i, n, name)) if case_progress else None,
            resume=resume,
        )
        return _flatten_ch21_result(cid, result, analytical, perf_counter() - started)
    except Exception as exc:
        return {
            "CaseID": cid, "Status": "FAILED", "Backend": BACKEND_CH21,
            "Error": str(exc), "RunDir": str(batch_dir / _safe_name(f"Case_{cid}", "Case")),
            "Runtime_s": float(perf_counter() - started),
        }


def _run_structural_case(
    row: dict,
    *,
    batch_dir: Path,
    abaqus_command: str,
    settings: AbaqusSettings,
    resume: bool,
) -> dict:
    cid = str(row.get("CaseID", "Case")); started = perf_counter()
    try:
        model, mats, lams = build_case_from_row(row); model.name = cid
        files = prepare_abaqus_run(batch_dir, model, mats, lams, settings, run_name=f"Case_{cid}", resume=resume)
        result = run_abaqus(files, abaqus_command, settings, resume=resume)
        return _flatten_structural_result(cid, result, perf_counter() - started)
    except Exception as exc:
        return {
            "CaseID": cid, "Status": "FAILED", "Backend": BACKEND_STRUCTURAL,
            "Error": str(exc), "RunDir": str(batch_dir / _safe_name(f"Case_{cid}", "Case")),
            "Runtime_s": float(perf_counter() - started),
        }


def _checkpoint(rows_by_id: dict[str, dict], order: list[str], output_file: Path) -> pd.DataFrame:
    rows = [rows_by_id[cid] for cid in order if cid in rows_by_id]
    df = pd.DataFrame(rows)
    if not df.empty:
        _write_table(df, output_file, sheet="Results")
    return df


def run_abaqus_batch(
    input_file: str | Path,
    output_file: str | Path,
    *,
    backend: str = BACKEND_CH21,
    abaqus_command: str = "abaqus",
    run_root: str | Path | None = None,
    periodic_settings: PeriodicHomogenizationSettings | None = None,
    structural_settings: AbaqusSettings | None = None,
    max_concurrent_cases: int = 1,
    case_limit: int | None = None,
    case_ids: Iterable[str] | None = None,
    resume: bool = True,
    compare_analytical: bool = True,
    progress: Callable[[int, int, str, str], None] | None = None,
) -> pd.DataFrame:
    """Run an Abaqus-backed database batch with restartable per-case checkpoints.

    For CH21, every RVE gets a deterministic case folder. ``resume=True`` reuses
    completed macro-strain states inside that folder, so an interrupted 21-state
    homogenisation continues from the first incomplete state rather than restarting.
    """
    if backend not in {BACKEND_STRUCTURAL, BACKEND_CH21}:
        raise ValueError(f"Unsupported Abaqus batch backend: {backend}")
    if int(max_concurrent_cases) < 1:
        raise ValueError("max_concurrent_cases must be at least 1")
    inp = Path(input_file).expanduser().resolve(); out = Path(output_file).expanduser().resolve(); out.parent.mkdir(parents=True, exist_ok=True)
    if not inp.exists():
        raise FileNotFoundError(inp)
    df = _read_table(inp); active = _active_rows(df, case_limit, case_ids)
    rows = active.to_dict("records")
    if not rows:
        empty = pd.DataFrame(columns=["CaseID", "Status", "Backend", "Error"]); _write_table(empty, out, sheet="Results"); return empty
    root = Path(run_root or (out.parent / "Abaqus_Batch_Runs")).expanduser().resolve(); root.mkdir(parents=True, exist_ok=True)
    suffix = "CH21" if backend == BACKEND_CH21 else "Structural"
    batch_dir = root / _safe_name(f"{inp.stem}_{suffix}", "AbaqusBatch"); batch_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "input_file": str(inp), "output_file": str(out), "backend": backend,
        "abaqus_command": str(abaqus_command), "resume": bool(resume),
        "max_concurrent_cases": int(max_concurrent_cases), "case_ids": [str(r.get("CaseID", "")) for r in rows],
        "periodic_settings": asdict(periodic_settings or PeriodicHomogenizationSettings()) if backend == BACKEND_CH21 else None,
        "structural_settings": asdict(structural_settings or AbaqusSettings()) if backend == BACKEND_STRUCTURAL else None,
        "compare_analytical": bool(compare_analytical),
    }
    (batch_dir / "batch_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    order = [str(r.get("CaseID", f"Case_{i+1}")) for i, r in enumerate(rows)]
    rows_by_id: dict[str, dict] = {}
    if resume and out.exists():
        try:
            old = _read_table(out)
            for rec in old.to_dict("records"):
                cid = str(rec.get("CaseID", ""))
                if cid in order and str(rec.get("Status", "")).upper() == "OK" and str(rec.get("Backend", "")) == backend:
                    rows_by_id[cid] = rec
        except Exception:
            pass

    pending = [r for r in rows if str(r.get("CaseID", "")) not in rows_by_id]
    total = len(rows); completed_count = len(rows_by_id)
    if progress:
        for cid in order:
            if cid in rows_by_id:
                progress(completed_count, total, cid, "SKIPPED - already complete")

    ps = periodic_settings or PeriodicHomogenizationSettings()
    ss = structural_settings or AbaqusSettings()

    def worker(r: dict) -> dict:
        if backend == BACKEND_CH21:
            internal = None
            # Detailed 1/21 progress is most useful/clean when only one Abaqus case runs at once.
            if progress and int(max_concurrent_cases) == 1:
                internal = lambda cid, i, n, name: progress(completed_count, total, cid, f"CH21 {i}/{n}: {name}")
            return _run_periodic_case(r, batch_dir=batch_dir, abaqus_command=abaqus_command, settings=ps,
                                      resume=resume, compare_analytical=compare_analytical, case_progress=internal)
        return _run_structural_case(r, batch_dir=batch_dir, abaqus_command=abaqus_command, settings=ss, resume=resume)

    if int(max_concurrent_cases) == 1:
        for r in pending:
            cid = str(r.get("CaseID", ""))
            if progress: progress(completed_count, total, cid, "STARTED")
            rec = worker(r); rows_by_id[cid] = rec; completed_count += 1
            _checkpoint(rows_by_id, order, out)
            if progress: progress(completed_count, total, cid, str(rec.get("Status", "")))
    else:
        with ThreadPoolExecutor(max_workers=int(max_concurrent_cases)) as ex:
            futs = {ex.submit(worker, r): str(r.get("CaseID", "")) for r in pending}
            for fut in as_completed(futs):
                cid = futs[fut]
                try: rec = fut.result()
                except Exception as exc: rec = {"CaseID": cid, "Status": "FAILED", "Backend": backend, "Error": str(exc)}
                rows_by_id[cid] = rec; completed_count += 1
                _checkpoint(rows_by_id, order, out)
                if progress: progress(completed_count, total, cid, str(rec.get("Status", "")))
    return _checkpoint(rows_by_id, order, out)
