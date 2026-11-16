from __future__ import annotations

from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
import json
import math
import os
import hashlib
import re
import shutil
import subprocess
from typing import Iterable

import numpy as np

from .models import RVEModel, Material, Laminate


LOAD_CASES = {
    "XZ shear": {"component": 1, "sign": 1.0, "description": "Top-edge translation in global X; bottom edge fixed."},
    "YZ shear": {"component": 2, "sign": 1.0, "description": "Top-edge translation in global Y; bottom edge fixed."},
    "Z compression": {"component": 3, "sign": -1.0, "description": "Top-edge translation in global Z; bottom edge fixed."},
}


@dataclass(slots=True)
class AbaqusSettings:
    element_type: str = "S4"
    mesh_size: float = 5.0
    displacement: float = 1.0
    load_case: str = "XZ shear"
    cpus: int = 1
    save_odb: bool = True

    def validate(self) -> None:
        if self.element_type.upper() != "S4":
            raise ValueError("This adapter currently supports the Abaqus S4 shell element only.")
        if not math.isfinite(self.mesh_size) or self.mesh_size <= 0:
            raise ValueError("Mesh size must be positive and finite.")
        if not math.isfinite(self.displacement) or self.displacement == 0:
            raise ValueError("Applied displacement must be finite and non-zero.")
        if self.load_case not in LOAD_CASES:
            raise ValueError(f"Unknown Abaqus load case: {self.load_case}")
        if int(self.cpus) < 1:
            raise ValueError("CPU count must be at least 1.")


@dataclass(slots=True)
class AbaqusRunFiles:
    run_dir: Path
    job_name: str
    input_file: Path
    extractor_file: Path
    snapshot_file: Path
    results_file: Path
    log_file: Path


def _safe_name(text: str, prefix: str = "X") -> str:
    s = re.sub(r"[^A-Za-z0-9_]", "_", str(text)).strip("_")
    if not s:
        s = prefix
    if not s[0].isalpha():
        s = f"{prefix}_{s}"
    return s[:70]


def detect_abaqus_executable() -> str | None:
    """Return a likely Abaqus command if one is on PATH."""
    candidates = ["abaqus", "abq2026", "abq2025", "abq2024", "abq2023", "abq2022", "abq2021", "abq2020"]
    for name in candidates:
        p = shutil.which(name)
        if p:
            return p
    return None


def _ids_lines(ids: Iterable[int], width: int = 16) -> str:
    vals = [str(int(x)) for x in ids]
    return "\n".join(", ".join(vals[i:i + width]) for i in range(0, len(vals), width))


def _abaqus_material_blocks(materials: dict[str, Material]) -> tuple[str, dict[str, str]]:
    chunks: list[str] = []
    names: dict[str, str] = {}
    used: set[str] = set()
    for i, (key, mat) in enumerate(materials.items(), 1):
        base = _safe_name(mat.name or key, "MAT")
        name = base
        k = 2
        while name.upper() in used:
            name = f"{base}_{k}"; k += 1
        used.add(name.upper()); names[key] = name
        g13 = float(mat.Gth if mat.Gth is not None else mat.G12)
        g23 = g13
        chunks.extend([
            f"*MATERIAL, NAME={name}",
            "*ELASTIC, TYPE=LAMINA",
            f"{mat.E1:.12g}, {mat.E2:.12g}, {mat.nu12:.12g}, {mat.G12:.12g}, {g13:.12g}, {g23:.12g}",
        ])
        if math.isfinite(float(mat.density_kg_m3)) and float(mat.density_kg_m3) > 0:
            # N-mm-s consistent mass unit is tonne; 1 kg/m^3 = 1e-12 tonne/mm^3.
            chunks.extend(["*DENSITY", f"{float(mat.density_kg_m3) * 1e-12:.12g}"])
    return "\n".join(chunks), names


def _divisions(length: float, mesh_size: float) -> int:
    """Number of elements along a length, robust to round-off.

    A 50 mm wall that comes out of a geometric cut as 50.00000000000003 mm must get
    the same 5 elements (10 mm mesh) as its exact periodic copy, or the copies on
    opposite RVE faces have unmatched nodes.
    """
    return max(1, int(math.ceil(float(length) / float(mesh_size) - 1e-9)))


def _build_shell_mesh(model: RVEModel, mesh_size: float):
    """Build a conforming S4 strip mesh for the vertical cell walls.

    Original 2-D RVE nodes become shared vertical node columns. Interior nodes are
    wall-local. This makes shell walls meeting at an RVE node share the complete
    vertical junction line.
    """
    nodes: dict[int, tuple[float, float, float]] = {}
    wall_elements: dict[int, list[tuple[int, tuple[int, int, int, int]]]] = {}
    top_ids: set[int] = set(); bottom_ids: set[int] = set()
    node_id = 1; elem_id = 1
    nz = _divisions(float(model.height), mesh_size)

    # Only create columns for nodes that are actually referenced by a wall.
    # Some analytical geometry definitions carry auxiliary lattice-corner nodes
    # that are intentionally not physical wall nodes; exporting those as free FE
    # nodes creates singular/disconnected Abaqus models.
    used_model_nodes = sorted({int(e.node_i) for e in model.elements} | {int(e.node_j) for e in model.elements})
    shared_columns: dict[int, list[int]] = {}
    for idx in used_model_nodes:
        x, y = np.asarray(model.nodes[idx], float)
        col = []
        for kz in range(nz + 1):
            z = float(model.height) * kz / nz
            nodes[node_id] = (float(x), float(y), z); col.append(node_id)
            if kz == 0: bottom_ids.add(node_id)
            if kz == nz: top_ids.add(node_id)
            node_id += 1
        shared_columns[idx] = col

    for wi, e in enumerate(model.elements, 1):
        p0 = np.asarray(model.nodes[e.node_i], float); p1 = np.asarray(model.nodes[e.node_j], float)
        length = float(np.linalg.norm(p1 - p0))
        if length <= np.finfo(float).eps:
            raise ValueError(f"Wall W{wi} has zero length.")
        ns = _divisions(length, mesh_size)
        columns: list[list[int]] = [shared_columns[e.node_i]]
        for ks in range(1, ns):
            t = ks / ns; xy = (1.0 - t) * p0 + t * p1; col = []
            for kz in range(nz + 1):
                z = float(model.height) * kz / nz
                nodes[node_id] = (float(xy[0]), float(xy[1]), z); col.append(node_id)
                if kz == 0: bottom_ids.add(node_id)
                if kz == nz: top_ids.add(node_id)
                node_id += 1
            columns.append(col)
        columns.append(shared_columns[e.node_j])
        elems: list[tuple[int, tuple[int, int, int, int]]] = []
        for ks in range(ns):
            for kz in range(nz):
                n1 = columns[ks][kz]
                n2 = columns[ks + 1][kz]
                n3 = columns[ks + 1][kz + 1]
                n4 = columns[ks][kz + 1]
                elems.append((elem_id, (n1, n2, n3, n4))); elem_id += 1
        wall_elements[wi] = elems
    return nodes, wall_elements, sorted(bottom_ids), sorted(top_ids)


def generate_abaqus_input(model: RVEModel, materials: dict[str, Material], laminates: list[Laminate], settings: AbaqusSettings, *, job_name: str = "RVE_S4") -> str:
    """Generate an Abaqus/Standard input deck using conforming S4 shell walls.

    The current FE adapter performs a displacement-controlled structural check
    (XZ shear, YZ shear, or Z compression). It intentionally does not claim to
    reproduce the analytical periodic homogenisation boundary conditions.
    """
    settings.validate()
    if not laminates:
        raise ValueError("At least one laminate stack is required for Abaqus export.")
    nodes, wall_elements, bottom_ids, top_ids = _build_shell_mesh(model, settings.mesh_size)
    material_text, mat_names = _abaqus_material_blocks(materials)

    out: list[str] = [
        "** RVE Studio Python - Abaqus S4 shell export",
        f"** Job: {job_name}",
        "** Units: N, mm, MPa (density converted to tonne/mm^3 when available)",
        "** 0-degree ply direction is the local wall direction.",
        "*HEADING",
        f"RVE Studio S4 model - {model.name}",
        "*NODE",
    ]
    for nid, xyz in sorted(nodes.items()):
        out.append(f"{nid}, {xyz[0]:.12g}, {xyz[1]:.12g}, {xyz[2]:.12g}")
    out.append("*ELEMENT, TYPE=S4")
    for wi in sorted(wall_elements):
        for eid, conn in wall_elements[wi]:
            out.append(f"{eid}, {conn[0]}, {conn[1]}, {conn[2]}, {conn[3]}")
    out.extend(["*NSET, NSET=BOTTOM", _ids_lines(bottom_ids), "*NSET, NSET=TOP", _ids_lines(top_ids)])

    for wi, elems in wall_elements.items():
        out.extend([f"*ELSET, ELSET=WALL_{wi}", _ids_lines(eid for eid, _ in elems)])

    out.append(material_text)

    for wi, e in enumerate(model.elements, 1):
        if e.laminate < 0 or e.laminate >= len(laminates):
            raise ValueError(f"Wall W{wi} references missing laminate index {e.laminate}.")
        lam = laminates[e.laminate]
        p0 = np.asarray(model.nodes[e.node_i], float); p1 = np.asarray(model.nodes[e.node_j], float)
        d = p1 - p0; d /= np.linalg.norm(d)
        # Local 1 axis = wall direction; local 2 axis = core-height/global Z.
        out.extend([
            f"*ORIENTATION, NAME=ORI_W{wi}, SYSTEM=RECTANGULAR",
            f"{d[0]:.12g}, {d[1]:.12g}, 0., 0., 0., 1.",
            "3, 0.",  # project the ply 1-axis about the wall normal (local 3): ply 1 = along the wall, not the height
            f"*SHELL SECTION, ELSET=WALL_{wi}, COMPOSITE, ORIENTATION=ORI_W{wi}",
        ])
        for pi, ply in enumerate(lam.plies, 1):
            if ply.material not in mat_names:
                raise ValueError(f"Wall W{wi} ply {pi} references material '{ply.material}', which is not in the material library.")
            out.append(f"{float(ply.thickness):.12g}, 3, {mat_names[ply.material]}, {float(ply.angle_deg):.12g}, PLY{pi}")

    lc = LOAD_CASES[settings.load_case]
    comp = int(lc["component"]); disp = float(settings.displacement) * float(lc["sign"])
    out.extend([
        "*STEP, NAME=LOAD, NLGEOM=NO",
        "*STATIC",
        "0.1, 1.0, 1e-06, 0.1",
        "*BOUNDARY",
        "BOTTOM, 1, 6, 0.",
        f"TOP, {comp}, {comp}, {disp:.12g}",
        "*OUTPUT, FIELD, VARIABLE=PRESELECT",
        "*NODE OUTPUT, NSET=TOP",
        "U, RF",
        "*ELEMENT OUTPUT",
        "S, E",
        "*OUTPUT, HISTORY, VARIABLE=PRESELECT",
        "*END STEP",
    ])
    return "\n".join(out) + "\n"


def _extractor_source() -> str:
    """Return an Abaqus-Python ODB extractor compatible with older releases.

    Abaqus ODB repositories are dict-like but older releases do not guarantee
    support for ``Repository.get``. Use indexed access and make every failure
    serialisable so the desktop application can report the real post-processing
    exception instead of a misleading missing-results message.
    """
    return r"""from __future__ import print_function
import json, os, sys, traceback
from odbAccess import openOdb

odb = None
out_path = 'results.json'
if len(sys.argv) >= 3:
    odb_path = sys.argv[-2]
    out_path = sys.argv[-1]
elif len(sys.argv) == 2:
    odb_path = sys.argv[-1]
else:
    odbs = [x for x in os.listdir('.') if x.lower().endswith('.odb')]
    if len(odbs) != 1:
        raise RuntimeError('ODB extractor requires an ODB path (or exactly one .odb in the run directory)')
    odb_path = odbs[0]

def repo_get(repo, name):
    # Abaqus Repository objects are dict-like, but .get() is not portable across releases.
    try:
        return repo[name]
    except Exception:
        pass
    try:
        target = str(name).upper()
        for key in repo.keys():
            if str(key).upper() == target:
                return repo[key]
    except Exception:
        pass
    return None

def write_json(data):
    with open(out_path, 'w') as f:
        json.dump(data, f, indent=2)

try:
    odb = openOdb(path=odb_path, readOnly=True)
    step = odb.steps['LOAD']
    if not step.frames:
        raise RuntimeError('LOAD step contains no output frames')
    frame = step.frames[-1]
    assembly = odb.rootAssembly

    top = repo_get(assembly.nodeSets, 'TOP')
    if top is None:
        for inst in assembly.instances.values():
            top = repo_get(inst.nodeSets, 'TOP')
            if top is not None:
                break
    if top is None:
        try:
            available = list(assembly.nodeSets.keys())
        except Exception:
            available = []
        raise RuntimeError('TOP node set was not found in the ODB. Assembly node sets: %s' % available)

    if 'RF' not in frame.fieldOutputs.keys():
        raise RuntimeError('RF field output is missing from the final LOAD frame')
    if 'U' not in frame.fieldOutputs.keys():
        raise RuntimeError('U field output is missing from the final LOAD frame')

    rf = frame.fieldOutputs['RF'].getSubset(region=top)
    u = frame.fieldOutputs['U'].getSubset(region=top)
    reaction = [0.0, 0.0, 0.0]
    max_u = [0.0, 0.0, 0.0]
    for v in rf.values:
        data = v.data
        for i in range(min(3, len(data))):
            reaction[i] += float(data[i])
    for v in u.values:
        data = v.data
        for i in range(min(3, len(data))):
            val = abs(float(data[i]))
            if val > max_u[i]:
                max_u[i] = val

    result = {
        'status': 'COMPLETED',
        'odb': odb_path,
        'step': 'LOAD',
        'top_reaction_force': reaction,
        'max_abs_top_displacement': max_u,
        'frames': len(step.frames),
    }
    try:
        s = frame.fieldOutputs['S']
        mises = []
        for v in s.values:
            try:
                mises.append(float(v.mises))
            except Exception:
                pass
        if mises:
            result['max_mises'] = max(mises)
    except Exception:
        pass
    write_json(result)
    print(json.dumps(result, indent=2))
except Exception as exc:
    error = {
        'status': 'ERROR',
        'odb': locals().get('odb_path', ''),
        'error': str(exc),
        'traceback': traceback.format_exc(),
    }
    try:
        write_json(error)
    except Exception:
        pass
    print(json.dumps(error, indent=2))
    sys.exit(2)
finally:
    if odb is not None:
        try:
            odb.close()
        except Exception:
            pass
"""


def _structural_input_signature(model: RVEModel, materials: dict[str, Material], laminates: list[Laminate], settings: AbaqusSettings) -> str:
    payload = {
        "model": {"nodes": np.asarray(model.nodes, float).tolist(), "elements": [asdict(e) for e in model.elements],
                  "a1": np.asarray(model.a1, float).tolist(), "a2": np.asarray(model.a2, float).tolist(), "height": float(model.height)},
        "materials": {str(k): asdict(v) for k, v in sorted(materials.items(), key=lambda kv: str(kv[0]))},
        "laminates": [{"name": l.name, "plies": [asdict(p) for p in l.plies]} for l in laminates],
        "settings": asdict(settings),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=True).encode("utf-8")).hexdigest()


def _completed_result(path: Path) -> dict | None:
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
        return raw if str(raw.get("status", "")).upper() == "COMPLETED" else None
    except Exception:
        return None


def prepare_abaqus_run(run_root: str | Path, model: RVEModel, materials: dict[str, Material], laminates: list[Laminate], settings: AbaqusSettings, *, run_name: str | None = None, resume: bool = False) -> AbaqusRunFiles:
    settings.validate()
    root = Path(run_root).expanduser().resolve(); root.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    job = _safe_name(run_name if run_name else f"{model.name}_{settings.load_case}_{stamp}", "RVE")
    run_dir = root / job
    signature = _structural_input_signature(model, materials, laminates, settings)
    if run_dir.exists():
        if not resume:
            raise FileExistsError(f"Abaqus structural run folder already exists: {run_dir}")
        old_file = run_dir / "model_snapshot.json"
        if not old_file.exists():
            raise RuntimeError(f"Cannot resume {run_dir}: model_snapshot.json is missing.")
        old = json.loads(old_file.read_text(encoding="utf-8"))
        if old.get("input_signature") and old.get("input_signature") != signature:
            raise RuntimeError(f"Cannot resume {run_dir}: input model/settings changed.")
    else:
        run_dir.mkdir(parents=False, exist_ok=False)
    files = AbaqusRunFiles(run_dir, job, run_dir / f"{job}.inp", run_dir / "extract_odb.py", run_dir / "model_snapshot.json", run_dir / "results.json", run_dir / "abaqus_run.log")
    if not (resume and files.input_file.exists()): files.input_file.write_text(generate_abaqus_input(model, materials, laminates, settings, job_name=job), encoding="utf-8")
    if not (resume and files.extractor_file.exists()): files.extractor_file.write_text(_extractor_source(), encoding="utf-8")
    snapshot = {
        "input_signature": signature,
        "model": {
            "name": model.name, "nodes": model.nodes.tolist(), "elements": [asdict(e) for e in model.elements],
            "a1": model.a1.tolist(), "a2": model.a2.tolist(), "height": model.height, "wall_groups": list(model.wall_groups),
        },
        "materials": {k: asdict(v) for k, v in materials.items()},
        "laminates": [{"name": l.name, "plies": [asdict(p) for p in l.plies]} for l in laminates],
        "abaqus_settings": asdict(settings),
        "note": "Abaqus S4 structural verification model. This release does not claim periodic FE homogenisation equivalence with the analytical solver.",
    }
    files.snapshot_file.write_text(json.dumps(snapshot, indent=2), encoding="utf-8")
    return files


def _run_shell_command(parts: list[str], cwd: Path, log_handle) -> int:
    if os.name == "nt":
        cmd = subprocess.list2cmdline(parts)
        p = subprocess.run(cmd, cwd=str(cwd), stdout=log_handle, stderr=subprocess.STDOUT, shell=True)
    else:
        p = subprocess.run(parts, cwd=str(cwd), stdout=log_handle, stderr=subprocess.STDOUT)
    return int(p.returncode)


def run_abaqus(files: AbaqusRunFiles, abaqus_command: str, settings: AbaqusSettings, *, resume: bool = False) -> dict:
    """Submit the prepared Abaqus job and extract ODB results.

    `abaqus_command` can be `abaqus`, an abq20xx command, or the full path to the
    SIMULIA command/batch file.
    """
    exe = str(abaqus_command).strip()
    if not exe:
        raise ValueError("Abaqus executable/command is empty.")
    if resume:
        ready = _completed_result(files.results_file)
        if ready is not None:
            return ready
    with files.log_file.open("a" if resume else "w", encoding="utf-8", errors="replace") as log:
        log.write(f"RVE Studio Abaqus run\nJob: {files.job_name}\nCommand: {exe}\n\n")
        rc = _run_shell_command([exe, f"job={files.job_name}", f"input={files.input_file.name}", f"cpus={int(settings.cpus)}", "interactive"], files.run_dir, log)
        if rc != 0:
            raise RuntimeError(f"Abaqus job returned exit code {rc}. See {files.log_file}")
        odb = files.run_dir / f"{files.job_name}.odb"
        if not odb.exists():
            raise RuntimeError(f"Abaqus completed without producing {odb.name}. See {files.log_file}")
        log.write("\n--- ODB extraction ---\n")
        log.flush()
        rc = _run_shell_command([exe, "python", files.extractor_file.name, odb.name, files.results_file.name], files.run_dir, log)
        log.flush()
        if not files.results_file.exists():
            try:
                tail = files.log_file.read_text(encoding="utf-8", errors="replace")[-5000:]
            except Exception:
                tail = ""
            raise RuntimeError(
                f"ODB result extraction produced no results file (Abaqus wrapper exit code {rc}). "
                f"See {files.log_file}.\n\nLast extractor output:\n{tail}"
            )
    result = json.loads(files.results_file.read_text(encoding="utf-8"))
    if str(result.get("status", "")).upper() != "COMPLETED":
        detail = result.get("traceback") or result.get("error") or "Unknown ODB extraction error"
        raise RuntimeError(f"ODB result extraction failed:\n{detail}\nRun log: {files.log_file}")
    lc = LOAD_CASES[settings.load_case]
    comp = int(lc["component"]) - 1
    rf = result.get("top_reaction_force", [math.nan, math.nan, math.nan])
    result["load_case"] = settings.load_case
    result["applied_displacement"] = float(settings.displacement) * float(lc["sign"])
    try:
        result["reaction_component"] = float(rf[comp])
        result["structural_stiffness_N_per_mm"] = abs(float(rf[comp]) / float(settings.displacement))
    except Exception:
        pass
    result["job_name"] = files.job_name
    result["run_dir"] = str(files.run_dir)
    files.results_file.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result
