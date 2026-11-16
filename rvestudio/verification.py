"""FE verification core: what is compared with what, and the plain-language findings.

Every FE model uses one FE RVE (``fe_rve.select_fe_rve``) and one of three
boundary-value problems:

=====================  ==================  =========================  =========================================
Model                  Side faces          Top and bottom faces       Terms
=====================  ==================  =========================  =========================================
PBC (strain energy)    periodic            periodic: U3+ - U3- = e33 h full 6x6 C (21 terms), full Q matrix
                                                                      (sigma33 = 0, 15 terms), engineering
                                                                      constants E1-E3, G12, G13, G23, nu_ij
DBC, shell terms       thesis Table 3.4    U3 free (in-plane cases)   Q11 Q22 Q12 Q66 Q44 Q55   (sigma33 = 0)
DBC, 3D terms          the same            U3 = 0                     C1111 C2222 C1122 C1212   (eps33 = 0)
=====================  ==================  =========================  =========================================

PBC: six unit macro strains on the periodic RVE; the stiffness follows from the strain
energy, C_ij = u_i^T K u_j / V (built-in) or C_ii = 2 W_i / (V e^2), C_ij = (W_ij - W_i - W_j) / (V e^2)
from the 21 single and paired strain states (Abaqus ALLSE).  Q is C condensed to sigma33 = 0
(Schur complement), order [11, 22, 12, 13, 23] with the paper's names Q44 = 13, Q55 = 23;
the engineering constants come from S = C^-1.  Every FE value is compared with the same
quantity of the analytical model (its C^H); no FE values from the literature are used.

Results are long-format rows (one per project, solver, model, mesh and term, with a unit),
the same for the built-in shell solver and for Abaqus, so the GUI, the reports and the
CLI read one structure.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import copy
import math
import time
from typing import Callable, Iterable

import numpy as np

from .models import RVEModel, Material, Laminate

SHELL_TERMS = ("Q11", "Q22", "Q12", "Q66", "Q44", "Q55")
SOLID_TERMS = ("C1111", "C2222", "C1122", "C1212", "C1313", "C2323")
# Full 6x6 stiffness, Voigt [11, 22, 33, 23, 13, 12], named by tensor indices (C1133, C2313, ...):
# unambiguous, unlike Voigt C44 (= C2323), which the paper's Q44 (= C1313) would contradict.
VOIGT_PAIRS = ("11", "22", "33", "23", "13", "12")
FULL_C_TERMS = tuple(f"C{VOIGT_PAIRS[i]}{VOIGT_PAIRS[j]}" for i in range(6) for j in range(i, 6))
# Full Q (sigma33 = 0) in the paper's numbering 1 = 11, 2 = 22, 6 = 12, 4 = 13, 5 = 23; order [11, 22, 12, 13, 23].
Q_ORDER = ("1", "2", "6", "4", "5")
FULL_Q_TERMS = ("Q11", "Q22", "Q12", "Q66", "Q44", "Q55", "Q16", "Q26", "Q45", "Q14", "Q15", "Q24", "Q25", "Q46", "Q56")
ENG_TERMS = ("E1", "E2", "E3", "G12", "G13", "G23", "nu12", "nu13", "nu23", "nu21", "nu31", "nu32")
GROUPS = {"C": FULL_C_TERMS, "Q": FULL_Q_TERMS, "eng": ENG_TERMS, "shell": SHELL_TERMS, "3d": SOLID_TERMS}
GROUP_TITLES = {"C": "Full stiffness matrix C (6×6, strain energy)",
                "Q": "Q matrix (σ33 = 0, from C)",
                "eng": "Engineering constants (from S = C⁻¹)",
                "shell": "DBC shell terms (σ33 = 0)", "3d": "DBC 3D stiffness terms (ε33 = 0)"}
GROUP_OF = {t: g for g, ts in GROUPS.items() for t in ts}
MODEL_GROUPS = {"PBC": ("C", "Q", "eng"), "DBC": ("shell", "3d")}
# Key of each group in the analytical reference and in an FE result dictionary.
REF_KEY = {"C": "C21", "Q": "Q15", "eng": "eng", "shell": "Q", "3d": "C"}
NEGLIGIBLE = 1e-6  # |term| below this fraction of the largest matrix entry: "≈ 0", no % difference
SOLVERS = ("built-in", "Abaqus")
BOUNDARY_CONDITIONS = [
    ("PBC (strain energy)", "periodic", "periodic through the height: U3⁺ − U3⁻ = ε33 h",
     "full 6×6 C (21 terms), Q (σ33 = 0), engineering constants"),
    ("DBC, shell terms", "prescribed (thesis Table 3.4)", "U3 free (in-plane cases); U1 = U2 = 0 (Q44, Q55)",
     "Q11, Q22, Q12, Q66, Q44, Q55"),
    ("DBC, 3D terms", "the same", "U3 = 0 on both faces", "C1111, C2222, C1122, C1212 (C1313, C2323 = Q44, Q55 cases)"),
]


def term_unit(term: str) -> str:
    return "-" if term.startswith("nu") else "MPa"


def full_c_terms(C) -> dict[str, float]:
    """The 21 independent terms of a symmetric 6x6 Voigt [11,22,33,23,13,12] stiffness."""
    C = np.asarray(C, float)
    return {f"C{VOIGT_PAIRS[i]}{VOIGT_PAIRS[j]}": float(C[i, j]) for i in range(6) for j in range(i, 6)}


def q_matrix(C) -> np.ndarray:
    """C condensed to sigma33 = 0 (Schur complement), 5x5 in the order [11, 22, 12, 13, 23]."""
    from .abaqus_homogenization import plane_stress_shell_reduction
    return plane_stress_shell_reduction(np.asarray(C, float))[0]


def q_terms(Q5) -> dict[str, float]:
    """Named terms of the 5x5 Q in the paper's numbering (Q44 = 13 shear, Q55 = 23 shear)."""
    Q5 = np.asarray(Q5, float)
    out = {}
    for i in range(5):
        for j in range(i, 5):
            a, b = sorted((Q_ORDER[i], Q_ORDER[j]))
            out[f"Q{a}{b}"] = float(Q5[i, j])
    return {t: out[t] for t in FULL_Q_TERMS}


def engineering_terms(C) -> dict[str, float]:
    """E1-E3, G12, G13, G23 and the six Poisson ratios from S = C^-1."""
    from .abaqus_homogenization import compliance_and_engineering_constants
    _S, eng = compliance_and_engineering_constants(np.asarray(C, float))
    return {t: float(eng[t]) for t in ENG_TERMS}


def pbc_values(C) -> dict:
    """Everything the PBC (strain-energy) model reports, from its 6x6 stiffness."""
    C = 0.5 * (np.asarray(C, float) + np.asarray(C, float).T)
    Q5 = q_matrix(C)
    return {"CH": C.tolist(), "Qmatrix": Q5.tolist(), "C21": full_c_terms(C), "Q15": q_terms(Q5),
            "eng": engineering_terms(C), "C": solid_terms(C)}


def _scale(group: str, values: dict) -> float:
    vals = [abs(float(v)) for v in values.values() if math.isfinite(_num(v))]
    return max(vals) if vals else 0.0


def _num(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def negligible(group: str, term: str, analytical_values: dict, fe_value=None) -> bool:
    """True for a matrix term that is zero to round-off in the analytical model and the FE."""
    if group not in ("C", "Q"):
        return False
    sc = _scale(group, analytical_values)
    if sc <= 0:
        return False
    small = abs(_num(analytical_values.get(term))) < NEGLIGIBLE * sc
    if fe_value is not None and math.isfinite(_num(fe_value)):
        small = small and abs(_num(fe_value)) < 1e-3 * sc
    return small
WALL_OPTIONS = {
    "cylindrical": ("paper", "cylindrical"),
    "anisotropic strip": ("anisotropicStrip", "cylindrical"),
    "free moment": ("anisotropicStrip", "freeMoment"),
    "twist-corrected": ("anisotropicStrip", "twistRelease"),
}


def signed(x: float, decimals: int = 2) -> str:
    """Signed fixed-point text with round-off shown as +0.00 (never -0.00)."""
    x = float(x)
    if not math.isfinite(x):
        return "–"
    if abs(x) < 0.5 * 10.0 ** (-decimals):
        x = 0.0
    return f"{x:+.{decimals}f}"


def pct(a: float, b: float) -> float:
    """100 (a - b) / b, NaN when undefined."""
    try:
        a, b = float(a), float(b)
    except (TypeError, ValueError):
        return math.nan
    return 100.0 * (a - b) / b if b and math.isfinite(a) and math.isfinite(b) else math.nan


def solid_terms(C) -> dict[str, float]:
    """C1111, C2222, C1122, C1212, C1313, C2323 of a 6x6 Voigt [11,22,33,23,13,12] matrix."""
    C = np.asarray(C, float)
    return {"C1111": float(C[0, 0]), "C2222": float(C[1, 1]), "C1122": float(C[0, 1]),
            "C1212": float(C[5, 5]), "C1313": float(C[4, 4]), "C2323": float(C[3, 3])}


# ------------------------------------------------------------------------------------ references
def _active_option(model: RVEModel) -> str:
    et = str(model.wall_element_type).lower().replace("_", "")
    if et in ("anisotropicstrip", "anisotropic", "fullanisotropic"):
        mode = str(model.strip_bending_mode).lower().replace("_", "")
        if mode in ("twistrelease", "twistcorrected", "twist"):
            return "twist-corrected"
        return "free moment" if mode in ("freemoment", "free", "anisotropic") else "anisotropic strip"
    return "cylindrical"


def _beam(model: RVEModel) -> str:
    from .walls import is_timoshenko
    return "timoshenko" if is_timoshenko(getattr(model, "wall_shear", "euler-bernoulli")) else "euler-bernoulli"


def analytical_reference(model: RVEModel, materials: dict, laminates: list) -> dict:
    """Analytical Q and C for the model's wall option, the other two bending options and the other
    beam theory (Euler-Bernoulli = the original model; Timoshenko = with wall shear flexibility)."""
    from .homogenization import homogenize_rve
    from .validation import extract_six
    active = _active_option(model)
    beam = _beam(model)
    other = "timoshenko" if beam == "euler-bernoulli" else "euler-bernoulli"
    label = {"euler-bernoulli": "Euler\u2013Bernoulli", "timoshenko": "Timoshenko"}
    runs = [(name, etype, mode, beam) for name, (etype, mode) in WALL_OPTIONS.items()]
    etype, mode = WALL_OPTIONS[active]
    runs.append((f"{active}, {label[other]}", etype, mode, other))
    options: dict[str, dict] = {}
    base = None
    for name, etype, mode, shear in runs:
        m = copy.copy(model)
        m.wall_element_type = etype
        m.strip_bending_mode = mode
        m.wall_shear = shear
        try:
            res = homogenize_rve(m, materials, laminates)
        except Exception as exc:  # an option can be undefined for a stack (e.g. singular D)
            options[name] = {"error": str(exc)}
            continue
        CH = np.asarray(res.CH, float)
        options[name] = {"Q": dict(zip(SHELL_TERMS, (float(x) for x in extract_six(res)))),
                         **pbc_values(CH), "beam": shear}
        if name == active:
            base = res
    if base is None:
        base = homogenize_rve(model, materials, laminates)
        CH = np.asarray(base.CH, float)
        options[active] = {"Q": dict(zip(SHELL_TERMS, (float(x) for x in extract_six(base)))),
                           **pbc_values(CH), "beam": beam}
    act = options[active]
    return {"active": active, "beam": beam, "options": options, "Q": act["Q"], "C": act["C"], "CH": act["CH"],
            "C21": act["C21"], "Q15": act["Q15"], "Qmatrix": act["Qmatrix"], "eng": act["eng"],
            "relative_density": float(base.relative_density),
            "density_kg_m3": float(base.physical.average_density_kg_m3)}


def paper_configuration_of(model: RVEModel, materials: dict[str, Material], laminates: list[Laminate]) -> int | None:
    """Number (1-14) of the benchmark configuration this RVE is identical to, else None.

    Compared by geometry, lattice, height and ply properties (not by names), so a
    configuration loaded into the designer is still labelled as such.
    """
    from .examples_data import paper_hexagonal_case
    X = np.asarray(model.nodes, float)
    conn = sorted(tuple(sorted((e.node_i, e.node_j))) for e in model.elements)

    def plies(lams, mats):
        out = []
        for lam in lams:
            for p in lam.plies:
                m = mats.get(p.material)
                if m is None:
                    return None
                out.append((round(float(p.angle_deg), 6), round(float(p.thickness), 6), round(float(m.E1), 3),
                            round(float(m.E2), 3), round(float(m.G12), 3), round(float(m.nu12), 6)))
        return out
    mine = plies(laminates, materials)
    if mine is None:
        return None
    for n in range(1, 15):
        m2, mats2, lams2 = paper_hexagonal_case(n)
        if (np.asarray(m2.nodes).shape != X.shape or not np.allclose(m2.nodes, X, atol=1e-6)
                or not np.allclose(m2.a1, model.a1, atol=1e-6) or not np.allclose(m2.a2, model.a2, atol=1e-6)
                or abs(float(m2.height) - float(model.height)) > 1e-9
                or sorted(tuple(sorted((e.node_i, e.node_j))) for e in m2.elements) != conn):
            continue
        if plies(lams2, mats2) == mine:
            return n
    return None


# ------------------------------------------------------------------------------------ built-in FE
def builtin_check(model: RVEModel, materials: dict, laminates: list, mesh_size: float, *,
                  models: Iterable[str] = ("DBC", "PBC"), displacement: float = 1.0,
                  face_wall_section: str = "half-thickness", fields: bool = False) -> dict:
    """Built-in shell FE of the PBC (strain energy: full C, Q, engineering constants) and DBC models."""
    from .shell_fe import run_dbc_shell_check, run_pbc_shell_check
    if not math.isfinite(float(mesh_size)) or float(mesh_size) < 2.5:
        raise ValueError("The built-in shell check needs a mesh size of at least 2.5 mm.")
    models = [m.upper() for m in models]
    out: dict = {"mesh_mm": float(mesh_size)}
    t0 = time.time()
    if "DBC" in models:
        d = run_dbc_shell_check(model, materials, laminates, mesh_size=mesh_size, displacement=displacement,
                                face_wall_section=face_wall_section, return_fields=fields)
        out["DBC"] = {"Q": d.get("Q"), "C": d.get("C"), "energies_Nmm": d["energies_Nmm"],
                      "nodes": d["nodes"], "elements": d["elements"]}
        if fields:
            out.setdefault("fields", {})["DBC"] = d["fields"]
    if "PBC" in models:
        r = run_pbc_shell_check(model, materials, laminates, mesh_size=mesh_size, return_fields=fields)
        out["PBC"] = {**pbc_values(r["CH"]), "symmetry_error": r["symmetry_error"],
                      "nodes": r["nodes"], "elements": r["elements"]}
        if fields:
            out.setdefault("fields", {})["PBC"] = r["fields"]
    out["seconds"] = round(time.time() - t0, 2)
    return out


def result_rows(project: str, solver: str, mesh_mm: float, fe: dict, analytical: dict) -> list[dict]:
    """Long-format rows for one project, solver and mesh: every FE term next to the analytical value.

    ``fe`` maps "PBC" and/or "DBC" to dictionaries holding the groups under :data:`REF_KEY`
    ("C21", "Q15", "eng" for PBC; "Q", "C" for DBC).  Matrix terms that are zero to round-off
    in both are marked ``negligible`` and get no % difference.
    """
    rows = []
    for model, groups in MODEL_GROUPS.items():
        res = fe.get(model) or {}
        for g in groups:
            vals = res.get(REF_KEY[g])
            if not vals:
                continue
            ref = analytical[REF_KEY[g]]
            for t in GROUPS[g]:
                if t not in vals:
                    continue
                v = _num(vals[t])
                neg = negligible(g, t, ref, v)
                rows.append({"project": project, "solver": solver, "model": model, "mesh_mm": float(mesh_mm),
                             "group": g, "term": t, "unit": term_unit(t), "FE": v, "analytical": float(ref[t]),
                             "analytical_vs_FE_pct": math.nan if neg else pct(ref[t], v), "negligible": bool(neg)})
    return rows


# ------------------------------------------------------------------------------------ summaries
def _f(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def best_values(rows: list[dict], project: str, solver: str, model: str, group: str,
                convergence: list[dict] | None = None) -> dict[str, tuple[float, str]]:
    """Value to quote per term: the Richardson value when three meshes converge, else the finest mesh."""
    out = {}
    sel = [r for r in rows if r["project"] == project and r["solver"] == solver and r["model"] == model
           and r["group"] == group]
    for t in GROUPS[group]:
        pts = sorted(((_f(r["mesh_mm"]), _f(r["FE"])) for r in sel if r["term"] == t), key=lambda p: p[0])
        if not pts:
            continue
        val, how = pts[0][1], f"{pts[0][0]:g} mm"
        for c in convergence or []:
            if (c.get("project"), c.get("solver"), c.get("model"), c.get("term")) == (project, solver, model, t) \
                    and c.get("status") == "monotone" and math.isfinite(_f(c.get("extrapolated"))):
                val, how = _f(c["extrapolated"]), "extrapolated"
        out[t] = (val, how)
    return out


def _alt_option(analytical: dict):
    """The analytical result with the other wall beam theory (same wall option), if computed."""
    return next(((k, o) for k, o in analytical.get("options", {}).items()
                 if ", " in k and "Q" in o and k.startswith(analytical.get("active", ""))), None)


def summary_table(rows: list[dict], project: str, group: str, analytical: dict,
                  convergence: list[dict] | None = None, *, include_negligible: bool = False
                  ) -> tuple[list[str], list[list], dict[str, str]]:
    """Header, rows and value notes for one project's term table (GUI and reports).

    Columns: term, unit, analytical (both beam theories), then for every FE source present its
    value and the analytical-vs-FE difference in %.  Matrix terms that are zero to round-off
    are left out unless ``include_negligible``; the notes say how many.
    """
    sources = [(s, m) for s in SOLVERS for m in MODEL_GROUPS if group in MODEL_GROUPS[m]
               and any(r["project"] == project and r["solver"] == s and r["model"] == m and r["group"] == group
                       for r in rows)]
    ref = analytical[REF_KEY[group]]
    label = {"DBC": "DBC" if group == "shell" else "DBC (U3 = 0)", "PBC": "PBC"}
    beam_name = {"euler-bernoulli": "Euler\u2013Bernoulli", "timoshenko": "Timoshenko"}
    alt = _alt_option(analytical)
    cols = ["term", "unit", f"analytical ({beam_name.get(analytical.get('beam', ''), 'model')})"]
    if alt:
        cols.append(f"analytical ({beam_name.get(alt[1].get('beam'), alt[0])})")
    best, notes = {}, {}
    for s, m in sources:
        best[(s, m)] = best_values(rows, project, s, m, group, convergence)
        cols += [f"{s} {label[m]}", "analytical vs FE (%)"]
        notes[f"{s} {label[m]}"] = ", ".join(sorted({h for (_v, h) in best[(s, m)].values()}))
    table, skipped = [], []
    for t in GROUPS[group]:
        fe_vals = [best[k].get(t, (math.nan, ""))[0] for k in sources]
        neg = negligible(group, t, ref) and all(negligible(group, t, ref, v) for v in fe_vals if math.isfinite(v))
        if neg and not include_negligible:
            skipped.append(t)
            continue
        line = [t, term_unit(t), float(ref[t])]
        if alt:
            line.append(float(alt[1][REF_KEY[group]][t]))
        for v in fe_vals:
            line += [v, math.nan if neg else pct(ref[t], v)]
        table.append(line)
    if skipped:
        notes["≈ 0"] = f"{len(skipped)} terms zero to round-off not listed ({', '.join(skipped)})"
    return cols, table, notes


MATRIX_LABELS = {"C": ("11", "22", "33", "23", "13", "12"), "Q": ("11", "22", "12", "13", "23")}


def matrix_from_terms(group: str, terms: dict) -> np.ndarray:
    """Symmetric 6x6 C (Voigt [11,22,33,23,13,12]) or 5x5 Q ([11,22,12,13,23]) from named terms."""
    if group == "C":
        M = np.full((6, 6), math.nan)
        for i in range(6):
            for j in range(i, 6):
                M[i, j] = M[j, i] = _num(terms.get(f"C{VOIGT_PAIRS[i]}{VOIGT_PAIRS[j]}"))
        return M
    if group == "Q":
        M = np.full((5, 5), math.nan)
        for i in range(5):
            for j in range(i, 5):
                a, b = sorted((Q_ORDER[i], Q_ORDER[j]))
                M[i, j] = M[j, i] = _num(terms.get(f"Q{a}{b}"))
        return M
    raise ValueError("Matrices exist for groups 'C' and 'Q'.")


def matrix_comparison(rows: list[dict], project: str, group: str, analytical: dict,
                      convergence: list[dict] | None = None) -> dict | None:
    """Analytical and FE matrix (Abaqus PBC if available, else built-in) and the % difference.

    Difference entries are NaN where the term is zero to round-off.  None without PBC results.
    """
    ref = analytical[REF_KEY[group]]
    for solver in reversed(SOLVERS):  # Abaqus first
        b = best_values(rows, project, solver, "PBC", group, convergence)
        if not b:
            continue
        fe = {t: v for t, (v, _h) in b.items()}
        A, F = matrix_from_terms(group, ref), matrix_from_terms(group, fe)
        D = np.full_like(A, math.nan)
        names = MATRIX_LABELS[group]
        for i in range(len(names)):
            for j in range(len(names)):
                if group == "C":
                    t = f"C{VOIGT_PAIRS[min(i, j)]}{VOIGT_PAIRS[max(i, j)]}"
                else:
                    a, c = sorted((Q_ORDER[i], Q_ORDER[j]))
                    t = f"Q{a}{c}"
                if not negligible(group, t, ref, fe.get(t)):
                    D[i, j] = pct(ref[t], fe.get(t, math.nan))
                else:
                    A[i, j] = F[i, j] = 0.0  # round-off shown as 0 (the rows keep the raw values)
        how = ", ".join(sorted({h for (_v, h) in b.values()}))
        return {"solver": solver, "labels": names, "analytical": A, "FE": F, "difference_pct": D, "how": how}
    return None


def findings(rows: list[dict], project: str, analytical: dict, rve_summary: str | None = None,
             convergence: list[dict] | None = None) -> list[str]:
    """Plain-language conclusions for one project (FE against the analytical model only)."""
    out = [f"FE RVE: {rve_summary}"] if rve_summary else []
    sel = [r for r in rows if r["project"] == project]
    if not sel:
        return out + ["No FE results yet."]
    what = {"C": "Full C (strain energy)", "Q": "Q matrix (σ33 = 0)", "eng": "Engineering constants",
            "shell": "DBC shell terms (σ33 = 0)", "3d": "DBC 3D terms (ε33 = 0)"}

    def worst(solver, model, group):
        b = best_values(rows, project, solver, model, group, convergence)
        ref = analytical[REF_KEY[group]]
        d = [(abs(pct(ref[t], v)), t, pct(ref[t], v)) for t, (v, _h) in b.items()
             if not negligible(group, t, ref, v) and math.isfinite(pct(ref[t], v))]
        return max(d) if d else None

    for solver in SOLVERS:
        for model, groups in MODEL_GROUPS.items():
            for g in groups:
                w = worst(solver, model, g)
                if w is not None:
                    out.append(f"{what[g]}, analytical vs {solver} {model}: largest difference {signed(w[2])}% ({w[1]}).")
    for solver in reversed(SOLVERS):  # which analytical wall option the periodic FE supports (Abaqus first)
        b = best_values(rows, project, solver, "PBC", "C", convergence)
        if not b:
            continue
        cand = []
        for name, opt in analytical.get("options", {}).items():
            if "C21" in opt:
                errs = [abs(pct(opt["C21"][t], v)) for t, (v, _h) in b.items()
                        if not negligible("C", t, opt["C21"], v) and math.isfinite(pct(opt["C21"][t], v))]
                if errs:
                    cand.append((max(errs), name))
        if cand:
            e = min(cand)[0]
            same = [n for ee, n in sorted(cand) if abs(ee - e) <= 1e-9]
            label = " = ".join(same) + (" (identical for this stack)" if len(same) > 1 else "")
            out.append(f"Closest analytical wall option to the {solver} PBC: {label} (largest difference {e:.2f}%).")
        break
    finest = min(_f(r["mesh_mm"]) for r in sel)
    if finest > 5.0:
        out.append(f"Finest mesh {finest:g} mm: the wall-bending terms (C1212, Q66, G12) converge from above; "
                   "use 5 mm or finer for values to quote.")
    return out


# ------------------------------------------------------------------------------------ one design
@dataclass(slots=True)
class VerificationSettings:
    mesh_size: float = 10.0
    run_dbc: bool = True
    run_pbc: bool = True
    displacement: float = 1.0
    face_wall_section: str = "half-thickness"
    return_fields: bool = True

    def validate(self) -> None:
        if not math.isfinite(self.mesh_size) or self.mesh_size < 2.5:
            raise ValueError("The built-in shell check needs a mesh size of at least 2.5 mm.")
        if not (self.run_dbc or self.run_pbc):
            raise ValueError("Select at least one of the DBC and PBC models.")


def verify_model(model: RVEModel, materials: dict[str, Material], laminates: list[Laminate],
                 settings: VerificationSettings | None = None,
                 progress: Callable[[str], None] | None = None) -> dict:
    """One design, one mesh: analytical vs built-in FE, no files written."""
    from .fe_rve import select_fe_rve, fe_rve_summary
    s = settings or VerificationSettings()
    s.validate()
    t0 = time.time()
    say = progress or (lambda _m: None)
    say("Selecting the FE RVE")
    _block, rve = select_fe_rve(model, face_wall_section=s.face_wall_section)
    say("Analytical model")
    ana = analytical_reference(model, materials, laminates)
    config = paper_configuration_of(model, materials, laminates)
    say(f"Built-in shell FE, {s.mesh_size:g} mm mesh")
    models = [m for m, on in (("DBC", s.run_dbc), ("PBC", s.run_pbc)) if on]
    fe = builtin_check(model, materials, laminates, s.mesh_size, models=models, displacement=s.displacement,
                       face_wall_section=s.face_wall_section, fields=s.return_fields)
    name = str(model.name) or "design"
    rows = result_rows(name, "built-in", s.mesh_size, fe, ana)
    return {"project": name, "model": name if config is None else f"{name} (= benchmark configuration {config})",
            "settings": asdict(s), "rve": rve, "rve_summary": fe_rve_summary(rve), "analytical": ana,
            "benchmark_configuration": config, "fe": fe, "rows": rows,
            "findings": findings(rows, name, ana, fe_rve_summary(rve)),
            "seconds": round(time.time() - t0, 2)}
