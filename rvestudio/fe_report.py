"""Self-contained HTML report and Excel workbook of an FE verification study.

``write_study_report(root)`` -> ``FE_STUDY_REPORT.html`` and ``write_study_workbook(root)``
-> ``FE_STUDY_RESULTS.xlsx``, from the result rows the built-in check and the Abaqus
collection write into the study folder.  Every project gets the same block: its FE
RVE, the findings, the full stiffness matrix C and the Q matrix (sigma33 = 0) from the
strain-energy PBC model next to the analytical ones, tables of C, Q and the engineering
constants, the DBC shell and 3D terms, difference charts and the deformed shapes.  Every
FE value is compared with the analytical model only.  ``write_verification_report`` renders the same block for one
:func:`rvestudio.verification.verify_model` result.

The HTML embeds every figure as PNG, so a report is one file that can be e-mailed.
"""
from __future__ import annotations

import csv
import html
import json
import math
from datetime import datetime
from pathlib import Path

import numpy as np
from matplotlib.figure import Figure

from . import __version__
from .verification import (signed, summary_table, findings, GROUPS, GROUP_TITLES, BOUNDARY_CONDITIONS,
                           SHELL_TERMS, SOLID_TERMS, best_values, pct, SOLVERS, REF_KEY, matrix_comparison,
                           negligible)
from .fe_plots import figure_data_uri, plot_rve_plan, project_figure, deformed_figure, plot_deformed, CASE_FIELDS

VOIGT = ("11", "22", "33", "23", "13", "12")

CSS = """
body{font-family:system-ui,-apple-system,'Segoe UI',Roboto,Arial,sans-serif;margin:0;background:#f6f7f9;color:#1d232b}
main{max-width:1180px;margin:0 auto;padding:24px 20px 60px}
h1{font-size:24px;margin:0 0 4px} h2{font-size:18px;margin:32px 0 10px;border-bottom:1px solid #d8dde3;padding-bottom:4px}
h3{font-size:15px;margin:20px 0 8px} h4{font-size:13.5px;margin:16px 0 6px}
.meta{color:#5b6672;font-size:13px;margin-bottom:18px}
.card{background:#fff;border:1px solid #dfe3e8;border-radius:8px;padding:14px 16px;margin:12px 0}
ul.findings{margin:6px 0 0 18px;padding:0} ul.findings li{margin:4px 0}
table{border-collapse:collapse;font-size:13px;background:#fff;margin:6px 0 4px}
th,td{border:1px solid #dfe3e8;padding:4px 8px;text-align:right;white-space:nowrap}
th{background:#eef1f4;font-weight:600} td.l,th.l{text-align:left}
td.ok{background:#e5f4e4} td.warn{background:#f9eeda} td.bad{background:#fbebeb}
img{max-width:100%;height:auto;display:block;margin:6px 0}
.note{color:#5b6672;font-size:12.5px}
.grid{display:grid;grid-template-columns:minmax(260px,360px) 1fr;gap:16px;align-items:start}
@media (max-width:760px){.grid{grid-template-columns:1fr}}
.wide{overflow-x:auto}
code{background:#eef1f4;padding:1px 4px;border-radius:3px}
nav a{margin-right:12px;font-size:13px}
"""
DIFF_NOTE = ("Differences in %: 100 (first − second) / second. Green |difference| ≤ 1 %, amber ≤ 5 %, "
             "red > 5 %. FE values: Richardson-extrapolated when three meshes converge monotonically, else the "
             "finest mesh (noted under each table).")


def _e(x) -> str:
    return html.escape(str(x))


def _num(v, fmt=".4g") -> str:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return _e(v) if v not in (None, "") else "&ndash;"
    return format(f, fmt) if math.isfinite(f) else "&ndash;"


def _diff_cell(v) -> str:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return "<td>&ndash;</td>"
    if not math.isfinite(f):
        return "<td>&ndash;</td>"
    cls = "ok" if abs(f) <= 1.0 else ("warn" if abs(f) <= 5.0 else "bad")
    return f'<td class="{cls}">{signed(f)}</td>'


def _table(header: list[str], rows: list[list], diff_cols: set[int] = frozenset(), left_cols: set[int] = frozenset({0})) -> str:
    h = "".join(f'<th class="{"l" if i in left_cols else ""}">{_e(c)}</th>' for i, c in enumerate(header))
    body = []
    for r in rows:
        cells = []
        for i, v in enumerate(r):
            if i in diff_cols:
                cells.append(_diff_cell(v))
            elif i in left_cols:
                cells.append(f'<td class="l">{_e(v)}</td>')
            else:
                cells.append(f"<td>{_num(v)}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    return f"<div class='wide'><table><thead><tr>{h}</tr></thead><tbody>{''.join(body)}</tbody></table></div>"


def _page(title: str, meta: str, body: str) -> str:
    return (f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' "
            f"content='width=device-width,initial-scale=1'><title>{_e(title)}</title><style>{CSS}</style></head>"
            f"<body><main><h1>{_e(title)}</h1><div class='meta'>{meta}</div>{body}</main></body></html>")


def _img(fig: Figure, alt: str, dpi: int = 120) -> str:
    return f"<img alt='{_e(alt)}' src='{figure_data_uri(fig, dpi)}'>"


def _bc_table() -> str:
    return _table(["model", "side faces", "top and bottom faces", "terms"], [list(r) for r in BOUNDARY_CONDITIONS],
                  left_cols={0, 1, 2, 3})


def _term_tables(rows, pid, analytical, conv) -> dict:
    out = {}
    for g in GROUPS:
        cols, table, notes = summary_table(rows, pid, g, analytical, conv)
        if len(cols) > 3 and table:
            out[g] = (cols, table, notes)
    return out


def _matrix_html(mc: dict, group: str) -> str:
    """Analytical | FE | difference, side by side, for the C (6x6) or Q (5x5) matrix."""
    lab = list(mc["labels"])

    def grid(M, diff=False):
        h = "<th></th>" + "".join(f"<th>{_e(x)}</th>" for x in lab)
        body = []
        for i, r in enumerate(lab):
            cells = [f"<th>{_e(r)}</th>"]
            for j in range(len(lab)):
                v = float(M[i, j])
                if diff:
                    cells.append(_diff_cell(v) if math.isfinite(v) else "<td class='note'>≈0</td>")
                else:
                    cells.append(f"<td>{_num(v, '.5g')}</td>")
            body.append("<tr>" + "".join(cells) + "</tr>")
        return f"<table><thead><tr>{h}</tr></thead><tbody>{''.join(body)}</tbody></table>"
    name = "C (6×6, Voigt 11, 22, 33, 23, 13, 12)" if group == "C" else "Q (σ33 = 0; 11, 22, 12, 13, 23)"
    return (f"<h4>{_e(name)}, MPa</h4><div class='wide' style='display:flex;gap:14px;flex-wrap:wrap'>"
            f"<div><b>Analytical</b>{grid(mc['analytical'])}</div>"
            f"<div><b>FE: {_e(mc['solver'])} PBC</b> <span class='note'>({_e(mc['how'])})</span>{grid(mc['FE'])}</div>"
            f"<div><b>Analytical vs FE (%)</b>{grid(mc['difference_pct'], True)}</div></div>")


def _project_block(pid: str, title: str, rve: dict | None, rve_text: str | None, analytical: dict,
                   rows: list[dict], conv: list[dict],
                   fields: dict | None = None, abaqus_fig: Figure | None = None) -> str:
    tables = _term_tables(rows, pid, analytical, conv)
    found = findings(rows, pid, analytical, None, conv)
    block = [f"<div class='card' id='{_e(pid)}'><h3>{_e(title)}</h3><div class='grid'><div>"]
    if rve:
        fig = Figure(figsize=(4.6, 4.4))
        plot_rve_plan(fig.add_subplot(111), rve, title="FE RVE")
        block.append(_img(fig, f"FE RVE {pid}", dpi=90))
    block.append("</div><div>")
    if rve_text:
        block.append(f"<p><b>FE RVE.</b> {_e(rve_text)}</p>")
    block.append("<b>Findings</b><ul class='findings'>" + "".join(f"<li>{_e(f)}</li>" for f in found) + "</ul>")
    block.append("</div></div>")
    for g in ("C", "Q"):
        mc = matrix_comparison(rows, pid, g, analytical, conv)
        if mc is not None:
            block.append(_matrix_html(mc, g))
    for g, (cols, table, notes) in tables.items():
        diff = {i for i, c in enumerate(cols) if c.endswith("(%)")}
        block.append(f"<h4>{_e(GROUP_TITLES[g])}</h4>" + _table(cols, table, diff_cols=diff, left_cols={0, 1}))
        if notes:
            block.append("<p class='note'>FE values: " + "; ".join(f"{_e(k)} {_e(v)}" for k, v in notes.items()) + ".</p>")
    if tables:
        block.append(_img(project_figure(None, {g: t[:2] for g, t in tables.items()}, Figure(figsize=(11, 4.2))),
                          "differences", dpi=100))
    for model in ("DBC", "PBC"):
        f = (fields or {}).get(model)
        if f is not None and len(f.get("u", {})):
            mesh = (fields or {}).get("mesh_mm")
            block.append(f"<h4>Built-in {model}: deformed shapes{f', {mesh:g} mm mesh' if mesh else ''}</h4>")
            block.append(_img(deformed_figure(f, model, Figure(figsize=(12, 5.2 if model == 'DBC' else 2.9))),
                              f"{model} deformed", dpi=80))
    if abaqus_fig is not None:
        block.append("<h4>Abaqus DBC: deformed shapes</h4>" + _img(abaqus_fig, "Abaqus deformed", dpi=80))
    block.append("</div>")
    return "".join(block)


# ----------------------------------------------------------------------------- one design
def verification_html(v: dict) -> str:
    """Report of one :func:`verification.verify_model` result (built-in solver)."""
    s = v["settings"]
    meta = (f"Model <b>{_e(v['model'])}</b> &middot; built-in shell solver, {s['mesh_size']:g} mm mesh &middot; "
            f"RVE Studio {_e(__version__)} &middot; {datetime.now():%Y-%m-%d %H:%M} &middot; {v.get('seconds', 0):.1f} s")
    body = ("<h2>Boundary-value problems</h2>" + _bc_table() + f"<p class='note'>{_e(DIFF_NOTE)}</p>" +
            _project_block(v["project"], v["model"], v["rve"], v["rve_summary"], v["analytical"],
                           v["rows"], [], (v.get("fe") or {}).get("fields")))
    return _page("RVE verification report", meta, body)


def write_verification_report(v: dict, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(verification_html(v), encoding="utf-8")
    return path


def write_verification_workbook(v: dict, path: str | Path) -> Path:
    import pandas as pd
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        pd.DataFrame({"finding": [v["rve_summary"]] + v["findings"]}).to_excel(xw, sheet_name="Summary", index=False)
        for g in GROUPS:
            cols, table, _n = summary_table(v["rows"], v["project"], g, v["analytical"], include_negligible=True)
            if table and len(cols) > 3:
                pd.DataFrame(table, columns=cols).to_excel(xw, sheet_name=SHEET_NAMES[g], index=False)
        pd.DataFrame(v["rows"]).to_excel(xw, sheet_name="All results", index=False)
    return path


# ----------------------------------------------------------------------------- FE study
def _read_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


def read_inp_mesh(path: str | Path) -> tuple[dict[int, tuple[float, float, float]], list[list[int]]]:
    """Nodes and S4 connectivity of an RVE Studio Abaqus deck."""
    nodes: dict[int, tuple[float, float, float]] = {}
    elems: list[list[int]] = []
    mode = None
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("**"):
            continue
        if line.startswith("*"):
            u = line.upper()
            mode = "n" if u.startswith("*NODE") and "OUTPUT" not in u and "FILE" not in u else (
                "e" if u.startswith("*ELEMENT") and "OUTPUT" not in u else None)
            continue
        if not line.strip():
            continue
        try:
            if mode == "n":
                p = line.split(",")
                nodes[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]) if len(p) > 3 else 0.0)
            elif mode == "e":
                p = [int(x) for x in line.split(",") if x.strip()]
                if len(p) >= 5:
                    elems.append(p[1:5])
        except ValueError:
            continue
    return nodes, elems


def abaqus_dbc_fields(run_dir: Path) -> dict | None:
    """{'coords', 'elements', 'u': {case: N x 3}} from the nodal_displacements.csv files of a DBC run."""
    cases = sorted(d for d in Path(run_dir).iterdir() if d.is_dir() and (d / "nodal_displacements.csv").exists())
    if not cases:
        return None
    inp = next(iter(cases[0].glob("*.inp")), None)
    if inp is None:
        return None
    nodes, elems = read_inp_mesh(inp)
    ids = sorted(nodes)
    idx = {n: i for i, n in enumerate(ids)}
    out = {"coords": np.array([nodes[n] for n in ids], float),
           "elements": np.array([[idx[n] for n in e] for e in elems if all(n in idx for n in e)], int), "u": {}}
    for cdir in cases:
        u_by = {int(r["label"]): (float(r["u1"]), float(r["u2"]), float(r["u3"])) for r in _read_csv(cdir / "nodal_displacements.csv")}
        out["u"][cdir.name.split("_", 1)[-1]] = np.array([u_by.get(n, (0.0, 0.0, 0.0)) for n in ids], float)
    return out


def project_rve(root: Path, study: dict, pid: str) -> dict | None:
    for run in study.get("runs", []):
        if run.get("project") == pid and run.get("check") == "dbc" and run.get("run_dir"):
            snap = root / run["run_dir"] / "thesis_dbc_snapshot.json"
            if snap.exists():
                dom = json.loads(snap.read_text(encoding="utf-8")).get("domain", {})
                if dom.get("wall_segments"):
                    return dom
    try:
        from .io import load_project
        from .fe_rve import select_fe_rve
        model, _m, _l = load_project(root / study["projects"][pid]["project"])
        return select_fe_rve(model)[1]
    except Exception:
        return None


SHEET_NAMES = {"C": "C matrix", "Q": "Q matrix", "eng": "Engineering constants", "shell": "DBC shell terms",
               "3d": "DBC 3D terms"}


def _overview(rows, study, conv) -> str:
    parts = []
    for g, model in (("Q", "PBC"), ("eng", "PBC"), ("shell", "DBC")):
        projects = list(study["projects"])
        grid, used = {}, set()
        for pid in projects:
            ana = study["projects"][pid]["analytical"][REF_KEY[g]]
            for solver in reversed(SOLVERS):  # Abaqus first
                b = best_values(rows, pid, solver, model, g, conv)
                if b:
                    used.add(solver)
                    for t, (v, _h) in b.items():
                        if not negligible(g, t, ana, v):
                            grid[(pid, t)] = pct(ana[t], v)
                    break
        if not grid:
            continue
        terms = [t for t in GROUPS[g] if any((p, t) in grid for p in projects)]
        M = np.array([[grid.get((p, t), math.nan) for t in terms] for p in projects], float)
        fig = Figure(figsize=(8.5, 0.3 * len(projects) + 1.7))
        ax = fig.add_subplot(111)
        lim = max(2.0, float(np.nanmax(np.abs(M))) if np.isfinite(M).any() else 2.0)
        im = ax.imshow(M, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
        for (i, j), val in np.ndenumerate(M):
            if math.isfinite(val):
                ax.text(j, i, signed(val, 1), ha="center", va="center", fontsize=7,
                        fontweight="bold" if abs(val) > 1 else "normal", color="k" if abs(val) < 0.6 * lim else "w")
        ax.set_xticks(range(len(terms))); ax.set_xticklabels(terms, fontsize=8)
        ax.set_yticks(range(len(projects))); ax.set_yticklabels(projects, fontsize=7)
        cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label("analytical − FE (%)", fontsize=8)
        ax.set_title(f"{GROUP_TITLES[g]}: analytical vs {' / '.join(sorted(used, reverse=True))} {model}", fontsize=9)
        fig.tight_layout()
        parts.append(_img(fig, f"overview {g}"))
    return "".join(parts)


def study_html(root: str | Path) -> str:
    from .fe_workflow import load_results, load_fields
    from .fe_rve import fe_rve_summary
    R = load_results(root)
    root, study, rows, conv = R["root"], R["study"], R["rows"], R["convergence"]
    counts = R["status"]["counts"]
    s = study["settings"]
    solvers = sorted({r["solver"] for r in rows}, reverse=True)
    meta = (f"Study <code>{_e(root.name)}</code> &middot; {len(study['projects'])} projects &middot; meshes "
            f"{', '.join(f'{h:g}' for h in s['meshes'])} mm &middot; models {', '.join(c.upper() for c in s['checks'])} &middot; "
            f"RVE Studio {_e(__version__)} &middot; {datetime.now():%Y-%m-%d %H:%M}")
    parts = ["<div class='card'><b>Results</b>: " + (", ".join(solvers) if solvers else "none yet") +
             ". <b>Abaqus jobs</b>: " + (", ".join(f"{_e(k)} {v}" for k, v in sorted(counts.items())) or "none prepared") +
             ". Only finished results are shown; nothing pending is estimated.</div>"]
    parts.append("<nav>" + "".join(f"<a href='#{_e(p)}'>{_e(p)}</a>" for p in study["projects"]) + "</nav>")
    parts.append("<h2>Boundary-value problems</h2>" + _bc_table() + f"<p class='note'>{_e(DIFF_NOTE)}</p>")
    ov = _overview(rows, study, conv)
    if ov:
        parts.append("<h2>Overview</h2>" + ov)
    parts.append("<h2>Projects</h2>")
    for pid, info in study["projects"].items():
        rve = project_rve(root, study, pid)
        ana = info["analytical"]
        title = pid + (f" (benchmark configuration {ana['paper_config']})" if ana.get("paper_config") else "")
        afig = None
        for run in study.get("runs", []):
            if run.get("project") == pid and run.get("check") == "dbc" and run.get("run_dir"):
                f = abaqus_dbc_fields(root / run["run_dir"])
                if f is not None:
                    afig = deformed_figure(f, "DBC", Figure(figsize=(12, 5.2)), solver="Abaqus")
                    break
        parts.append(_project_block(pid, title, rve, fe_rve_summary(rve) if rve else None, ana,
                                    rows, conv, load_fields(root / "results" / "builtin_fields" / f"{pid}.npz"), afig))
    multi = [c for c in conv if int(c["meshes"]) >= 2]
    if multi:
        parts.append("<h2>Mesh convergence</h2>")
        parts.append(_table(["project", "solver", "model", "term", "unit", "meshes", "finest (mm)", "FE finest",
                             "FE reference", "analytical", "analytical vs FE (%)", "convergence", "order p",
                             "GCI fine (%)", "asymptotic ratio", "mesh verdict"],
                            [[c["project"], c["solver"], c["model"], c["term"], c.get("unit", ""), c["meshes"],
                              c["finest_mesh_mm"], c["FE_finest"], c["FE_reference"], c["analytical"],
                              c["analytical_vs_FE_pct"], c.get("status", ""), c.get("p", ""), c.get("gci_fine_pct", ""),
                              c.get("asymptotic_ratio", ""), c.get("verdict", "")]
                             for c in multi],
                            diff_cols={10}, left_cols={0, 1, 2, 3, 4, 11, 15}))
        parts.append("<p class='note'>A term is <b>converged</b> when the three meshes agree to round-off, or when they "
                     "converge monotonically with a plausible observed order (0.5 &le; p &le; 4), lie in the asymptotic range "
                     "(asymptotic ratio GCI<sub>32</sub>/(r<sup>p</sup> GCI<sub>21</sub>) within 1 &plusmn; 0.15) and the "
                     "fine-mesh GCI is at most 1 %. Otherwise add a finer mesh (halve the finest size) and collect again.</p>")
    if R["pbc_audit"]:
        parts.append("<h3>Abaqus PBC audit (strain energy vs macro reactions)</h3>")
        parts.append(_table(["project", "mesh (mm)", "load cases", "energy vs reaction C", "reaction reciprocity",
                             "positive definite"],
                            [[a["project"], a["mesh_mm"], a.get("cases", ""), a.get("energy_vs_reaction_error", ""),
                              a.get("reaction_reciprocity_error", ""), a.get("positive_definite", "")]
                             for a in R["pbc_audit"]], left_cols={0, 5}))
    bad = [j for j in R["status"]["jobs"] if j.get("state") in ("FAILED", "NOT PREPARED")]
    if bad:
        parts.append("<h2>Problems</h2>")
        parts.append(_table(["model", "project", "mesh", "case", "state", "note"],
                            [[j["check"].upper(), j["project"], j["mesh"], j["case"], j["state"], j.get("note", "")]
                             for j in bad[:300]], left_cols={0, 1, 3, 4, 5}))
    return _page("FE verification study", meta, "".join(parts))


def write_study_report(root: str | Path) -> Path:
    root = Path(root)
    out = root / "FE_STUDY_REPORT.html"
    out.write_text(study_html(root), encoding="utf-8")
    return out


def write_study_workbook(root: str | Path) -> Path:
    import pandas as pd
    from .fe_workflow import load_results
    R = load_results(root)
    root, study, rows, conv = R["root"], R["study"], R["rows"], R["convergence"]
    out = root / "FE_STUDY_RESULTS.xlsx"
    with pd.ExcelWriter(out, engine="openpyxl") as xw:
        summ = []
        for pid, info in study["projects"].items():
            ana = info["analytical"]
            for f in findings(rows, pid, ana, None, conv):
                summ.append({"project": pid, "finding": f})
        pd.DataFrame(summ or [{"project": "", "finding": "No results yet."}]).to_excel(xw, sheet_name="Summary", index=False)
        for g, name in SHEET_NAMES.items():
            frames = []
            for pid, info in study["projects"].items():
                ana = info["analytical"]
                cols, table, _n = summary_table(rows, pid, g, ana, conv, include_negligible=True)
                if len(cols) > 3 and table:
                    df = pd.DataFrame(table, columns=cols)
                    df.insert(0, "project", pid)
                    frames.append(df)
            if frames:
                pd.concat(frames, ignore_index=True).to_excel(xw, sheet_name=name, index=False)
        if rows:
            pd.DataFrame(rows).to_excel(xw, sheet_name="All results", index=False)
        if conv:
            pd.DataFrame(conv).to_excel(xw, sheet_name="Convergence", index=False)
        if R["pbc_audit"]:
            pd.DataFrame(R["pbc_audit"]).to_excel(xw, sheet_name="PBC audit", index=False)
        if R["status"]["jobs"]:
            pd.DataFrame(R["status"]["jobs"]).to_excel(xw, sheet_name="Abaqus jobs", index=False)
        pt = root / "results" / "paper_table.csv"
        if pt.exists() and pt.stat().st_size > 1:
            try:
                pd.read_csv(pt).to_excel(xw, sheet_name="Paper table", index=False)
            except Exception:
                pass
    return out
