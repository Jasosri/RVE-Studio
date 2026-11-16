"""Matplotlib drawings of FE RVEs, boundary conditions, deformed shapes and comparisons.

Every function draws on a given Axes or Figure, so the same picture is used in the
GUI (FigureCanvasQTAgg), in the HTML report (PNG) and in saved figures.  Nothing
here needs Qt or Abaqus.
"""
from __future__ import annotations

import io
import base64
import math

import numpy as np
from matplotlib.figure import Figure
from matplotlib.patches import Polygon, Rectangle, FancyArrowPatch, Circle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from .verification import signed

FACE_COLOURS = {"X1_0": "#0072B2", "X1_P": "#D55E00", "X2_0": "#009E73", "X2_Q": "#CC79A7"}
FACE_WALL = "#E69F00"
WALL = "0.25"
CELL_FILL = "#DCE8F7"
CASE_FIELDS = {"Q11": ((0,), "u1"), "Q22": ((1,), "u2"), "Q12": ((0, 1), "(u1+u2)/2"),
               "Q66": ((1,), "u2"), "Q44": ((2,), "u3"), "Q55": ((2,), "u3"),
               "C1111": ((0,), "u1"), "C2222": ((1,), "u2"), "C1122": ((0, 1), "(u1+u2)/2"), "C1212": ((1,), "u2")}
PBC_FIELDS = {"11": (0, "u1"), "22": (1, "u2"), "33": (2, "u3"), "23": (2, "u3"), "13": (2, "u3"), "12": (1, "u2")}
PBC_TITLES = {"11": "unit E11", "22": "unit E22", "33": "unit E33", "23": "unit G23", "13": "unit G13",
              "12": "unit G12"}


# ----------------------------------------------------------------------------- plan of the FE RVE
def plot_rve_plan(ax, rve: dict, *, title: str | None = None, dims: bool = True, faces: bool = True,
                  legend: bool = True, lw: float = 1.8) -> None:
    """Plan of the FE RVE from its diagnostics (``fe_rve.select_fe_rve`` diag)."""
    p, q = float(rve["p"]), float(rve["q"])
    poly = rve.get("central_cell_polygon_rve")
    if poly:
        ax.add_patch(Polygon(np.asarray(poly), closed=True, fc=CELL_FILL, ec="none", zorder=0))
    face = set(int(i) for i in rve.get("face_walls", []))
    for k, (x0, y0, x1, y1) in enumerate(rve.get("wall_segments", []), 1):
        if k in face:
            ax.plot([x0, x1], [y0, y1], color=FACE_WALL, lw=lw * 1.5, solid_capstyle="butt", zorder=3)
        else:
            ax.plot([x0, x1], [y0, y1], color=WALL, lw=lw, solid_capstyle="round", zorder=2)
    ax.add_patch(Rectangle((0, 0), p, q, fill=False, ec="0.55", lw=0.7, ls=(0, (4, 3)), zorder=1))
    s = max(p, q)
    if faces:
        off = 0.05 * s
        ax.text(-off, q / 2, "x1 = 0", rotation=90, ha="center", va="center", fontsize=7, color=FACE_COLOURS["X1_0"])
        ax.text(p + off, q / 2, "x1 = p", rotation=90, ha="center", va="center", fontsize=7, color=FACE_COLOURS["X1_P"])
        ax.text(p / 2, -off, "x2 = 0", ha="center", va="center", fontsize=7, color=FACE_COLOURS["X2_0"])
        ax.text(p / 2, q + off, "x2 = q", ha="center", va="center", fontsize=7, color=FACE_COLOURS["X2_Q"])
    if dims:
        y = -0.13 * s
        ax.annotate("", xy=(p, y), xytext=(0, y), arrowprops=dict(arrowstyle="<->", lw=0.7, color="0.3"))
        ax.text(p / 2, y - 0.03 * s, f"p = {p:.1f} mm", ha="center", va="top", fontsize=7)
        x = -0.13 * s
        ax.annotate("", xy=(x, q), xytext=(x, 0), arrowprops=dict(arrowstyle="<->", lw=0.7, color="0.3"))
        ax.text(x - 0.03 * s, q / 2, f"q = {q:.1f} mm", rotation=90, ha="right", va="center", fontsize=7)
    if legend:
        ax.plot([], [], color=WALL, lw=lw, label="wall")
        if face:
            ax.plot([], [], color=FACE_WALL, lw=lw * 1.5, label="wall in a face (half laminate)")
        ax.add_patch(Rectangle((0, 0), 0, 0, fc=CELL_FILL, label=f"complete {rve.get('central_cell', {}).get('name', 'cell')}"))
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, 0.02), fontsize=7, frameon=False,
                  ncol=1, handlelength=1.5, borderaxespad=0.0)
    pad = 0.22 * s
    ax.set_xlim(-pad, p + 0.12 * s)
    ax.set_ylim(-pad - (0.28 * s if legend else 0.0), q + 0.12 * s)
    ax.set_aspect("equal")
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=9)


# ----------------------------------------------------------------------------- replicated design + RVE choice
RVE_CHOSEN = "#0072B2"
RVE_AUTO = "0.45"


def plot_rve_layout(ax, layout: dict, *, rect: tuple[float, float, float, float] | None = None,
                    rve: dict | None = None, drag: tuple[float, float, float, float] | None = None,
                    grid: bool = True, centres: bool = True, title: str | None = None) -> None:
    """The design replicated in the FE frame (x1 along a1) with the automatic and chosen FE RVE.

    ``layout`` comes from :func:`rvestudio.fe_rve.fe_rve_layout`; ``rect`` = (x0, y0, w, h)
    of the chosen RVE (:func:`rvestudio.fe_rve.spec_rectangle`) and ``rve`` its diagnostics,
    used to highlight the walls it contains; ``drag`` is a rectangle being drawn.
    """
    x_lo, x_hi, y_lo, y_hi = layout["extent"]
    S = np.asarray(layout["segments"], float).reshape(-1, 4)
    for x0, y0, x1, y1 in S:
        ax.plot([x0, x1], [y0, y1], color="0.72", lw=1.0, solid_capstyle="round", zorder=1)
    if grid:
        # Rectangular periodic cells: any RVE whose sides are whole multiples is exactly periodic.
        rw, rh = layout["rect_width"], layout["rect_height"]
        ax0, ay0 = layout["auto"]["window_origin_local"]
        for k in range(math.floor((x_lo - ax0) / rw), math.ceil((x_hi - ax0) / rw) + 1):
            ax.axvline(ax0 + k * rw, color="0.88", lw=0.6, zorder=0)
        for k in range(math.floor((y_lo - ay0) / rh), math.ceil((y_hi - ay0) / rh) + 1):
            ax.axhline(ay0 + k * rh, color="0.88", lw=0.6, zorder=0)
    if centres and layout.get("cell_centres"):
        C = np.asarray(layout["cell_centres"], float)
        ax.plot(C[:, 0], C[:, 1], "+", color="0.6", ms=5, mew=0.8, zorder=1)
    a = layout["auto"]
    ax.add_patch(Rectangle(tuple(a["window_origin_local"]), a["p"], a["q"], fill=False, ec=RVE_AUTO,
                           lw=1.1, ls=(0, (5, 3)), zorder=3, label="automatic RVE"))
    if rect is not None:
        x0, y0, w, h = rect
        if rve is not None and rve.get("wall_segments"):
            face = set(int(i) for i in rve.get("face_walls", []))
            for k, (u0, v0, u1, v1) in enumerate(rve["wall_segments"], 1):
                ax.plot([x0 + u0, x0 + u1], [y0 + v0, y0 + v1], color=FACE_WALL if k in face else WALL,
                        lw=3.0 if k in face else 1.7, solid_capstyle="round", zorder=6 if k in face else 4)
            poly = rve.get("central_cell_polygon_rve")
            if poly:
                P = np.asarray(poly, float) + np.array([x0, y0])
                ax.add_patch(Polygon(P, closed=True, fc=CELL_FILL, ec="none", zorder=0.5))
        ax.add_patch(Rectangle((x0, y0), w, h, fill=False, ec=RVE_CHOSEN, lw=2.0, zorder=5, label="FE RVE"))
        ax.plot([x0 + w / 2], [y0 + h / 2], "x", color=RVE_CHOSEN, ms=7, mew=1.6, zorder=6)
    if drag is not None:
        dx0, dy0, dx1, dy1 = drag
        ax.add_patch(Rectangle((min(dx0, dx1), min(dy0, dy1)), abs(dx1 - dx0), abs(dy1 - dy0),
                               fc=RVE_CHOSEN, alpha=0.08, ec=RVE_CHOSEN, lw=1.0, ls=":", zorder=6))
    ax.plot([], [], color=WALL, lw=1.7, label="walls in the FE RVE")
    if rve is not None and rve.get("face_walls"):
        ax.plot([], [], color=FACE_WALL, lw=2.4, label="wall in a face (half laminate)")
    if grid:
        ax.plot([], [], color="0.85", lw=0.8, label="rectangular periodic cells")
    if centres:
        ax.plot([], [], "+", color="0.6", ms=5, label="cell centres (snap points)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=7, frameon=False, borderaxespad=0.0)
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    ax.set_aspect("equal")
    ax.set_xlabel("x1 along a1 (mm)", fontsize=8)
    ax.set_ylabel("x2 (mm)", fontsize=8)
    ax.tick_params(labelsize=7)
    for sp in ax.spines.values():
        sp.set_color("0.8")
    if title:
        ax.set_title(title, fontsize=9)


# ----------------------------------------------------------------------------- boundary conditions
def plot_bc_panel(ax, rve: dict, case: str, displacement: float = 1.0) -> None:
    """Plan with the Table 3.4 / Table 3 boundary conditions of one DBC case."""
    from .abaqus_thesis_dbc import _case_constraints
    p, q = float(rve["p"]), float(rve["q"])
    s = max(p, q)
    plot_rve_plan(ax, rve, dims=False, faces=False, legend=False, lw=1.2)
    off = 0.035 * s
    geom = {"X1_0": ((-off, 0), (-off, q), np.array([-1.0, 0.0])), "X1_P": ((p + off, 0), (p + off, q), np.array([1.0, 0.0])),
            "X2_0": ((0, -off), (p, -off), np.array([0.0, -1.0])), "X2_Q": ((0, q + off), (p, q + off), np.array([0.0, 1.0]))}
    per_face: dict[str, list] = {}
    for f, dof, val in _case_constraints(case, displacement):
        if f in geom:
            per_face.setdefault(f, []).append((dof, val))
    for f, items in per_face.items():
        (x0, y0), (x1, y1), nrm = geom[f]
        ax.plot([x0, x1], [y0, y1], color=FACE_COLOURS[f], lw=2.2, solid_capstyle="butt")
        labels = []
        for dof, val in items:
            if val != 0 and dof in (1, 2):
                v = np.array([1.0, 0.0]) if dof == 1 else np.array([0.0, 1.0])
                for frac in (0.25, 0.75):
                    base = np.array([x0, y0]) + frac * (np.array([x1, y1]) - np.array([x0, y0])) + 0.012 * s * nrm
                    ax.add_patch(FancyArrowPatch(base, base + 0.1 * s * v, arrowstyle="-|>", mutation_scale=8,
                                                 color="k", lw=1.0, zorder=5))
                labels.append(f"U{dof}=d")
            elif val != 0 and dof == 3:
                for frac in (0.25, 0.75):
                    base = np.array([x0, y0]) + frac * (np.array([x1, y1]) - np.array([x0, y0])) + 0.04 * s * nrm
                    ax.add_patch(Circle(base, 0.018 * s, color="k", zorder=6))
                labels.append("U3=d")
            else:
                labels.append(f"U{dof}=0")
        if abs(nrm[0]) > 0:
            mid = np.array([x0, q / 2]) + nrm * 0.16 * s
            ax.text(mid[0], mid[1], "\n".join(labels), fontsize=6.5, ha="center", va="center", color=FACE_COLOURS[f])
        else:
            mid = np.array([p / 2, y0]) + nrm * 0.1 * s
            ax.text(mid[0], mid[1], "  ".join(labels), fontsize=6.5, ha="center", va="center", color=FACE_COLOURS[f])
    note = ("top/bottom edges: U1 = U2 = 0" if case in ("Q44", "Q55") else
            "top/bottom edges: U3 = 0" if case.startswith("C") else "top/bottom free (+ one U3 gauge)")
    ax.text(p / 2, -0.25 * s, note, ha="center", fontsize=6.5, color="0.3")
    ax.set_xlim(-0.32 * s, p + 0.32 * s)
    ax.set_ylim(-0.3 * s, q + 0.25 * s)


# ----------------------------------------------------------------------------- 3D shapes
def _clean3d(ax, elev=28, azim=-62):
    ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()
    try:
        ax.set_proj_type("ortho")
    except Exception:
        pass


def plot_deformed(ax, coords, elements, u, *, component=(0,), label: str = "", scale: float | None = None,
                  cmap: str = "turbo", vmin: float | None = None, vmax: float | None = None) -> float:
    """Deformed shell mesh coloured by a displacement component; returns the magnification used."""
    X = np.asarray(coords, float)
    U = np.asarray(u, float)
    E = np.asarray(elements, int)
    span = float(np.max(np.ptp(X, axis=0)))
    umax = float(np.max(np.abs(U))) if U.size else 0.0
    if scale is None:
        scale = 0.08 * span / umax if umax > 0 else 1.0
        scale = float(f"{scale:.2g}")
    val = U[:, list(component)].mean(axis=1)
    lo = float(val.min()) if vmin is None else vmin
    hi = float(val.max()) if vmax is None else vmax
    if hi - lo < 1e-15:
        hi = lo + 1.0
    Xd = X + scale * U
    import matplotlib
    cm = matplotlib.colormaps[cmap]
    fv = val[E].mean(axis=1)
    polys = Xd[E]
    ax.add_collection3d(Poly3DCollection(polys, facecolors=cm((fv - lo) / (hi - lo)), edgecolor=(0, 0, 0, 0.08), linewidth=0.1))
    ax.set_xlim(Xd[:, 0].min(), Xd[:, 0].max())
    ax.set_ylim(Xd[:, 1].min(), Xd[:, 1].max())
    ax.set_zlim(Xd[:, 2].min(), Xd[:, 2].max())
    ax.set_box_aspect((np.ptp(Xd[:, 0]), np.ptp(Xd[:, 1]), max(np.ptp(Xd[:, 2]), 1e-9)))
    _clean3d(ax)
    if label:
        ax.set_title(f"{label}, x{scale:g}", fontsize=8)
    return scale


def plot_mesh(ax, coords, elements, *, colour="#56B4E9", lw=0.15) -> None:
    X = np.asarray(coords, float)
    E = np.asarray(elements, int)
    ax.add_collection3d(Poly3DCollection(X[E], facecolor=colour, edgecolor="0.25", linewidth=lw, alpha=0.9))
    ax.set_xlim(X[:, 0].min(), X[:, 0].max()); ax.set_ylim(X[:, 1].min(), X[:, 1].max()); ax.set_zlim(X[:, 2].min(), X[:, 2].max())
    ax.set_box_aspect((np.ptp(X[:, 0]), np.ptp(X[:, 1]), np.ptp(X[:, 2])))
    _clean3d(ax)


# ----------------------------------------------------------------------------- comparisons
def plot_summary_differences(ax, cols: list[str], table: list[list], *, title: str = "") -> None:
    """Grouped bars of every '(%)' column of a ``verification.summary_table``."""
    palette = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00"]
    terms = [r[0] for r in table]
    series = []
    for k, c in enumerate(cols):
        if c.endswith("(%)"):
            name = cols[k - 1] if c == "analytical vs FE (%)" else c[:-4]
            label = f"analytical vs {name}" if c == "analytical vs FE (%)" else name
            series.append((k, label))
    y = np.arange(len(terms))
    h = 0.8 / max(len(series), 1)
    vals_all = []
    for n, (k, lab) in enumerate(series):
        vals = [float(r[k]) if isinstance(r[k], (int, float)) else math.nan for r in table]
        off = (n - (len(series) - 1) / 2) * h
        ax.barh(y + off, vals, height=h * 0.92, color=palette[n % len(palette)], label=lab)
        for yy, v in zip(y + off, vals):
            if math.isfinite(v):
                vals_all.append(abs(v))
                ax.text(v, yy, " " + signed(v, 1) + " ", va="center", ha="left" if v >= 0 else "right", fontsize=6.5)
    ax.axvline(0, color="0.3", lw=0.8)
    for x in (-1, 1):
        ax.axvline(x, color="0.75", lw=0.6, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(terms)
    ax.invert_yaxis()
    ax.set_xlabel("difference (%) = 100 (first − second) / second;  dashed: ±1 %", fontsize=7.5)
    if series:
        leg = ax.legend(fontsize=7, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=1,
                        title=title or None, title_fontsize=8.5)
        leg._legend_box.align = "center"
    elif title:
        ax.set_title(title, fontsize=8.5)
    lim = max(vals_all + [2.0])
    ax.set_xlim(-1.35 * lim, 1.35 * lim)
    ax.tick_params(labelsize=7.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def plot_ch_difference(ax, C_ref, C_fe, *, title: str = "") -> None:
    """6x6 map of (C_ref - C_fe) / |C_fe| in % on entries that are not negligible (NaN = not identified)."""
    A = np.asarray(C_ref, float)
    C = np.asarray(C_fe, float)
    scale = max(float(np.nanmax(np.abs(C))), 1e-30)
    D = np.full((6, 6), np.nan)
    for i in range(6):
        for j in range(6):
            if math.isfinite(C[i, j]) and max(abs(A[i, j]), abs(C[i, j])) >= 1e-6 * scale and abs(C[i, j]) > 1e-12 * scale:
                D[i, j] = 100.0 * (A[i, j] - C[i, j]) / abs(C[i, j])
    lim = max(2.0, float(np.nanmax(np.abs(D))) if np.isfinite(D).any() else 2.0)
    im = ax.imshow(D, cmap="RdBu_r", vmin=-lim, vmax=lim)
    lab = ["11", "22", "33", "23", "13", "12"]
    for (i, j), v in np.ndenumerate(D):
        if math.isfinite(v):
            ax.text(j, i, signed(v, 1), ha="center", va="center", fontsize=6.5, color="k" if abs(v) < 0.6 * lim else "w")
        else:
            ax.text(j, i, "·", ha="center", va="center", fontsize=7, color="0.6")
    ax.set_xticks(range(6)); ax.set_xticklabels(lab, fontsize=7)
    ax.set_yticks(range(6)); ax.set_yticklabels(lab, fontsize=7)
    for sp in ax.spines.values():
        sp.set_visible(False)
    if title:
        ax.set_title(title, fontsize=8)
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.ax.tick_params(labelsize=6)
    cb.set_label("difference (%)", fontsize=7)


# ----------------------------------------------------------------------------- composite figures
def project_figure(rve: dict | None, tables: dict, fig: Figure | None = None) -> Figure:
    """FE RVE plan and the difference bars of the Q-matrix, engineering-constant and DBC tables (up to three)."""
    fig = fig or Figure(figsize=(12, 4.8))
    fig.clear()
    groups = [g for g in ("Q", "eng", "shell", "3d") if g in tables and tables[g][1]][:3]
    n = (1 if rve else 0) + len(groups)
    if n == 0:
        return fig
    widths = ([0.8] if rve else []) + [1.15] * len(groups)
    gs = fig.add_gridspec(1, n, width_ratios=widths)
    k = 0
    if rve:
        plot_rve_plan(fig.add_subplot(gs[0, 0]), rve, title="FE RVE")
        k = 1
    from .verification import GROUP_TITLES
    for g in groups:
        cols, table = tables[g][0], tables[g][1]
        plot_summary_differences(fig.add_subplot(gs[0, k]), cols, table, title=GROUP_TITLES[g])
        k += 1
    fig.subplots_adjust(left=0.03, right=0.98, top=0.72, bottom=0.14, wspace=0.42)
    return fig


def deformed_figure(fields: dict | None, model: str = "DBC", fig: Figure | None = None,
                    cases: list[str] | None = None, solver: str = "built-in") -> Figure:
    """Deformed shapes of one model: the ten DBC cases or the six PBC unit strains.

    ``fields`` = {'coords', 'elements', 'u': {case: N x 3}} (``fe_workflow.load_fields``
    gives one per model; ``verification.builtin_check(..., fields=True)['fields']``
    too).
    """
    fig = fig or Figure(figsize=(12, 5.5))
    fig.clear()
    if not fields or not fields.get("u"):
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "No displacement fields yet (run the check)", ha="center", va="center")
        ax.axis("off")
        return fig
    keys = cases or list(fields["u"].keys())
    ncol = min(6 if model == "PBC" else 5, len(keys))
    nrow = int(math.ceil(len(keys) / ncol))
    for k, key in enumerate(keys):
        ax = fig.add_subplot(nrow, ncol, k + 1, projection="3d")
        if model == "DBC":
            comp, lab = CASE_FIELDS.get(key, ((0,), "u"))
            title = f"{key}: {lab}"
        else:
            c, lab = PBC_FIELDS.get(key, (0, "u"))
            comp = (c,)
            title = f"{PBC_TITLES.get(key, key)}: {lab}"
        plot_deformed(ax, fields["coords"], fields["elements"], fields["u"][key], component=comp, label=title)
    fig.subplots_adjust(left=0.0, right=1.0, top=0.94, bottom=0.02, wspace=0.0, hspace=0.12)
    return fig


def figure_png(fig: Figure, dpi: int = 130) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    return buf.getvalue()


def figure_data_uri(fig: Figure, dpi: int = 130) -> str:
    return "data:image/png;base64," + base64.b64encode(figure_png(fig, dpi)).decode("ascii")
