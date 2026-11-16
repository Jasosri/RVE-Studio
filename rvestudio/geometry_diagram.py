"""Parameter diagrams for the Geometry Library.

A diagram shows one geometry family at its reference parameter values: the primitive
cell's walls in bold, each parameter's walls in their own colour, a dimension line on
one representative wall per parameter, neighbouring cells faded for context, the
unit cell (dashed) and the lattice vectors a1, a2.  It is drawn from the family's own
generator and from the same parameter -> wall evidence the Geometry Inspector reports
(``geometry_analysis``), so it always matches what the parameters actually control:

* built-in families (HEX, SQUARE, MIXED, TRI): wall-length parameters on their walls;
  the ratio R is drawn on the inclined walls as Ld = R·Lh;
* wizard-created W/H families: W and H as the drawing width and height;
* any other plug-in: the geometry at its reference values with the parameter values listed.

Nothing here needs Qt: the GUI embeds the Figure, tests save it.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Mapping

import re
import numpy as np
from matplotlib.figure import Figure
from matplotlib.patches import Polygon

from .geometry import Geometry, BUILTIN_GEOMETRIES
from .models import Element

PALETTE = ("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
FADED = "0.82"
OTHER_WALL = "0.35"
CELL_FILL = "#F2F6FB"

# How each built-in parameter is drawn: label on the dimension line and a sentence.
PARAMETER_TEXT = {
    ("HEX", "L1"): ("L1", "length of the wall normal to a1"),
    ("HEX", "L2"): ("L2", "length of the inclined walls"),
    ("SQUARE", "Lx"): ("Lx", "length of each wall segment parallel to a1 (junction to cell edge)"),
    ("SQUARE", "Ly"): ("Ly", "length of each wall segment normal to a1 (junction to cell edge)"),
    ("MIXED", "Lh"): ("Lh", "length of the reference walls parallel to a1"),
    ("MIXED", "R"): ("Ld = R·Lh", "R = Ld/Lh, the inclined-wall length Ld over Lh"),
    ("TRI", "Lh"): ("Lh", "length of the reference walls parallel to a1"),
    ("TRI", "R"): ("Ld = R·Lh", "R = Ld/Lh, the inclined-wall length Ld over Lh"),
}


# ----------------------------------------------------------------------------- data
def reference_values(entry: Mapping) -> dict[str, float]:
    """Reference value of every library parameter: its 'reference', else the middle of its range."""
    out = {}
    for p in entry.get("parameters", []) or []:
        lo, hi = float(p.get("min", 0.0)), float(p.get("max", 0.0))
        ref = p.get("reference")
        out[str(p["name"])] = float(ref) if ref not in (None, "") else 0.5 * (lo + hi)
    return out


def library_geometry(entry: Mapping, registry: Mapping | None = None,
                     values: Mapping[str, float] | None = None) -> tuple[Geometry, dict[str, float]]:
    """Generate a library entry's geometry at ``values`` (default: reference values)."""
    reg = dict(registry) if registry is not None else dict(BUILTIN_GEOMETRIES)
    gid = str(entry.get("function") or entry.get("id", "")).upper()
    if gid not in reg:
        gid = str(entry.get("id", "")).upper()
    if gid not in reg:
        raise ValueError(f"No generator for geometry {entry.get('id')!r} (is its plug-in in the plugins folder?).")
    vals = dict(reference_values(entry))
    if values:
        vals.update({k: float(v) for k, v in values.items()})
    return reg[gid](vals), vals


def parameter_walls(geom: Geometry, values: Mapping[str, float]) -> tuple[str, dict[str, list[int]], str]:
    """(kind, {parameter: 0-based wall indices}, family id) for the diagram.

    kind = "walls" (built-in families, from the Geometry Inspector's evidence),
           "box"   (wizard W/H families: drawing width and height),
           "none"  (other plug-ins).
    """
    from .geometry_analysis import analyse_geometry
    elements = [Element(int(i), int(j), 0, float(w)) for (i, j), w in zip(geom.connectivity, geom.wall_weights)]
    try:
        s = analyse_geometry(np.asarray(geom.nodes, float), elements, geom.a1, geom.a2)
    except Exception:
        s = None
    names = list(values)
    if s is not None and s.detected_geometry_id and s.parameter_sources and \
            set(s.parameter_sources) >= set(names) and names:
        return "walls", {k: [int(i) - 1 for i in s.parameter_sources[k]] for k in names}, s.detected_geometry_id
    if "W/H" in (geom.note or "") and len(names) >= 2:
        return "box", {}, str(geom.id)
    # Plug-ins: a parameter controls the walls whose group name contains it as a token ('inclined-L2-A' -> L2).
    groups = [str(g) for g in (geom.wall_groups or [])]
    if names and len(groups) == len(geom.connectivity):
        by = {n: [k for k, g in enumerate(groups) if n in re.split(r"[-_ ]+", g)] for n in names}
        if all(by.values()):
            return "walls", by, str(geom.id)
    return "none", {}, str(geom.id)


# ----------------------------------------------------------------------------- drawing helpers
def _walls(geom: Geometry) -> list[tuple[np.ndarray, np.ndarray]]:
    P = np.asarray(geom.nodes, float)
    return [(P[int(i)], P[int(j)]) for i, j in geom.connectivity]


def _seg_dist(x: np.ndarray, p: np.ndarray, q: np.ndarray) -> float:
    d = q - p
    dd = float(d @ d)
    t = 0.0 if dd <= 0 else float(np.clip((x - p) @ d / dd, 0.0, 1.0))
    return float(np.linalg.norm(x - (p + t * d)))


def _outward(p: np.ndarray, q: np.ndarray, centre: np.ndarray) -> np.ndarray:
    d = q - p
    n = np.array([-d[1], d[0]]) / max(float(np.linalg.norm(d)), 1e-300)
    return n if float(((p + q) / 2 - centre) @ n) >= 0 else -n


def _clearance(segments, k: int, off: float, n: np.ndarray, avoid=()) -> float:
    """Distance from the dimension line of wall k (offset along n) to every other wall and label."""
    p, q = segments[k]
    pts = [p + t * (q - p) + n * off for t in (0.15, 0.5, 0.85)] + [(p + q) / 2 + n * 1.55 * off]
    others = [s for i, s in enumerate(segments) if i != k]
    d = min((_seg_dist(x, a, b) for x in pts for a, b in others), default=1e9)
    for c in avoid:
        d = min(d, min(float(np.linalg.norm(x - c)) for x in pts))
    return d


def _best_dimension(segments, walls, off: float, centre: np.ndarray, avoid=()):
    """(wall, normal) whose dimension line sits furthest from the other walls; outward wins ties."""
    best = None
    for k in walls:
        p, q = segments[k]
        out = _outward(p, q, centre)
        for n, bonus in ((out, 1e-6), (-out, 0.0)):
            c = _clearance(segments, k, off, n, avoid) + bonus
            if best is None or c > best[0]:
                best = (c, k, n)
    return best[1], best[2]


def _dimension(ax, p, q, label: str, colour: str, off: float, n: np.ndarray, *, fontsize: float = 10) -> np.ndarray:
    a, b = p + n * off, q + n * off
    ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="<|-|>", color=colour, lw=1.1,
                                                    shrinkA=0, shrinkB=0, mutation_scale=9), zorder=6)
    for e, s in ((p, a), (q, b)):
        ax.plot([e[0], s[0] + n[0] * 0.25 * off], [e[1], s[1] + n[1] * 0.25 * off], color=colour, lw=0.6,
                ls=(0, (2, 2)), zorder=5)
    m = (a + b) / 2 + n * 0.55 * off
    ax.text(m[0], m[1], label, color=colour, ha="center", va="center", fontsize=fontsize, fontweight="bold",
            zorder=7, bbox=dict(fc="white", ec="none", alpha=0.9, pad=0.6))
    return m


def _lattice_key(fig_or_ax, ax, geom: Geometry) -> None:
    ins = ax.inset_axes([0.0, 0.0, 0.22, 0.22])
    a1, a2 = np.asarray(geom.a1, float), np.asarray(geom.a2, float)
    s = max(np.linalg.norm(a1), np.linalg.norm(a2))
    for v, name in ((a1, "a1"), (a2, "a2")):
        u = v / s
        ins.annotate("", xy=u, xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color="k", lw=1.1,
                                                                 shrinkA=0, shrinkB=0, mutation_scale=9))
        ins.text(*(u * 1.22), name, fontsize=8, style="italic", ha="center", va="center")
    ins.patch.set_facecolor("white")
    ins.patch.set_alpha(0.85)
    ins.set_xlim(-0.9, 1.4)
    ins.set_ylim(-0.3, 1.4)
    ins.set_aspect("equal")
    ins.axis("off")


# ----------------------------------------------------------------------------- public
def draw_geometry_diagram(ax, geom: Geometry, values: Mapping[str, float], *, title: str | None = None,
                          descriptions: Mapping[str, str] | None = None, legend: bool = True,
                          fontsize: float = 10) -> dict:
    """Draw a parameter diagram on ``ax``.  Returns {"kind", "family", "parameters": [(name, text)]}."""
    kind, walls_of, fam = parameter_walls(geom, values)
    segs = _walls(geom)
    a1, a2 = np.asarray(geom.a1, float), np.asarray(geom.a2, float)
    P = np.asarray(geom.nodes, float)
    if segs:
        # Geometry centre: length-weighted centroid of the walls; the unit cell is drawn about it.
        Lw = np.array([max(float(np.linalg.norm(q - p)), 1e-300) for p, q in segs])
        centre = (np.array([(p + q) / 2 for p, q in segs]) * Lw[:, None]).sum(axis=0) / Lw.sum()
    else:
        centre = P.mean(axis=0)
    colour = {k: PALETTE[i % len(PALETTE)] for i, k in enumerate(walls_of)}
    owner = {w: k for k, ws in walls_of.items() for w in ws}

    # Neighbours (faded) then the primitive cell (bold, coloured by parameter).
    for i in (-2, -1, 0, 1, 2):          # two rings, so the clipped view is filled for sheared lattices too
        for j in (-2, -1, 0, 1, 2):
            if i == j == 0:
                continue
            t = i * a1 + j * a2
            for p, q in segs:
                ax.plot([p[0] + t[0], q[0] + t[0]], [p[1] + t[1], q[1] + t[1]], color=FADED, lw=1.1,
                        solid_capstyle="round", zorder=1)
    o = centre - 0.5 * (a1 + a2)
    ax.add_patch(Polygon(np.array([o, o + a1, o + a1 + a2, o + a2]), closed=True, fc=CELL_FILL,
                         ec="0.55", lw=0.8, ls=(0, (4, 3)), zorder=0))
    for k, (p, q) in enumerate(segs):
        c = colour.get(owner.get(k), OTHER_WALL)
        ax.plot([p[0], q[0]], [p[1], q[1]], color=c, lw=3.0, solid_capstyle="round", zorder=3)
        # wall number at one third from Node I: keeps it clear of the dimension label at the wall middle
        ax.text(*(p + (q - p) / 3), f"W{k + 1}", fontsize=fontsize * 0.62, color="0.35", ha="center", va="center",
                zorder=4, bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.3))

    scale = max(np.linalg.norm(a1), np.linalg.norm(a2))
    off = 0.09 * scale
    described: list[tuple[str, str]] = []
    desc = dict(descriptions or {})
    if kind == "walls":
        placed: list[np.ndarray] = []
        for name, ws in walls_of.items():
            if not ws:
                continue
            k, n = _best_dimension(segs, ws, off, centre, placed)
            label, text = PARAMETER_TEXT.get((fam, name), (name, desc.get(name, "")))
            v = float(values.get(name, math.nan))
            placed.append(_dimension(ax, segs[k][0], segs[k][1], label, colour[name], off, n, fontsize=fontsize))
            described.append((name, f"{text} — walls {', '.join(f'W{w + 1}' for w in ws)}; reference {v:.4g}"))
    elif kind == "box":
        names = list(values)
        lo, hi = P.min(axis=0), P.max(axis=0)
        # Width / height parameters are matched to the drawing size by value (wizard families).
        w_name = min(names, key=lambda n: abs(float(values[n]) - float(hi[0] - lo[0])))
        h_name = min([n for n in names if n != w_name] or names, key=lambda n: abs(float(values[n]) - float(hi[1] - lo[1])))
        for name, (p, q), c in ((w_name, (np.array([lo[0], lo[1]]), np.array([hi[0], lo[1]])), PALETTE[0]),
                                (h_name, (np.array([hi[0], lo[1]]), np.array([hi[0], hi[1]])), PALETTE[1])):
            _dimension(ax, p, q, name, c, off * 1.6, _outward(p, q, (lo + hi) / 2), fontsize=fontsize)
        described += [(w_name, f"drawing width (all x-coordinates and a1, a2 scale with it); reference {values[w_name]:.4g}"),
                      (h_name, f"drawing height (all y-coordinates and a1, a2 scale with it); reference {values[h_name]:.4g}")]
        for n in names:
            if n not in (w_name, h_name):
                described.append((n, f"{desc.get(n, 'plug-in parameter')}; reference {float(values[n]):.4g}"))
    else:
        for n, v in values.items():
            described.append((n, f"{desc.get(n, 'plug-in parameter (see the plug-in)')}; reference {float(v):.4g}"))

    if legend:
        for name, ws in walls_of.items():
            label = PARAMETER_TEXT.get((fam, name), (name,))[0]
            ax.plot([], [], color=colour[name], lw=3, label=f"{label}: " + ", ".join(f"W{w + 1}" for w in ws))
        if kind == "walls" and any(k not in owner for k in range(len(segs))):
            ax.plot([], [], color=OTHER_WALL, lw=3, label="other walls")
        ax.plot([], [], color=FADED, lw=1.1, label="neighbouring cells")
        ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=fontsize * 0.75, frameon=False,
                  borderaxespad=0.0, handlelength=1.6)
    _lattice_key(None, ax, geom)
    cell = np.array([o, o + a1, o + a1 + a2, o + a2] + [x for s in segs for x in s])
    c0 = 0.5 * (cell.min(axis=0) + cell.max(axis=0))
    span = 0.5 * float(np.ptp(cell, axis=0).max()) + 0.32 * scale
    ax.set_xlim(c0[0] - span, c0[0] + span)
    ax.set_ylim(c0[1] - span, c0[1] + span)
    ax.set_aspect("equal")
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=fontsize * 1.1)
    return {"kind": kind, "family": fam, "parameters": described}


def entry_diagram(entry: Mapping, registry: Mapping | None = None, fig: Figure | None = None,
                  values: Mapping[str, float] | None = None, *, fontsize: float = 10) -> tuple[Figure, dict]:
    """Figure with the parameter diagram of one Geometry Library entry."""
    fig = fig or Figure(figsize=(6.4, 4.6))
    fig.clear()
    ax = fig.add_subplot(111)
    geom, vals = library_geometry(entry, registry, values)
    desc = {}
    text = str(entry.get("description", ""))
    for p in entry.get("parameters", []) or []:
        if p.get("description"):
            desc[str(p["name"])] = str(p["description"])
    info = draw_geometry_diagram(ax, geom, vals, title=f"{entry.get('id', '')} — {entry.get('name', '')}",
                                 descriptions=desc, fontsize=fontsize)
    info["description"] = text
    info["ranges"] = {str(p["name"]): (float(p.get("min", math.nan)), float(p.get("max", math.nan)))
                      for p in entry.get("parameters", []) or []}
    fig.subplots_adjust(left=0.02, right=0.70, top=0.92, bottom=0.03)
    return fig, info


def caption_html(info: Mapping) -> str:
    """Rich-text caption for a diagram: one line per parameter with its meaning and range."""
    rows = []
    for name, text in info.get("parameters", []):
        lo, hi = info.get("ranges", {}).get(name, (math.nan, math.nan))
        rng = f" &nbsp;<span style='color:#5B6B82'>range {lo:.4g} … {hi:.4g}</span>" if math.isfinite(lo) else ""
        rows.append(f"<b>{name}</b>: {text}{rng}")
    return "<br>".join(rows) if rows else "Fixed geometry: no independent parameters."


def builtin_guide_figure(fig: Figure | None = None) -> Figure:
    """The four built-in families side by side (Parameter guide)."""
    from .geometry import DEFAULT_GEOMETRY_LIBRARY
    fig = fig or Figure(figsize=(12, 9))
    fig.clear()
    for k, entry in enumerate(DEFAULT_GEOMETRY_LIBRARY[:4]):
        ax = fig.add_subplot(2, 2, k + 1)
        geom, vals = library_geometry(entry)
        draw_geometry_diagram(ax, geom, vals, title=f"{entry['id']} — {entry['name']}", fontsize=9)
    fig.subplots_adjust(left=0.01, right=0.86, top=0.95, bottom=0.02, wspace=0.55, hspace=0.12)
    return fig


def save_library_diagrams(entries, out_dir: str | Path, registry: Mapping | None = None) -> list[Path]:
    """PNG diagram of every library entry (for documentation and reports)."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = []
    for e in entries:
        fig, _ = entry_diagram(e, registry)
        p = out / f"geometry_{str(e.get('id', 'X')).lower()}.png"
        fig.savefig(p, dpi=130)
        paths.append(p)
    return paths
