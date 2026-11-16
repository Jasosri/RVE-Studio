"""Wall-local axes: what Node I -> Node J means, and a drawing of it.

Every wall runs from its **Node I** (start) to its **Node J** (end).  That order fixes
the wall's local frame, used by the laminate, the strain-energy wall element, the wall
load recovery and the buckling screen:

* 1̄ (1bar) -- along the wall, positive from Node I towards Node J;
* 3̄ (3bar) -- along the core height (x3);
* 2̄ (2bar) -- the wall normal, completing the right-handed frame (2̄ = 3̄ x 1̄, i.e. 1̄
  turned +90 deg in the x1-x2 plane).

Ply angles are measured from 1̄.  Reversing a wall (J -> I) reverses 1̄ and 2̄, which
mirrors +/-theta plies and turns an unsymmetric stack over: for off-axis or unsymmetric
laminates the effective properties change (up to about 40 % of the largest C^H entry in
the 14 benchmark configurations).  The benchmark honeycomb is defined with every wall
pointing *into* the three-wall junction; :func:`default_honeycomb_orientation` reads it
from the benchmark model.

Nothing here needs Qt: the GUI embeds the drawing, tests save it.
"""
from __future__ import annotations

import math
from typing import Iterable, Sequence

import numpy as np
from matplotlib.patches import FancyArrowPatch, Polygon

AXIS_1 = "1̄"
AXIS_2 = "2̄"
AXIS_3 = "3̄"
# Matplotlib (mathtext) and Qt rich-text forms of the same symbols.
MT = {1: r"$\bar{1}$", 2: r"$\bar{2}$", 3: r"$\bar{3}$"}
HTML = {k: f'<span style="text-decoration:overline">{k}</span>' for k in (1, 2, 3)}


def html_note() -> str:
    """Rich-text version of the convention for Qt labels."""
    b1, b2, b3 = HTML[1], HTML[2], HTML[3]
    return (f"<b>Node I → Node J</b> sets each wall's positive local <b>{b1}</b> axis (along the wall); "
            f"{b3} is the core height and {b2} the wall normal. Ply angles are measured from {b1}, so reversing I "
            f"and J changes off-axis and unsymmetric laminates. " + default_orientation_text())

SHORT_NOTE = ("Node I → Node J sets each wall's positive local 1̄ axis (along the wall). 3̄ is the core height and "
              "2̄ the wall normal. Ply angles are measured from 1̄, so reversing I and J changes off-axis and "
              "unsymmetric laminates.")

COLUMN_TIPS = {
    "Node I": "Start node of the wall. The local 1̄ axis points from Node I to Node J.",
    "Node J": "End node of the wall. The local 1̄ axis points from Node I to Node J.",
    "Angle": "Angle of the local 1̄ axis (Node I → Node J) from the x1 axis, counter-clockwise, degrees.",
    "Length": "Wall length between Node I and Node J.",
    "Stack ID": "Laminate stack of the wall. Ply angles are measured from the wall's 1̄ axis.",
    "Weight": "Wall weight (1 = a full physical wall; 0.5 = a wall shared with the neighbouring cell).",
    "Wall group": "Free label used to group walls (e.g. in results and the Geometry Inspector).",
}


def default_honeycomb_orientation() -> list[tuple[int, int]]:
    """(Node I, Node J), 1-based, of the benchmark honeycomb walls (benchmark configuration 1)."""
    from .examples_data import paper_hexagonal_case
    m, _mats, _lams = paper_hexagonal_case(1)
    return [(int(e.node_i) + 1, int(e.node_j) + 1) for e in m.elements]


def default_orientation_text() -> str:
    pairs = default_honeycomb_orientation()
    junction = {j for _i, j in pairs}
    tail = (f" (every wall points into junction node {next(iter(junction))})" if len(junction) == 1 else "")
    return ("Default honeycomb (benchmark configurations 1–14): "
            + ", ".join(f"W{k} = N{i} → N{j}" for k, (i, j) in enumerate(pairs, 1)) + tail + ".")


def long_note() -> str:
    """Plain-text explanation for help panels and dialogs."""
    return (f"{SHORT_NOTE}\n\n{default_orientation_text()} This orientation reproduces all 84 reference "
            "stiffness values within 0.3 %; reversing a wall changes the results of the off-axis configurations.")


def wall_frames(nodes: Sequence[Sequence[float]], walls: Iterable[Sequence]) -> list[dict]:
    """Local frame of every wall.

    ``walls`` rows are (node_i, node_j, ...) **1-based** (as in the Wall table).
    Each entry: id, node_i, node_j, length, angle_deg (of 1̄ from x1), e1 (1̄), e2 (2̄),
    mid (mid-point), and ``ok`` (False for a missing node or zero length).
    """
    P = np.asarray(nodes, float).reshape(-1, 2)
    out = []
    for k, row in enumerate(walls, 1):
        try:
            i, j = int(row[0]), int(row[1])
            a, b = P[i - 1], P[j - 1]
            if i < 1 or j < 1:
                raise IndexError
        except (IndexError, TypeError, ValueError):
            out.append({"id": k, "node_i": row[0] if row else None, "node_j": row[1] if len(row) > 1 else None,
                        "ok": False, "length": math.nan, "angle_deg": math.nan})
            continue
        d = b - a
        L = float(np.hypot(*d))
        if L <= 1e-12:
            out.append({"id": k, "node_i": i, "node_j": j, "ok": False, "length": 0.0, "angle_deg": math.nan})
            continue
        e1 = d / L
        out.append({"id": k, "node_i": i, "node_j": j, "ok": True, "length": L,
                    "angle_deg": math.degrees(math.atan2(d[1], d[0])), "e1": e1,
                    "e2": np.array([-e1[1], e1[0]]), "mid": 0.5 * (a + b), "a": a, "b": b})
    return out


def reverse_walls(walls: list[Sequence], rows: Iterable[int]) -> list[list]:
    """Copy of ``walls`` (1-based rows) with Node I and Node J swapped in the given 0-based rows."""
    rows = set(int(r) for r in rows)
    out = []
    for k, w in enumerate(walls):
        w = list(w)
        if k in rows and len(w) >= 2:
            w[0], w[1] = w[1], w[0]
        out.append(w)
    return out


def draw_wall_axes(ax, nodes: Sequence[Sequence[float]], walls: Iterable[Sequence], a1=None, a2=None, *,
                   highlight: Iterable[int] = (), title: str | None = None, fontsize: float = 9,
                   neighbours: bool = True) -> list[dict]:
    """Draw the walls with an arrow Node I → Node J and the local 1̄, 2̄ axes at mid-wall.

    ``walls`` rows are 1-based (node_i, node_j, ...); ``highlight`` holds 0-based wall rows.
    Returns :func:`wall_frames`.
    """
    P = np.asarray(nodes, float).reshape(-1, 2)
    frames = wall_frames(P, walls)
    good = [f for f in frames if f["ok"]]
    hl = set(int(h) for h in highlight)
    pts = np.vstack([P] + [np.array([f["a"], f["b"]]) for f in good]) if len(P) else np.zeros((1, 2))
    scale = float(np.ptp(pts, axis=0).max()) if len(pts) > 1 else 1.0
    scale = scale if scale > 0 else 1.0
    if neighbours and a1 is not None and a2 is not None:
        A1, A2 = np.asarray(a1, float), np.asarray(a2, float)
        for s1 in (-1, 0, 1):
            for s2 in (-1, 0, 1):
                if s1 == s2 == 0:
                    continue
                t = s1 * A1 + s2 * A2
                for f in good:
                    ax.plot([f["a"][0] + t[0], f["b"][0] + t[0]], [f["a"][1] + t[1], f["b"][1] + t[1]],
                            color="0.86", lw=1.0, zorder=1)
        o = P.mean(axis=0) - 0.5 * (A1 + A2)
        ax.add_patch(Polygon(np.array([o, o + A1, o + A1 + A2, o + A2]), closed=True, fill=False, ec="0.7",
                             lw=0.7, ls=(0, (4, 3)), zorder=0))
    arrow = 0.10 * scale
    for f in good:
        k = f["id"] - 1
        col = "#D55E00" if k in hl else "#172033"
        lw = 2.6 if k in hl else 1.8
        # The wall arrow itself is the local 1̄ axis: Node I -> Node J.
        ax.add_patch(FancyArrowPatch(f["a"], f["b"], arrowstyle="-|>", mutation_scale=14 if k in hl else 12,
                                     color=col, lw=lw, shrinkA=5, shrinkB=10, zorder=3))
        m, e1, e2, L = f["mid"], f["e1"], f["e2"], f["length"]
        # 1̄ label beside the arrow, three quarters along; 2̄ normal arrow at mid-wall.
        q1 = f["a"] + 0.72 * L * e1 + 0.045 * scale * e2
        ax.text(q1[0], q1[1], MT[1], color="#0072B2", fontsize=fontsize * 1.05, ha="center", va="center",
                zorder=6, fontweight="bold")
        ax.add_patch(FancyArrowPatch(m, m + arrow * e2, arrowstyle="-|>", mutation_scale=9, color="#009E73",
                                     lw=1.3, shrinkA=0, shrinkB=0, zorder=5))
        q2 = m + 1.3 * arrow * e2
        ax.text(q2[0], q2[1], MT[2], color="#009E73", fontsize=fontsize * 1.05, ha="center", va="center",
                zorder=6, fontweight="bold")
        lab = m - 0.075 * scale * e2
        ax.text(lab[0], lab[1], f"W{f['id']}: N{f['node_i']}→N{f['node_j']}", fontsize=fontsize * 0.85, color=col,
                ha="center", va="center", zorder=6, rotation=_readable(f["angle_deg"]), rotation_mode="anchor",
                bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.4))
    for n, (x, y) in enumerate(P, 1):
        ax.plot([x], [y], "o", ms=5, mfc="white", mec="#172033", zorder=7)
        ax.text(x, y, f"  N{n}", fontsize=fontsize * 0.9, color="#172033", ha="left", va="bottom", zorder=8)
    ax.plot([], [], color="#172033", lw=1.8, marker=">", ms=5,
            label=f"wall arrow = positive local {MT[1]} axis, Node I → Node J")
    ax.plot([], [], color="#009E73", lw=1.3, label=f"{MT[2]}: wall normal; {MT[3]} = core height (out of page)")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.02), fontsize=fontsize * 0.85, frameon=False, ncol=1)
    c0 = 0.5 * (pts.min(axis=0) + pts.max(axis=0))
    span = 0.5 * scale + 0.25 * scale
    ax.set_xlim(c0[0] - span, c0[0] + span)
    ax.set_ylim(c0[1] - span, c0[1] + span)
    ax.set_aspect("equal")
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=fontsize * 1.1)
    return frames


def _readable(angle_deg: float) -> float:
    """Text rotation along a wall that never reads upside down."""
    a = (angle_deg + 180.0) % 360.0 - 180.0
    if a > 90:
        a -= 180
    elif a < -90:
        a += 180
    return a
