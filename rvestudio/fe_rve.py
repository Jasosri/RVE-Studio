"""Finite-height FE RVE for the source-study prescribed-displacement (DBC) check.

The thesis (Fig. 3.6, Tables 3.3-3.4) and the paper (Fig. 8, Table 3) do not use a
one-primitive-cell window.  For the honeycomb they use a rectangular block that
holds one complete hexagon in the middle:

* x1 runs along the zig-zag (a1) direction and x2 along the straight walls;
* the block is two rectangular periodic cells long in x1 and one in x2,
  i.e. p = 2|a1| and q = |a1 + 2 a2| (173.2 x 150 mm for l = 50 mm walls),
  four primitive cells in all;
* the x1 = 0 and x1 = p faces pass through junction columns, so the straight
  walls there lie *in* the faces; the x2 = 0 and x2 = q faces cut the straight
  walls at mid-length.

A wall lying in a face is shared with the neighbouring block.  It is modelled
with half of its laminate (each ply at half thickness, by default) so that the
block holds exactly the material of its primitive cells and the relative
density is preserved.  The reported Q44/Q55 values of the source study require exactly this: with
half face walls the antiplane-shear field of a regular honeycomb is exactly affine
and Q44 = Q55, as reported for the cross-ply cores; with full face walls Q55
comes out about 10 % too high.

The same construction is applied automatically to any periodic wall network
(the FE RVE used by every DBC and PBC model, the FE study and the built-in checks):

1. x1 is aligned with a1 and the smallest rectangular periodic cell is found;
2. the distinct open cells of the network are found and ranked by inscribed
   circle (largest first) -- the block is centred on a complete cell;
3. the block is the smallest of 1x1, 2x1, 1x2, 2x2, 3x1, 1x3, 3x2, 2x3 or 3x3
   rectangular cells that holds that cell strictly inside;
4. a block whose faces would leave a wall piece shorter than 5 % of the median
   wall is passed over for the next candidate (no slivers for the mesher);
5. the wall measure is checked: shared-face-weighted wall length = cells x
   primitive wall length, area = cells x primitive area.

For the built-in hexagonal family this reproduces the source-study RVE exactly; the
square, triangular and mixed families get 2x2, 2x1 and 2x2 blocks of four
primitive cells.  ``select_fe_rve`` is the public entry point; the older name
``build_published_fe_rve`` is kept.
"""
from __future__ import annotations

import math
from typing import Iterable

import numpy as np
from matplotlib.path import Path as _MplPath

from .models import RVEModel, Element

FACE_WALL_SECTIONS = ("half-thickness", "half-stiffness", "full")
CELL_CANDIDATES = ((1, 1), (2, 1), (1, 2), (2, 2), (3, 1), (1, 3), (3, 2), (2, 3), (3, 3))
SLIVER_RATIO = 0.05
# How the FE RVE of a design is chosen (RVEModel.fe_rve):
#   auto   -- the rule of the source study below (default; stored as None)
#   cells  -- a block of n1 x n2 rectangular periodic cells, centred on a complete cell
#   window -- a rectangle drawn over the replicated design: n1 x n2 rectangular cells
#             (so it stays exactly periodic) placed with its centre at a chosen point
FE_RVE_MODES = ("auto", "cells", "window")
MAX_RVE_CELLS = 12
RVE_RULE = ("smallest of 1x1 ... 3x3 rectangular periodic cells, centred on the largest complete cell, "
            "holding that cell strictly inside, with no wall piece shorter than 5 % of the median wall "
            "(thesis Fig. 3.6b / paper Fig. 8b)")
CELL_NAMES = {3: "triangle", 4: "quadrilateral", 5: "pentagon", 6: "hexagon", 8: "octagon"}


# ----------------------------------------------------------------------------------
# Small geometric helpers
# ----------------------------------------------------------------------------------

def _rotation_to_a1(model: RVEModel) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Proper rotation taking a1 to +x1; a2 is flipped (same lattice) if needed."""
    a1 = np.asarray(model.a1, float).reshape(2)
    a2 = np.asarray(model.a2, float).reshape(2)
    n1 = float(np.linalg.norm(a1))
    if not math.isfinite(n1) or n1 <= np.finfo(float).eps:
        raise ValueError("Lattice vector a1 must be non-zero.")
    e1 = a1 / n1
    e2 = np.array([-e1[1], e1[0]], float)
    R = np.vstack((e1, e2))
    a1r, a2r = R @ a1, R @ a2
    if a2r[1] < 0.0:
        a2r = -a2r  # same lattice, keeps the frame right-handed (no mirror image)
    if abs(float(a2r[1])) <= 1e-12 * max(n1, 1.0):
        raise ValueError("Lattice vectors a1 and a2 are parallel.")
    return R, a1r, a2r


def rectangular_cell(a1: np.ndarray, a2: np.ndarray, search_limit: int = 8, tol: float = 1e-8) -> dict:
    """Smallest rectangular periodic cell with b1 along +x1 (a1 must lie on +x1)."""
    A = np.column_stack((a1, a2))
    best = None
    lim = int(search_limit)
    for m in range(1, lim + 1):
        for n in range(-lim, lim + 1):
            b1 = A @ np.array([m, n], float)
            if abs(float(b1[1])) > tol * max(float(np.linalg.norm(b1)), 1.0) or b1[0] <= 0:
                continue
            for k in range(-lim, lim + 1):
                for l in range(1, lim + 1):
                    b2 = A @ np.array([k, l], float)
                    if abs(float(b2[0])) > tol * max(float(np.linalg.norm(b2)), 1.0) or b2[1] <= 0:
                        continue
                    count = abs(m * l - n * k)
                    if count == 0:
                        continue
                    score = (count, float(b1[0]) + float(b2[1]), (m, n, k, l))
                    if best is None or score < best[0]:
                        best = (score, (m, n), (k, l), b1, b2, count)
    if best is None:
        raise ValueError("No rectangular periodic cell aligned with a1 was found; this lattice "
                         "has no rectangular block, so the source-study DBC RVE cannot be built.")
    _, c1, c2, b1, b2, count = best
    return {"b1_coefficients": list(c1), "b2_coefficients": list(c2),
            "width": float(b1[0]), "height": float(b2[1]), "primitive_cells": int(count)}


def _tile_segments(nodes: np.ndarray, elements: Iterable[Element], a1: np.ndarray, a2: np.ndarray,
                   centre: np.ndarray, reach: float) -> list[tuple[np.ndarray, np.ndarray, int, int]]:
    """Copies of the primitive walls whose bounding box comes within ``reach`` of ``centre``."""
    A = np.column_stack((a1, a2))
    elems = list(elements)
    used = sorted({e.node_i for e in elems} | {e.node_j for e in elems})
    frac = np.linalg.solve(A, np.asarray(nodes, float)[used].T).T
    cf = np.linalg.solve(A, np.asarray(centre, float))
    span = np.abs(np.linalg.solve(A, np.eye(2))) @ np.array([reach, reach])
    i0 = math.floor(cf[0] - span[0] - frac[:, 0].max()) - 1
    i1 = math.ceil(cf[0] + span[0] - frac[:, 0].min()) + 1
    j0 = math.floor(cf[1] - span[1] - frac[:, 1].max()) - 1
    j1 = math.ceil(cf[1] + span[1] - frac[:, 1].min()) + 1
    out = []
    for ii in range(i0, i1 + 1):
        for jj in range(j0, j1 + 1):
            shift = A @ np.array([ii, jj], float)
            for wi, e in enumerate(elems):
                p0 = np.asarray(nodes[e.node_i], float) + shift
                p1 = np.asarray(nodes[e.node_j], float) + shift
                lo = np.minimum(p0, p1) - reach
                hi = np.maximum(p0, p1) + reach
                if np.all(centre >= lo) and np.all(centre <= hi):
                    out.append((p0, p1, wi, int(e.laminate)))
    return out


def _point_segment_distance(P: np.ndarray, S0: np.ndarray, S1: np.ndarray) -> np.ndarray:
    """Distances from points P (n,2) to segments S0-S1 (m,2) -> (n,m)."""
    d = S1 - S0
    L2 = np.maximum(np.einsum('ij,ij->i', d, d), np.finfo(float).tiny)
    t = np.clip(((P[:, None, :] - S0[None, :, :]) * d[None, :, :]).sum(-1) / L2[None, :], 0.0, 1.0)
    proj = S0[None, :, :] + t[..., None] * d[None, :, :]
    return np.linalg.norm(P[:, None, :] - proj, axis=-1)


def _face_polygon(segments: list[tuple[np.ndarray, np.ndarray]], point: np.ndarray, scale: float) -> np.ndarray:
    """Vertices (counter-clockwise) of the face of a planar wall network containing ``point``."""
    tol = 1e-9 * max(scale, 1.0)
    grid = 1e-7 * max(scale, 1.0)
    key = lambda xy: (int(round(float(xy[0]) / grid)), int(round(float(xy[1]) / grid)))
    ids: dict[tuple[int, int], int] = {}
    xy: list[np.ndarray] = []
    adj: dict[int, set[int]] = {}

    def nid(p):
        k = key(p)
        if k not in ids:
            ids[k] = len(xy); xy.append(np.asarray(p, float)); adj[ids[k]] = set()
        return ids[k]

    for p0, p1 in segments:
        a, b = nid(p0), nid(p1)
        if a != b:
            adj[a].add(b); adj[b].add(a)
    X = np.asarray(xy)
    # Cast a slightly tilted ray to find the first wall hit; that wall bounds the face.
    ang = 0.1234567
    r = np.array([math.cos(ang), math.sin(ang)])
    best = None
    for a in adj:
        for b in adj[a]:
            if b < a:
                continue
            p, q = X[a], X[b]
            d = q - p
            M = np.column_stack((r, -d))
            if abs(np.linalg.det(M)) < 1e-14:
                continue
            t, s = np.linalg.solve(M, p - point)
            if t > tol and -1e-12 <= s <= 1 + 1e-12 and (best is None or t < best[0]):
                best = (t, a, b)
    if best is None:
        raise ValueError("Could not locate the wall network around the cell centre.")
    _, a, b = best
    # Orient the first edge so that the face lies on its left.
    if (X[b, 0] - X[a, 0]) * (point[1] - X[a, 1]) - (X[b, 1] - X[a, 1]) * (point[0] - X[a, 0]) < 0:
        a, b = b, a
    start = (a, b)
    poly = [a]
    u, v = a, b
    for _ in range(10000):
        poly.append(v)
        back = math.atan2(X[u, 1] - X[v, 1], X[u, 0] - X[v, 0])
        best_w, best_turn = None, None
        for w in adj[v]:
            if w == u and len(adj[v]) > 1:
                continue
            ang_w = math.atan2(X[w, 1] - X[v, 1], X[w, 0] - X[v, 0])
            turn = (back - ang_w) % (2 * math.pi)  # clockwise sweep from the way back
            if turn <= 1e-12:
                turn += 2 * math.pi
            if best_turn is None or turn < best_turn:
                best_w, best_turn = w, turn
        u, v = v, best_w
        if (u, v) == start:
            break
    else:
        raise ValueError("The cell containing the RVE centre is not closed.")
    P = X[poly[:-1]]
    area = 0.5 * float(np.sum(P[:, 0] * np.roll(P[:, 1], -1) - np.roll(P[:, 0], -1) * P[:, 1]))
    if area <= 0:
        raise ValueError("Cell boundary walk failed (non-positive area).")
    return P


def _polygon_centroid(P: np.ndarray) -> np.ndarray:
    x, y = P[:, 0], P[:, 1]
    xn, yn = np.roll(x, -1), np.roll(y, -1)
    cr = x * yn - xn * y
    A = 0.5 * cr.sum()
    return np.array([((x + xn) * cr).sum() / (6 * A), ((y + yn) * cr).sum() / (6 * A)])


def _polygon_area(P: np.ndarray) -> float:
    x, y = P[:, 0], P[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def _simplify_polygon(P: np.ndarray, tol: float) -> np.ndarray:
    """Drop collinear vertices (walls split at mid-points) so the corner count is the cell's."""
    keep = []
    n = len(P)
    for i in range(n):
        a, b, c = P[i - 1], P[i], P[(i + 1) % n]
        cross = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])
        if abs(cross) > tol * max(np.linalg.norm(b - a) * np.linalg.norm(c - b), tol):
            keep.append(i)
    return P[keep] if len(keep) >= 3 else P


def _cell_name(P: np.ndarray) -> str:
    n = len(P)
    if n == 4:
        sides = np.linalg.norm(np.roll(P, -1, axis=0) - P, axis=1)
        e = np.roll(P, -1, axis=0) - P
        right = all(abs(float(e[i] @ e[(i + 1) % 4])) <= 1e-6 * sides[i] * sides[(i + 1) % 4] for i in range(4))
        equal = float(np.ptp(sides)) <= 1e-6 * float(sides.max())
        return "square" if (right and equal) else "rectangle" if right else "rhombus" if equal else "quadrilateral"
    return CELL_NAMES.get(n, f"{n}-gon")


def find_cells(nodes: np.ndarray, elements: list[Element], a1: np.ndarray, a2: np.ndarray,
               samples: int = 64, max_cells: int = 8) -> list[dict]:
    """Distinct open cells of a periodic wall network, largest inscribed circle first.

    Each entry: centroid, polygon (counter-clockwise), area, inradius (sampled),
    corners (after merging collinear vertices) and a name ("hexagon", ...).
    """
    A = np.column_stack((a1, a2))
    scale = max(float(np.linalg.norm(a1)), float(np.linalg.norm(a2)))
    origin = A @ np.array([0.5, 0.5])
    segs = _tile_segments(nodes, elements, a1, a2, origin, 2.5 * scale)
    if not segs:
        raise ValueError("The wall network has no walls near the cell.")
    S0 = np.array([s[0] for s in segs]); S1 = np.array([s[1] for s in segs])
    g = (np.arange(samples) + 0.5) / samples
    U, V = np.meshgrid(g, g, indexing='ij')
    P = (A @ np.vstack((U.ravel(), V.ravel()))).T
    dist = _point_segment_distance(P, S0, S1).min(axis=1)
    covered = dist <= 1e-9 * scale
    shifts = [A @ np.array([i, j], float) for i in (-2, -1, 0, 1, 2) for j in (-2, -1, 0, 1, 2)]
    cells: list[dict] = []
    for k in np.argsort(-dist):
        if covered[k] or dist[k] <= 1e-9 * scale:
            continue
        try:
            poly = _face_polygon([(s[0], s[1]) for s in segs], P[k], scale)
        except ValueError:
            covered[k] = True
            continue
        path = _MplPath(poly)
        hit = np.zeros(len(P), bool)
        for sh in shifts:
            hit |= path.contains_points(P - sh)
        hit[k] = True
        covered |= hit
        cen = _polygon_centroid(poly)
        frac = np.linalg.solve(A, cen)
        if any(np.linalg.norm(((np.linalg.solve(A, c["centroid"]) - frac + 0.5) % 1.0) - 0.5, ord=np.inf) < 1e-6
               for c in cells):
            continue
        simple = _simplify_polygon(poly, 1e-9)
        cells.append({"centroid": cen, "polygon": poly, "area": _polygon_area(poly), "inradius": float(dist[k]),
                      "corners": int(len(simple)), "name": _cell_name(simple)})
        if len(cells) >= max_cells or covered.all():
            break
    if not cells:
        raise ValueError("No open cell was found in the wall network.")
    return cells


def cell_centre(nodes: np.ndarray, elements: list[Element], a1: np.ndarray, a2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Centroid of the largest cell (max inscribed circle) and that cell's polygon."""
    c = find_cells(nodes, elements, a1, a2)[0]
    return c["centroid"], c["polygon"]


def _clip(p0: np.ndarray, p1: np.ndarray, L: float, W: float, tol: float):
    """Clip segment p0-p1 to the box [0,L]x[0,W]."""
    d = p1 - p0
    t0, t1 = 0.0, 1.0
    for ax, hi in ((0, L), (1, W)):
        if abs(d[ax]) <= tol:
            if p0[ax] < -tol or p0[ax] > hi + tol:
                return None
            continue
        a, b = (0.0 - p0[ax]) / d[ax], (hi - p0[ax]) / d[ax]
        t0, t1 = max(t0, min(a, b)), min(t1, max(a, b))
        if t0 > t1:
            return None
    q0 = np.clip(p0 + t0 * d, [0.0, 0.0], [L, W])
    q1 = np.clip(p0 + t1 * d, [0.0, 0.0], [L, W])
    if np.linalg.norm(q1 - q0) <= tol:
        return None
    return q0, q1


def _snap(v: float, targets: tuple[float, ...], tol: float) -> float:
    for t in targets:
        if abs(v - t) <= tol:
            return float(t)
    return float(v)


# ----------------------------------------------------------------------------------
# User-defined FE RVE (stored on the design as RVEModel.fe_rve)
# ----------------------------------------------------------------------------------

def normalise_fe_rve_spec(spec) -> dict | None:
    """Validate an FE RVE choice; ``None`` (or mode "auto") means the automatic rule.

    Accepted forms::

        None / {"mode": "auto"}
        {"mode": "cells",  "cells": [n1, n2], "centre_cell": k}      # centre_cell optional
        {"mode": "window", "cells": [n1, n2], "centre": [x1, x2]}    # FE frame, mm
    """
    if spec is None:
        return None
    if not isinstance(spec, dict):
        raise ValueError("FE RVE specification must be a dictionary or None.")
    mode = str(spec.get("mode", "auto")).strip().lower()
    if mode not in FE_RVE_MODES:
        raise ValueError(f"FE RVE mode must be one of {FE_RVE_MODES}, not {mode!r}.")
    if mode == "auto":
        return None
    try:
        n1, n2 = (int(round(float(c))) for c in spec.get("cells", ()))
    except (TypeError, ValueError):
        raise ValueError("FE RVE 'cells' must be two integers (rectangular cells along x1 and x2).") from None
    if min(n1, n2) < 1 or max(n1, n2) > MAX_RVE_CELLS:
        raise ValueError(f"FE RVE cell counts must be between 1 and {MAX_RVE_CELLS}; got {n1} x {n2}.")
    out: dict = {"mode": mode, "cells": [n1, n2]}
    if mode == "cells":
        k = spec.get("centre_cell")
        if k is not None:
            if int(k) < 0:
                raise ValueError("FE RVE 'centre_cell' must be 0 or more.")
            out["centre_cell"] = int(k)
    else:
        try:
            cx, cy = (float(c) for c in spec.get("centre", ()))
        except (TypeError, ValueError):
            raise ValueError("A drawn FE RVE needs 'centre': [x1, x2] in the FE frame.") from None
        if not (math.isfinite(cx) and math.isfinite(cy)):
            raise ValueError("FE RVE centre must be finite.")
        out["centre"] = [cx, cy]
    return out


def describe_fe_rve_spec(spec) -> str:
    """Short label of a design's FE RVE choice, e.g. for lists and file names."""
    s = normalise_fe_rve_spec(spec)
    if s is None:
        return "automatic"
    n1, n2 = s["cells"]
    if s["mode"] == "cells":
        return f"{n1} x {n2} cells" + (f", cell type {s['centre_cell']}" if "centre_cell" in s else "")
    return f"drawn {n1} x {n2} cells at ({s['centre'][0]:.4g}, {s['centre'][1]:.4g})"


def _cell_at(cell_list: list[dict], A: np.ndarray, point: np.ndarray) -> tuple[int, np.ndarray] | None:
    """The open cell (type index, polygon image) that contains ``point``, if any."""
    Ainv = np.linalg.inv(A)
    for ci, c in enumerate(cell_list):
        base = np.round(Ainv @ (point - c["centroid"]))
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                shift = A @ (base + np.array([di, dj], float))
                poly = c["polygon"] + shift
                if _MplPath(poly).contains_point(point):
                    return ci, poly
    return None


def fe_rve_layout(model: RVEModel, *, reach: float = 1.4) -> dict:
    """Everything needed to draw the replicated design in the FE frame and place an RVE on it.

    Returns the tiled walls (``segments``, rows x0, y0, x1, y1 in the FE frame, x1 along a1),
    the rectangular periodic cell size, the automatic RVE (centre, size, diagnostics),
    the centres of all open cells in view (for snapping) and the view extent.
    """
    R, a1, a2 = _rotation_to_a1(model)
    nodes = (R @ np.asarray(model.nodes, float).T).T
    elems = list(model.elements)
    A = np.column_stack((a1, a2))
    rect = rectangular_cell(a1, a2)
    cell_list = find_cells(nodes, elems, a1, a2)
    _block, auto = select_fe_rve(model, use_model_spec=False)
    auto_centre = np.asarray(auto["window_origin_local"], float) + 0.5 * np.array([auto["p"], auto["q"]])
    # A square view about the automatic RVE, ``reach`` times its larger side either way:
    # room to draw blocks of about three times the automatic one.
    r = reach * max(auto["p"], auto["q"], 2.0 * rect["width"], 2.0 * rect["height"])
    half = np.array([r, r])
    lo, hi = auto_centre - half, auto_centre + half
    segs = []
    for p0, p1, wi, lam in _tile_segments(nodes, elems, a1, a2, auto_centre, float(np.linalg.norm(half)) + 1.0):
        c = _clip(p0 - lo, p1 - lo, *(hi - lo), 1e-9)
        if c is not None:
            segs.append([*(c[0] + lo), *(c[1] + lo), wi, lam])
    Ainv = np.linalg.inv(A)
    centres = []
    for ci, c in enumerate(cell_list):
        f_lo = np.floor(Ainv @ (np.array([lo, [lo[0], hi[1]], hi, [hi[0], lo[1]]]) - c["centroid"]).T).min(axis=1) - 1
        f_hi = np.ceil(Ainv @ (np.array([lo, [lo[0], hi[1]], hi, [hi[0], lo[1]]]) - c["centroid"]).T).max(axis=1) + 1
        for i in range(int(f_lo[0]), int(f_hi[0]) + 1):
            for j in range(int(f_lo[1]), int(f_hi[1]) + 1):
                pt = c["centroid"] + A @ np.array([i, j], float)
                if np.all(pt >= lo) and np.all(pt <= hi):
                    centres.append([float(pt[0]), float(pt[1]), ci])
    S = np.asarray(segs, float).reshape(-1, 6)
    # Candidate geometry centres: cell centres first, then junctions and wall mid-points
    # (in triangular networks the mirror centres are junctions and wall mid-points).
    cand = [[x, y, 0] for x, y, _k in centres]
    seen = {(round(x, 6), round(y, 6)) for x, y, _k in cand}
    if len(S):
        pts = [(1, S[:, 0:2]), (1, S[:, 2:4]), (2, 0.5 * (S[:, 0:2] + S[:, 2:4]))]
        for kind, P in pts:
            for x, y in P:
                key = (round(float(x), 6), round(float(y), 6))
                if key not in seen and lo[0] < x < hi[0] and lo[1] < y < hi[1]:
                    seen.add(key)
                    cand.append([float(x), float(y), kind])
    return {"segments": S[:, :4], "centre_candidates": cand, "wall_index": S[:, 4].astype(int).tolist(), "laminate": S[:, 5].astype(int).tolist(),
            "rect_width": rect["width"], "rect_height": rect["height"], "rectangular_cell": rect,
            "auto": auto, "auto_centre": auto_centre.tolist(), "auto_cells": list(auto["cells_x1_x2"]),
            "cell_centres": centres, "cell_types": auto["cell_types"],
            "cell_base_centroids": [[float(c["centroid"][0]), float(c["centroid"][1])] for c in cell_list],
            "extent": [float(lo[0]), float(hi[0]), float(lo[1]), float(hi[1])],
            "rotation_to_local": R.tolist()}


def snap_to_cell(layout: dict, point) -> list[float]:
    """Nearest open-cell centre to ``point`` (FE frame)."""
    C = np.asarray(layout["cell_centres"], float).reshape(-1, 3)
    if not len(C):
        return [float(point[0]), float(point[1])]
    k = int(np.argmin(np.linalg.norm(C[:, :2] - np.asarray(point, float)[None, :], axis=1)))
    return [float(C[k, 0]), float(C[k, 1])]


def window_from_drag(layout: dict, x0: float, y0: float, x1: float, y1: float, *,
                     snap_centre: bool = True) -> dict:
    """FE RVE spec from a rectangle dragged over the replicated design.

    The size is rounded to whole rectangular periodic cells (so the RVE stays exactly
    periodic and holds a whole number of primitive cells).  The centre is moved to the
    nearest cell centre -- the geometry centre -- by default, so the RVE faces cut the
    walls at matching points as in the automatic RVE; ``snap_centre=False`` keeps the
    centre where it was drawn.
    """
    w, h = abs(float(x1) - float(x0)), abs(float(y1) - float(y0))
    n1 = int(min(MAX_RVE_CELLS, max(1, round(w / layout["rect_width"]))))
    n2 = int(min(MAX_RVE_CELLS, max(1, round(h / layout["rect_height"]))))
    centre = [0.5 * (float(x0) + float(x1)), 0.5 * (float(y0) + float(y1))]
    if snap_centre:
        centre = snap_to_cell(layout, centre)  # quick rule; centred_window() also checks mirror symmetry
    return {"mode": "window", "cells": [n1, n2], "centre": centre}


def centred_window(model: RVEModel, layout: dict, cells, centre, *, snap: bool = True,
                   search: int = 40) -> tuple[dict, str]:
    """Window spec centred on the geometry wherever possible.

    With ``snap`` the centre moves to the nearest *centre of the geometry* that gives a
    mirror-symmetric block without slivers: cell centres are tried first, then junctions
    and wall mid-points (the mirror centres of triangular networks), nearest first among
    the ``search`` closest candidates.  If none qualifies the nearest cell centre is used,
    and if that block cannot be built the centre stays where it was given.
    Returns (spec, note) -- the note says what was done, for the user.
    """
    raw = normalise_fe_rve_spec({"mode": "window", "cells": list(cells), "centre": [float(c) for c in centre]})
    if not snap:
        return raw, "Centre kept where drawn (not centred on the geometry)."
    n1, n2 = raw["cells"]
    R, a1, a2 = _rotation_to_a1(model)
    nodes = (R @ np.asarray(model.nodes, float).T).T
    elems = list(model.elements)
    prim_area = abs(float(np.linalg.det(np.column_stack((a1, a2)))))
    prim_len = float(sum(np.linalg.norm(nodes[e.node_j] - nodes[e.node_i]) for e in elems))
    rect = layout["rectangular_cell"]
    C = np.asarray(layout.get("centre_candidates") or [[x, y, 0] for x, y, _k in layout["cell_centres"]], float)
    if not len(C):
        return raw, "No geometry centre in view; centre kept where drawn."
    d = np.linalg.norm(C[:, :2] - np.asarray(raw["centre"])[None, :], axis=1)
    order = np.argsort(d + 1e-9 * C[:, 2])[:max(1, int(search))]
    names = {0: "cell centre", 1: "junction", 2: "wall mid-point"}
    first_cell = None
    for k in order:
        pt = C[k, :2]
        try:
            b = _cut_block(model, nodes, elems, a1, a2, rect, pt, n1, n2, prim_len, prim_area)
        except ValueError:
            continue
        if C[k, 2] == 0 and first_cell is None:
            first_cell = pt
        if b["mirror"]["about_x1_mid"] and b["mirror"]["about_x2_mid"] and b["min_fragment_ratio"] >= SLIVER_RATIO:
            return (dict(raw, centre=[float(pt[0]), float(pt[1])]),
                    f"Centred on the geometry: nearest {names[int(C[k, 2])]} giving a mirror-symmetric RVE.")
    if first_cell is not None:
        return (dict(raw, centre=[float(first_cell[0]), float(first_cell[1])]),
                "Centred on the nearest cell centre (no mirror-symmetric placement of this size nearby).")
    return raw, "Could not centre the RVE on the geometry; centre kept where drawn."


def spec_rectangle(layout: dict, spec) -> tuple[float, float, float, float]:
    """(x0, y0, width, height) of an FE RVE spec in the FE frame of ``layout``."""
    s = normalise_fe_rve_spec(spec)
    if s is None:
        a = layout["auto"]
        return (*a["window_origin_local"], a["p"], a["q"])
    n1, n2 = s["cells"]
    L, W = n1 * layout["rect_width"], n2 * layout["rect_height"]
    if s["mode"] == "window":
        c = s["centre"]
    else:
        # select_fe_rve centres an explicit block on the base image of the cell type
        bases = layout["cell_base_centroids"]
        c = bases[min(int(s.get("centre_cell", 0)), len(bases) - 1)]
    return (c[0] - 0.5 * L, c[1] - 0.5 * W, L, W)


# ----------------------------------------------------------------------------------
# Public builder
# ----------------------------------------------------------------------------------

def _holds(poly: np.ndarray, centre: np.ndarray, L: float, W: float, margin: float) -> bool:
    rel = poly - (centre - 0.5 * np.array([L, W]))
    return bool(rel[:, 0].min() > margin and rel[:, 0].max() < L - margin and
                rel[:, 1].min() > margin and rel[:, 1].max() < W - margin)


def _mirror_flags(X: np.ndarray, elements: list[Element], L: float, W: float, tol: float) -> dict:
    """Is the clipped network symmetric about x1 = p/2 and about x2 = q/2?"""
    def key(a, b):
        pa = (round(float(a[0]) / tol), round(float(a[1]) / tol)); pb = (round(float(b[0]) / tol), round(float(b[1]) / tol))
        return tuple(sorted((pa, pb)))
    segs = {key(X[e.node_i], X[e.node_j]) for e in elements}
    def mirrored(fx):
        return all(key(fx(X[e.node_i]), fx(X[e.node_j])) in segs for e in elements)
    return {"about_x1_mid": mirrored(lambda v: np.array([L - v[0], v[1]])),
            "about_x2_mid": mirrored(lambda v: np.array([v[0], W - v[1]]))}


def _cut_block(model: RVEModel, nodes: np.ndarray, elems: list[Element], a1: np.ndarray, a2: np.ndarray,
               rect: dict, centre: np.ndarray, n1: int, n2: int, prim_len: float, prim_area: float) -> dict:
    """Clip the tiled network to the n1 x n2 block centred on ``centre``; check its wall measure."""
    scale = max(rect["width"], rect["height"], 1.0)
    tol = max(1e-9 * scale, 1e-9)
    L, W = n1 * rect["width"], n2 * rect["height"]
    origin = centre - 0.5 * np.array([L, W])
    segs = _tile_segments(nodes, elems, a1, a2, centre, 0.75 * (L + W))
    groups_in = list(model.wall_groups) if model.wall_groups else []
    node_lookup: dict[tuple[int, int], int] = {}
    new_nodes: list[np.ndarray] = []
    new_elements: list[Element] = []
    groups: list[str] = []
    seen: set[tuple[int, int]] = set()
    face_of: list[str | None] = []
    source_wall: list[int] = []
    grid = 1e-7 * scale

    def node_id(pt: np.ndarray) -> int:
        k = (int(round(float(pt[0]) / grid)), int(round(float(pt[1]) / grid)))
        if k not in node_lookup:
            node_lookup[k] = len(new_nodes); new_nodes.append(np.asarray(pt, float))
        return node_lookup[k]

    for p0, p1, wi, lam in segs:
        c = _clip(p0 - origin, p1 - origin, L, W, tol)
        if c is None:
            continue
        q0, q1 = c
        q0 = np.array([_snap(q0[0], (0.0, L), 1e-7 * scale), _snap(q0[1], (0.0, W), 1e-7 * scale)])
        q1 = np.array([_snap(q1[0], (0.0, L), 1e-7 * scale), _snap(q1[1], (0.0, W), 1e-7 * scale)])
        a, b = node_id(q0), node_id(q1)
        if a == b:
            continue
        k = (min(a, b), max(a, b))
        if k in seen:
            continue
        seen.add(k)
        face = None
        for name, ax, val in (("X1_0", 0, 0.0), ("X1_P", 0, L), ("X2_0", 1, 0.0), ("X2_Q", 1, W)):
            if q0[ax] == val and q1[ax] == val:
                face = name
        new_elements.append(Element(a, b, lam, 1.0))
        face_of.append(face)
        source_wall.append(int(wi) + 1)
        base = groups_in[wi] if wi < len(groups_in) else f"W{wi + 1}"
        groups.append(base + (f" | face {face}" if face else ""))
    if not new_elements:
        raise ValueError("No walls fall inside the FE RVE.")
    X = np.asarray(new_nodes, float)
    lengths = np.array([np.linalg.norm(X[e.node_j] - X[e.node_i]) for e in new_elements])
    factors = np.array([0.5 if f else 1.0 for f in face_of])
    count = n1 * n2 * rect["primitive_cells"]
    expected = count * prim_len
    effective = float(np.sum(lengths * factors))
    err = abs(effective - expected) / max(expected, np.finfo(float).eps)
    if err > 1e-7:
        raise ValueError(f"FE RVE material check failed: shared-wall-weighted wall length {effective:.6g} mm "
                         f"vs {expected:.6g} mm for {count} cells (relative error {err:.2e}).")
    area = L * W
    if abs(area - count * prim_area) > 1e-9 * area:
        raise ValueError("FE RVE area is not a whole number of primitive cells.")
    prim_walls = [float(np.linalg.norm(nodes[e.node_j] - nodes[e.node_i])) for e in elems]
    med = float(np.median(prim_walls)) if prim_walls else 1.0
    return {"X": X, "elements": new_elements, "face_of": face_of, "groups": groups, "lengths": lengths,
            "factors": factors, "count": count, "expected": expected, "effective": effective, "err": err,
            "L": L, "W": W, "origin": origin, "source_wall": source_wall,
            "min_fragment": float(lengths.min()), "min_fragment_ratio": float(lengths.min() / max(med, 1e-12)),
            "mirror": _mirror_flags(X, new_elements, L, W, 1e-6 * scale)}


def select_fe_rve(model: RVEModel, *, cells: tuple[int, int] | None = None, centre_cell: int | None = None,
                  window: tuple[float, float] | None = None, face_wall_section: str = "half-thickness",
                  use_model_spec: bool = True) -> tuple[RVEModel, dict]:
    """FE RVE of any periodic wall network: automatic (source-study RVE rule, see module docstring)
    or user-defined.

    ``cells`` forces the block size (n1 x n2 rectangular cells) and ``centre_cell``
    the cell type to centre on (index into :func:`find_cells`, 0 = largest).
    ``window`` (with ``cells``) places the block with its centre at that point of the
    FE frame -- a rectangle drawn over the replicated design.  Without any of these the
    design's own choice (``model.fe_rve``, see :func:`normalise_fe_rve_spec`) is used,
    unless ``use_model_spec`` is False; with none, the automatic rule.
    Returns the clipped wall network in the x1-x2 frame of the block (origin at a
    corner, 0 <= x1 <= p, 0 <= x2 <= q) and a diagnostics dictionary.  Walls that
    lie in a face are listed in ``face_walls`` (1-based) and carry the section
    factor in ``wall_section_factors``.
    """
    if face_wall_section not in FACE_WALL_SECTIONS:
        raise ValueError(f"face_wall_section must be one of {FACE_WALL_SECTIONS}")
    if not model.elements:
        raise ValueError("The RVE contains no walls.")
    for wi, e in enumerate(model.elements, 1):
        if abs(float(e.weight) - 1.0) > 1e-12:
            raise ValueError(f"Wall W{wi} has weight {e.weight:g}. The FE RVE needs explicit physical walls with weight 1.")
    explicit = cells is not None or centre_cell is not None or window is not None
    spec = normalise_fe_rve_spec(getattr(model, "fe_rve", None)) if (use_model_spec and not explicit) else None
    if spec is not None:
        cells = tuple(spec["cells"])
        if spec["mode"] == "cells":
            centre_cell = spec.get("centre_cell")
        else:
            window = tuple(spec["centre"])
    if cells is not None and (len(cells) != 2 or min(int(c) for c in cells) < 1):
        raise ValueError("RVE cell counts must be at least 1 x 1.")
    if window is not None and cells is None:
        raise ValueError("A drawn FE RVE needs its size in rectangular cells (cells=(n1, n2)).")
    mode = "window" if window is not None else "cells" if (cells is not None or centre_cell is not None) else "auto"

    R, a1, a2 = _rotation_to_a1(model)
    nodes = (R @ np.asarray(model.nodes, float).T).T
    elems = list(model.elements)
    prim_area = abs(float(np.linalg.det(np.column_stack((a1, a2)))))
    prim_len = float(sum(np.linalg.norm(nodes[e.node_j] - nodes[e.node_i]) for e in elems))
    rect = rectangular_cell(a1, a2)
    cell_list = find_cells(nodes, elems, a1, a2)
    scale = max(rect["width"], rect["height"], 1.0)
    margin = 1e-6 * scale
    if centre_cell is not None:
        if not 0 <= int(centre_cell) < len(cell_list):
            raise ValueError(f"centre_cell must be 0..{len(cell_list) - 1} ({len(cell_list)} cell types found).")
        order = [int(centre_cell)]
    else:
        order = list(range(len(cell_list)))
    blocks = [(int(cells[0]), int(cells[1]))] if cells is not None else list(CELL_CANDIDATES)

    tried: list[dict] = []
    chosen = None
    fallback = None
    warnings: list[str] = []
    if window is not None:
        n1, n2 = int(cells[0]), int(cells[1])
        centre = np.asarray(window, float).reshape(2)
        built = _cut_block(model, nodes, elems, a1, a2, rect, centre, n1, n2, prim_len, prim_area)
        A = np.column_stack((a1, a2))
        hit = _cell_at(cell_list, A, centre)
        ci = hit[0] if hit else 0
        poly = hit[1] if hit else None
        L, W = n1 * rect["width"], n2 * rect["height"]
        inside = poly is not None and _holds(poly, centre, L, W, -margin)
        cell_rec = dict(cell_list[ci])
        cell_rec["centroid"] = _polygon_centroid(poly) if poly is not None else centre
        cell_rec["polygon"] = poly if inside else np.zeros((0, 2))
        cell_list_view = list(cell_list)
        cell_list_view[ci] = cell_rec
        chosen = (ci, n1, n2, built)
        tried.append({"cell": ci, "cells_x1_x2": [n1, n2], "result": "user-drawn rectangle"})
        if built["min_fragment_ratio"] < SLIVER_RATIO:
            warnings.append(f"The drawn rectangle leaves a wall piece of {built['min_fragment']:.3g} mm "
                            f"({100 * built['min_fragment_ratio']:.1f} % of the median wall); move it slightly "
                            f"or snap it to a cell centre to avoid a sliver for the mesher.")
        if not inside:
            warnings.append("No complete cell lies at the centre of the drawn rectangle.")
        if not (built["mirror"]["about_x1_mid"] and built["mirror"]["about_x2_mid"]):
            _, auto_d = select_fe_rve(model, use_model_spec=False, face_wall_section=face_wall_section)
            if auto_d["mirror_symmetry"]["about_x1_mid"] and auto_d["mirror_symmetry"]["about_x2_mid"]:
                warnings.append("The drawn rectangle is not mirror-symmetric although the automatic RVE is: its "
                                "faces cut the walls at unmatched points, so the DBC (finite-height) terms can differ "
                                "strongly from the periodic values. PBC terms are unaffected. Snap the centre to a "
                                "cell for a DBC comparison.")
        order = []
    for ci in order:
        c = cell_list[ci]
        for n1, n2 in blocks:
            L, W = n1 * rect["width"], n2 * rect["height"]
            if cells is None and not _holds(c["polygon"], c["centroid"], L, W, margin):
                tried.append({"cell": ci, "cells_x1_x2": [n1, n2], "result": "cell not inside"})
                continue
            try:
                built = _cut_block(model, nodes, elems, a1, a2, rect, c["centroid"], n1, n2, prim_len, prim_area)
            except ValueError as exc:
                tried.append({"cell": ci, "cells_x1_x2": [n1, n2], "result": f"rejected: {exc}"})
                continue
            if cells is not None or built["min_fragment_ratio"] >= SLIVER_RATIO:
                tried.append({"cell": ci, "cells_x1_x2": [n1, n2], "result": "selected"})
                chosen = (ci, n1, n2, built)
                break
            tried.append({"cell": ci, "cells_x1_x2": [n1, n2],
                          "result": f"sliver: shortest wall piece {built['min_fragment']:.3g}"})
            if fallback is None or built["min_fragment_ratio"] > fallback[3]["min_fragment_ratio"]:
                fallback = (ci, n1, n2, built)
        if chosen is not None:
            break
    if chosen is None and fallback is not None:
        chosen = fallback
        warnings.append(f"Every candidate block leaves a short wall piece; the best one "
                        f"({chosen[3]['min_fragment']:.3g} mm) was used. Refine the mesh near it.")
    if chosen is None:
        raise ValueError("No block of up to 3 x 3 rectangular cells holds a complete cell of this wall network; "
                         "give the cell counts explicitly (cells=(n1, n2)).")
    ci, n1, n2, b = chosen
    cell = cell_list_view[ci] if window is not None else cell_list[ci]
    face_of = b["face_of"]
    section_factors = np.ones_like(b["factors"]) if face_wall_section == "full" else b["factors"]
    L, W, origin, X = b["L"], b["W"], b["origin"], b["X"]
    dbc_model = RVEModel(
        X, b["elements"], np.array([L, 0.0]), np.array([0.0, W]), float(model.height),
        periodic_tol=max(float(model.periodic_tol), 1e-9), anchor_node=0, fix_anchor_rotation=False,
        wall_element_type=model.wall_element_type, strip_bending_mode=model.strip_bending_mode,
        length_unit=model.length_unit, gravity=model.gravity,
        name=f"{model.name} - FE RVE", wall_groups=b["groups"],
    )
    face_idx = [i + 1 for i, f in enumerate(face_of) if f]
    diag = {
        "mode": "thesis_table_3_4_DBC",
        "rve": {"auto": "source-study finite-height RVE (thesis Fig. 3.6b, Table 3.3; paper Fig. 8b), applied automatically",
                "cells": f"user-defined block of {n1} x {n2} rectangular periodic cells centred on a complete cell",
                "window": f"user-drawn rectangle of {n1} x {n2} rectangular periodic cells"}[mode],
        "rve_rule": RVE_RULE if mode == "auto" else
                    ("block size chosen by the user; centred on a complete cell" if mode == "cells" else
                     "drawn over the replicated design; size rounded to whole rectangular periodic cells, "
                     "centre where drawn"),
        "rve_mode": mode,
        "rve_spec": (None if mode == "auto" else
                     {"mode": mode, "cells": [n1, n2],
                      **({"centre": [float(x) for x in window]} if mode == "window" else {}),
                      **({"centre_cell": int(centre_cell)} if (mode == "cells" and centre_cell is not None) else {})}),
        "top_bottom_periodic": False,
        "finite_height_core_validation": True,
        "p": float(L), "q": float(W), "height": float(model.height),
        "rve_area": float(L * W), "rve_volume": float(L * W * model.height),
        "cells_x1_x2": [n1, n2],
        "rectangular_cell": rect,
        "primitive_cells": int(b["count"]),
        "primitive_area": float(prim_area),
        "primitive_wall_length": float(prim_len),
        "expected_wall_length": float(b["expected"]),
        "actual_wall_length": float(b["lengths"].sum()),
        "effective_wall_length": float(b["effective"]),
        "wall_measure_relative_error": float(b["err"]),
        "face_walls": face_idx,
        "face_wall_faces": {str(i + 1): f for i, f in enumerate(face_of) if f},
        "face_wall_length": float(b["lengths"][b["factors"] < 1].sum()),
        "face_wall_section": face_wall_section,
        "wall_section_factors": [float(x) for x in section_factors],
        "source_wall": list(b["source_wall"]),
        "wall_segments": [[float(X[e.node_i][0]), float(X[e.node_i][1]), float(X[e.node_j][0]), float(X[e.node_j][1])]
                          for e in b["elements"]],
        "cell_centroid_rve": [float(x) for x in (cell["centroid"] - origin)],
        "central_cell_polygon_rve": [[float(x) for x in (pt - origin)] for pt in cell["polygon"]],
        "central_cell": {"index": int(ci), "name": cell["name"], "corners": int(cell["corners"]),
                         "area": float(cell["area"]), "inradius": float(cell["inradius"])},
        "cell_types": [{"index": i, "name": c["name"], "corners": int(c["corners"]), "area": float(c["area"])}
                       for i, c in enumerate(cell_list)],
        "min_wall_piece": float(b["min_fragment"]),
        "min_wall_piece_ratio": float(b["min_fragment_ratio"]),
        "mirror_symmetry": b["mirror"],
        "selection": {"automatic": mode == "auto", "mode": mode, "user_defined": mode != "auto",
                      "tried": tried, "warnings": warnings},
        "window_origin_local": [float(x) for x in origin],
        "rotation_to_local": R.tolist(),
        "generated_nodes": int(len(X)), "generated_walls": int(len(b["elements"])),
        "coordinate_frame": "x1 along a1; x2 normal to a1 in the core plane; x3 = core height; origin at an RVE corner",
    }
    return dbc_model, diag


def build_published_fe_rve(model: RVEModel, *, cells: tuple[int, int] | None = None,
                           face_wall_section: str = "half-thickness",
                           centre_cell: int | None = None) -> tuple[RVEModel, dict]:
    """Backward-compatible name of :func:`select_fe_rve`."""
    return select_fe_rve(model, cells=cells, centre_cell=centre_cell, face_wall_section=face_wall_section)


def fe_rve_summary(diag: dict) -> str:
    """One-paragraph, plain-language description of a selected FE RVE."""
    p, q, h = diag.get("p", 0.0), diag.get("q", 0.0), diag.get("height", 0.0)
    n1, n2 = diag.get("cells_x1_x2", [1, 1])
    cell = diag.get("central_cell", {})
    sym = diag.get("mirror_symmetry", {})
    sym_txt = ("mirror-symmetric about both mid-planes" if sym.get("about_x1_mid") and sym.get("about_x2_mid") else
               "mirror-symmetric about x1 = p/2" if sym.get("about_x1_mid") else
               "mirror-symmetric about x2 = q/2" if sym.get("about_x2_mid") else "not mirror-symmetric")
    nf = len(diag.get("face_walls", []))
    face_txt = (f"{nf} wall piece{'s' if nf != 1 else ''} lie in a face and carry "
                f"{'the full laminate' if diag.get('face_wall_section') == 'full' else 'half the laminate'}"
                if nf else "no wall lies in a face")
    mode = diag.get("rve_mode", "auto")
    where = (f"centred on a complete {cell.get('name', 'cell')}" if mode != "window" else
             (f"drawn by the user, a complete {cell.get('name', 'cell')} at its centre"
              if diag.get("central_cell_polygon_rve") else "drawn by the user"))
    lead = {"auto": "", "cells": "User-defined FE RVE: ", "window": "User-drawn FE RVE: "}.get(mode, "")
    txt = (f"{lead}{p:.1f} x {q:.1f} x {h:.1f} mm block of {n1} x {n2} rectangular cells "
           f"({diag.get('primitive_cells', '?')} primitive cells), {where}; "
           f"{diag.get('generated_walls', '?')} wall pieces, {face_txt}; {sym_txt}; "
           f"shortest wall piece {diag.get('min_wall_piece', float('nan')):.1f} mm.")
    if mode != "auto":
        txt += (" The source study used this RVE; the DBC (finite-height) terms depend on "
                "the RVE, the PBC terms should not.")
    for w in diag.get("selection", {}).get("warnings", []):
        txt += " Warning: " + w
    return txt
