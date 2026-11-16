from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Iterable, Sequence

import numpy as np

from .models import Element


@dataclass(slots=True)
class GeometryIssue:
    severity: str  # "error" or "warning"
    message: str


@dataclass(slots=True)
class WallGeometry:
    wall_id: int
    node_i: int
    node_j: int
    dx: float
    dy: float
    length: float
    angle_deg: float
    midpoint: tuple[float, float]
    weight: float = 1.0


@dataclass(slots=True)
class GeometrySummary:
    n_nodes: int
    n_walls: int
    bounds: tuple[float, float, float, float]
    width: float
    height_2d: float
    centroid: tuple[float, float]
    cell_area: float
    lattice_length_1: float
    lattice_length_2: float
    lattice_angle_deg: float
    p: float
    q: float
    r_signed: float
    r: float
    total_wall_length: float
    weighted_wall_length: float
    node_degree: list[int]
    walls: list[WallGeometry] = field(default_factory=list)
    periodic_pairs: np.ndarray | None = None
    detected_geometry_id: str = ""
    detected_geometry_name: str = ""
    detected_parameters: dict[str, float] = field(default_factory=dict)
    detection_confidence: float = 0.0
    parameter_sources: dict[str, list[int]] = field(default_factory=dict)
    issues: list[GeometryIssue] = field(default_factory=list)

    @property
    def is_valid(self) -> bool:
        return not any(i.severity == "error" for i in self.issues)


def _as_elements(elements: Iterable[Element | Sequence[float]]) -> list[Element]:
    out: list[Element] = []
    for e in elements:
        if isinstance(e, Element):
            out.append(e)
        else:
            row = list(e)
            if len(row) < 2:
                raise ValueError("Each wall must provide node_i and node_j")
            weight = float(row[2]) if len(row) >= 3 else 1.0
            out.append(Element(int(row[0]), int(row[1]), 0, weight))
    return out


def characteristic_pqr(a1: Sequence[float], a2: Sequence[float]) -> tuple[float, float, float, float]:
    """Return original-model characteristic p, q, signed r and |r|.

    The lattice is expressed in a local frame whose first axis is aligned with a1:
    a1 = [p, 0], a2 = [r_signed, q_signed].  The reported q is |q_signed| and
    r is |r_signed|.  This reduces to the benchmark honeycomb p/q/r convention for the
    standard honeycomb RVEs while remaining meaningful for a rotated sketch.
    """
    u = np.asarray(a1, dtype=float).reshape(2)
    v = np.asarray(a2, dtype=float).reshape(2)
    p = float(np.linalg.norm(u))
    if not math.isfinite(p) or p <= np.finfo(float).eps:
        return math.nan, math.nan, math.nan, math.nan
    e1 = u / p
    e2 = np.array([-e1[1], e1[0]])
    r_signed = float(np.dot(v, e1))
    q_signed = float(np.dot(v, e2))
    return p, abs(q_signed), r_signed, abs(r_signed)


def detect_periodic_pairs(nodes: np.ndarray, a1: Sequence[float], a2: Sequence[float], tol: float = 1e-8) -> np.ndarray:
    """Preview periodically equivalent node groups without altering the solver.

    Returns rows [master, slave, n1, n2] using zero-based node indices and the
    same fractional-coordinate criterion as the periodic solver.
    """
    pts = np.asarray(nodes, dtype=float)
    cell = np.column_stack((np.asarray(a1, float), np.asarray(a2, float)))
    if pts.ndim != 2 or pts.shape[1] != 2 or not np.all(np.isfinite(pts)):
        return np.empty((0, 4), float)
    if abs(float(np.linalg.det(cell))) <= np.finfo(float).eps:
        return np.empty((0, 4), float)
    xi = np.linalg.solve(cell, pts.T).T
    assigned = np.zeros(len(pts), bool)
    pairs: list[list[float]] = []
    for i in range(len(pts)):
        if assigned[i]:
            continue
        group = [i]
        assigned[i] = True
        for j in range(i + 1, len(pts)):
            if assigned[j]:
                continue
            dd = xi[j] - xi[i]
            if np.linalg.norm(dd - np.rint(dd), ord=np.inf) <= tol:
                group.append(j)
                assigned[j] = True
        if len(group) > 1:
            master = group[0]
            for slave in group[1:]:
                shift = np.rint(xi[slave] - xi[master]).astype(int)
                pairs.append([master, slave, int(shift[0]), int(shift[1])])
    return np.asarray(pairs, float).reshape((-1, 4)) if pairs else np.empty((0, 4), float)



def _detect_builtin_geometry(pts: np.ndarray, walls: list[WallGeometry], degrees: list[int], a1: np.ndarray) -> tuple[str, str, dict[str, float], float, dict[str, list[int]]]:
    """Recognise built-in topology families and identify parameter evidence.

    Recognition is intentionally conservative.  Topology must match exactly; wall
    directions are then resolved relative to ``a1`` so rotation of the sketch does
    not change the family.  The confidence value reflects how consistently the
    walls assigned to each parameter share a common length.
    """
    n, nw = len(pts), len(walls)
    signature = (n, nw, tuple(sorted(degrees)))
    known = {
        (5, 3, (0, 1, 1, 1, 3)): ("HEX", "Hexagonal"),
        (7, 4, (0, 0, 1, 1, 1, 1, 4)): ("SQUARE", "Square-family"),
        (7, 6, (1, 1, 1, 1, 1, 1, 6)): ("MIXED", "Mixed rhombus-triangular"),
        (9, 6, (0, 0, 1, 1, 1, 1, 1, 1, 6)): ("TRI", "Triangular-family"),
    }
    if signature not in known or not walls:
        return "", "Custom / arbitrary periodic RVE", {}, 0.0, {}
    gid, name = known[signature]
    u = np.asarray(a1, float).reshape(2)
    Lu = float(np.linalg.norm(u))
    if not math.isfinite(Lu) or Lu <= np.finfo(float).eps:
        return gid, name, {}, 0.55, {}
    e1 = u / Lu
    parallel: list[WallGeometry] = []
    perpendicular: list[WallGeometry] = []
    other: list[WallGeometry] = []
    for w in walls:
        d = np.array([w.dx, w.dy], float) / w.length
        c = abs(float(np.dot(d, e1)))
        if c >= math.cos(math.radians(12.0)):
            parallel.append(w)
        elif c <= math.sin(math.radians(12.0)):
            perpendicular.append(w)
        else:
            other.append(w)

    def mean_group(group):
        return float(np.mean([w.length for w in group])) if group else math.nan

    def source_ids(group):
        return [int(w.wall_id) for w in group]

    def consistency(group):
        if not group:
            return 0.0
        vals=np.asarray([w.length for w in group],float)
        m=float(np.mean(vals))
        if len(vals)==1 or m <= np.finfo(float).eps:
            return 1.0
        # 0% spread -> 1.0; 10% coefficient of variation -> 0.0.
        cv=float(np.std(vals)/m)
        return max(0.0, min(1.0, 1.0-cv/0.10))

    if gid == "HEX":
        g1=perpendicular
        g2=other if other else [w for w in walls if w not in g1]
        params={"L1":mean_group(g1), "L2":mean_group(g2)}
        sources={"L1":source_ids(g1), "L2":source_ids(g2)}
        quality=(consistency(g1)+consistency(g2))/2
    elif gid == "SQUARE":
        params={"Lx":mean_group(parallel), "Ly":mean_group(perpendicular)}
        sources={"Lx":source_ids(parallel), "Ly":source_ids(perpendicular)}
        quality=(consistency(parallel)+consistency(perpendicular))/2
    else:
        gh=parallel
        gd=other if other else perpendicular
        Lh=mean_group(gh); Ld=mean_group(gd)
        params={"Lh":Lh, "R":(Ld/Lh if math.isfinite(Lh) and Lh>0 and math.isfinite(Ld) else math.nan)}
        sources={"Lh":source_ids(gh), "R":source_ids(gd)}
        quality=(consistency(gh)+consistency(gd))/2
    confidence=max(0.55,min(1.0,0.60+0.40*quality))
    return gid, name, params, confidence, sources

def analyse_geometry(
    nodes: Sequence[Sequence[float]],
    elements: Iterable[Element | Sequence[float]],
    a1: Sequence[float],
    a2: Sequence[float],
    *,
    periodic_tol: float = 1e-8,
    zero_tol: float = 1e-10,
    duplicate_tol: float = 1e-8,
) -> GeometrySummary:
    pts = np.asarray(nodes, dtype=float)
    if pts.size == 0:
        pts = np.empty((0, 2), dtype=float)
    elems = _as_elements(elements)
    issues: list[GeometryIssue] = []

    if pts.ndim != 2 or pts.shape[1] != 2:
        raise ValueError("nodes must be an N x 2 array")
    n = len(pts)
    if n == 0:
        bounds = (math.nan, math.nan, math.nan, math.nan)
        centroid = (math.nan, math.nan)
        width = height_2d = 0.0
        issues.append(GeometryIssue("error", "No nodes have been defined."))
    else:
        finite = np.all(np.isfinite(pts), axis=1)
        if not np.all(finite):
            issues.append(GeometryIssue("error", "All node coordinates must be finite."))
        safe = pts[finite] if np.any(finite) else np.zeros((1, 2))
        xmin, ymin = np.min(safe, axis=0)
        xmax, ymax = np.max(safe, axis=0)
        bounds = (float(xmin), float(ymin), float(xmax), float(ymax))
        width, height_2d = float(xmax - xmin), float(ymax - ymin)
        c = np.mean(safe, axis=0)
        centroid = (float(c[0]), float(c[1]))

    # Near-duplicate nodes make drawing/periodicity ambiguous.
    for i in range(n):
        for j in range(i + 1, n):
            if np.all(np.isfinite(pts[[i, j]])) and np.linalg.norm(pts[j] - pts[i]) <= duplicate_tol:
                issues.append(GeometryIssue("warning", f"Nodes {i+1} and {j+1} are coincident or nearly coincident."))

    degrees = [0] * n
    walls: list[WallGeometry] = []
    seen_edges: set[tuple[int, int]] = set()
    adjacency: dict[int, set[int]] = {i: set() for i in range(n)}
    total = weighted = 0.0
    for k, e in enumerate(elems, 1):
        ni, nj = int(e.node_i), int(e.node_j)
        if ni < 0 or nj < 0 or ni >= n or nj >= n:
            issues.append(GeometryIssue("error", f"Wall {k} references a node outside 1..{n}."))
            continue
        if ni == nj:
            issues.append(GeometryIssue("error", f"Wall {k} starts and ends at node {ni+1}."))
            continue
        key = tuple(sorted((ni, nj)))
        if key in seen_edges:
            issues.append(GeometryIssue("warning", f"Wall {k} duplicates an existing connection between nodes {key[0]+1} and {key[1]+1}."))
        seen_edges.add(key)
        d = pts[nj] - pts[ni]
        L = float(np.linalg.norm(d))
        if not math.isfinite(L) or L <= zero_tol:
            issues.append(GeometryIssue("error", f"Wall {k} has zero or invalid length."))
            continue
        ang = math.degrees(math.atan2(float(d[1]), float(d[0])))
        midpoint = tuple((0.5 * (pts[ni] + pts[nj])).tolist())
        walls.append(WallGeometry(k, ni, nj, float(d[0]), float(d[1]), L, ang, midpoint, float(e.weight)))
        total += L
        weighted += float(e.weight) * L
        degrees[ni] += 1
        degrees[nj] += 1
        adjacency[ni].add(nj)
        adjacency[nj].add(ni)

    if not elems:
        issues.append(GeometryIssue("error", "No walls have been defined."))
    if n:
        active = [i for i, d in enumerate(degrees) if d > 0]
        if active:
            visited = {active[0]}
            stack = [active[0]]
            while stack:
                cur = stack.pop()
                for nxt in adjacency[cur]:
                    if nxt not in visited:
                        visited.add(nxt)
                        stack.append(nxt)
            if any(i not in visited for i in active):
                issues.append(GeometryIssue("error", "The wall network contains more than one disconnected component."))

    u = np.asarray(a1, dtype=float).reshape(2)
    v = np.asarray(a2, dtype=float).reshape(2)
    det = float(np.linalg.det(np.column_stack((u, v))))
    cell_area = abs(det)
    L1 = float(np.linalg.norm(u))
    L2 = float(np.linalg.norm(v))
    if not np.all(np.isfinite(np.r_[u, v])):
        issues.append(GeometryIssue("error", "Periodic lattice vectors must be finite."))
    if L1 <= zero_tol or L2 <= zero_tol or cell_area <= zero_tol:
        issues.append(GeometryIssue("error", "Periodic lattice vectors must be non-zero and linearly independent."))
        angle = math.nan
    else:
        cosang = float(np.dot(u, v) / (L1 * L2))
        angle = math.degrees(math.acos(max(-1.0, min(1.0, cosang))))

    p, q, r_signed, r = characteristic_pqr(u, v)
    pairs = detect_periodic_pairs(pts, u, v, periodic_tol) if n else np.empty((0, 4), float)
    periodic_nodes = set()
    for row in pairs:
        periodic_nodes.add(int(row[0])); periodic_nodes.add(int(row[1]))
    isolated = [i + 1 for i, d in enumerate(degrees) if d == 0 and i not in periodic_nodes]
    if isolated:
        issues.append(GeometryIssue("warning", "Disconnected non-periodic node(s): " + ", ".join(map(str, isolated))))
    if n and cell_area > zero_tol and len(pairs) == 0:
        issues.append(GeometryIssue("warning", "No periodic duplicate-node pairs were detected. This is acceptable only if explicit pairs will be supplied."))

    detected_id, detected_name, detected_parameters, detection_confidence, parameter_sources = _detect_builtin_geometry(pts, walls, degrees, u)

    return GeometrySummary(
        n_nodes=n,
        n_walls=len(elems),
        bounds=bounds,
        width=width,
        height_2d=height_2d,
        centroid=centroid,
        cell_area=cell_area,
        lattice_length_1=L1,
        lattice_length_2=L2,
        lattice_angle_deg=angle,
        p=p,
        q=q,
        r_signed=r_signed,
        r=r,
        total_wall_length=total,
        weighted_wall_length=weighted,
        node_degree=degrees,
        walls=walls,
        periodic_pairs=pairs,
        detected_geometry_id=detected_id,
        detected_geometry_name=detected_name,
        detected_parameters=detected_parameters,
        detection_confidence=detection_confidence,
        parameter_sources=parameter_sources,
        issues=issues,
    )


def summary_text(summary: GeometrySummary, core_height: float | None = None) -> str:
    def f(v: float) -> str:
        return "—" if not math.isfinite(v) else f"{v:.6g}"

    lines = [
        "GEOMETRY SUMMARY",
        f"Nodes: {summary.n_nodes}",
        f"Walls: {summary.n_walls}",
        f"Drawing width × height: {f(summary.width)} × {f(summary.height_2d)}",
        f"Bounds [xmin,ymin,xmax,ymax]: {[f(v) for v in summary.bounds]}",
        f"Centroid: ({f(summary.centroid[0])}, {f(summary.centroid[1])})",
        f"Node coordination: {summary.node_degree}",
        f"Cell area |a1 × a2|: {f(summary.cell_area)}",
        f"|a1|: {f(summary.lattice_length_1)}",
        f"|a2|: {f(summary.lattice_length_2)}",
        "",
        "Characteristic cell dimensions",
        f"p: {f(summary.p)}",
        f"q: {f(summary.q)}",
        f"r offset: {f(summary.r)}  (signed {f(summary.r_signed)})",
        f"Core height h: {f(float(core_height)) if core_height is not None else '—'}",
        f"RVE volume: {f(summary.cell_area*float(core_height)) if core_height is not None else '—'}",
        "",
        f"Total wall length: {f(summary.total_wall_length)}",
        f"Weighted wall length: {f(summary.weighted_wall_length)}",
        f"Detected periodic pairs: {0 if summary.periodic_pairs is None else len(summary.periodic_pairs)}",
        "",
        "GEOMETRY LIBRARY MATCH",
        f"Family: {summary.detected_geometry_name}",
    ]
    if summary.detected_geometry_id:
        lines.append(f"Library ID: {summary.detected_geometry_id}")
        lines.append(f"Match confidence: {100.0*summary.detection_confidence:.1f}%")
        for key, value in summary.detected_parameters.items():
            src=summary.parameter_sources.get(key,[])
            source_text=("  from " + ", ".join(f"W{i}" for i in src)) if src else ""
            lines.append(f"{key}: {f(value)}{source_text}")
    else:
        lines.append("No built-in family matched. The generic periodic solver can still use this geometry.")
    lines += [
        "",
        "WALL PARAMETERS (wall geometry of the source study: lᵢ, θᵢ)",
    ]
    for w in summary.walls:
        lines.append(f"W{w.wall_id}: N{w.node_i+1} → N{w.node_j+1}   l={f(w.length)}   θ={f(w.angle_deg)}°   weight={f(w.weight)}")
    lines.append("")
    lines.append("VALIDATION")
    if not summary.issues:
        lines.append("✓ Geometry checks passed.")
    else:
        for issue in summary.issues:
            icon = "✕" if issue.severity == "error" else "⚠"
            lines.append(f"{icon} {issue.message}")
    return "\n".join(lines)
