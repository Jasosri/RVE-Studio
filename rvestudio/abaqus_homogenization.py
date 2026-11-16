from __future__ import annotations

from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
import json
import math
import os
import shutil
import hashlib
from typing import Callable

import numpy as np

from .models import RVEModel, Material, Laminate, Element
from .abaqus import (
    _safe_name, _ids_lines, _abaqus_material_blocks, _build_shell_mesh,
    _run_shell_command,
)

# Engineering-strain Voigt order used by RVE Studio.
VOIGT_LABELS = ("E11", "E22", "E33", "G23", "G13", "G12")
STRESS_LABELS = ("S11", "S22", "S33", "T23", "T13", "T12")

PBC_DOMAIN_PRIMITIVE = "primitive"
PBC_DOMAIN_LATTICE_SUPERCELL = "lattice_supercell"
PBC_DOMAIN_RECTANGULAR_SUPERCELL = "rectangular_supercell"
PBC_DOMAIN_SYMMETRIC_RECTANGULAR = "symmetric_rectangular"
PBC_DOMAIN_PUBLISHED = "published_rve"  # the thesis/paper DBC block, used periodically

CH_METHOD_DIRECT = "direct_reaction6"
CH_METHOD_ENERGY = "energy21"
# Five unit strains without E33: U3+ - U3- = 0 between top and bottom in every state, so the
# periodic model gives the 3D stiffness terms at eps33 = 0 (the FE verification default).
CH_METHOD_DIRECT5 = "direct_reaction5"
PBC5_COMPONENTS = (0, 1, 3, 4, 5)


def _strain_tensor_from_engineering(v: np.ndarray) -> np.ndarray:
    """Convert [e11,e22,e33,g23,g13,g12] to the symmetric small-strain tensor."""
    e = np.asarray(v, float).reshape(6)
    return np.array([
        [e[0], e[5] / 2.0, e[4] / 2.0],
        [e[5] / 2.0, e[1], e[3] / 2.0],
        [e[4] / 2.0, e[3] / 2.0, e[2]],
    ], float)


def homogenization_cases(amplitude: float = 1.0) -> list[dict]:
    """Return the 21 strain states needed to identify a symmetric 6x6 stiffness."""
    a = float(amplitude)
    cases: list[dict] = []
    for i, name in enumerate(VOIGT_LABELS):
        v = np.zeros(6); v[i] = a
        cases.append({"index": len(cases) + 1, "name": name, "components": [i], "strain": v})
    for i in range(6):
        for j in range(i + 1, 6):
            v = np.zeros(6); v[i] = a; v[j] = a
            cases.append({"index": len(cases) + 1, "name": f"{VOIGT_LABELS[i]}_{VOIGT_LABELS[j]}", "components": [i, j], "strain": v})
    return cases


@dataclass(slots=True)
class PeriodicHomogenizationSettings:
    element_type: str = "S4"
    mesh_size: float = 5.0
    strain_amplitude: float = 1.0e-3
    cpus: int = 1
    save_odb: bool = True
    stress_cross_check: bool = True
    identification_method: str = CH_METHOD_DIRECT  # six PBC macro-reaction columns; optional energy21
    significance_tolerance: float = 1.0e-4
    reference_length: float = 1.0  # mm; converts dimensionless macro strain to control-node displacement
    domain_mode: str = PBC_DOMAIN_PUBLISHED  # automatic FE RVE (same block as the DBC check)
    supercell_repeat_x: int = 2
    supercell_repeat_y: int = 2
    supercell_search_limit: int = 6
    supercell_orthogonality_tolerance: float = 1.0e-8
    strict_pbc_audit: bool = True

    def validate(self) -> None:
        if self.element_type.upper() != "S4":
            raise ValueError("Periodic homogenisation currently supports S4 only.")
        if not math.isfinite(self.mesh_size) or self.mesh_size <= 0:
            raise ValueError("Mesh size must be positive and finite.")
        if not math.isfinite(self.strain_amplitude) or self.strain_amplitude <= 0:
            raise ValueError("Macro strain amplitude must be positive and finite.")
        if int(self.cpus) < 1:
            raise ValueError("CPU count must be at least 1.")
        if not math.isfinite(self.significance_tolerance) or self.significance_tolerance < 0:
            raise ValueError("Significance tolerance must be non-negative and finite.")
        if not math.isfinite(self.reference_length) or self.reference_length <= 0:
            raise ValueError("Reference length must be positive and finite.")
        if self.identification_method not in {CH_METHOD_DIRECT, CH_METHOD_ENERGY, CH_METHOD_DIRECT5}:
            raise ValueError("Identification method must be direct_reaction5, direct_reaction6 or energy21.")
        if self.domain_mode not in {PBC_DOMAIN_PRIMITIVE, PBC_DOMAIN_LATTICE_SUPERCELL, PBC_DOMAIN_RECTANGULAR_SUPERCELL, PBC_DOMAIN_SYMMETRIC_RECTANGULAR, PBC_DOMAIN_PUBLISHED}:
            raise ValueError(f"Unknown periodic domain mode: {self.domain_mode}")
        if int(self.supercell_repeat_x) < 1 or int(self.supercell_repeat_y) < 1:
            raise ValueError("Supercell repeats must be at least 1 in each direction.")
        if int(self.supercell_search_limit) < 1:
            raise ValueError("Supercell integer search limit must be at least 1.")
        if not math.isfinite(self.supercell_orthogonality_tolerance) or self.supercell_orthogonality_tolerance <= 0:
            raise ValueError("Supercell orthogonality tolerance must be positive and finite.")


@dataclass(slots=True)
class HomogenizationCaseFiles:
    index: int
    name: str
    strain: list[float]
    case_dir: Path
    job_name: str
    input_file: Path
    extractor_file: Path
    results_file: Path
    log_file: Path


@dataclass(slots=True)
class PeriodicHomogenizationFiles:
    run_dir: Path
    snapshot_file: Path
    summary_file: Path
    log_file: Path
    cases: list[HomogenizationCaseFiles]


def _validate_periodic_fe_model(model: RVEModel) -> None:
    if not model.elements:
        raise ValueError("The RVE contains no walls.")
    used = set()
    for wi, e in enumerate(model.elements, 1):
        if e.node_i == e.node_j:
            raise ValueError(f"Wall W{wi} has identical end nodes.")
        if not (0 <= e.node_i < len(model.nodes) and 0 <= e.node_j < len(model.nodes)):
            raise ValueError(f"Wall W{wi} contains an invalid node index.")
        if abs(float(e.weight) - 1.0) > 1e-12:
            raise ValueError(
                f"Wall W{wi} has multiplicity/weight {e.weight:g}. Periodic Abaqus FE requires the explicit physical wall geometry; "
                "weighted analytical walls must be expanded into explicit walls first."
            )
        used.update((e.node_i, e.node_j))
    if not used:
        raise ValueError("No connected RVE nodes were found.")
    area = abs(float(np.linalg.det(np.column_stack((model.a1, model.a2)))))
    if area <= np.finfo(float).eps or model.height <= 0:
        raise ValueError("The RVE must have a positive lattice area and core height.")



def find_rectangular_supercell_basis(
    model: RVEModel,
    *,
    search_limit: int = 6,
    orthogonality_tolerance: float = 1.0e-8,
) -> dict:
    """Find the smallest commensurate orthogonal in-plane lattice basis.

    The returned vectors are integer combinations of the primitive ``a1`` and
    ``a2``.  Therefore the rectangular domain is a true periodic supercell and
    its area must be an integer multiple of the primitive-cell area.
    """
    A = np.column_stack((np.asarray(model.a1, float), np.asarray(model.a2, float)))
    best = None
    lim = int(search_limit)
    tol = float(orthogonality_tolerance)
    for m in range(-lim, lim + 1):
        for n in range(-lim, lim + 1):
            if m == 0 and n == 0:
                continue
            c1 = np.array([m, n], dtype=int)
            b1 = A @ c1
            l1 = float(np.linalg.norm(b1))
            if l1 <= np.finfo(float).eps:
                continue
            # Deterministic orientation: the first vector points into +X where
            # possible (or +Y if it is vertical).
            if b1[0] < -1.0e-12 or (abs(float(b1[0])) <= 1.0e-12 and b1[1] < 0.0):
                continue
            for p in range(-lim, lim + 1):
                for q in range(-lim, lim + 1):
                    if p == 0 and q == 0:
                        continue
                    c2 = np.array([p, q], dtype=int)
                    coeff = np.column_stack((c1, c2))
                    det_i = int(round(float(np.linalg.det(coeff))))
                    if det_i == 0:
                        continue
                    b2 = A @ c2
                    l2 = float(np.linalg.norm(b2))
                    if l2 <= np.finfo(float).eps:
                        continue
                    physical_det = float(np.linalg.det(np.column_stack((b1, b2))))
                    if physical_det <= 1.0e-12:
                        continue
                    orth_error = abs(float(np.dot(b1, b2))) / (l1 * l2)
                    if orth_error > tol:
                        continue
                    # Minimise primitive-cell count first.  Then prefer a first
                    # vector close to global +X and a reasonably compact rectangle.
                    angle_x = abs(math.atan2(float(b1[1]), float(b1[0])))
                    aspect_penalty = abs(math.log(l1 / l2))
                    score = (
                        abs(det_i), angle_x, aspect_penalty, l1 + l2,
                        tuple(int(x) for x in c1), tuple(int(x) for x in c2),
                    )
                    if best is None or score < best[0]:
                        best = (score, c1, c2, b1, b2, orth_error, abs(det_i))
    if best is None:
        raise ValueError(
            "No commensurate rectangular supercell was found within the integer "
            f"search limit ±{lim}. Increase the search limit or use the primitive-cell PBC mode."
        )
    _, c1, c2, b1, b2, orth_error, count = best
    return {
        "b1_coefficients": c1.astype(int).tolist(),
        "b2_coefficients": c2.astype(int).tolist(),
        "b1": b1.tolist(),
        "b2": b2.tolist(),
        "primitive_cells_in_minimum_rectangle": int(count),
        "orthogonality_relative_error": float(orth_error),
    }


def _clip_segment_to_unit_square(f0: np.ndarray, f1: np.ndarray, tol: float = 1.0e-10):
    """Clip a line segment in fractional coordinates to 0<=u,v<=1."""
    a = np.asarray(f0, float).reshape(2)
    b = np.asarray(f1, float).reshape(2)
    d = b - a
    t0, t1 = 0.0, 1.0
    for axis in range(2):
        if abs(float(d[axis])) <= tol:
            if a[axis] < -tol or a[axis] > 1.0 + tol:
                return None
            continue
        x0 = -a[axis] / d[axis]
        x1 = (1.0 - a[axis]) / d[axis]
        lo, hi = min(x0, x1), max(x0, x1)
        t0, t1 = max(t0, lo), min(t1, hi)
        if t0 > t1 + tol:
            return None
    p0 = np.clip(a + t0 * d, 0.0, 1.0)
    p1 = np.clip(a + t1 * d, 0.0, 1.0)
    if float(np.linalg.norm(p1 - p0)) <= tol:
        return None
    return p0, p1


def build_lattice_supercell(
    model: RVEModel,
    *,
    repeats_a1: int = 2,
    repeats_a2: int = 2,
) -> tuple[RVEModel, dict]:
    """Build an exact N x M lattice-tiled supercell by translating the primitive RVE.

    Unlike the rectangular-cut mode, this routine never clips, stretches, rotates,
    or remaps an RVE wall.  Every copy is the original primitive wall network
    translated by ``i*a1 + j*a2``.  Coincident nodes and genuinely coincident walls
    are merged only to create a conforming FE model.  The periodic boundary is the
    outer parallelogram spanned by ``N*a1`` and ``M*a2``.
    """
    _validate_periodic_fe_model(model)
    na, nb = int(repeats_a1), int(repeats_a2)
    if na < 1 or nb < 1:
        raise ValueError("Lattice-supercell repeats must be at least 1 x 1.")

    A1 = np.asarray(model.a1, float)
    A2 = np.asarray(model.a2, float)
    B1 = na * A1
    B2 = nb * A2
    primitive_area = abs(float(np.linalg.det(np.column_stack((A1, A2)))))
    primitive_count = na * nb
    super_area = abs(float(np.linalg.det(np.column_stack((B1, B2)))))
    expected_area = primitive_count * primitive_area
    area_error = abs(super_area - expected_area) / max(expected_area, np.finfo(float).eps)
    if area_error > 1.0e-12:
        raise ValueError(f"Lattice-supercell area consistency failed: relative error {area_error:.3e}.")

    scale = max(float(np.linalg.norm(B1)), float(np.linalg.norm(B2)), 1.0)
    coord_tol = max(scale * 1.0e-10, 1.0e-9)
    node_lookup: dict[tuple[int, int], int] = {}
    new_nodes: list[np.ndarray] = []
    new_elements: list[Element] = []
    new_groups: list[str] = []
    wall_lookup: dict[tuple[int, int], tuple[int, str]] = {}

    def node_id(xy: np.ndarray) -> int:
        key = tuple(int(x) for x in np.rint(np.asarray(xy, float) / coord_tol))
        if key in node_lookup:
            return node_lookup[key]
        idx = len(new_nodes)
        node_lookup[key] = idx
        new_nodes.append(np.asarray(xy, float).copy())
        return idx

    for ia in range(na):
        for ib in range(nb):
            shift = ia * A1 + ib * A2
            for wi, e in enumerate(model.elements):
                p0 = np.asarray(model.nodes[e.node_i], float) + shift
                p1 = np.asarray(model.nodes[e.node_j], float) + shift
                n0, n1 = node_id(p0), node_id(p1)
                if n0 == n1:
                    raise ValueError("A translated supercell wall collapsed to zero length.")
                key = tuple(sorted((n0, n1)))
                group = model.wall_groups[wi] if wi < len(model.wall_groups) else ""
                if key in wall_lookup:
                    old_lam, old_group = wall_lookup[key]
                    if old_lam != int(e.laminate):
                        raise ValueError(
                            "Coincident translated walls have different laminate assignments; "
                            "the primitive RVE is not materially periodic."
                        )
                    continue
                wall_lookup[key] = (int(e.laminate), group)
                new_elements.append(Element(n0, n1, int(e.laminate), 1.0))
                new_groups.append(group)

    nodes_array = np.asarray(new_nodes, float)
    primitive_wall_length = float(sum(
        np.linalg.norm(np.asarray(model.nodes[e.node_j], float) - np.asarray(model.nodes[e.node_i], float))
        for e in model.elements
    ))
    actual_wall_length = float(sum(
        np.linalg.norm(nodes_array[e.node_j] - nodes_array[e.node_i]) for e in new_elements
    ))
    expected_wall_length = primitive_count * primitive_wall_length
    wall_error = abs(actual_wall_length - expected_wall_length) / max(expected_wall_length, np.finfo(float).eps)
    if wall_error > 1.0e-10:
        raise ValueError(
            "Exact lattice tiling changed the primitive-cell wall measure. "
            f"Expected {expected_wall_length:.12g}, got {actual_wall_length:.12g}."
        )

    super_model = RVEModel(
        nodes_array, new_elements, B1, B2, float(model.height),
        periodic_tol=max(float(model.periodic_tol), 1.0e-9), anchor_node=0,
        fix_anchor_rotation=False, wall_element_type=model.wall_element_type,
        strip_bending_mode=model.strip_bending_mode, length_unit=model.length_unit,
        gravity=model.gravity, name=f"{model.name} - lattice supercell {na}x{nb}",
        wall_groups=new_groups,
    )
    diag = {
        "domain_mode": PBC_DOMAIN_LATTICE_SUPERCELL,
        "construction": "exact primitive-RVE translations; no clipping or geometric remapping",
        "primitive_a1": A1.tolist(), "primitive_a2": A2.tolist(),
        "repeat_a1": na, "repeat_a2": nb,
        "repeat_x": na, "repeat_y": nb,
        "primitive_cell_count": primitive_count,
        "supercell_b1": B1.tolist(), "supercell_b2": B2.tolist(),
        "primitive_area": primitive_area, "supercell_area": super_area,
        "expected_area_from_cell_count": expected_area,
        "area_consistency_relative_error": area_error,
        "primitive_wall_length": primitive_wall_length,
        "expected_supercell_wall_length": expected_wall_length,
        "actual_supercell_wall_length": actual_wall_length,
        "wall_measure_relative_error": wall_error,
        "generated_nodes": int(len(nodes_array)), "generated_walls": int(len(new_elements)),
    }
    return super_model, diag


def build_rectangular_supercell(
    model: RVEModel,
    *,
    repeats_x: int = 2,
    repeats_y: int = 2,
    search_limit: int = 6,
    orthogonality_tolerance: float = 1.0e-8,
    prefer_mirror_symmetry: bool = False,
) -> tuple[RVEModel, dict]:
    """Repeat a primitive periodic network and cut a true rectangular supercell.

    A small phase shift is deliberately used before clipping.  This avoids putting
    a complete shell wall exactly on both opposite supercell boundaries, which
    would double-count that wall's strain energy in a zero-thickness shell model.
    The generated wall-length measure is checked against the exact integer number
    of primitive cells before the model is accepted.
    """
    _validate_periodic_fe_model(model)
    rx, ry = int(repeats_x), int(repeats_y)
    if rx < 1 or ry < 1:
        raise ValueError("Rectangular supercell repeats must be at least 1 x 1.")
    basis = find_rectangular_supercell_basis(
        model, search_limit=search_limit, orthogonality_tolerance=orthogonality_tolerance
    )
    A = np.column_stack((np.asarray(model.a1, float), np.asarray(model.a2, float)))
    c1 = np.asarray(basis["b1_coefficients"], int)
    c2 = np.asarray(basis["b2_coefficients"], int)
    b1_min = A @ c1
    b2_min = A @ c2
    B1 = rx * b1_min
    B2 = ry * b2_min
    B = np.column_stack((B1, B2))
    primitive_area = abs(float(np.linalg.det(A)))
    base_count = int(basis["primitive_cells_in_minimum_rectangle"])
    primitive_count = int(base_count * rx * ry)
    super_area = abs(float(np.linalg.det(B)))
    expected_area = primitive_count * primitive_area
    area_error = abs(super_area - expected_area) / max(expected_area, np.finfo(float).eps)
    if area_error > 1.0e-10:
        raise ValueError(f"Supercell area consistency failed: relative error {area_error:.3e}.")

    used = sorted({int(e.node_i) for e in model.elements} | {int(e.node_j) for e in model.elements})
    primitive_frac = np.linalg.solve(A, np.asarray(model.nodes, float).T).T
    used_frac = primitive_frac[used]
    cu, cv = rx * c1, ry * c2
    corners = np.array([[0.0, 0.0], cu, cv, cu + cv], float)
    min_corner, max_corner = corners.min(axis=0), corners.max(axis=0)
    min_node, max_node = used_frac.min(axis=0), used_frac.max(axis=0)
    i0 = math.floor(float(min_corner[0] - max_node[0])) - 2
    i1 = math.ceil(float(max_corner[0] - min_node[0])) + 2
    j0 = math.floor(float(min_corner[1] - max_node[1])) - 2
    j1 = math.ceil(float(max_corner[1] - min_node[1])) + 2

    primitive_wall_length = float(sum(
        np.linalg.norm(np.asarray(model.nodes[e.node_j], float) - np.asarray(model.nodes[e.node_i], float))
        for e in model.elements
    ))
    expected_wall_length = primitive_count * primitive_wall_length
    scale = max(float(np.linalg.norm(B1)), float(np.linalg.norm(B2)), 1.0)
    coord_tol = max(scale * 1.0e-10, 1.0e-9)

    # Generic offsets work better than exact half-cell shifts for square/triangular
    # networks where a half shift can still place complete walls on the boundary.
    phase_candidates = [
        (0.371, 0.229), (0.173, 0.417), (0.293, 0.613),
        (0.5, 0.5), (0.25, 0.25), (0.17, 0.43),
    ]
    if prefer_mirror_symmetry:
        # Rational phase shifts place candidate mirrors on cell walls, nodes,
        # or cell centres. Only accept those meeting all geometric audits.
        phase_candidates = [(i / 8.0, j / 8.0) for i in range(8) for j in range(8)] + phase_candidates


    def build_for_phase(phase):
        origin = A @ np.asarray(phase, float)
        node_lookup: dict[tuple[int, int], int] = {}
        new_nodes: list[np.ndarray] = []
        new_elements: list[Element] = []
        new_groups: list[str] = []
        wall_lookup: dict[tuple[int, int], tuple[int, str]] = {}

        def node_id(xy: np.ndarray) -> int:
            key = tuple(int(x) for x in np.rint(np.asarray(xy, float) / coord_tol))
            if key in node_lookup:
                return node_lookup[key]
            idx = len(new_nodes)
            node_lookup[key] = idx
            new_nodes.append(np.asarray(xy, float))
            return idx

        for ii in range(i0, i1 + 1):
            for jj in range(j0, j1 + 1):
                shift = A @ np.array([ii, jj], float) - origin
                for wi, e in enumerate(model.elements):
                    p0 = np.asarray(model.nodes[e.node_i], float) + shift
                    p1 = np.asarray(model.nodes[e.node_j], float) + shift
                    clipped = _clip_segment_to_unit_square(
                        np.linalg.solve(B, p0), np.linalg.solve(B, p1)
                    )
                    if clipped is None:
                        continue
                    q0 = B @ clipped[0]
                    q1 = B @ clipped[1]
                    n0, n1 = node_id(q0), node_id(q1)
                    if n0 == n1:
                        continue
                    key = tuple(sorted((n0, n1)))
                    group = model.wall_groups[wi] if wi < len(model.wall_groups) else ""
                    if key in wall_lookup:
                        old_lam, old_group = wall_lookup[key]
                        if old_lam != int(e.laminate):
                            raise ValueError(
                                "A supercell boundary produced coincident walls with different laminate assignments. "
                                "The primitive RVE is not materially periodic."
                            )
                        continue
                    wall_lookup[key] = (int(e.laminate), group)
                    new_elements.append(Element(n0, n1, int(e.laminate), 1.0))
                    new_groups.append(group)
        if not new_elements:
            raise ValueError("Rectangular supercell clipping produced no physical walls.")
        nodes_array = np.asarray(new_nodes, float)
        total_length = float(sum(
            np.linalg.norm(nodes_array[e.node_j] - nodes_array[e.node_i]) for e in new_elements
        ))
        rel_length_error = abs(total_length - expected_wall_length) / max(expected_wall_length, np.finfo(float).eps)
        return nodes_array, new_elements, new_groups, total_length, rel_length_error, origin

    built = []
    for phase in phase_candidates:
        try:
            candidate = build_for_phase(phase)
            if prefer_mirror_symmetry and candidate[4] <= 1.0e-7:
                # Local import avoids introducing a module import cycle.
                from .rectangular_selection import _reflection_error
                tmp_nodes,tmp_elements,tmp_groups,_,_,_=candidate
                tmp = RVEModel(tmp_nodes, tmp_elements, np.asarray(B1,float),np.asarray(B2,float),float(model.height))
                err = _reflection_error(tmp)
                symmetry_penalty = max(err['x_geometry_relative'],err['y_geometry_relative'])
                built.append((candidate[4], symmetry_penalty, phase, candidate))
            elif not prefer_mirror_symmetry:
                built.append((candidate[4], 0.0, phase, candidate))
        except Exception:
            continue
    if not built:
        raise ValueError("Unable to construct a clipped rectangular periodic supercell from the primitive network.")
    # For a symmetric selection, reflection error dominates, provided the
    # required physical wall measure and periodic rectangle are correct.
    if prefer_mirror_symmetry:
        built.sort(key=lambda row: (row[1],row[0]))
    else:
        built.sort(key=lambda row: row[0])
    length_error, symmetry_error, phase, candidate = built[0]
    nodes_array, elements, groups, total_length, _, origin = candidate
    if length_error > 1.0e-7:
        raise ValueError(
            "Rectangular supercell wall-measure check failed. The clipped domain does not contain exactly the expected "
            f"periodic material amount (relative wall-length error {length_error:.3e})."
        )

    super_model = RVEModel(
        nodes_array, elements, np.asarray(B1, float), np.asarray(B2, float), float(model.height),
        periodic_tol=max(float(model.periodic_tol), 1.0e-9), anchor_node=0,
        fix_anchor_rotation=False, wall_element_type=model.wall_element_type,
        strip_bending_mode=model.strip_bending_mode, length_unit=model.length_unit,
        gravity=model.gravity, name=f"{model.name} - rectangular supercell {rx}x{ry}",
        wall_groups=groups,
    )
    diag = {
        "domain_mode": PBC_DOMAIN_RECTANGULAR_SUPERCELL,
        "primitive_a1": np.asarray(model.a1, float).tolist(),
        "primitive_a2": np.asarray(model.a2, float).tolist(),
        "minimum_rectangle_b1_coefficients": basis["b1_coefficients"],
        "minimum_rectangle_b2_coefficients": basis["b2_coefficients"],
        "minimum_rectangle_primitive_cells": base_count,
        "repeat_x": rx, "repeat_y": ry,
        "primitive_cell_count": primitive_count,
        "supercell_b1": np.asarray(B1, float).tolist(),
        "supercell_b2": np.asarray(B2, float).tolist(),
        "rectangle_width": float(np.linalg.norm(B1)),
        "rectangle_height": float(np.linalg.norm(B2)),
        "orthogonality_relative_error": abs(float(np.dot(B1, B2))) / (float(np.linalg.norm(B1)) * float(np.linalg.norm(B2))),
        "primitive_area": primitive_area,
        "supercell_area": super_area,
        "expected_area_from_cell_count": expected_area,
        "area_consistency_relative_error": area_error,
        "primitive_wall_length": primitive_wall_length,
        "expected_supercell_wall_length": expected_wall_length,
        "actual_supercell_wall_length": total_length,
        "wall_measure_relative_error": float(length_error),
        "phase_shift_primitive_fraction": [float(phase[0]), float(phase[1])],
        "phase_shift_cartesian": np.asarray(origin, float).tolist(),
        "mirror_geometry_relative_error": float(symmetry_error) if prefer_mirror_symmetry else None,
        "generated_nodes": int(len(nodes_array)),
        "generated_walls": int(len(elements)),
    }
    return super_model, diag


def build_published_periodic_domain(model: RVEModel) -> tuple[RVEModel, dict]:
    """The source-study DBC block (thesis Fig. 3.6b / paper Fig. 8b) used as a periodic cell.

    Walls lying in a face appear twice (x1 = 0 and x1 = p, or x2 = 0 and x2 = q) and
    are exact periodic images of each other.  Both copies are kept so that every
    face node has its periodic mate, and each copy is given half of the stiffness
    (all moduli halved, same plies), so that together they act exactly as one full wall.
    """
    from .fe_rve import build_published_fe_rve
    block, d = build_published_fe_rve(model, face_wall_section="half-stiffness")
    p, q = float(d["p"]), float(d["q"])
    prim_area = float(d["primitive_area"])
    count = int(d["primitive_cells"])
    super_model = RVEModel(
        np.asarray(block.nodes, float), list(block.elements), np.array([p, 0.0]), np.array([0.0, q]),
        float(model.height), periodic_tol=max(float(model.periodic_tol), 1.0e-9), anchor_node=0,
        fix_anchor_rotation=False, wall_element_type=model.wall_element_type,
        strip_bending_mode=model.strip_bending_mode, length_unit=model.length_unit,
        gravity=model.gravity, name=f"{model.name} - source-study RVE (periodic)", wall_groups=list(block.wall_groups),
    )
    diag = {
        "domain_mode": PBC_DOMAIN_PUBLISHED,
        "primitive_a1": np.asarray(model.a1, float).tolist(),
        "primitive_a2": np.asarray(model.a2, float).tolist(),
        "primitive_cell_count": count,
        "primitive_area": prim_area,
        "supercell_area": p * q,
        "expected_area_from_cell_count": count * prim_area,
        "area_consistency_relative_error": abs(p * q - count * prim_area) / (count * prim_area),
        "repeat_x": int(d["cells_x1_x2"][0]), "repeat_y": int(d["cells_x1_x2"][1]),
        "supercell_b1": [p, 0.0], "supercell_b2": [0.0, q],
        "rectangle_width": p, "rectangle_height": q,
        "primitive_wall_length": float(d["primitive_wall_length"]),
        "expected_supercell_wall_length": float(d["expected_wall_length"]),
        "actual_supercell_wall_length": float(d["actual_wall_length"]),
        "effective_wall_length": float(d["effective_wall_length"]),
        "wall_measure_relative_error": float(d["wall_measure_relative_error"]),
        "face_walls": d["face_walls"], "face_wall_faces": d["face_wall_faces"],
        "face_wall_section": "half-stiffness", "wall_section_factors": d["wall_section_factors"],
        "rve": d["rve"], "rve_rule": d["rve_rule"],
        "central_cell_polygon_rve": d["central_cell_polygon_rve"],
        "generated_nodes": int(len(block.nodes)), "generated_walls": int(len(block.elements)),
    }
    return super_model, diag


def resolve_periodic_domain(model: RVEModel, settings: PeriodicHomogenizationSettings) -> tuple[RVEModel, dict]:
    """Return the physical FE domain requested for periodic homogenisation."""
    settings.validate()
    if settings.domain_mode == PBC_DOMAIN_PUBLISHED:
        try:
            return build_published_periodic_domain(model)
        except ValueError as exc:
            # No rectangular block exists for this lattice (or no complete cell fits in 3 x 3 cells):
            # fall back to exact copies of the primitive cell, which is always periodic.
            dom, diag = build_lattice_supercell(model, repeats_a1=settings.supercell_repeat_x,
                                                repeats_a2=settings.supercell_repeat_y)
            diag["fallback_from"] = PBC_DOMAIN_PUBLISHED
            diag["fallback_reason"] = str(exc)
            return dom, diag
    if settings.domain_mode == PBC_DOMAIN_PRIMITIVE:
        area = abs(float(np.linalg.det(np.column_stack((model.a1, model.a2)))))
        return model, {
            "domain_mode": PBC_DOMAIN_PRIMITIVE,
            "primitive_cell_count": 1,
            "primitive_area": area,
            "supercell_area": area,
            "area_consistency_relative_error": 0.0,
            "repeat_x": 1, "repeat_y": 1,
        }
    if settings.domain_mode == PBC_DOMAIN_LATTICE_SUPERCELL:
        return build_lattice_supercell(
            model,
            repeats_a1=settings.supercell_repeat_x,
            repeats_a2=settings.supercell_repeat_y,
        )
    if settings.domain_mode == PBC_DOMAIN_SYMMETRIC_RECTANGULAR:
        from .rectangular_selection import inspect_rectangular_rve
        selected, diag = inspect_rectangular_rve(
            model, repeats_x=settings.supercell_repeat_x,
            repeats_y=settings.supercell_repeat_y,
            search_limit=settings.supercell_search_limit,
            symmetry_tolerance=1e-6,
        )
        if not diag['mirror_geometry_symmetric']:
            raise ValueError('Requested a symmetric rectangular FE domain, but neither phase-shifted '
                             'rectangle passed the two-axis mirror audit. Select another geometry '
                             '(e.g. symmetric SQUARE or MIXED) or use the unrestricted rectangular cut.')
        diag['domain_mode']=PBC_DOMAIN_SYMMETRIC_RECTANGULAR
        return selected, diag
    return build_rectangular_supercell(
        model,
        repeats_x=settings.supercell_repeat_x,
        repeats_y=settings.supercell_repeat_y,
        search_limit=settings.supercell_search_limit,
        orthogonality_tolerance=settings.supercell_orthogonality_tolerance,
    )


def periodic_bc_audit(
    nodes: dict[int, tuple[float, float, float]],
    wall_elements: dict[int, list[tuple[int, tuple[int, int, int, int]]]],
    model: RVEModel,
    pairs: list[tuple[int, int, np.ndarray]],
) -> dict:
    """Audit PBC pairing, connectivity and rigid-body rank before Abaqus runs."""
    ids = sorted(nodes)
    xyz = np.array([nodes[i] for i in ids], float)
    id_to_row = {nid: k for k, nid in enumerate(ids)}
    frac = _fractional_coordinates(xyz, model)
    ftol = max(1.0e-7, float(model.periodic_tol) * 10.0)

    parent = {nid: nid for nid in ids}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[rb] = ra
    for master, slave, _ in pairs:
        union(int(master), int(slave))
    groups: dict[int, set[int]] = {}
    for nid in ids:
        groups.setdefault(find(nid), set()).add(nid)

    coverage = {}
    unpaired: set[int] = set()
    for axis, label in enumerate(("a1", "a2", "a3")):
        lower = [ids[k] for k in range(len(ids)) if abs(float(frac[k, axis])) <= ftol]
        upper = [ids[k] for k in range(len(ids)) if abs(float(frac[k, axis] - 1.0)) <= ftol]
        bad = []
        for nid in lower:
            g = groups.get(find(nid), {nid})
            if not any(abs(float(frac[id_to_row[x], axis] - 1.0)) <= ftol for x in g):
                bad.append(nid)
        for nid in upper:
            g = groups.get(find(nid), {nid})
            if not any(abs(float(frac[id_to_row[x], axis])) <= ftol for x in g):
                bad.append(nid)
        unpaired.update(bad)
        coverage[label] = {
            "lower_boundary_nodes": len(lower),
            "upper_boundary_nodes": len(upper),
            "unpaired_nodes": len(set(bad)),
            "complete": len(set(bad)) == 0 and len(lower) > 0 and len(upper) > 0,
        }

    deltas = np.array([np.asarray(d, float) for _, _, d in pairs], float) if pairs else np.empty((0, 3))
    translation_rank = int(np.linalg.matrix_rank(deltas, tol=1.0e-9)) if len(deltas) else 0

    # Rigid-body matrix q=[tx,ty,tz,wx,wy,wz].  Fixed anchor translations plus
    # homogeneous periodic displacement jumps must have rank six.
    anchor_id = min(ids)
    xa = np.asarray(nodes[anchor_id], float)
    rb_rows = []
    for dof in range(3):
        row = np.zeros(6); row[dof] = 1.0
        # omega x x = -[x]_x omega
        cross_coeff = np.array([[0.0, xa[2], -xa[1]], [-xa[2], 0.0, xa[0]], [xa[1], -xa[0], 0.0]])
        row[3:] = cross_coeff[dof]
        rb_rows.append(row)
    for _, _, delta in pairs:
        dx = np.asarray(delta, float)
        cross_coeff = np.array([[0.0, dx[2], -dx[1]], [-dx[2], 0.0, dx[0]], [dx[1], -dx[0], 0.0]])
        for dof in range(3):
            row = np.zeros(6); row[3:] = cross_coeff[dof]; rb_rows.append(row)
    rigid_rank = int(np.linalg.matrix_rank(np.asarray(rb_rows, float), tol=1.0e-9))

    # FE connectivity using S4 element edges.
    adjacency = {nid: set() for nid in ids}
    used_mesh_nodes = set()
    for elems in wall_elements.values():
        for _, conn in elems:
            cc = list(conn); used_mesh_nodes.update(cc)
            for a, b in zip(cc, cc[1:] + cc[:1]):
                adjacency[a].add(b); adjacency[b].add(a)
    def component_count(graph):
        count = 0; seen = set()
        for start in sorted(used_mesh_nodes):
            if start in seen: continue
            count += 1; stack = [start]; seen.add(start)
            while stack:
                cur = stack.pop()
                for nxt in graph[cur]:
                    if nxt not in seen: seen.add(nxt); stack.append(nxt)
        return count

    physical_components = component_count(adjacency)
    periodic_graph = {nid: set(nei) for nid, nei in adjacency.items()}
    for master, slave, _ in pairs:
        periodic_graph[int(master)].add(int(slave)); periodic_graph[int(slave)].add(int(master))
    constrained_components = component_count(periodic_graph)

    test_eps = np.array([0.13, -0.07, 0.05, 0.03, -0.02, 0.11])
    E = _strain_tensor_from_engineering(test_eps)
    affine_residual = 0.0
    for _, _, delta in pairs:
        exact = E @ np.asarray(delta, float)
        via_controls = _macro_coefficients(np.asarray(delta, float), 1.0) @ test_eps
        affine_residual = max(affine_residual, float(np.linalg.norm(exact - via_controls, ord=np.inf)))

    pair_keys = [(min(int(a), int(b)), max(int(a), int(b))) for a, b, _ in pairs]
    duplicate_pairs = len(pair_keys) - len(set(pair_keys))
    passed = (
        constrained_components == 1 and translation_rank == 3 and rigid_rank == 6 and
        not unpaired and duplicate_pairs == 0 and all(v["complete"] for v in coverage.values()) and
        affine_residual <= 1.0e-10
    )
    return {
        "passed": bool(passed),
        "mesh_nodes": len(nodes),
        "mesh_elements": int(sum(len(x) for x in wall_elements.values())),
        "periodic_pair_count": len(pairs),
        "physical_connected_components": int(physical_components),
        "connected_components": int(constrained_components),
        "periodic_translation_rank": translation_rank,
        "rigid_body_constraint_rank": rigid_rank,
        "unpaired_boundary_nodes": len(unpaired),
        "duplicate_periodic_pairs": int(duplicate_pairs),
        "affine_field_max_residual": float(affine_residual),
        "boundary_coverage": coverage,
    }

def _fractional_coordinates(xyz: np.ndarray, model: RVEModel) -> np.ndarray:
    A = np.array([
        [float(model.a1[0]), float(model.a2[0]), 0.0],
        [float(model.a1[1]), float(model.a2[1]), 0.0],
        [0.0, 0.0, float(model.height)],
    ], float)
    return np.linalg.solve(A, np.asarray(xyz, float).T).T


def _mod1(v: np.ndarray, tol: float) -> np.ndarray:
    out = np.mod(np.asarray(v, float), 1.0)
    out[np.abs(out) <= tol] = 0.0
    out[np.abs(out - 1.0) <= tol] = 0.0
    return out


def _wrapped_inf_distance(a: np.ndarray, b: np.ndarray) -> float:
    d = np.abs(np.asarray(a) - np.asarray(b))
    d = np.minimum(d, 1.0 - d)
    return float(np.max(d))


def periodic_mesh_pairs(nodes: dict[int, tuple[float, float, float]], model: RVEModel, tol: float | None = None) -> list[tuple[int, int, np.ndarray]]:
    """Find periodic mesh mates using a periodic 3-D fractional-coordinate KD-tree.

    This is near-linear in normal shell meshes, instead of checking every pair of
    mesh nodes (quadratic time). Every candidate is checked against integer
    lattice translation to prevent false periodic ties from the search radius.
    """
    from scipy.spatial import cKDTree
    if tol is None:
        tol = max(float(model.periodic_tol), 1.0e-8)
    tol = float(tol)
    ids = sorted(nodes)
    xyz = np.array([nodes[i] for i in ids], dtype=float)
    frac = _fractional_coordinates(xyz, model)
    frac_mod = _mod1(frac, tol)
    # Use the wrapped (toroidal) metric to recognise mates on 0/1 boundaries.
    tree = cKDTree(frac_mod, boxsize=1.0)
    candidate_pairs = tree.query_pairs(tol * np.sqrt(3.0) * (1.0 + 1.0e-7), output_type='ndarray')
    parent = np.arange(len(ids), dtype=np.int64)
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = int(parent[x])
        return x
    for i, j in candidate_pairs:
        i = int(i); j = int(j)
        if _wrapped_inf_distance(frac_mod[j], frac_mod[i]) > tol:
            continue
        df = frac[j] - frac[i]
        if np.linalg.norm(df - np.rint(df), ord=np.inf) > max(tol * 5.0, 1e-7):
            continue
        ri = find(i); rj = find(j)
        if ri != rj:
            parent[rj] = ri
    groups: dict[int, list[int]] = {}
    for i in range(len(ids)):
        groups.setdefault(find(i), []).append(i)
    pairs = []
    for group in sorted(groups.values(), key=lambda g: min(g)):
        if len(group) < 2:
            continue
        group.sort()
        master = group[0]
        for slave in group[1:]:
            delta = xyz[slave] - xyz[master]
            if np.linalg.norm(delta) > 1e-12:
                pairs.append((ids[master], ids[slave], delta))
    return pairs


def _equation_lines(terms: list[tuple[int, int, float]]) -> list[str]:
    terms = [(int(n), int(d), float(c)) for n, d, c in terms if abs(float(c)) > 1e-14]
    if len(terms) < 2:
        raise ValueError("A periodic equation must contain at least two non-zero terms.")
    lines = ["*EQUATION", str(len(terms))]
    # Abaqus accepts multiple node,dof,coefficient triplets; wrap for readability/old readers.
    for k in range(0, len(terms), 4):
        chunk = terms[k:k + 4]
        lines.append(", ".join(f"{n}, {d}, {c:.16g}" for n, d, c in chunk))
    return lines


def _macro_coefficients(delta: np.ndarray, reference_length: float) -> np.ndarray:
    """Return 3x6 coefficients mapping macro control displacements to displacement jump.

    Each macro control node U1 equals engineering_strain * reference_length.
    """
    dx = np.asarray(delta, float).reshape(3)
    B = np.zeros((3, 6), float)
    for k in range(6):
        unit = np.zeros(6); unit[k] = 1.0
        B[:, k] = _strain_tensor_from_engineering(unit) @ dx / float(reference_length)
    return B


def generate_periodic_homogenization_input(
    model: RVEModel,
    materials: dict[str, Material],
    laminates: list[Laminate],
    settings: PeriodicHomogenizationSettings,
    strain: np.ndarray,
    *,
    job_name: str = "RVE_HOM",
    resolved_domain: bool = False,
    domain_diagnostics: dict | None = None,
) -> tuple[str, dict]:
    """Generate one Abaqus/Standard S4 periodic macro-strain case.

    The cell is periodic in the two in-plane lattice vectors AND in the height
    vector a3=(0,0,h). This describes a 3-D repeating material surrogate,
    not a finite-height honeycomb core with physical top/bottom surfaces. Six scalar control nodes represent the engineering
    strain vector [e11,e22,e33,g23,g13,g12].
    """
    settings.validate()
    if resolved_domain:
        # The caller has already constructed the exact physical FE domain.  Do
        # not re-create it from the primitive model: the same resolved geometry
        # must be written to every one of the 21 Abaqus decks.
        domain_diag = dict(domain_diagnostics or {})
        if not domain_diag:
            area = abs(float(np.linalg.det(np.column_stack((model.a1, model.a2)))))
            domain_diag = {
                "domain_mode": PBC_DOMAIN_PRIMITIVE,
                "primitive_cell_count": 1,
                "primitive_area": area,
                "supercell_area": area,
            }
    else:
        source_model = model
        model, domain_diag = resolve_periodic_domain(source_model, settings)
    _validate_periodic_fe_model(model)
    strain = np.asarray(strain, float).reshape(6)
    if not laminates:
        raise ValueError("At least one laminate stack is required for Abaqus homogenisation.")
    if settings.domain_mode == PBC_DOMAIN_SYMMETRIC_RECTANGULAR:
        from .rectangular_selection import _reflection_error, _signature
        mirror = _reflection_error(model, [_signature(l) for l in laminates])
        if any(mirror[axis+'_same_stack_relative'] > 1e-6 for axis in ('x','y')):
            raise ValueError('Rectangular FE wall graph is geometrically symmetric but the material/ply stacks '
                             'do not match on mirrored wall segments.')
        # A 0/90 ply axis is invariant under in-plane mirror reversal; off-axis
        # angles require an orientation-aware mirrored stack, not yet certified.
        if any(abs(np.sin(np.deg2rad(2.0 * float(p.angle_deg)))) > 1e-8
               for lay in laminates for p in lay.plies):
            raise ValueError('Symmetry-certified rectangular FE mode currently requires 0/90-degree '
                             'ply angles. For off-axis stacks select the unrestricted rectangular '
                             'periodic cut; it remains valid as a PBC cell but is not reflection-certified.')
        domain_diag.update({'mirror_stack_symmetric':True,'mirror_ply_axis_symmetric':True,
                            'mirror_geometry_symmetric':True})

    nodes, wall_elements, _, _ = _build_shell_mesh(model, settings.mesh_size)
    if not nodes:
        raise ValueError("The shell mesh contains no nodes.")
    pairs = periodic_mesh_pairs(nodes, model)
    if not pairs:
        raise ValueError(
            "No periodic mesh-node pairs were found. Check a1/a2, the cell boundary and whether the drawn RVE contains translated boundary copies."
        )
    audit = periodic_bc_audit(nodes, wall_elements, model, pairs)
    if settings.strict_pbc_audit and not audit["passed"]:
        raise ValueError(
            "Periodic boundary-condition audit failed before Abaqus input generation: "
            f"components={audit['connected_components']}, translation rank={audit['periodic_translation_rank']}/3, "
            f"rigid-body rank={audit['rigid_body_constraint_rank']}/6, "
            f"unpaired boundary nodes={audit['unpaired_boundary_nodes']}, "
            f"duplicate pairs={audit['duplicate_periodic_pairs']}."
        )

    max_nid = max(nodes)
    macro_ids = [max_nid + i for i in range(1, 7)]
    anchor_id = min(nodes)
    factors = [float(x) for x in (domain_diag.get("wall_section_factors") or [1.0] * len(model.elements))]
    if len(factors) != len(model.elements):
        raise ValueError("Wall section factors do not match the periodic wall list.")
    mats_all = dict(materials)
    if any(f < 1.0 for f in factors):
        # Periodic images of a face wall: two copies, each with half of the moduli.
        for key, mat in materials.items():
            dens = float(mat.density_kg_m3)
            mats_all[key + "__HALF"] = Material(
                name=f"{mat.name}_HALF", E1=0.5 * mat.E1, E2=0.5 * mat.E2, G12=0.5 * mat.G12, nu12=mat.nu12,
                Gth=0.5 * float(mat.Gth if mat.Gth is not None else mat.G12),
                density_kg_m3=0.5 * dens if math.isfinite(dens) else math.nan)
    material_text, mat_names = _abaqus_material_blocks(mats_all)

    out: list[str] = [
        "** RVE Studio Python - full periodic energy homogenisation",
        f"** Job: {job_name}",
        "** Engineering strain order: E11,E22,E33,G23,G13,G12",
        "** CAUTION: 3-D PERIODIC MATERIAL SURROGATE. Opposite z=0,h faces are artificially paired.",
        "** Source-matched finite-height core FE uses thesis Table 3.4 displacement BC, not this model.",
        "** Shear components use the symmetric engineering-strain convention (gamma=2*epsilon_ij).",
        "** PBC: u+ - u- = E_macro . (x+ - x-); shell rotations are periodic. Only anchor translations are fixed; anchor rotations remain free.",
        f"** Periodic domain: {domain_diag.get('domain_mode', PBC_DOMAIN_PRIMITIVE)}; primitive-cell count: {domain_diag.get('primitive_cell_count', 1)}",
        "*HEADING",
        f"RVE Studio periodic S4 homogenisation - {model.name}",
        "*NODE",
    ]
    for nid, xyz in sorted(nodes.items()):
        out.append(f"{nid}, {xyz[0]:.12g}, {xyz[1]:.12g}, {xyz[2]:.12g}")
    # Six independent scalar control nodes, all using DOF 1.
    for k, nid in enumerate(macro_ids):
        out.append(f"{nid}, 0., 0., 0.")

    out.append("*ELEMENT, TYPE=S4")
    for wi in sorted(wall_elements):
        for eid, conn in wall_elements[wi]:
            out.append(f"{eid}, {conn[0]}, {conn[1]}, {conn[2]}, {conn[3]}")
    for wi, elems in wall_elements.items():
        out.extend([f"*ELSET, ELSET=WALL_{wi}", _ids_lines(eid for eid, _ in elems)])
    out.extend(["*NSET, NSET=PBC_ANCHOR", str(anchor_id)])
    for k, nid in enumerate(macro_ids):
        out.extend([f"*NSET, NSET=MACRO_{k + 1}", str(nid)])
    out.extend(["*NSET, NSET=MACRO_ALL", _ids_lines(macro_ids)])
    out.append(material_text)

    for wi, e in enumerate(model.elements, 1):
        if e.laminate < 0 or e.laminate >= len(laminates):
            raise ValueError(f"Wall W{wi} references missing laminate index {e.laminate}.")
        lam = laminates[e.laminate]
        p0 = np.asarray(model.nodes[e.node_i], float); p1 = np.asarray(model.nodes[e.node_j], float)
        d = p1 - p0; d /= np.linalg.norm(d)
        out.extend([
            f"*ORIENTATION, NAME=ORI_W{wi}, SYSTEM=RECTANGULAR",
            f"{d[0]:.12g}, {d[1]:.12g}, 0., 0., 0., 1.",
            "3, 0.",  # project the ply 1-axis about the wall normal (local 3): ply 1 = along the wall, not the height
            f"*SHELL SECTION, ELSET=WALL_{wi}, COMPOSITE, ORIENTATION=ORI_W{wi}",
        ])
        half = factors[wi - 1] < 1.0
        if half:
            out.insert(len(out) - 4, f"** Wall W{wi} lies in an RVE face; its periodic image is also present. Each copy carries half of the moduli.")
        for pi, ply in enumerate(lam.plies, 1):
            if ply.material not in mat_names:
                raise ValueError(f"Wall W{wi} ply {pi} references material '{ply.material}', which is not in the material library.")
            mname = mat_names[ply.material + "__HALF"] if half else mat_names[ply.material]
            out.append(f"{float(ply.thickness):.12g}, 3, {mname}, {float(ply.angle_deg):.12g}, PLY{pi}")

    # Periodic translations and rotations. The macro control coefficients provide the affine jump.
    for master, slave, delta in pairs:
        B = _macro_coefficients(delta, settings.reference_length)
        for dof in range(1, 4):
            terms = [(slave, dof, 1.0), (master, dof, -1.0)]
            for k, mid in enumerate(macro_ids):
                c = -B[dof - 1, k]
                if abs(c) > 1e-14:
                    terms.append((mid, 1, c))
            out.extend(_equation_lines(terms))
        for dof in range(4, 7):
            out.extend(_equation_lines([(slave, dof, 1.0), (master, dof, -1.0)]))

    amp = float(settings.strain_amplitude)
    out.extend([
        "*STEP, NAME=HOM, NLGEOM=NO",
        "*STATIC",
        "0.1, 1.0, 1e-08, 0.1",
        "*BOUNDARY",
        "PBC_ANCHOR, 1, 3, 0.",
    ])
    for k in range(6):
        value = float(strain[k]) * float(settings.reference_length)
        out.append(f"MACRO_{k + 1}, 1, 1, {value:.16g}")
        # Stabilise unused translational directions of the scalar control node.
        out.append(f"MACRO_{k + 1}, 2, 3, 0.")
    out.extend([
        "*OUTPUT, FIELD",
        "*NODE OUTPUT, NSET=MACRO_ALL",
        "U, RF",
        "*ELEMENT OUTPUT",
        "S, E",
        "*OUTPUT, HISTORY",
        "*ENERGY OUTPUT",
        "ALLSE",
        "*END STEP",
    ])

    meta = {
        "macro_node_ids": macro_ids,
        "anchor_node_id": anchor_id,
        "periodic_pair_count": len(pairs),
        "boundary_physics": "3d_periodic_material_surrogate",
        "top_bottom_periodic": True,
        "finite_height_core_validation": False,
        "mesh_node_count": len(nodes),
        "mesh_element_count": sum(len(v) for v in wall_elements.values()),
        "strain": strain.tolist(),
        "strain_amplitude": amp,
        "reference_length": float(settings.reference_length),
        "periodic_domain": domain_diag,
        "pbc_audit": audit,
        "analysis_a1": np.asarray(model.a1, float).tolist(),
        "analysis_a2": np.asarray(model.a2, float).tolist(),
        "analysis_geometry_node_count": int(len(model.nodes)),
        "analysis_geometry_wall_count": int(len(model.elements)),
        "analysis_geometry_bounds": {
            "xmin": float(np.min(np.asarray(model.nodes, float)[:, 0])),
            "xmax": float(np.max(np.asarray(model.nodes, float)[:, 0])),
            "ymin": float(np.min(np.asarray(model.nodes, float)[:, 1])),
            "ymax": float(np.max(np.asarray(model.nodes, float)[:, 1])),
        },
    }
    return "\n".join(out) + "\n", meta


def _homogenization_extractor_source() -> str:
    """Abaqus-Python 2021-compatible extractor for ALLSE and macro reactions."""
    return r'''from __future__ import print_function
import json, os, sys, traceback
from odbAccess import openOdb

odb = None
out_path = 'case_results.json'
if len(sys.argv) >= 3:
    odb_path = sys.argv[-2]
    out_path = sys.argv[-1]
elif len(sys.argv) == 2:
    odb_path = sys.argv[-1]
else:
    odbs = [x for x in os.listdir('.') if x.lower().endswith('.odb')]
    if len(odbs) != 1:
        raise RuntimeError('Periodic homogenisation extractor requires an ODB path (or exactly one .odb)')
    odb_path = odbs[0]

def repo_get(repo, name):
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

def find_allse(step):
    found = []
    try:
        keys = list(step.historyRegions.keys())
    except Exception:
        keys = []
    for rkey in keys:
        try:
            region = step.historyRegions[rkey]
            hkeys = list(region.historyOutputs.keys())
        except Exception:
            continue
        for hkey in hkeys:
            if str(hkey).upper() == 'ALLSE':
                data = region.historyOutputs[hkey].data
                if data:
                    found.append((str(rkey), float(data[-1][1])))
    if not found:
        raise RuntimeError('ALLSE history output was not found. Available history regions: %s' % keys)
    # A whole-model request should yield one value. If Abaqus exposes more than one,
    # select the largest-magnitude total and preserve all candidates for diagnostics.
    found.sort(key=lambda x: abs(x[1]), reverse=True)
    return found[0][1], found

try:
    odb = openOdb(path=odb_path, readOnly=True)
    step = repo_get(odb.steps, 'HOM')
    if step is None:
        raise RuntimeError('HOM step was not found in ODB')
    if not step.frames:
        raise RuntimeError('HOM step contains no output frames')
    frame = step.frames[-1]
    assembly = odb.rootAssembly
    macro = repo_get(assembly.nodeSets, 'MACRO_ALL')
    if macro is None:
        for ikey in assembly.instances.keys():
            macro = repo_get(assembly.instances[ikey].nodeSets, 'MACRO_ALL')
            if macro is not None:
                break
    if macro is None:
        raise RuntimeError('MACRO_ALL node set was not found in the ODB')
    if 'RF' not in frame.fieldOutputs.keys():
        raise RuntimeError('RF field output is missing from final HOM frame')
    if 'U' not in frame.fieldOutputs.keys():
        raise RuntimeError('U field output is missing from final HOM frame')

    rf = frame.fieldOutputs['RF'].getSubset(region=macro)
    uu = frame.fieldOutputs['U'].getSubset(region=macro)
    reactions = {}
    displacements = {}
    for v in rf.values:
        data = v.data
        reactions[str(int(v.nodeLabel))] = float(data[0]) if len(data) else 0.0
    for v in uu.values:
        data = v.data
        displacements[str(int(v.nodeLabel))] = float(data[0]) if len(data) else 0.0
    allse, candidates = find_allse(step)

    result = {
        'status': 'COMPLETED',
        'odb': odb_path,
        'step': 'HOM',
        'strain_energy_ALLSE': float(allse),
        'allse_candidates': candidates,
        'macro_reaction_u1': reactions,
        'macro_displacement_u1': displacements,
        'frames': len(step.frames),
    }
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
'''


def reconstruct_ch_from_energies(case_results: list[dict], volume: float, amplitude: float) -> np.ndarray:
    """Recover all 21 independent terms from six pure + fifteen paired energies."""
    if len(case_results) != 21:
        raise ValueError(f"Expected 21 homogenisation case results, received {len(case_results)}.")
    V = float(volume); a = float(amplitude)
    if V <= 0 or a <= 0:
        raise ValueError("Volume and strain amplitude must be positive.")
    energies = {}
    for r in case_results:
        comps = tuple(int(x) for x in r["components"])
        energies[comps] = float(r["strain_energy"])
    C = np.zeros((6, 6), float)
    for i in range(6):
        W = energies[(i,)]
        C[i, i] = 2.0 * W / (V * a * a)
    for i in range(6):
        for j in range(i + 1, 6):
            W = energies[(i, j)]
            cij = (W - energies[(i,)] - energies[(j,)]) / (V * a * a)
            C[i, j] = C[j, i] = cij
    return C


def reconstruct_ch_from_macro_reactions(case_results: list[dict], volume: float, amplitude: float, reference_length: float) -> tuple[np.ndarray, float]:
    """Recover C columns from generalized macro reactions in the six pure cases.

    Abaqus reaction sign is normalised case-by-case so the active diagonal stiffness
    is positive. The unsymmetrised matrix and its reciprocity error are returned.
    """
    pure = {tuple(r["components"]): r for r in case_results if len(r["components"]) == 1}
    if len(pure) != 6:
        raise ValueError("Six pure strain cases are required for the macro-reaction cross-check.")
    C = np.zeros((6, 6), float)
    scale = float(reference_length) / (float(volume) * float(amplitude))
    for j in range(6):
        r = pure[(j,)]
        reactions = np.asarray(r.get("macro_reactions", [math.nan] * 6), float)
        if reactions.shape != (6,) or not np.all(np.isfinite(reactions)):
            raise ValueError(f"Macro reactions are incomplete for pure case {VOIGT_LABELS[j]}.")
        sign = 1.0 if reactions[j] >= 0 else -1.0
        C[:, j] = sign * reactions * scale
    denom = max(float(np.linalg.norm(C)), np.finfo(float).eps)
    symmetry_error = float(np.linalg.norm(C - C.T) / denom)
    return C, symmetry_error


def reconstruct_c_direct5(case_results: list[dict], volume: float, amplitude: float,
                          reference_length: float) -> tuple[np.ndarray, dict]:
    """3D stiffness at eps33 = 0 from the five periodic states E11, E22, G23, G13, G12.

    Column j of C is the macro-reaction vector of state j (one global sign, fixed by
    the energy-work balance).  Returns a 6x6 matrix whose 5x5 block over
    [11, 22, 23, 13, 12] is symmetrised; row and column 33 are NaN (no E33 state).
    """
    V, a, L = float(volume), float(amplitude), float(reference_length)
    if min(V, a, L) <= 0 or not all(map(math.isfinite, (V, a, L))):
        raise ValueError("Volume, amplitude and reference length must be positive and finite")
    keyed = {}
    for row in case_results:
        comps = tuple(row.get('components', ()))
        if len(comps) != 1 or comps[0] not in PBC5_COMPONENTS or comps in keyed:
            raise ValueError('Five unique pure-strain states (E11, E22, G23, G13, G12) are required')
        rf = np.asarray(row.get('macro_reactions', []), dtype=float)
        if rf.shape != (6,) or not np.all(np.isfinite(rf)):
            raise ValueError(f'Missing/non-finite six-component macro reactions: {comps}')
        W = float(row.get('strain_energy', float('nan')))
        if not math.isfinite(W) or W <= 0:
            raise ValueError(f'Non-positive or missing strain energy: {comps}')
        keyed[comps] = (rf, W)
    if len(keyed) != len(PBC5_COMPONENTS):
        raise ValueError(f'Expected five unit-strain results, got {len(keyed)}')
    idx = list(PBC5_COMPONENTS)
    R = np.column_stack([keyed[(j,)][0] for j in idx])          # 6 x 5 reactions
    energies = np.asarray([keyed[(j,)][1] for j in idx], float)
    raw_works = np.array([R[j, k] for k, j in enumerate(idx)]) * L * a
    sign = 1.0 if np.sum(raw_works) >= 0 else -1.0
    works = sign * raw_works
    if np.any(works <= 0):
        raise ValueError('Inconsistent macro-reaction signs across the five PBC load cases')
    work_errors = np.abs(works - 2.0 * energies) / np.maximum(2.0 * energies, np.finfo(float).tiny)
    C = np.full((6, 6), np.nan)
    blk = sign * R[idx, :] * L / (V * a)                          # 5 x 5 over [11,22,23,13,12]
    symmetry_error = float(np.linalg.norm(blk - blk.T) / max(np.linalg.norm(blk), np.finfo(float).eps))
    C[np.ix_(idx, idx)] = 0.5 * (blk + blk.T)
    return C, {
        'reaction_sign': int(sign),
        'work_error_by_case': work_errors.tolist(),
        'max_energy_reaction_work_error': float(np.max(work_errors)),
        'macro_reaction_symmetry_error': symmetry_error,
        'states': [VOIGT_LABELS[j] for j in idx],
        'note': 'No E33 state: U3+ - U3- = 0 between top and bottom in every state; C33 row/column not identified.',
    }


def solid_terms(C) -> dict:
    """C1111, C2222, C1122, C1212, C1313, C2323 (MPa) of a 6x6 Voigt [11,22,33,23,13,12] matrix."""
    C = np.asarray(C, float)
    return {"C1111": float(C[0, 0]), "C2222": float(C[1, 1]), "C1122": float(C[0, 1]),
            "C1212": float(C[5, 5]), "C1313": float(C[4, 4]), "C2323": float(C[3, 3])}


def reconstruct_ch_direct_from_reactions(case_results: list[dict], volume: float,
                                         amplitude: float, reference_length: float) -> tuple[np.ndarray, dict]:
    """Direct PBC homogenisation from six independent strain/reaction states.

    Macro node U1 = Lref * engineering strain; conjugate RF has units N.
    The sign is established ONCE using the recorded strain energies. It must
    agree for all six columns; per-column sign correction would conceal an
    inconsistent Abaqus model. Returns the *unsymmetrised* six-column tensor
    and quantitative work/reciprocity audits.
    """
    V, a, L = float(volume), float(amplitude), float(reference_length)
    if min(V, a, L) <= 0 or not all(map(math.isfinite, (V, a, L))):
        raise ValueError("Volume, amplitude and reference length must be positive and finite")
    if len(case_results) != 6:
        raise ValueError(f"Expected exactly six unit-strain results, got {len(case_results)}")
    keyed = {}
    for row in case_results:
        comps = tuple(row.get('components', ()))
        if len(comps) != 1 or comps[0] not in range(6) or comps in keyed:
            raise ValueError('Six unique pure-strain columns are required')
        rf = np.asarray(row.get('macro_reactions', []), dtype=float)
        if rf.shape != (6,) or not np.all(np.isfinite(rf)):
            raise ValueError(f'Missing/non-finite six-component macro reactions: {comps}')
        W = float(row.get('strain_energy', float('nan')))
        if not math.isfinite(W) or W <= 0:
            raise ValueError(f'Non-positive or missing strain energy: {comps}')
        keyed[comps] = (rf, W)
    if len(keyed) != 6:
        raise ValueError('Incomplete pure-strain cases')
    reactions = np.column_stack([keyed[(j,)][0] for j in range(6)])
    energies = np.asarray([keyed[(j,)][1] for j in range(6)], float)
    raw_works = np.diag(reactions) * L * a
    # A uniform global sign can account for the RF sign convention; independent
    # sign changes for each loading would be physically unjustified.
    sign = 1.0 if np.sum(raw_works) >= 0 else -1.0
    works = sign * raw_works
    work_errors = np.abs(works - 2.0*energies) / np.maximum(2.0*energies, np.finfo(float).tiny)
    if np.any(works <= 0):
        raise ValueError('Inconsistent macro-reaction signs across the six PBC load cases')
    C = sign * reactions * L / (V*a)
    symmetry_error = float(np.linalg.norm(C-C.T) / max(np.linalg.norm(C), np.finfo(float).eps))
    return C, {
        'reaction_sign': int(sign),
        'work_error_by_case': work_errors.tolist(),
        'max_energy_reaction_work_error': float(np.max(work_errors)),
        'macro_reaction_symmetry_error': symmetry_error,
        'direct_tensor_is_symmetrized_for_output': True,
        'sign_note': 'Single global RF sign chosen from the six energy-work checks',
    }


def rve_area_volume_diagnostics(model: RVEModel) -> dict:
    """Return the lattice-area/volume quantities used for energy homogenisation.

    The Chapter-3 hexagonal formulation uses V_RVE = p*q*h.  For a general
    oblique periodic cell, p is taken along a1 and q is the perpendicular height
    of a2 relative to a1, so p*q == |det([a1,a2])| exactly.  r is the signed
    projection of a2 along a1.  This makes the area definition explicit instead
    of accidentally using a drawing bounding box or solid-wall area.
    """
    a1 = np.asarray(model.a1, float).reshape(2)
    a2 = np.asarray(model.a2, float).reshape(2)
    p = float(np.linalg.norm(a1))
    if p <= np.finfo(float).eps:
        raise ValueError("a1 has zero length; RVE area is undefined.")
    e1 = a1 / p
    r = float(np.dot(a2, e1))
    area = abs(float(np.linalg.det(np.column_stack((a1, a2)))))
    q = area / p
    h = float(model.height)
    return {
        "a1_length_p": p,
        "a2_parallel_offset_r": r,
        "a2_perpendicular_height_q": q,
        "lattice_area_det": area,
        "chapter3_area_p_times_q": p * q,
        "area_consistency_relative_error": abs(area - p * q) / max(area, np.finfo(float).eps),
        "height": h,
        "rve_volume": area * h,
        "area_definition": "|det([a1,a2])| = p*q; volume = area*height",
    }


def plane_stress_shell_reduction(C: np.ndarray) -> tuple[np.ndarray, dict]:
    """Eliminate epsilon_33 under sigma_33=0 from a general 6x6 engineering C.

    Input order is [E11,E22,E33,G23,G13,G12].  Returned 5x5 order is
    [E11,E22,G12,G13,G23].  The Schur complement is valid even when the
    general anisotropic normal/shear coupling terms are non-zero.
    """
    C = np.asarray(C, float).reshape(6, 6)
    C = 0.5 * (C + C.T)
    if abs(float(C[2, 2])) <= np.finfo(float).eps:
        raise ValueError("C33 is zero; plane-stress reduction cannot eliminate epsilon_33.")
    keep = [0, 1, 5, 4, 3]  # 11,22,12,13,23
    cr = C[np.ix_(keep, [2])]
    Q = C[np.ix_(keep, keep)] - (cr @ cr.T) / float(C[2, 2])
    Q = 0.5 * (Q + Q.T)
    legacy = {
        "Q11": float(Q[0, 0]),
        "Q22": float(Q[1, 1]),
        "Q12": float(Q[0, 1]),
        "Q66": float(Q[2, 2]),  # 12 shear
        "Q44": float(Q[3, 3]),  # 13 shear in thesis notation
        "Q55": float(Q[4, 4]),  # 23 shear in thesis notation
    }
    return Q, legacy


# FE values reported in thesis Tables 3.5-3.10 (prescribed-displacement FE).  Used only by the
# Validation & Benchmarks page (paper comparison); FE verification does not use them (1.13).
_THESIS_FE_Q = {
    "Q11": [322.63,591.66,477.88,398.85,399.82,226.96,380.83,1785.94,1760.67,1662.92,445.37,435.58,1715.09,356.44],
    "Q22": [322.61,591.64,477.85,398.84,399.84,226.96,380.92,1785.82,1760.67,1662.83,445.36,435.60,1714.91,356.38],
    "Q12": [319.03,584.02,474.61,392.19,395.67,223.49,372.85,1772.24,1732.95,1653.49,438.49,430.08,1694.91,351.42],
    "Q66": [2.47,5.37,2.25,4.75,2.93,2.39,5.70,9.63,20.25,6.72,4.83,3.90,14.47,3.46],
    "Q44": [335.65,149.54,149.53,229.92,229.91,385.91,756.43,170.41,170.47,170.47,249.05,398.00,170.47,543.17],
    "Q55": [337.46,149.44,149.44,228.69,228.02,385.24,735.73,169.84,170.36,170.36,246.76,371.94,170.36,518.36],
}


def thesis_fe_q_reference(model: RVEModel) -> dict | None:
    """Return Table 3.5-3.10 FE Q values for built-in 'Hexagonal config N'."""
    import re
    m = re.search(r"config(?:uration)?\s*(\d+)", str(model.name), re.I)
    if not m:
        return None
    n = int(m.group(1))
    if not (1 <= n <= 14):
        return None
    return {k: float(v[n - 1]) for k, v in _THESIS_FE_Q.items()}


def compliance_and_engineering_constants(C: np.ndarray) -> tuple[np.ndarray, dict]:
    C = np.asarray(C, float).reshape(6, 6)
    S = np.linalg.inv(C)
    S = 0.5 * (S + S.T)
    vals = {
        "E1": 1.0 / S[0, 0], "E2": 1.0 / S[1, 1], "E3": 1.0 / S[2, 2],
        "G23": 1.0 / S[3, 3], "G13": 1.0 / S[4, 4], "G12": 1.0 / S[5, 5],
        "nu12": -S[1, 0] / S[0, 0], "nu21": -S[0, 1] / S[1, 1],
        "nu13": -S[2, 0] / S[0, 0], "nu31": -S[0, 2] / S[2, 2],
        "nu23": -S[2, 1] / S[1, 1], "nu32": -S[1, 2] / S[2, 2],
    }
    return S, {k: float(v) for k, v in vals.items()}


def independent_terms(C: np.ndarray, significance_tolerance: float = 1e-4) -> list[dict]:
    C = np.asarray(C, float).reshape(6, 6)
    rows = []
    for i in range(6):
        for j in range(i, 6):
            denom = math.sqrt(abs(float(C[i, i] * C[j, j]))) if C[i, i] * C[j, j] != 0 else math.nan
            eta = abs(float(C[i, j])) / denom if math.isfinite(denom) and denom > 0 else math.nan
            rows.append({
                "i": i + 1, "j": j + 1, "term": f"C{i + 1}{j + 1}",
                "value": float(C[i, j]), "relative_coupling": float(eta) if math.isfinite(eta) else math.nan,
                "significant": bool(i == j or (math.isfinite(eta) and eta >= significance_tolerance)),
            })
    return rows



def _write_periodic_domain_svg(path: Path, model: RVEModel, domain_diag: dict) -> None:
    """Write a lightweight 2-D preview of the exact wall network sent to Abaqus."""
    pts=np.asarray(model.nodes,float)
    if pts.size == 0:
        return
    xmin,xmax=float(np.min(pts[:,0])),float(np.max(pts[:,0]))
    ymin,ymax=float(np.min(pts[:,1])),float(np.max(pts[:,1]))
    w=max(xmax-xmin,1.0); h=max(ymax-ymin,1.0)
    pad=0.04*max(w,h)
    xmin-=pad; xmax+=pad; ymin-=pad; ymax+=pad
    width=900.0; height=max(420.0,min(1100.0,width*(ymax-ymin)/(xmax-xmin)))
    def xy(p):
        x=(float(p[0])-xmin)/(xmax-xmin)*width
        y=height-(float(p[1])-ymin)/(ymax-ymin)*height
        return x,y
    lines=[
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height:.0f}" viewBox="0 0 {width:.0f} {height:.0f}">',
        '<rect width="100%" height="100%" fill="white"/>',
    ]
    # For exact lattice supercells, draw the primitive parallelogram tiling first.
    if domain_diag.get("domain_mode") == PBC_DOMAIN_LATTICE_SUPERCELL:
        pa1=np.asarray(domain_diag.get("primitive_a1",[0.,0.]),float)
        pa2=np.asarray(domain_diag.get("primitive_a2",[0.,0.]),float)
        na=int(domain_diag.get("repeat_a1",domain_diag.get("repeat_x",1)))
        nb=int(domain_diag.get("repeat_a2",domain_diag.get("repeat_y",1)))
        # Internal primitive-cell boundaries are visual guides only.
        for ia in range(na+1):
            p0=ia*pa1; p1=ia*pa1+nb*pa2
            x1,y1=xy(p0); x2,y2=xy(p1)
            lines.append(f'<line x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}" stroke="#9aa4b2" stroke-width="1" stroke-dasharray="5,5"/>')
        for ib in range(nb+1):
            p0=ib*pa2; p1=na*pa1+ib*pa2
            x1,y1=xy(p0); x2,y2=xy(p1)
            lines.append(f'<line x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}" stroke="#9aa4b2" stroke-width="1" stroke-dasharray="5,5"/>')
    for e in model.elements:
        x1,y1=xy(pts[e.node_i]); x2,y2=xy(pts[e.node_j])
        lines.append(f'<line x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}" stroke="black" stroke-width="1.6"/>')
    # Highlight the actual outer periodic boundary.
    if domain_diag.get("domain_mode") in {PBC_DOMAIN_LATTICE_SUPERCELL, PBC_DOMAIN_RECTANGULAR_SUPERCELL, PBC_DOMAIN_PUBLISHED}:
        b1=np.asarray(model.a1,float); b2=np.asarray(model.a2,float)
        boundary=[np.zeros(2),b1,b1+b2,b2,np.zeros(2)]
        for p0,p1 in zip(boundary[:-1],boundary[1:]):
            x1,y1=xy(p0); x2,y2=xy(p1)
            lines.append(f'<line x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}" stroke="#d62728" stroke-width="2.4"/>')
    for pxy in pts:
        x,y=xy(pxy); lines.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="2.2" fill="black"/>')
    label=(f"{domain_diag.get('domain_mode','domain')} | primitive cells: {domain_diag.get('primitive_cell_count',1)} | "
           f"geometry nodes: {len(model.nodes)} | walls: {len(model.elements)}")
    lines.append(f'<text x="12" y="22" font-family="Arial" font-size="14">{label}</text>')
    lines.append('</svg>')
    path.write_text('\n'.join(lines),encoding='utf-8')


def _assert_deck_uses_resolved_domain(meta: dict, analysis_model: RVEModel, domain_diag: dict) -> None:
    """Guard against accidentally handing the primitive geometry to Abaqus."""
    if int(meta.get('analysis_geometry_node_count',-1)) != len(analysis_model.nodes):
        raise RuntimeError('Internal Abaqus handoff error: deck node geometry does not match the resolved periodic domain.')
    if int(meta.get('analysis_geometry_wall_count',-1)) != len(analysis_model.elements):
        raise RuntimeError('Internal Abaqus handoff error: deck wall geometry does not match the resolved periodic domain.')
    if not np.allclose(np.asarray(meta.get('analysis_a1'),float),np.asarray(analysis_model.a1,float),rtol=0,atol=1e-10):
        raise RuntimeError('Internal Abaqus handoff error: deck a1 differs from the resolved periodic domain.')
    if not np.allclose(np.asarray(meta.get('analysis_a2'),float),np.asarray(analysis_model.a2,float),rtol=0,atol=1e-10):
        raise RuntimeError('Internal Abaqus handoff error: deck a2 differs from the resolved periodic domain.')
    if domain_diag.get('domain_mode') in {PBC_DOMAIN_LATTICE_SUPERCELL, PBC_DOMAIN_RECTANGULAR_SUPERCELL, PBC_DOMAIN_PUBLISHED}:
        if int(domain_diag.get('primitive_cell_count',1)) <= 1:
            raise RuntimeError('Repeated-supercell mode did not create more than one primitive cell.')
        area=abs(float(np.linalg.det(np.column_stack((analysis_model.a1,analysis_model.a2)))))
        expected=float(domain_diag.get('supercell_area',area))
        if not math.isclose(area,expected,rel_tol=1e-10,abs_tol=1e-9):
            raise RuntimeError('Internal Abaqus handoff error: resolved supercell area differs from the requested supercell area.')

def _periodic_input_signature(
    model: RVEModel,
    materials: dict[str, Material],
    laminates: list[Laminate],
    settings: PeriodicHomogenizationSettings,
) -> str:
    """Stable signature used to protect resumable batch folders from stale inputs."""
    payload = {
        "model": {
            "nodes": np.asarray(model.nodes, float).tolist(),
            "elements": [asdict(e) for e in model.elements],
            "a1": np.asarray(model.a1, float).tolist(),
            "a2": np.asarray(model.a2, float).tolist(),
            "height": float(model.height),
            "wall_groups": list(model.wall_groups),
        },
        "materials": {str(k): asdict(v) for k, v in sorted(materials.items(), key=lambda kv: str(kv[0]))},
        "laminates": [{"name": l.name, "plies": [asdict(p) for p in l.plies]} for l in laminates],
        "settings": asdict(settings),
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=True).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _periodic_files_from_directory(run_dir: Path, snapshot: dict) -> PeriodicHomogenizationFiles:
    """Reconstruct the file bundle for a prepared periodic FE run."""
    cases: list[HomogenizationCaseFiles] = []
    for case in snapshot.get("cases", []):
        case_name = str(case.get("directory") or f"{int(case['index']):02d}_{case['name']}")
        cdir = run_dir / case_name
        job = _safe_name(f"HOM_{case_name}", "HOM")
        cases.append(HomogenizationCaseFiles(
            int(case["index"]), str(case["name"]), list(case["strain"]), cdir, job,
            cdir / f"{job}.inp", cdir / "extract_hom_odb.py", cdir / "case_results.json", cdir / "abaqus_case.log",
        ))
    return PeriodicHomogenizationFiles(
        run_dir, run_dir / "homogenization_snapshot.json", run_dir / "homogenized_results.json",
        run_dir / "homogenization_run.log", cases,
    )


def _completed_json(path: Path) -> dict | None:
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
        return raw if str(raw.get("status", "")).upper() == "COMPLETED" else None
    except Exception:
        return None


def prepare_periodic_homogenization(
    run_root: str | Path,
    model: RVEModel,
    materials: dict[str, Material],
    laminates: list[Laminate],
    settings: PeriodicHomogenizationSettings,
    *,
    run_name: str | None = None,
    resume: bool = False,
) -> PeriodicHomogenizationFiles:
    settings.validate(); _validate_periodic_fe_model(model)
    analysis_model, domain_diag = resolve_periodic_domain(model, settings)
    _validate_periodic_fe_model(analysis_model)
    root = Path(run_root).expanduser().resolve(); root.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    _tag = {CH_METHOD_DIRECT: 'CH6', CH_METHOD_DIRECT5: 'C5'}.get(settings.identification_method, 'CH21')
    batch_name = _safe_name(run_name if run_name else f"{model.name}_PBC_{_tag}_{stamp}", "RVE_HOM")
    run_dir = root / batch_name
    desired_signature = _periodic_input_signature(model, materials, laminates, settings)
    snapshot = run_dir / "homogenization_snapshot.json"
    summary = run_dir / "homogenized_results.json"
    log_file = run_dir / "homogenization_run.log"
    if run_dir.exists():
        if not resume:
            raise FileExistsError(f"Abaqus homogenisation folder already exists: {run_dir}")
        if not snapshot.exists():
            raise RuntimeError(f"Cannot resume {run_dir}: homogenization_snapshot.json is missing.")
        old = json.loads(snapshot.read_text(encoding="utf-8"))
        old_sig = str(old.get("input_signature", ""))
        if old_sig and old_sig != desired_signature:
            raise RuntimeError(
                f"Cannot resume {run_dir}: the RVE/material/laminate/settings signature has changed. "
                "Use a new batch folder or remove the stale case folder."
            )
    else:
        run_dir.mkdir(parents=False, exist_ok=False)
    extractor = _homogenization_extractor_source()
    cases: list[HomogenizationCaseFiles] = []
    case_metadata = []
    # Reuse the same audited S4 mesh, PBC equations, shell sections and material
    # cards for every macro-strain state; only control-node displacements change.
    # This avoids repeated construction of the same fine 2x2 supercell mesh.
    template_text = None
    template_meta = None
    selected_cases = homogenization_cases(settings.strain_amplitude)
    if settings.identification_method == CH_METHOD_DIRECT:
        selected_cases = selected_cases[:6]
    elif settings.identification_method == CH_METHOD_DIRECT5:
        selected_cases = [c for c in selected_cases[:6] if c["components"][0] in PBC5_COMPONENTS]
    for case in selected_cases:
        case_name = f"{case['index']:02d}_{case['name']}"
        cdir = run_dir / case_name; cdir.mkdir(exist_ok=True)
        job = _safe_name(f"HOM_{case_name}", "HOM")
        inp = cdir / f"{job}.inp"; ex = cdir / "extract_hom_odb.py"; res = cdir / "case_results.json"; clog = cdir / "abaqus_case.log"
        if template_text is None:
            text, meta = generate_periodic_homogenization_input(
                analysis_model, materials, laminates, settings, case["strain"], job_name=job,
                resolved_domain=True, domain_diagnostics=domain_diag,
            )
            template_text = text
            template_meta = meta
            original_job = job
        else:
            import re as _re
            text = template_text.replace(f"** Job: {original_job}", f"** Job: {job}", 1)
            for k in range(6):
                new_value = float(case['strain'][k]) * float(settings.reference_length)
                old_line = _re.compile(rf"(?m)^MACRO_{k + 1}, 1, 1, [^\r\n]*$")
                replacement = f"MACRO_{k + 1}, 1, 1, {new_value:.16g}"
                text, count = old_line.subn(replacement, text)
                if count != 1:
                    raise RuntimeError(f"Expected one boundary line for MACRO_{k + 1}; found {count}")
            meta = {**template_meta, 'strain': np.asarray(case['strain'], float).tolist()}
        _assert_deck_uses_resolved_domain(meta, analysis_model, domain_diag)
        if not (resume and inp.exists()): inp.write_text(text, encoding="utf-8")
        if not (resume and ex.exists()): ex.write_text(extractor, encoding="utf-8")
        cases.append(HomogenizationCaseFiles(case["index"], case["name"], np.asarray(case["strain"], float).tolist(), cdir, job, inp, ex, res, clog))
        case_metadata.append({**case, "strain": np.asarray(case["strain"], float).tolist(), "deck_meta": meta, "directory": case_name})
    area_diag = rve_area_volume_diagnostics(analysis_model)
    area = float(area_diag["lattice_area_det"])
    snap = {
        "mode": {CH_METHOD_DIRECT: "general_anisotropic_periodic_direct_CH",
                 CH_METHOD_DIRECT5: "periodic_3d_terms_eps33_zero"}.get(settings.identification_method,
                                                                         "general_anisotropic_periodic_energy_homogenization"),
        "input_signature": desired_signature,
        "strain_convention": "engineering [E11,E22,E33,G23,G13,G12]; symmetric 3D small strain",
        "identification": {CH_METHOD_DIRECT: "6 pure engineering-strain cases; all 6 macro reaction forces per case form the C^H columns",
                           CH_METHOD_DIRECT5: "5 pure engineering-strain cases (no E33: U3+ - U3- = 0 top to bottom); macro reactions form the C columns"}.get(
                               settings.identification_method, "21 energy cases (6 pure + 15 paired); no constitutive symmetry beyond major elasticity symmetry is imposed"),
        "model": {"name": model.name, "nodes": model.nodes.tolist(), "elements": [asdict(e) for e in model.elements], "a1": model.a1.tolist(), "a2": model.a2.tolist(), "height": model.height, "wall_groups": list(model.wall_groups)},
        "analysis_model": {"name": analysis_model.name, "nodes": analysis_model.nodes.tolist(), "elements": [asdict(e) for e in analysis_model.elements], "a1": analysis_model.a1.tolist(), "a2": analysis_model.a2.tolist(), "height": analysis_model.height, "wall_groups": list(analysis_model.wall_groups)},
        "periodic_domain": domain_diag,
        "materials": {k: asdict(v) for k, v in materials.items()},
        "laminates": [{"name": l.name, "plies": [asdict(p) for p in l.plies]} for l in laminates],
        "settings": asdict(settings),
        "rve_area": area, "rve_volume": area * float(analysis_model.height),
        "rve_geometry_diagnostics": area_diag,
        "cases": case_metadata,
        "theory_note": ("Direct 6-column PBC macro-reaction homogenisation; elastic engineering constants are derived from inverse CH, not used to construct CH." if settings.identification_method == CH_METHOD_DIRECT else "Full 21-case PBC energy identification of a general symmetric 6x6 stiffness matrix."),
    }
    snapshot.write_text(json.dumps(snap, indent=2), encoding="utf-8")
    (run_dir / "periodic_domain.json").write_text(json.dumps(domain_diag, indent=2), encoding="utf-8")
    if case_metadata:
        (run_dir / "pbc_audit.json").write_text(json.dumps(case_metadata[0].get("deck_meta", {}).get("pbc_audit", {}), indent=2), encoding="utf-8")
    _write_periodic_domain_svg(run_dir / "periodic_domain_preview.svg", analysis_model, domain_diag)
    if case_metadata:
        first_meta=case_metadata[0].get("deck_meta", {})
        (run_dir / "deck_geometry_check.json").write_text(json.dumps({
            "domain_mode": domain_diag.get("domain_mode"),
            "primitive_cell_count": domain_diag.get("primitive_cell_count", 1),
            "geometry_nodes_sent_to_deck": first_meta.get("analysis_geometry_node_count"),
            "geometry_walls_sent_to_deck": first_meta.get("analysis_geometry_wall_count"),
            "mesh_nodes_sent_to_deck": first_meta.get("mesh_node_count"),
            "mesh_elements_sent_to_deck": first_meta.get("mesh_element_count"),
            "geometry_bounds": first_meta.get("analysis_geometry_bounds"),
            "analysis_a1": first_meta.get("analysis_a1"),
            "analysis_a2": first_meta.get("analysis_a2"),
            "supercell_area": domain_diag.get("supercell_area"),
        }, indent=2), encoding="utf-8")
    if domain_diag.get("domain_mode") in {PBC_DOMAIN_LATTICE_SUPERCELL, PBC_DOMAIN_RECTANGULAR_SUPERCELL, PBC_DOMAIN_PUBLISHED}:
        with (run_dir / "supercell_nodes.csv").open("w", encoding="utf-8") as f:
            f.write("node_id,x_mm,y_mm\n")
            for i, xy in enumerate(np.asarray(analysis_model.nodes, float), 1):
                f.write(f"{i},{xy[0]:.16g},{xy[1]:.16g}\n")
        with (run_dir / "supercell_walls.csv").open("w", encoding="utf-8") as f:
            f.write("wall_id,node_i,node_j,laminate_index,group\n")
            for i, e in enumerate(analysis_model.elements, 1):
                group = analysis_model.wall_groups[i-1] if i-1 < len(analysis_model.wall_groups) else ""
                f.write(f"{i},{e.node_i+1},{e.node_j+1},{e.laminate},\"{str(group).replace(chr(34), chr(34)*2)}\"\n")
    return PeriodicHomogenizationFiles(run_dir, snapshot, summary, log_file, cases)


def _read_completed_case(case: HomogenizationCaseFiles, macro_ids: list[int] | None = None) -> dict:
    raw = json.loads(case.results_file.read_text(encoding="utf-8"))
    if str(raw.get("status", "")).upper() != "COMPLETED":
        raise RuntimeError(raw.get("traceback") or raw.get("error") or f"Case {case.name} ODB extraction failed")
    result = {
        "index": case.index,
        "name": case.name,
        "strain": list(case.strain),
        "components": [i for i, x in enumerate(case.strain) if abs(float(x)) > 0],
        "strain_energy": float(raw["strain_energy_ALLSE"]),
        "odb": raw.get("odb", ""),
    }
    # Macro nodes are numbered consecutively at the end of the shell mesh; infer order from displacement magnitudes/labels
    # via the snapshot deck metadata in the caller when available. Here preserve the keyed mapping too.
    result["macro_reactions_by_node"] = raw.get("macro_reaction_u1", {})
    result["macro_displacements_by_node"] = raw.get("macro_displacement_u1", {})
    if macro_ids:
        result["macro_reactions"] = [float(raw.get("macro_reaction_u1", {}).get(str(int(n)), math.nan)) for n in macro_ids]
    return result


def run_periodic_homogenization(
    files: PeriodicHomogenizationFiles,
    abaqus_command: str,
    settings: PeriodicHomogenizationSettings,
    progress: Callable[[int, int, str], None] | None = None,
    *,
    resume: bool = False,
) -> dict:
    """Run six reaction cases (default) or 21 energy cases to identify C^H.

    With ``resume=True``, any load case with a completed ``case_results.json`` is
    reused. A completed homogenized summary is returned immediately. This makes
    long FE batches restartable at the individual macro-strain-state level.
    """
    settings.validate()
    exe = str(abaqus_command).strip()
    if not exe:
        raise ValueError("Abaqus executable/command is empty.")
    if resume:
        ready = _completed_json(files.summary_file)
        if ready is not None:
            if progress:
                progress(len(files.cases), len(files.cases), "already completed")
            return ready
    snap = json.loads(files.snapshot_file.read_text(encoding="utf-8"))
    V = float(snap["rve_volume"])
    deck_cases = snap["cases"]
    completed: list[dict] = []
    with files.log_file.open("a" if resume else "w", encoding="utf-8", errors="replace") as master_log:
        master_log.write("RVE Studio full periodic Abaqus homogenisation\n")
        master_log.write(f"Cases: {len(files.cases)}\nCommand: {exe}\n\n")
        for idx, case in enumerate(files.cases, 1):
            macro_ids = list(deck_cases[idx - 1]["deck_meta"]["macro_node_ids"])
            if resume and _completed_json(case.results_file) is not None:
                if progress:
                    progress(idx, len(files.cases), case.name + " (resumed)")
                master_log.write(f"[{idx:02d}/{len(files.cases):02d}] {case.name} -- reused completed result\n"); master_log.flush()
                completed.append(_read_completed_case(case, macro_ids))
                continue
            if progress:
                progress(idx, len(files.cases), case.name)
            master_log.write(f"[{idx:02d}/{len(files.cases):02d}] {case.name}\n"); master_log.flush()
            with case.log_file.open("a" if resume else "w", encoding="utf-8", errors="replace") as log:
                rc = _run_shell_command([exe, f"job={case.job_name}", f"input={case.input_file.name}", f"cpus={int(settings.cpus)}", "interactive"], case.case_dir, log)
                if rc != 0:
                    raise RuntimeError(f"Abaqus case {case.name} returned exit code {rc}. See {case.log_file}")
                odb = case.case_dir / f"{case.job_name}.odb"
                if not odb.exists():
                    raise RuntimeError(f"Abaqus case {case.name} produced no ODB. See {case.log_file}")
                log.write("\n--- ODB homogenisation extraction ---\n"); log.flush()
                rc2 = _run_shell_command([exe, "python", case.extractor_file.name, odb.name, case.results_file.name], case.case_dir, log)
                if not case.results_file.exists():
                    raise RuntimeError(f"Case {case.name} ODB extraction produced no results file (wrapper code {rc2}). See {case.log_file}")
            row = _read_completed_case(case, macro_ids)
            completed.append(row)
            if not settings.save_odb:
                try:
                    (case.case_dir / f"{case.job_name}.odb").unlink()
                except Exception:
                    pass

    method = str(snap.get('settings', {}).get('identification_method', CH_METHOD_ENERGY))
    if method == CH_METHOD_DIRECT5:
        C5, audit5 = reconstruct_c_direct5(completed, V, settings.strain_amplitude, settings.reference_length)
        result = {
            "status": "COMPLETED", "mode": "Periodic 3D stiffness terms, eps33 = 0 (5 macro-reaction states)",
            "identification_method": method, "CH": [[None if not math.isfinite(x) else float(x) for x in row] for row in C5],
            "C_terms": solid_terms(C5), "case_count": len(completed), "direct_reaction_audit": audit5,
            "macro_reaction_symmetry_error": audit5["macro_reaction_symmetry_error"],
            "energy_vs_reaction_relative_error": audit5["max_energy_reaction_work_error"],
            "strain_convention": "engineering [E11,E22,E33,G23,G13,G12]; E33 not applied",
            "boundary_condition_note": "Periodic side faces; top and bottom tied with U3+ - U3- = 0 (no E33 state); periodic shell rotations; anchor translations fixed.",
            "rve_area": float(snap.get("rve_area", math.nan)), "rve_volume": V,
            "periodic_domain": snap.get("periodic_domain", {}),
            "pbc_audit": (deck_cases[0].get("deck_meta", {}).get("pbc_audit", {}) if deck_cases else {}),
            "strain_amplitude": float(settings.strain_amplitude), "cases": completed, "run_dir": str(files.run_dir),
        }
        files.summary_file.write_text(json.dumps(result, indent=2), encoding="utf-8")
        with (files.run_dir / "C_terms_MPa.csv").open("w", encoding="utf-8") as f:
            f.write("term,value_MPa\n")
            for k, v in result["C_terms"].items():
                f.write(f"{k},{v:.16g}\n")
        return result
    direct_audit = {}
    C_raw_direct = None
    if method == CH_METHOD_DIRECT:
        C_raw_direct, direct_audit = reconstruct_ch_direct_from_reactions(
            completed, V, settings.strain_amplitude, settings.reference_length)
        # Preserve the raw reciprocity error before reporting major-symmetric CH.
        C_primary = 0.5 * (C_raw_direct + C_raw_direct.T)
        source = 'six_macro_reaction_columns'
    else:
        C_primary = reconstruct_ch_from_energies(completed, V, settings.strain_amplitude)
        source = 'twenty_one_strain_energies'
    eig = np.linalg.eigvalsh(C_primary)
    cond = float(np.linalg.cond(C_primary))
    compliance_error = None
    try:
        S, eng = compliance_and_engineering_constants(C_primary)
        compliance_method = "inverse"
    except Exception as exc:
        # Preserve C^H even if the homogenised cell has a mechanism or is numerically singular.
        S = np.linalg.pinv(C_primary)
        S = 0.5 * (S + S.T)
        eng = {}
        compliance_method = "pseudo-inverse"
        compliance_error = str(exc)
    terms = independent_terms(C_primary, settings.significance_tolerance)
    Q5, q_legacy = plane_stress_shell_reduction(C_primary)
    result = {
        "status": "COMPLETED",
        "mode": "Direct PBC C^H (6 macro-reaction columns)" if method == CH_METHOD_DIRECT else "PBC C^H (21 energy terms)",
        "identification_method": method,
        "CH": C_primary.tolist(),
        "CH_source": source,
        "case_count": len(completed),
        "strain_convention": "engineering [E11,E22,E33,G23,G13,G12]",
        **({"CH_energy": C_primary.tolist()} if method == CH_METHOD_ENERGY else {"CH_macro_reaction": C_primary.tolist(), "CH_macro_reaction_raw": C_raw_direct.tolist(), "direct_reaction_audit": direct_audit}),
        "compliance": S.tolist(),
        "compliance_method": compliance_method,
        "engineering_constants": eng,
        "independent_terms": terms,
        "Q_plane_stress_general": Q5.tolist(),
        "Q_plane_stress_labels": ["E11", "E22", "G12", "G13", "G23"],
        "Q_thesis_terms": q_legacy,
        "stiffness_eigenvalues": eig.tolist(),
        "stiffness_condition_number": cond,
        "positive_definite": bool(np.all(eig > 0.0)),
        "rve_area": float(snap.get("rve_area", math.nan)),
        "rve_volume": V,
        "rve_geometry_diagnostics": snap.get("rve_geometry_diagnostics", {}),
        "periodic_domain": snap.get("periodic_domain", {}),
        "pbc_audit": (deck_cases[0].get("deck_meta", {}).get("pbc_audit", {}) if deck_cases else {}),
        "boundary_condition_note": (
            "Periodic affine displacement PBC on a rectangular repeated supercell (b1,b2,a3); periodic shell rotations; only anchor translations U1-U3 fixed. "
            "Supercell area is an integer multiple of the primitive-cell area. This is not the prescribed-displacement FE boundary condition of thesis Table 3.4."
            if snap.get("periodic_domain", {}).get("domain_mode") in {PBC_DOMAIN_LATTICE_SUPERCELL, PBC_DOMAIN_RECTANGULAR_SUPERCELL, PBC_DOMAIN_PUBLISHED} else
            "Periodic affine displacement PBC on primitive a1,a2,a3; periodic shell rotations; only anchor translations U1-U3 fixed. This is not the prescribed-displacement FE boundary condition of thesis Table 3.4."
        ),
        "strain_amplitude": float(settings.strain_amplitude),
        "cases": completed,
        "run_dir": str(files.run_dir),
    }
    if compliance_error:
        result["compliance_warning"] = compliance_error
    if settings.stress_cross_check and method == CH_METHOD_ENERGY:
        try:
            C_rf, symerr = reconstruct_ch_from_macro_reactions(completed, V, settings.strain_amplitude, settings.reference_length)
            C_rf_sym = 0.5 * (C_rf + C_rf.T)
            denom = max(float(np.linalg.norm(C_primary)), np.finfo(float).eps)
            result["CH_macro_reaction_raw"] = C_rf.tolist()
            result["CH_macro_reaction_symmetric"] = C_rf_sym.tolist()
            result["macro_reaction_symmetry_error"] = symerr
            result["energy_vs_reaction_relative_error"] = float(np.linalg.norm(C_primary - C_rf_sym) / denom)
        except Exception as exc:
            result["macro_reaction_cross_check_error"] = str(exc)
    if method == CH_METHOD_DIRECT:
        result['macro_reaction_symmetry_error'] = direct_audit['macro_reaction_symmetry_error']
        result['energy_vs_reaction_relative_error'] = direct_audit['max_energy_reaction_work_error']
    files.summary_file.write_text(json.dumps(result, indent=2), encoding="utf-8")
    # Research-friendly flat exports in addition to the complete JSON record.
    np.savetxt(files.run_dir / "CH_FE_MPa.csv", C_primary, delimiter=",", header=",".join(VOIGT_LABELS), comments="")
    if method == CH_METHOD_ENERGY:
        np.savetxt(files.run_dir / "CH_energy_MPa.csv", C_primary, delimiter=",", header=",".join(VOIGT_LABELS), comments="")
    else:
        np.savetxt(files.run_dir / "CH_direct_MPa.csv", C_primary, delimiter=",", header=",".join(VOIGT_LABELS), comments="")
        np.savetxt(files.run_dir / "CH_reaction_unsymmetrised_MPa.csv", C_raw_direct, delimiter=",", header=",".join(VOIGT_LABELS), comments="")
    np.savetxt(files.run_dir / "compliance_1_per_MPa.csv", S, delimiter=",", header=",".join(VOIGT_LABELS), comments="")
    np.savetxt(files.run_dir / "Q_plane_stress_general_MPa.csv", Q5, delimiter=",", header=",".join(["E11","E22","G12","G13","G23"]), comments="")
    with (files.run_dir / "Q_thesis_terms_MPa.csv").open("w", encoding="utf-8") as f:
        f.write("term,value_MPa\n")
        for _k in ["Q11","Q22","Q12","Q66","Q44","Q55"]:
            f.write(f"{_k},{q_legacy[_k]:.16g}\n")
    with (files.run_dir / "independent_21_terms.csv").open("w", encoding="utf-8") as f:
        f.write("term,value_MPa,relative_coupling,significant\n")
        for t in terms:
            f.write(f"{t['term']},{t['value']:.16g},{t['relative_coupling']:.16g},{int(t['significant'])}\n")
    with (files.run_dir / "case_energies.csv").open("w", encoding="utf-8") as f:
        f.write("index,name,strain_energy_Nmm,strain_vector\n")
        for r in completed:
            f.write(f"{r['index']},{r['name']},{r['strain_energy']:.16g},\"{' '.join(format(float(x),'.16g') for x in r['strain'])}\"\n")
    return result
