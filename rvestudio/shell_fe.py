"""Built-in linear shell finite elements for checking the source-study DBC RVE.

This is a small, self-contained cross-check of the Abaqus decks written by
:mod:`rvestudio.abaqus_thesis_dbc`.  It meshes exactly the same RVE with exactly
the same node sets and prescribed displacements, and solves the six Table 3.4
cases with a flat four-node laminated shell:

* membrane: bilinear quadrilateral with Wilson incompatible modes (no in-plane
  shear locking; exact for rectangles, which is what the wall mesher produces);
* bending and transverse shear: Reissner-Mindlin plate with MITC4 assumed
  transverse shear strains (no shear locking), 2x2 Gauss integration;
* laminate: full A, B, D from classical lamination theory in the Abaqus wall
  orientation (local 1 along the wall, local 2 vertical, normal = 1 x 2,
  ply 1 at the bottom) and transverse shear 5/6 * sum(G t);
* drilling rotations: a negligible spring (1e-6 of the bending stiffness).

It is not a replacement for Abaqus S4; it lets the RVE, the boundary sets and the
energy-to-stiffness reconstruction be verified before a licence is needed, and
gives an independent number to compare the Abaqus results with.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

from .laminate import qbar_orthotropic
from .models import Laminate, Material, RVEModel
from .abaqus import _build_shell_mesh

_G = (-1.0 / math.sqrt(3.0), 1.0 / math.sqrt(3.0))
_XI = np.array([-1.0, 1.0, 1.0, -1.0])
_ETA = np.array([-1.0, -1.0, 1.0, 1.0])


def _factor(A):
    """Sparse LU of a symmetric positive-definite shell stiffness matrix.

    Symmetric mode (minimum-degree ordering of A + A^T, diagonal pivots) is 4-7 times faster than the general
    COLAMD ordering for these matrices and gives the same solution to round-off; the general ordering is kept
    as the fallback for a matrix that needs off-diagonal pivoting."""
    try:
        lu = splu(A, permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0, options=dict(SymmetricMode=True))
        d = np.abs(lu.U.diagonal())
        if np.all(np.isfinite(d)) and d.min() > 1e-13 * d.max():
            return lu
    except RuntimeError:
        pass
    return splu(A, permc_spec="COLAMD")


def shell_section(lam: Laminate, materials: dict[str, Material], *, thickness_factor: float = 1.0,
                  modulus_factor: float = 1.0, shear_correction: float = 5.0 / 6.0):
    """A, B, D (3x3) and transverse shear As (2x2) of a laminate in its own 1-2 frame."""
    t = np.array([p.thickness for p in lam.plies], float) * float(thickness_factor)
    tt = float(t.sum())
    z = np.concatenate(([-tt / 2], -tt / 2 + np.cumsum(t)))
    A = np.zeros((3, 3)); B = np.zeros((3, 3)); D = np.zeros((3, 3)); As = np.zeros((2, 2))
    for k, ply in enumerate(lam.plies):
        m = materials[ply.material]
        th = math.radians(float(ply.angle_deg))
        Q = qbar_orthotropic(m, th) * float(modulus_factor)
        z0, z1 = z[k], z[k + 1]
        A += Q * (z1 - z0)
        B += Q * (z1 ** 2 - z0 ** 2) / 2.0
        D += Q * (z1 ** 3 - z0 ** 3) / 3.0
        g = float(m.Gth if m.Gth is not None else m.G12) * float(modulus_factor)  # deck: G13 = G23 = Gth
        As += np.eye(2) * g * (z1 - z0)
    return A, B, D, As * float(shear_correction)


def element_stiffness(a: float, b: float, A: np.ndarray, B: np.ndarray, D: np.ndarray, As: np.ndarray,
                      drill: float | None = None) -> np.ndarray:
    """24x24 local stiffness of an a x b rectangle; DOFs per node [u, v, w, rx, ry, rz]."""
    n = 28
    K = np.zeros((n, n))
    Cg = np.block([[A, B], [B, D]])
    detJ = a * b / 4.0
    # MITC4 tying vectors (rectangle): gamma_xz on the edges eta = +-1, gamma_yz on xi = +-1.
    gA = np.zeros(n); gC = np.zeros(n); gB = np.zeros(n); gD = np.zeros(n)
    w = lambda i: 6 * i + 2
    rx = lambda i: 6 * i + 3
    ry = lambda i: 6 * i + 4
    gA[w(2)] += 1 / a; gA[w(3)] -= 1 / a; gA[ry(2)] += 0.5; gA[ry(3)] += 0.5
    gC[w(1)] += 1 / a; gC[w(0)] -= 1 / a; gC[ry(0)] += 0.5; gC[ry(1)] += 0.5
    gB[w(3)] += 1 / b; gB[w(0)] -= 1 / b; gB[rx(0)] -= 0.5; gB[rx(3)] -= 0.5
    gD[w(2)] += 1 / b; gD[w(1)] -= 1 / b; gD[rx(1)] -= 0.5; gD[rx(2)] -= 0.5
    for xi in _G:
        for eta in _G:
            dNx = 0.25 * _XI * (1 + _ETA * eta) * 2.0 / a
            dNy = 0.25 * _ETA * (1 + _XI * xi) * 2.0 / b
            Bg = np.zeros((6, n))
            for i in range(4):
                Bg[0, 6 * i] = dNx[i]
                Bg[1, 6 * i + 1] = dNy[i]
                Bg[2, 6 * i] = dNy[i]; Bg[2, 6 * i + 1] = dNx[i]
                Bg[3, 6 * i + 4] = dNx[i]
                Bg[4, 6 * i + 3] = -dNy[i]
                Bg[5, 6 * i + 4] = dNy[i]; Bg[5, 6 * i + 3] = -dNx[i]
            dP1x = -2.0 * xi * 2.0 / a
            dP2y = -2.0 * eta * 2.0 / b
            Bg[0, 24] = dP1x            # u: alpha1 (1 - xi^2)
            Bg[2, 25] = dP2y            # u: alpha2 (1 - eta^2)
            Bg[2, 26] = dP1x            # v: alpha3 (1 - xi^2)
            Bg[1, 27] = dP2y            # v: alpha4 (1 - eta^2)
            K += Bg.T @ Cg @ Bg * detJ
            Bs = np.vstack((0.5 * (1 + eta) * gA + 0.5 * (1 - eta) * gC,
                            0.5 * (1 + xi) * gD + 0.5 * (1 - xi) * gB))
            K += Bs.T @ As @ Bs * detJ
    Kdd = K[:24, :24]; Kda = K[:24, 24:]; Kaa = K[24:, 24:]
    Ke = Kdd - Kda @ np.linalg.solve(Kaa, Kda.T)
    if drill is None:
        drill = 1e-6 * float(np.mean(np.diag(D))) if np.mean(np.diag(D)) > 0 else 1e-9
    for i in range(4):
        Ke[6 * i + 5, 6 * i + 5] += drill
    return 0.5 * (Ke + Ke.T)


@dataclass
class ShellModel:
    node_ids: list[int]
    coords: np.ndarray
    K: object
    ndof: int
    index: dict[int, int] = field(default_factory=dict)
    wall_elements: dict = field(default_factory=dict)


def assemble(model: RVEModel, materials: dict[str, Material], laminates: list[Laminate], mesh_size: float,
             wall_section_factors: list[float] | None = None, face_wall_section: str = "half-thickness"):
    """Mesh the wall network exactly as the Abaqus deck does and assemble K (6 DOF/node)."""
    nodes, wall_elements, _, _ = _build_shell_mesh(model, mesh_size)
    ids = sorted(nodes)
    index = {nid: k for k, nid in enumerate(ids)}
    X = np.array([nodes[i] for i in ids], float)
    ndof = 6 * len(ids)
    factors = list(wall_section_factors or [1.0] * len(model.elements))
    rows, cols, vals = [], [], []
    for wi, elems in wall_elements.items():
        e = model.elements[wi - 1]
        f = float(factors[wi - 1])
        if face_wall_section == "half-stiffness":
            sec = shell_section(laminates[e.laminate], materials, modulus_factor=f)
        elif face_wall_section == "half-thickness":
            sec = shell_section(laminates[e.laminate], materials, thickness_factor=f)
        else:
            sec = shell_section(laminates[e.laminate], materials)
        conn = np.array([[index[n] for n in c] for _, c in elems], int)
        x1, x2, x4 = X[conn[0, 0]], X[conn[0, 1]], X[conn[0, 3]]
        a = float(np.linalg.norm(x2 - x1)); b = float(np.linalg.norm(x4 - x1))
        ex = (x2 - x1) / a; ey = (x4 - x1) / b; ez = np.cross(ex, ey)
        if abs(float(ex @ ey)) > 1e-9:
            raise ValueError(f"Wall W{wi}: shell elements are not rectangular.")
        R = np.vstack((ex, ey, ez))
        T = np.kron(np.eye(8), R)
        kl = element_stiffness(a, b, *sec)
        kg = T.T @ kl @ T
        dofs = (6 * conn[:, :, None] + np.arange(6)[None, None, :]).reshape(len(conn), 24)
        rows.append(np.repeat(dofs, 24, axis=1).ravel())
        cols.append(np.tile(dofs, (1, 24)).ravel())
        vals.append(np.tile(kg.ravel(), len(conn)))
    K = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(ndof, ndof)).tocsc()
    return nodes, wall_elements, ShellModel(ids, X, K, ndof, index, wall_elements)


def solve_prescribed(sm: ShellModel, bc_map: dict[tuple[int, int], float],
                     cache: dict | None = None) -> tuple[float, np.ndarray]:
    """Solve K u = 0 with prescribed translations; returns strain energy and u.

    ``cache`` (a dict kept by the caller) reuses the factorisation for load cases
    that prescribe the same degrees of freedom with different values.
    """
    pres = np.array(sorted(6 * sm.index[n] + (d - 1) for (n, d) in bc_map), int)
    val = np.zeros(sm.ndof)
    for (n, d), v in bc_map.items():
        val[6 * sm.index[n] + (d - 1)] = v
    K = sm.K
    key = pres.tobytes()
    if cache is not None and key in cache:
        free, lu, Kfp = cache[key]
    else:
        mask = np.ones(sm.ndof, bool); mask[pres] = False
        free = np.nonzero(mask)[0]
        Kf = K[free]
        lu = _factor(Kf[:, free].tocsc())
        Kfp = Kf[:, pres]
        if cache is not None:
            cache[key] = (free, lu, Kfp)
    u = val.copy()
    u[free] = lu.solve(-(Kfp @ val[pres]))
    energy = 0.5 * float(u @ (K @ u))
    return energy, u


def run_dbc_shell_check(model: RVEModel, materials: dict[str, Material], laminates: list[Laminate], *,
                        mesh_size: float = 10.0, displacement: float = 1.0,
                        face_wall_section: str = "half-thickness", cells: tuple[int, int] | None = None,
                        cases: tuple[str, ...] | None = None, return_fields: bool = False) -> dict:
    """Finite-height DBC cases on the automatic FE RVE with the built-in shell element.

    Shell terms (``out['Q']``): thesis Table 3.4, top/bottom U3 free (sigma33 = 0).
    3D terms (``out['C']``): the same side faces with U3 = 0 on the top and bottom
    faces (eps33 = 0); C1313 and C2323 from the Q44 and Q55 cases.  All ten cases
    share one assembly.  With ``return_fields`` the mesh (node coordinates, element
    connectivity as indices) and the translation field of every case are returned.
    """
    from .abaqus_thesis_dbc import (DBC_CASES, build_thesis_dbc_rectangle, _boundary_node_sets,
                                    _case_constraints, _add_minimal_rbm_gauge, reconstruct_dbc_terms)
    dbc, diag = build_thesis_dbc_rectangle(model, cells=cells, face_wall_section=face_wall_section)
    nodes, wall_elements, sm = assemble(dbc, materials, laminates, mesh_size,
                                        diag.get("wall_section_factors"), face_wall_section)
    p, q, h = float(diag["p"]), float(diag["q"]), float(diag["height"])
    bsets = _boundary_node_sets(nodes, p, q, h)
    energies: dict[str, float] = {}
    max_u: dict[str, float] = {}
    fields: dict[str, list] = {}
    cache: dict = {}
    for case in (cases or DBC_CASES):
        bc: dict[tuple[int, int], float] = {}
        for set_name, dof, v in _case_constraints(case, float(displacement)):
            for nid in bsets[set_name]:
                bc[(int(nid), int(dof))] = float(v)
        _add_minimal_rbm_gauge(nodes, bc, bsets.get("X3_0", []))
        W, u = solve_prescribed(sm, bc, cache)
        energies[case] = W
        max_u[case] = float(np.max(np.abs(u.reshape(-1, 6)[:, :3])))
        if return_fields:
            fields[case] = u.reshape(-1, 6)[:, :3].tolist()
    out = {"energies_Nmm": energies, "max_translation_mm": max_u, "domain": diag,
           "mesh_size": float(mesh_size), "displacement": float(displacement),
           "nodes": len(nodes), "elements": int(sum(len(v) for v in wall_elements.values())),
           "dofs": int(sm.ndof), "face_wall_section": face_wall_section,
           "element": "flat 4-node: Wilson incompatible-mode membrane + MITC4 plate, full CLT A/B/D"}
    out.update(reconstruct_dbc_terms(energies, p, q, h, float(displacement)))
    if return_fields:
        out["fields"] = {"node_ids": list(sm.node_ids), "coords": sm.coords.tolist(),
                         "elements": [[sm.index[nn] for nn in conn] for w in wall_elements.values() for _, conn in w],
                         "u": fields}
    return out


# ----------------------------------------------------------------------------------
# Built-in 3D-periodic check: full C^H on the same FE RVE, no Abaqus
# ----------------------------------------------------------------------------------

VOIGT = ("11", "22", "33", "23", "13", "12")
# Periodic strain states: all six unit strains.  The top and bottom faces are periodic through
# a3 = (0, 0, h), so U3+ - U3- = E33 h; the full 6x6 C follows from the strain energies.
PBC_STATES = (0, 1, 2, 3, 4, 5)
SOLID_TERMS = ("C1111", "C2222", "C1122", "C1212", "C1313", "C2323")


def solid_terms(C) -> dict[str, float]:
    """The six 3D stiffness terms (MPa) of a 6x6 Voigt [11,22,33,23,13,12] matrix."""
    C = np.asarray(C, float)
    return {"C1111": float(C[0, 0]), "C2222": float(C[1, 1]), "C1122": float(C[0, 1]),
            "C1212": float(C[5, 5]), "C1313": float(C[4, 4]), "C2323": float(C[3, 3])}


def _json_matrix(M) -> list:
    return [[None if not math.isfinite(float(x)) else float(x) for x in row] for row in np.asarray(M, float)]


def _unit_strain(k: int) -> np.ndarray:
    """Tensor strain of the k-th unit engineering strain, Voigt [11, 22, 33, 23, 13, 12]."""
    e = np.zeros((3, 3))
    if k < 3:
        e[k, k] = 1.0
    else:
        i, j = {3: (1, 2), 4: (0, 2), 5: (0, 1)}[k]
        e[i, j] = e[j, i] = 0.5
    return e


def run_pbc_shell_check(model: RVEModel, materials: dict[str, Material], laminates: list[Laminate], *,
                        mesh_size: float = 10.0, return_fields: bool = False) -> dict:
    """Full 6x6 periodic (PBC) stiffness from strain energy, built-in shell element, no Abaqus.

    Same FE RVE and mesh as the Abaqus periodic decks (a wall lying in a face is kept
    on both faces at half moduli so that the two periodic images act as one wall).
    Every boundary node is tied to its periodic image, u_s - u_m = E (x_s - x_m) for
    the translations and equal rotations, the top and bottom faces through
    a3 = (0, 0, h).  The six unit strains (E11, E22, E33, G23, G13, G12) are solved with
    one factorisation and C_ij = u_i^T K u_j / V: the diagonal is twice the strain energy
    per volume of a unit strain, the off-diagonal the mutual energy -- the same as
    (W_ij - W_i - W_j) / V from the paired states that Abaqus runs.
    """
    import time
    from scipy.sparse import csr_matrix
    from .abaqus_homogenization import build_published_periodic_domain, periodic_mesh_pairs
    t0 = time.time()
    pm, d = build_published_periodic_domain(model)
    nodes, wall_elements, sm = assemble(pm, materials, laminates, mesh_size, d.get("wall_section_factors"),
                                        "half-stiffness")
    pairs = periodic_mesh_pairs(nodes, pm)
    n = len(sm.node_ids)
    rep = np.arange(n)
    off = np.zeros((n, 3))
    for m_id, s_id, delta in pairs:
        i, j = sm.index[m_id], sm.index[s_id]
        rep[j] = i
        off[j] = np.asarray(delta, float)
    if np.any(rep[rep] != rep):
        raise ValueError("Periodic pairing is not one level deep; cannot build the periodic map.")
    reps = np.unique(rep)
    rmap = -np.ones(n, int)
    rmap[reps] = np.arange(len(reps))
    nred = 6 * len(reps)
    rows = np.arange(6 * n)
    cols = (6 * rmap[rep][:, None] + np.arange(6)[None, :]).ravel()
    T = csr_matrix((np.ones(6 * n), (rows, cols)), shape=(6 * n, nred))
    K = sm.K.tocsr()
    Kr = (T.T @ K @ T).tocsc()
    anchor = int(rmap[reps[0]])
    fixed = {6 * anchor, 6 * anchor + 1, 6 * anchor + 2}
    free = np.array([k for k in range(nred) if k not in fixed], int)
    lu = _factor(Kr[free][:, free].tocsc())
    p, q, h = float(d["rectangle_width"]), float(d["rectangle_height"]), float(model.height)
    V = p * q * h
    U = {}
    for k in PBC_STATES:
        g = np.zeros(6 * n)
        g.reshape(-1, 6)[:, :3] = off @ _unit_strain(k).T
        rhs = -(T.T @ (K @ g))
        ur = np.zeros(nred)
        ur[free] = lu.solve(rhs[free])
        U[k] = T @ ur + g
    KU = {k: K @ u for k, u in U.items()}
    C = np.full((6, 6), np.nan)
    for i in PBC_STATES:
        for j in PBC_STATES:
            C[i, j] = float(U[i] @ KU[j]) / V
    blk = np.ix_(PBC_STATES, PBC_STATES)
    asym = float(np.max(np.abs(C[blk] - C[blk].T)) / max(np.max(np.abs(C[blk])), 1e-30))
    C[blk] = 0.5 * (C[blk] + C[blk].T)
    out = {"CH": _json_matrix(C), "C": solid_terms(C), "symmetry_error": asym, "volume": V, "domain": d,
           "mesh_size": float(mesh_size), "nodes": n, "elements": int(sum(len(v) for v in wall_elements.values())),
           "periodic_pairs": len(pairs), "dofs_reduced": int(len(free)), "seconds": round(time.time() - t0, 3),
           "element": "flat 4-node: Wilson incompatible-mode membrane + MITC4 plate, full CLT A/B/D",
           "boundary_conditions": ("periodic side faces; top and bottom periodic through a3 = (0, 0, h): "
                                   "U3+ - U3- = E33 h, U1, U2 and rotations periodic"),
           "method": "strain energy, C_ij = u_i^T K u_j / V, six unit strains"}
    if return_fields:
        out["fields"] = {"node_ids": list(sm.node_ids), "coords": sm.coords.tolist(),
                         "elements": [[sm.index[nn] for nn in conn] for w in wall_elements.values() for _, conn in w],
                         "u": {VOIGT[k]: U[k].reshape(-1, 6)[:, :3].tolist() for k in PBC_STATES}}
    return out
