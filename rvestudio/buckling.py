from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np
from scipy.linalg import eig

from .models import RVEModel, HomogenizationResult, LaminateSection
from .load_recovery import WallLoadRecoveryResult


@dataclass(slots=True)
class PlateBucklingResult:
    """Linear CLPT buckling result for one rectangular laminated wall plate.

    The plate occupies ``0<=x<=a`` and ``0<=y<=b`` and is simply supported on
    all four edges for transverse displacement.  The prescribed pre-buckling
    membrane resultants use the laminate convention [Nx, Ny, Nxy] in N/mm;
    compression is negative.  ``lambda_cr`` scales the complete combined load
    vector [Nx, Ny, Nxy] to the first bifurcation state.

    For unsymmetric laminates, ``use_abd_condensation=True`` uses the standard
    free-membrane condensed bending matrix D* = D - B A^-1 B.  This is a
    reduced-order CLPT screening model; exact unsymmetric-laminate buckling can
    depend on the in-plane edge restraints and requires a coupled plate model.
    """

    lambda_cr: float
    eigenvalues: np.ndarray
    mode_vector: np.ndarray
    dominant_m: int | None
    dominant_n: int | None
    D_effective: np.ndarray
    K_bending: np.ndarray
    K_geometric: np.ndarray
    modes: list[tuple[int, int]]


@dataclass(slots=True)
class WallBucklingSummary:
    walls: list[dict]
    minimum_lambda: float
    critical_wall: int | None
    safety_factor: float
    passed: bool
    settings: dict


def condensed_bending_matrix(section: LaminateSection, use_abd_condensation: bool = True) -> np.ndarray:
    """Return the CLPT bending matrix used by the wall-plate buckling model."""
    D = np.asarray(section.Doriginal, dtype=float)
    if use_abd_condensation:
        A = np.asarray(section.A, dtype=float)
        B = np.asarray(section.B, dtype=float)
        if np.linalg.cond(A) > 1e13:
            raise np.linalg.LinAlgError("Laminate A matrix is singular/ill-conditioned for ABD condensation")
        D = D - B @ np.linalg.solve(A, B)
    D = 0.5 * (D + D.T)
    vals = np.linalg.eigvalsh(D)
    if not np.all(np.isfinite(vals)) or float(np.min(vals)) <= 0.0:
        raise ValueError("Effective laminate bending matrix is not positive definite")
    return D


def _sc_integral(m: int, p: int, L: float) -> float:
    """Integral_0^L sin(m*pi*x/L) cos(p*pi*x/L) dx."""
    def term(k: int) -> float:
        if k == 0:
            return 0.0
        return (1.0 - (-1.0) ** k) / float(k)
    return float(L / (2.0 * math.pi) * (term(m + p) + term(m - p)))


def _ritz_matrices(D: np.ndarray, a: float, b: float, Nx: float, Ny: float, Nxy: float,
                   m_max: int, n_max: int) -> tuple[np.ndarray, np.ndarray, list[tuple[int, int]]]:
    if a <= 0 or b <= 0:
        raise ValueError("Plate dimensions a and b must be positive")
    if m_max < 1 or n_max < 1:
        raise ValueError("m_max and n_max must be >= 1")
    modes = [(m, n) for m in range(1, m_max + 1) for n in range(1, n_max + 1)]
    nm = len(modes)
    Kb = np.zeros((nm, nm), dtype=float)
    Kg = np.zeros((nm, nm), dtype=float)
    ab4 = a * b / 4.0

    alphas = np.array([m * math.pi / a for m, _ in modes], dtype=float)
    betas = np.array([n * math.pi / b for _, n in modes], dtype=float)

    for i, (mi, ni) in enumerate(modes):
        ai = alphas[i]; bi = betas[i]
        for j in range(i, nm):
            mj, nj = modes[j]; aj = alphas[j]; bj = betas[j]
            same = (mi == mj and ni == nj)
            kij = 0.0
            gij = 0.0
            if same:
                # [kappa_x, kappa_y, kappa_xy] =
                # [-a^2 sin sin, -b^2 sin sin, 2ab cos cos].
                kij += ab4 * (
                    D[0, 0] * ai**2 * aj**2
                    + D[0, 1] * (ai**2 * bj**2 + bi**2 * aj**2)
                    + D[1, 1] * bi**2 * bj**2
                    + D[2, 2] * (2.0 * ai * bi) * (2.0 * aj * bj)
                )
                gij += ab4 * (Nx * ai * aj + Ny * bi * bj)

            # D16/D26 and membrane-shear Nxy couple different sine modes.
            sx_ij = _sc_integral(mi, mj, a)   # sin(i) cos(j)
            sx_ji = _sc_integral(mj, mi, a)   # sin(j) cos(i)
            sy_ij = _sc_integral(ni, nj, b)
            sy_ji = _sc_integral(nj, ni, b)

            kx_i_ks_j = (-ai**2) * (2.0 * aj * bj) * sx_ij * sy_ij
            ks_i_kx_j = (2.0 * ai * bi) * (-aj**2) * sx_ji * sy_ji
            ky_i_ks_j = (-bi**2) * (2.0 * aj * bj) * sx_ij * sy_ij
            ks_i_ky_j = (2.0 * ai * bi) * (-bj**2) * sx_ji * sy_ji
            kij += D[0, 2] * (kx_i_ks_j + ks_i_kx_j)
            kij += D[1, 2] * (ky_i_ks_j + ks_i_ky_j)

            # phi_x_i * phi_y_j and its transpose contribution.
            # cos(mi)sin(mj) = SC(mj,mi); sin(ni)cos(nj) = SC(ni,nj).
            gx_i_gy_j = ai * bj * sx_ji * sy_ij
            gy_i_gx_j = bi * aj * sx_ij * sy_ji
            gij += Nxy * (gx_i_gy_j + gy_i_gx_j)

            Kb[i, j] = Kb[j, i] = kij
            Kg[i, j] = Kg[j, i] = gij

    return Kb, Kg, modes


def laminated_plate_buckling(section: LaminateSection, a: float, b: float,
                              Nx: float, Ny: float, Nxy: float, *,
                              m_max: int = 8, n_max: int = 8,
                              use_abd_condensation: bool = True) -> PlateBucklingResult:
    """Solve the simply-supported multi-mode CLPT linear buckling eigenproblem.

    ``Nx``, ``Ny`` and ``Nxy`` are the *current* combined pre-buckling membrane
    resultants.  The smallest positive eigenvalue scales the entire load vector
    to bifurcation.  Thus lambda_cr=1 is buckling at the current load and
    lambda_cr=1.5 gives a 1.5 linear buckling load factor.
    """
    D = condensed_bending_matrix(section, use_abd_condensation)
    Kb, Kg, modes = _ritz_matrices(D, float(a), float(b), float(Nx), float(Ny), float(Nxy), int(m_max), int(n_max))
    if np.linalg.norm(Kg, ord="fro") <= np.finfo(float).eps:
        return PlateBucklingResult(float("inf"), np.array([], float), np.zeros(len(modes)), None, None, D, Kb, Kg, modes)

    vals, vecs = eig(Kb, -Kg)
    vals = np.asarray(vals)
    scale = np.maximum(1.0, np.abs(vals.real))
    valid = np.isfinite(vals.real) & np.isfinite(vals.imag) & (np.abs(vals.imag) <= 1e-7 * scale) & (vals.real > 1e-10)
    if not np.any(valid):
        return PlateBucklingResult(float("inf"), np.array([], float), np.zeros(len(modes)), None, None, D, Kb, Kg, modes)

    real_vals = vals.real[valid]
    valid_indices = np.flatnonzero(valid)
    order = np.argsort(real_vals)
    eigvals = real_vals[order]
    idx = int(valid_indices[order[0]])
    v = np.asarray(vecs[:, idx].real, dtype=float)
    vmax = float(np.max(np.abs(v))) if v.size else 0.0
    if vmax > 0:
        v = v / vmax
    dom = int(np.argmax(np.abs(v))) if v.size else -1
    dm, dn = modes[dom] if dom >= 0 else (None, None)
    return PlateBucklingResult(float(eigvals[0]), eigvals, v, dm, dn, D, Kb, Kg, modes)


def evaluate_wall_buckling(model: RVEModel, result: HomogenizationResult, recovery: WallLoadRecoveryResult, *,
                           m_max: int = 8, n_max: int = 8, safety_factor: float = 1.5,
                           use_abd_condensation: bool = True, convergence_check: bool = True) -> WallBucklingSummary:
    """Evaluate all RVE walls under the recovered combined biaxial+shear load."""
    if safety_factor <= 0:
        raise ValueError("safety_factor must be positive")
    if len(recovery.walls) != len(model.elements):
        raise ValueError("Wall load recovery does not match the current RVE element count")
    rows: list[dict] = []
    min_lambda = float("inf")
    critical_wall: int | None = None

    lower_m = max(1, int(m_max) - 2)
    lower_n = max(1, int(n_max) - 2)
    for idx, (element, loadrow, edata) in enumerate(zip(model.elements, recovery.walls, result.element_data)):
        section = result.sections[element.laminate]
        a = float(edata["L"])
        b = float(model.height)
        Nx = float(loadrow["N11_N_per_mm"])
        Ny = float(loadrow["N33_N_per_mm"])
        Nxy = float(loadrow["N13_N_per_mm"])
        br = laminated_plate_buckling(section, a, b, Nx, Ny, Nxy,
                                      m_max=m_max, n_max=n_max,
                                      use_abd_condensation=use_abd_condensation)
        conv = float("nan")
        lower_lambda = float("nan")
        if convergence_check and (lower_m != m_max or lower_n != n_max):
            low = laminated_plate_buckling(section, a, b, Nx, Ny, Nxy,
                                           m_max=lower_m, n_max=lower_n,
                                           use_abd_condensation=use_abd_condensation)
            lower_lambda = low.lambda_cr
            if math.isfinite(br.lambda_cr) and math.isfinite(lower_lambda):
                conv = abs(br.lambda_cr - lower_lambda) / max(abs(br.lambda_cr), 1e-14)
            elif br.lambda_cr == lower_lambda:
                conv = 0.0
        lam = float(br.lambda_cr)
        if lam < min_lambda:
            min_lambda = lam
            critical_wall = idx + 1
        passed = bool(lam >= safety_factor)
        rows.append({
            "wall": idx + 1,
            "group": model.wall_groups[idx] if idx < len(model.wall_groups) else "",
            "laminate": int(element.laminate) + 1,
            "a_wall_length_mm": a,
            "b_core_height_mm": b,
            "N11_N_per_mm": Nx,
            "N33_N_per_mm": Ny,
            "N13_N_per_mm": Nxy,
            "lambda_cr": lam,
            "required_safety_factor": float(safety_factor),
            "buckling_margin": lam / float(safety_factor) if math.isfinite(lam) else float("inf"),
            "pass": passed,
            "dominant_m": br.dominant_m if br.dominant_m is not None else "",
            "dominant_n": br.dominant_n if br.dominant_n is not None else "",
            "mode_terms": int(m_max) * int(n_max),
            "lower_mode_lambda": lower_lambda,
            "mode_convergence_relative_error": conv,
            "D11_eff_Nmm": float(br.D_effective[0, 0]),
            "D22_eff_Nmm": float(br.D_effective[1, 1]),
            "D12_eff_Nmm": float(br.D_effective[0, 1]),
            "D16_eff_Nmm": float(br.D_effective[0, 2]),
            "D26_eff_Nmm": float(br.D_effective[1, 2]),
            "D66_eff_Nmm": float(br.D_effective[2, 2]),
            "laminate_symmetry_indicator": float(section.symmetry_indicator),
        })

    return WallBucklingSummary(
        walls=rows,
        minimum_lambda=float(min_lambda),
        critical_wall=critical_wall,
        safety_factor=float(safety_factor),
        passed=bool(min_lambda >= safety_factor),
        settings={
            "boundary_condition": "SSSS",
            "theory": "CLPT multi-mode Ritz/Navier",
            "m_max": int(m_max),
            "n_max": int(n_max),
            "use_abd_condensation": bool(use_abd_condensation),
            "convergence_check": bool(convergence_check),
            "note": "D*=D-B A^-1 B is a reduced-order treatment of membrane-bending coupling; exact unsymmetric-laminate buckling depends on in-plane edge restraint.",
        },
    )
