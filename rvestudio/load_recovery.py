from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable
import numpy as np

from .models import RVEModel, HomogenizationResult
from .walls import wall_transform


MACRO_LABELS = ("11", "22", "33", "23", "13", "12")


@dataclass(slots=True)
class WallLoadRecoveryResult:
    """Recovered wall membrane resultants for one macroscopic load state.

    Sign convention follows the homogenization solver. Normal compression is
    negative.  Wall resultants are per unit length (N/mm when the model uses
    MPa and mm): [N_11bar, N_33bar, N_13bar].
    """

    control_mode: str
    macro_strain: np.ndarray
    macro_stress: np.ndarray
    displacement: np.ndarray
    lagrange_multipliers: np.ndarray
    walls: list[dict]
    validation: dict

    @property
    def minimum_normal_resultant(self) -> float:
        if not self.walls:
            return float("nan")
        return float(min(min(w["N11_N_per_mm"], w["N33_N_per_mm"]) for w in self.walls))


def _as_six(values: Iterable[float], name: str) -> np.ndarray:
    arr = np.asarray(list(values), dtype=float).reshape(-1)
    if arr.size != 6 or not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain six finite values in [11,22,33,23,13,12] order")
    return arr


def _relative_error(a: float, b: float, floor: float = 1e-14) -> float:
    return float(abs(a - b) / max(abs(a), abs(b), floor))


def recover_wall_loads(
    model: RVEModel,
    result: HomogenizationResult,
    *,
    macro_strain: Iterable[float] | None = None,
    macro_stress: Iterable[float] | None = None,
) -> WallLoadRecoveryResult:
    """Recover wall-local membrane resultants for an arbitrary macro load state.

    Exactly one of ``macro_strain`` or ``macro_stress`` must be supplied.
    Engineering macro vectors use the same convention as the homogenization
    solver: [eps11, eps22, eps33, gamma23, gamma13, gamma12].  For stress
    control, [sigma11, sigma22, sigma33, tau23, tau13, tau12] is converted to
    strain using the already-computed homogenized compliance.

    The wall-local plate strain ordering is [eps_11bar, eps_33bar,
    gamma_13bar], and the corresponding laminate membrane resultants are
    [N_11bar, N_33bar, N_13bar] = A [eps_11bar, eps_33bar, gamma_13bar].
    """

    if (macro_strain is None) == (macro_stress is None):
        raise ValueError("Provide exactly one of macro_strain or macro_stress")

    if macro_strain is not None:
        E = _as_six(macro_strain, "macro_strain")
        control_mode = "strain"
    else:
        target_sigma = _as_six(macro_stress, "macro_stress")
        E = np.asarray(result.S, float) @ target_sigma
        control_mode = "stress"

    sigma = np.asarray(result.CH, float) @ E
    u = np.asarray(result.Umacro, float) @ E
    lag = np.asarray(result.Lambda, float) @ E

    # Global consistency checks use exactly the assembled matrices solved by
    # homogenize_rve, so they test the recovery against the strain-energy system.
    internal = np.asarray(result.K, float) @ u + np.asarray(result.H, float) @ E
    constraint_reaction = np.asarray(result.constraints_C, float).T @ lag
    equilibrium_residual = internal + constraint_reaction
    eq_scale = max(float(np.linalg.norm(internal)), float(np.linalg.norm(constraint_reaction)), 1.0)
    equilibrium_relative = float(np.linalg.norm(equilibrium_residual) / eq_scale)

    micro_energy = float(
        0.5 * u @ np.asarray(result.K, float) @ u
        + u @ np.asarray(result.H, float) @ E
        + 0.5 * E @ np.asarray(result.R, float) @ E
    )
    macro_energy = float(0.5 * result.volume * E @ np.asarray(result.CH, float) @ E)
    macro_work_half = float(0.5 * result.volume * E @ sigma)

    wall_rows: list[dict] = []
    element_energy_sum = 0.0
    max_membrane_force_error = 0.0
    max_n33_conjugate_error = 0.0

    for idx, (element, edata) in enumerate(zip(model.elements, result.element_data)):
        L = float(edata["L"])
        alpha = float(edata["alpha"])
        edofs = np.asarray(edata["edofs"], dtype=int)
        weight = float(element.weight)
        dglobal = u[edofs]
        T = np.asarray(edata.get("T", wall_transform(alpha)), float)
        dlocal = T.T @ dglobal

        eps11 = float((dlocal[4] - dlocal[0]) / L)
        eps33 = float(E[2])
        gamma13 = float((dlocal[6] - dlocal[2]) / L)
        eps_wall = np.array([eps11, eps33, gamma13], dtype=float)

        section = result.sections[element.laminate]
        A = np.asarray(section.A, float)
        N = A @ eps_wall
        N11, N33, N13 = (float(N[0]), float(N[1]), float(N[2]))

        membrane_energy = float(0.5 * weight * L * model.height * (eps_wall @ N))

        # Use the exact local element matrices when available. Older stored
        # result objects can still be processed, but new v1.3.5 results carry
        # these matrices so the identities below are independently checked.
        Kloc = edata.get("Klocal")
        gloc = edata.get("glocal")
        rloc = edata.get("rlocal")
        if Kloc is not None and gloc is not None and rloc is not None:
            Kloc = np.asarray(Kloc, float)
            gloc = np.asarray(gloc, float)
            rloc = float(rloc)
            qlocal = weight * (Kloc @ dlocal + gloc * eps33)
            expected_membrane_q = weight * model.height * np.array(
                [-N11, 0.0, -N13, 0.0, N11, 0.0, N13, 0.0], dtype=float
            )
            membrane_idx = np.array([0, 2, 4, 6], dtype=int)
            ferr = float(np.linalg.norm(qlocal[membrane_idx] - expected_membrane_q[membrane_idx]))
            fscale = max(float(np.linalg.norm(expected_membrane_q[membrane_idx])), 1.0)
            membrane_force_error = ferr / fscale

            n33_generalized = weight * float(gloc @ dlocal + rloc * eps33)
            n33_expected = weight * model.height * L * N33
            # Scale by the wall's largest resultant: a wall that carries (almost) no N33 would otherwise
            # compare two round-off numbers and report a spurious error.
            n33_scale = weight * model.height * L * max(abs(N11), abs(N33), abs(N13))
            n33_conjugate_error = _relative_error(n33_generalized, n33_expected, floor=max(n33_scale, 1e-300))

            total_element_energy = weight * float(
                0.5 * dlocal @ Kloc @ dlocal + dlocal @ gloc * eps33 + 0.5 * rloc * eps33 * eps33
            )
        else:
            membrane_force_error = float("nan")
            n33_conjugate_error = float("nan")
            total_element_energy = float("nan")

        if np.isfinite(total_element_energy):
            element_energy_sum += total_element_energy
        if np.isfinite(membrane_force_error):
            max_membrane_force_error = max(max_membrane_force_error, membrane_force_error)
        if np.isfinite(n33_conjugate_error):
            max_n33_conjugate_error = max(max_n33_conjugate_error, n33_conjugate_error)

        group = model.wall_groups[idx] if idx < len(model.wall_groups) else ""
        wall_rows.append({
            "wall": idx + 1,
            "group": group,
            "laminate": int(element.laminate) + 1,
            "weight": weight,
            "length_mm": L,
            "angle_deg": float(np.degrees(alpha)),
            "eps11_local": eps11,
            "eps33_local": eps33,
            "gamma13_local": gamma13,
            "N11_N_per_mm": N11,
            "N33_N_per_mm": N33,
            "N13_N_per_mm": N13,
            "membrane_energy_Nmm": membrane_energy,
            "total_element_energy_Nmm": total_element_energy,
            "membrane_force_identity_error": membrane_force_error,
            "N33_conjugate_identity_error": n33_conjugate_error,
        })

    constraint_lhs = np.asarray(result.constraints_C, float) @ u
    constraint_rhs = np.asarray(result.constraints_Dmacro, float) @ E
    constraint_norm = float(np.linalg.norm(constraint_lhs - constraint_rhs))
    constraint_scale = max(float(np.linalg.norm(constraint_lhs)), float(np.linalg.norm(constraint_rhs)), 1.0)

    validation = {
        "equilibrium_relative_residual": equilibrium_relative,
        "energy_relative_error": _relative_error(micro_energy, macro_energy),
        "macro_work_relative_error": _relative_error(micro_energy, macro_work_half),
        "element_energy_sum_relative_error": _relative_error(element_energy_sum, micro_energy),
        "max_membrane_force_identity_error": max_membrane_force_error,
        "max_N33_conjugate_identity_error": max_n33_conjugate_error,
        "micro_energy_Nmm": micro_energy,
        "macro_energy_Nmm": macro_energy,
        "element_energy_sum_Nmm": element_energy_sum,
        "constraint_residual_norm": constraint_norm,
        "constraint_relative_residual": float(constraint_norm / constraint_scale),
    }
    if control_mode == "stress":
        target_sigma = _as_six(macro_stress, "macro_stress")
        denom = max(float(np.linalg.norm(target_sigma)), 1.0)
        validation["stress_control_relative_error"] = float(np.linalg.norm(sigma - target_sigma) / denom)

    return WallLoadRecoveryResult(
        control_mode=control_mode,
        macro_strain=E,
        macro_stress=sigma,
        displacement=u,
        lagrange_multipliers=lag,
        walls=wall_rows,
        validation=validation,
    )
