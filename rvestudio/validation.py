"""Reproduce the 14 laminated-wall cases from Sriharan et al. (2024).

Reference: Engineering Structures 299, 117106,
https://doi.org/10.1016/j.engstruct.2023.117106

Paper Table 3–9 proposed strain-energy values are reported in MPa.  The
publication defines its transverse shear matrix in the order [13, 23]
(S44=C1313, S55=C2323), unlike the standard internal 3D Voigt [23,13].
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from .examples_data import paper_hexagonal_case
from .homogenization import homogenize_rve

PAPER_DOI = '10.1016/j.engstruct.2023.117106'
PAPER_URL = 'https://www.sciencedirect.com/science/article/pii/S0141029623015213'
PROPERTY_NAMES = ['Q11','Q22','Q12','Q66','Q44','Q55']

# Reference values of the strain-energy model of Sriharan et al. (2024), NOT FE results.
# Paper Table 3–9; rounded values transcribed from the article. Q44/Q55 use
# the paper's Eq. (5) assignment, rather than generic 3D Voigt order.
PAPER_PUBLISHED = np.array([
[324.71,324.71,319.85,2.43,341.61,341.61],
[595.17,595.17,584.25,5.46,149.39,149.39],
[482.53,482.53,478.29,2.12,149.39,149.39],
[403.27,403.27,393.65,4.81,231.45,231.45],
[401.50,401.50,395.42,3.04,231.42,231.42],
[228.00,228.00,223.19,2.40,384.70,384.70],
[385.23,385.23,370.15,7.54,745.01,745.01],
[1793.80,1793.80,1773.10,10.31,170.32,170.32],
[1807.10,1807.10,1759.80,23.65,170.32,170.32],
[1681.60,1681.60,1667.60,6.97,170.32,170.32],
[448.45,448.45,438.80,4.83,243.00,243.00],
[439.45,439.45,428.81,5.32,429.05,429.05],
[1799.70,1799.70,1767.20,16.26,170.32,170.32],
[360.14,360.14,350.77,4.68,581.06,581.06],
], dtype=float)


def extract_six(result) -> np.ndarray:
    """Return [Q11,Q22,Q12,Q66,Q44(13),Q55(23)] in paper notation."""
    q = result.props.Qplane
    s = result.props.Qshear
    return np.array([q[0,0], q[1,1], q[0,1], q[2,2], s[0,0], s[1,1]], dtype=float)


BENDING_OPTIONS = {
    "cylindrical": ("paper", "cylindrical"),
    "twist-corrected": ("anisotropicStrip", "twistRelease"),
    "anisotropic strip": ("anisotropicStrip", "cylindrical"),
    "free moment": ("anisotropicStrip", "freeMoment"),
}


def bending_option(name: str) -> tuple[str, str]:
    """(wall_element_type, strip_bending_mode) of a bending option name."""
    key = str(name or "cylindrical").strip().lower().replace("_", " ").replace("-", " ")
    for k, v in BENDING_OPTIONS.items():
        if key == k.replace("-", " ") or key.startswith(k.split()[0][:5]):
            return v
    raise ValueError(f"Unknown bending option {name!r}; use one of {list(BENDING_OPTIONS)}")


def run_all14_benchmark(tolerance_pct: float = 1.0, wall_shear: str = 'euler-bernoulli',
                        bending: str = 'cylindrical', twist_release=None) -> dict:
    """The 14 benchmark designs of Sriharan et al. (2024) against the reference values.

    ``passed`` means all 84 entries agree within the chosen relative tolerance. The reference values are those
    of the original model (Euler-Bernoulli walls, cylindrical bending), which the software reproduces; with any
    other wall beam (``wall_shear``) or bending option (``bending``: cylindrical, twist-corrected, anisotropic
    strip, free moment) the differences show how far that option departs from the original model.
    ``twist_release`` = (a, b, n) overrides the calibrated twist-release coefficients of the twist-corrected wall.
    """
    from .walls import is_timoshenko
    beam = 'timoshenko' if is_timoshenko(wall_shear) else 'euler-bernoulli'
    etype, mode = bending_option(bending)
    if not np.isfinite(tolerance_pct) or tolerance_pct <= 0:
        raise ValueError('tolerance_pct must be a finite positive percentage')
    calc = np.empty_like(PAPER_PUBLISHED)
    for i in range(14):
        model, mats, lams = paper_hexagonal_case(i+1)
        model.wall_shear = beam; model.wall_element_type = etype; model.strip_bending_mode = mode
        model.twist_release = None if twist_release is None else tuple(twist_release)
        calc[i] = extract_six(homogenize_rve(model, mats, lams))
    error = 100.0 * (calc - PAPER_PUBLISHED) / PAPER_PUBLISHED
    abs_error = np.abs(error)
    detail = pd.DataFrame({
        'Config': np.repeat(np.arange(1, 15), 6),
        'Property': PROPERTY_NAMES * 14,
        'Calculated_MPa': calc.reshape(-1),
        'Paper_MPa': PAPER_PUBLISHED.reshape(-1),
        'Error_pct': error.reshape(-1),
        'Within_tolerance': (abs_error <= tolerance_pct).reshape(-1),
    })
    summary = pd.DataFrame({
        'Property': PROPERTY_NAMES,
        'MAE_pct': np.mean(abs_error, axis=0),
        'MaxAbs_pct': np.max(abs_error, axis=0),
        'Within_tolerance': np.sum(abs_error <= tolerance_pct, axis=0),
        'Cases': [14] * 6,
    })
    return {
        'reference': PAPER_URL,
        'doi': PAPER_DOI,
        'tolerance_pct': float(tolerance_pct),
        'wall_shear': beam,
        'bending': [k for k, v in BENDING_OPTIONS.items() if v == (etype, mode)][0],
        'calculated': calc,
        'paper_published': PAPER_PUBLISHED.copy(),
        'error_vs_paper_pct': error,
        'detail': detail,
        'summary': summary,
        'max_abs_paper_pct': float(np.max(abs_error)),
        'within_tolerance': int(np.count_nonzero(abs_error <= tolerance_pct)),
        'total': int(abs_error.size),
        'passed': bool(np.all(abs_error <= tolerance_pct)),
    }
