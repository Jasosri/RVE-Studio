"""Stable, GUI-independent API for RVE Studio.

All stress/stiffness values use the unit of the supplied material moduli;
geometry and ply thicknesses share the model length unit. Density is kg/m3.
The 3D CH matrix retains Voigt [11,22,33,23,13,12]. The *paper-aligned*
2x2 transverse Qshear is ordered [13,23], so Q44=G13 and Q55=G23.
"""
from __future__ import annotations
from pathlib import Path
from typing import Mapping
import numpy as np
from .models import Material, Laminate, RVEModel, HomogenizationResult
from .examples_data import paper_hexagonal_case
from .homogenization import homogenize_rve
from .io import load_project
from .validation import extract_six, run_all14_benchmark, PROPERTY_NAMES, PAPER_URL

__all__ = ['solve', 'solve_example', 'solve_project', 'summary', 'benchmark', 'PAPER_URL']


def solve(model: RVEModel, materials: Mapping[str, Material], laminates: list[Laminate]) -> HomogenizationResult:
    """Run periodic strain-energy homogenisation for a specified RVE."""
    if not isinstance(model, RVEModel):
        raise TypeError('model must be an RVEModel')
    if not materials:
        raise ValueError('at least one material is required')
    if not laminates:
        raise ValueError('at least one laminate is required')
    return homogenize_rve(model, materials, laminates)


def solve_example(config: int = 1) -> HomogenizationResult:
    """Run one of the 14 wall configurations in the referenced journal paper."""
    model, materials, laminates = paper_hexagonal_case(config)
    return solve(model, materials, laminates)


def solve_project(path: str | Path) -> HomogenizationResult:
    """Load an .rvejson project and solve it in one call."""
    return solve(*load_project(path))


from .reliability import assess


def summary(result: HomogenizationResult, *, include_ch: bool = False) -> dict:
    """Return a JSON-serialisable summary with explicit shear axis labels."""
    p = result.props
    values = dict(zip(PROPERTY_NAMES, map(float, extract_six(result))))
    output = {
        'software': 'RVE Studio',
        'reference': PAPER_URL,
        'stress_unit': 'same as input material moduli',
        'voigt_order_3d': ['11','22','33','23','13','12'],
        'transverse_shear_order': ['13','23'],
        'stiffness': values,
        'engineering_constants': {
            name: float(getattr(p, name)) for name in
            ('E1','E2','E3','G12','G13','G23','nu12','nu21','nu13','nu31','nu23','nu32')
        },
        'rve': {'area': float(result.area), 'volume': float(result.volume)},
        'reliability': assess(result).as_dict(),
    }
    if include_ch:
        output['CH'] = np.asarray(result.CH).tolist()
        output['Qplane'] = np.asarray(p.Qplane).tolist()
        output['Qshear'] = np.asarray(p.Qshear).tolist()
    return output


def benchmark(*, tolerance_pct: float = 1.0, wall_shear: str = 'euler-bernoulli', bending: str = 'cylindrical',
              twist_release=None) -> dict:
    """The 14 benchmark designs against the reference values (pandas tables and numeric matrices).

    wall_shear: 'euler-bernoulli' (original model, default) or 'timoshenko'; bending: 'cylindrical' (default),
    'twist-corrected', 'anisotropic strip' or 'free moment'; twist_release: (a, b, n) of the twist-corrected wall."""
    return run_all14_benchmark(tolerance_pct=tolerance_pct, wall_shear=wall_shear, bending=bending, twist_release=twist_release)
