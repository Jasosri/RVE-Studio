from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence
import math
import numpy as np


@dataclass(slots=True)
class Material:
    """Orthotropic lamina material. Moduli use one consistent stress unit."""
    name: str
    E1: float
    E2: float
    G12: float
    nu12: float
    Gth: Optional[float] = None
    density_kg_m3: float = math.nan
    source: str = ""
    source_url: str = ""
    cost_per_kg: float = math.nan

    def __post_init__(self) -> None:
        for name in ("E1", "E2", "G12"):
            if not math.isfinite(float(getattr(self, name))) or float(getattr(self, name)) <= 0:
                raise ValueError(f"{self.name}: {name} must be positive and finite")
        if not math.isfinite(float(self.nu12)):
            raise ValueError(f"{self.name}: nu12 must be finite")
        if self.Gth is None:
            self.Gth = float(self.G12)
        if not math.isfinite(float(self.Gth)) or float(self.Gth) <= 0:
            raise ValueError(f"{self.name}: Gth must be positive and finite")
        if math.isfinite(float(self.cost_per_kg)) and float(self.cost_per_kg) < 0:
            raise ValueError(f"{self.name}: cost_per_kg must be non-negative when provided")


@dataclass(slots=True)
class Ply:
    material: str
    angle_deg: float
    thickness: float

    def __post_init__(self) -> None:
        if not math.isfinite(float(self.thickness)) or float(self.thickness) <= 0:
            raise ValueError("Ply thickness must be positive and finite")


@dataclass(slots=True)
class Laminate:
    plies: list[Ply]
    name: str = "Laminate"
    section_override: Optional["LaminateSection"] = None

    @property
    def thickness(self) -> float:
        return float(sum(p.thickness for p in self.plies))


@dataclass(slots=True)
class LaminateSection:
    A: np.ndarray
    B: np.ndarray
    Bneutral: np.ndarray
    Doriginal: np.ndarray
    D: np.ndarray
    Seq: float
    A44: float
    kappa_shear: float
    d_neutral: float
    thickness: float
    z: np.ndarray = field(repr=False)
    Qbar: np.ndarray = field(repr=False)
    symmetry_indicator: float = 0.0


@dataclass(slots=True)
class ReducedBending:
    Dfull: np.ndarray
    Dsection: np.ndarray
    Dstrip: float
    curvature_ratio: np.ndarray
    used_B_reduction: bool
    symmetry_measure: float
    mode: str


@dataclass(slots=True)
class Element:
    node_i: int
    node_j: int
    laminate: int = 0
    weight: float = 1.0


@dataclass(slots=True)
class RVEModel:
    nodes: np.ndarray
    elements: list[Element]
    a1: np.ndarray
    a2: np.ndarray
    height: float
    periodic_pairs: Optional[np.ndarray] = None
    periodic_tol: float = 1e-8
    anchor_node: int = 0
    fix_anchor_rotation: bool = False
    wall_element_type: str = "paper"
    strip_bending_mode: str = "cylindrical"
    length_unit: str = "mm"
    gravity: float = 9.80665
    name: str = "RVE"
    wall_groups: list[str] = field(default_factory=list)
    # Wall beam theory: "euler-bernoulli" (the original model [Sriharan et al. 2024]) or
    # "timoshenko" (adds wall shear flexibility through beta = 12 D11 / (L^2 S_eq)).
    wall_shear: str = "euler-bernoulli"
    # FE RVE used for FE verification (rvestudio.fe_rve).  None = automatic rule of the source study;
    # {"mode": "cells", "cells": [n1, n2], "centre_cell": k} = a chosen block of rectangular
    # periodic cells; {"mode": "window", "cells": [n1, n2], "centre": [x1, x2]} = a rectangle
    # drawn over the replicated design (centre in the FE frame, x1 along a1).
    fe_rve: Optional[dict] = None
    # Twist-corrected bending option: coefficients (a, b, n) of eta = (a + b psi^2) lambda^n; None = the calibrated
    # defaults of rvestudio.laminate (2.0, 2.15, 0.73).
    twist_release: Optional[tuple] = None

    def __post_init__(self) -> None:
        self.nodes = np.asarray(self.nodes, dtype=float)
        self.a1 = np.asarray(self.a1, dtype=float).reshape(2)
        self.a2 = np.asarray(self.a2, dtype=float).reshape(2)
        if self.nodes.ndim != 2 or self.nodes.shape[1] != 2:
            raise ValueError("nodes must be an N x 2 array")
        if self.height <= 0:
            raise ValueError("height must be positive")
        if abs(np.linalg.det(np.column_stack((self.a1, self.a2)))) <= np.finfo(float).eps:
            raise ValueError("Periodic lattice vectors are linearly dependent")
        if not self.wall_groups:
            self.wall_groups = ["" for _ in self.elements]
        elif len(self.wall_groups) != len(self.elements):
            raise ValueError("wall_groups must be empty or have one entry per element")


@dataclass(slots=True)
class EffectiveProperties:
    E1: float
    E2: float
    E3: float
    G23: float
    G13: float
    G12: float
    nu12: float
    nu21: float
    nu13: float
    nu31: float
    nu23: float
    nu32: float
    Qplane: np.ndarray
    Qshear: np.ndarray  # Paper order [13, 23]: Q44=G13, Q55=G23
    normal_shear_coupling_norm: float


@dataclass(slots=True)
class PhysicalProperties:
    length_unit: str
    cell_area: float
    rve_volume: float
    solid_volume: float
    void_volume: float
    relative_density: float
    total_mass_kg: float
    weight_n: float
    average_density_kg_m3: float
    effective_solid_density_kg_m3: float
    relative_density_from_mass: float
    has_complete_density: bool
    missing_material_names: list[str]
    material_breakdown: list[dict]
    wall_breakdown: list[dict]
    total_material_cost: float = math.nan
    effective_cost_per_kg: float = math.nan
    cost_per_rve_m3: float = math.nan
    has_complete_cost: bool = False
    missing_cost_material_names: list[str] = field(default_factory=list)


@dataclass(slots=True)
class HomogenizationResult:
    CH: np.ndarray
    S: np.ndarray
    props: EffectiveProperties
    Umacro: np.ndarray
    Lambda: np.ndarray
    K: np.ndarray
    H: np.ndarray
    R: np.ndarray
    constraints_C: np.ndarray
    constraints_Dmacro: np.ndarray
    constraint_info: dict
    sections: list[LaminateSection]
    element_data: list[dict]
    area: float
    volume: float
    physical: PhysicalProperties

    @property
    def relative_density(self) -> float:
        return self.physical.relative_density

    @property
    def mass_kg(self) -> float:
        return self.physical.total_mass_kg
