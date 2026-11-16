"""RVE Studio geometry plug-in: Kagome (trihexagonal) cell.

Corner-sharing triangles around hexagonal openings; every node joins four walls. L1 is the length of the walls
parallel to a1 and L2 the length of the inclined walls (L1 = L2: regular Kagome, in-plane isotropic and
stretch-dominated, E = rho* E_wall / 3 and nu12 = 1/3 for slender walls).
"""
import math
import numpy as np
from rvestudio.geometry import Geometry

GEOMETRY_ID = 'KAGOME'


def generate(p):
    L1, L2 = float(p.get('L1', 50.0)), float(p.get('L2', 50.0))
    if not 0 < L1 < 2 * L2:
        raise ValueError('Kagome cell needs 0 < L1 < 2 L2')
    h = math.sqrt(L2 * L2 - 0.25 * L1 * L1)          # height of each triangle
    nodes = np.array([[0, 0], [L1, 0], [0.5 * L1, h],  # A, B, C: one triangle
                      [2 * L1, 0], [L1, 2 * h], [0, 2 * h]], float)   # images A + a1, A + a2, B + a2 - a1
    walls = np.array([[0, 1], [1, 3],                  # walls along a1 (I -> J in +x1)
                      [0, 2], [2, 4],                  # walls along +60 deg
                      [1, 2], [2, 5]])                 # walls along +120 deg
    return Geometry(GEOMETRY_ID, 'Kagome (trihexagonal)', nodes, walls,
                    np.array([2 * L1, 0.0]), np.array([L1, 2 * h]),
                    ['horizontal-L1-A', 'horizontal-L1-B', 'inclined-L2-A', 'inclined-L2-B', 'inclined-L2-C', 'inclined-L2-D'],
                    np.ones(6), 'L1 = L2: regular Kagome', L1, L2)
