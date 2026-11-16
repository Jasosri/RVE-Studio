"""Automatic FE RVE (v1.10): one rule for every topology, used by every FE path."""
import math

import numpy as np
import pytest

from rvestudio.geometry import BUILTIN_GEOMETRIES, DEFAULT_GEOMETRY_LIBRARY
from rvestudio.models import RVEModel, Element
from rvestudio.fe_rve import select_fe_rve, fe_rve_summary, find_cells, SLIVER_RATIO
from rvestudio.examples_data import paper_hexagonal_case
from rvestudio.abaqus_homogenization import (PeriodicHomogenizationSettings, resolve_periodic_domain,
                                             PBC_DOMAIN_PUBLISHED, PBC_DOMAIN_LATTICE_SUPERCELL)


def _model(gid, par, h=150.0):
    g = BUILTIN_GEOMETRIES[gid](par)
    els = [Element(int(i), int(j), 0, float(w)) for (i, j), w in zip(g.connectivity, g.wall_weights)]
    return RVEModel(g.nodes, els, g.a1, g.a2, h, name=gid, wall_groups=g.wall_groups)


@pytest.mark.parametrize('gid,par,cells,cell_name', [
    ('HEX', {'L1': 50, 'L2': 50}, [2, 1], 'hexagon'),
    ('SQUARE', {'Lx': 50, 'Ly': 50}, [2, 2], 'square'),
    ('SQUARE', {'Lx': 40, 'Ly': 60}, [2, 2], 'rectangle'),
    ('TRI', {'Lh': 50, 'R': 1.0}, [2, 1], 'triangle'),
    ('MIXED', {'Lh': 40, 'R': 2.0}, [2, 2], 'rhombus'),
])
def test_every_family_gets_a_four_cell_block_centred_on_a_complete_cell(gid, par, cells, cell_name):
    dm, d = select_fe_rve(_model(gid, par))
    assert d['cells_x1_x2'] == cells
    assert d['primitive_cells'] == 4
    assert d['wall_measure_relative_error'] < 1e-7
    assert d['central_cell']['name'] == cell_name
    assert d['min_wall_piece_ratio'] >= SLIVER_RATIO
    assert d['selection']['automatic'] and not d['selection']['warnings']
    P = np.asarray(d['central_cell_polygon_rve'])
    assert P[:, 0].min() > 0 and P[:, 0].max() < d['p'] and P[:, 1].min() > 0 and P[:, 1].max() < d['q']
    assert len(d['wall_segments']) == len(dm.elements)
    assert f'complete {cell_name}' in fe_rve_summary(d)


def test_published_honeycomb_block_is_unchanged_and_symmetric():
    m, _, _ = paper_hexagonal_case(1)
    _, d = select_fe_rve(m)
    assert d['p'] == pytest.approx(2 * 86.60254037844386) and d['q'] == pytest.approx(150.0)
    assert d['face_walls'] == [1, 9, 10, 16]
    assert d['mirror_symmetry'] == {'about_x1_mid': True, 'about_x2_mid': True}
    assert d['effective_wall_length'] == pytest.approx(600.0)


def test_random_designs_of_every_family_all_get_an_rve():
    rng = np.random.default_rng(7)
    for gd in DEFAULT_GEOMETRY_LIBRARY:
        for _ in range(6):
            par = {p['name']: float(rng.uniform(p['min'], p['max'])) for p in gd['parameters']}
            _, d = select_fe_rve(_model(gd['id'], par))
            assert d['wall_measure_relative_error'] < 1e-7
            assert d['min_wall_piece_ratio'] >= SLIVER_RATIO


def test_overrides_block_size_and_centre_cell():
    m = _model('MIXED', {'Lh': 40, 'R': 2.0})
    _, auto = select_fe_rve(m)
    _, big = select_fe_rve(m, cells=(3, 3))
    assert big['cells_x1_x2'] == [3, 3] and big['primitive_cells'] == 9
    assert not big['selection']['automatic']
    names = [c['name'] for c in auto['cell_types']]
    k = names.index('triangle')
    _, tri = select_fe_rve(m, centre_cell=k)
    assert tri['central_cell']['name'] == 'triangle'
    with pytest.raises(ValueError):
        select_fe_rve(m, centre_cell=99)


def test_cell_types_are_found_once_each():
    m = _model('MIXED', {'Lh': 40, 'R': 2.0})
    from rvestudio.fe_rve import _rotation_to_a1
    R, a1, a2 = _rotation_to_a1(m)
    nodes = (R @ np.asarray(m.nodes, float).T).T
    cells = find_cells(nodes, list(m.elements), a1, a2)
    areas = sorted(c['area'] for c in cells)
    prim = abs(float(np.linalg.det(np.column_stack((a1, a2)))))
    assert sum(areas) == pytest.approx(prim, rel=1e-9)   # the cell types tile the primitive cell
    assert cells[0]['inradius'] >= cells[-1]['inradius']


def test_periodic_default_is_the_automatic_rve_with_lattice_fallback():
    assert PeriodicHomogenizationSettings().domain_mode == PBC_DOMAIN_PUBLISHED
    m, _, _ = paper_hexagonal_case(1)
    dom, d = resolve_periodic_domain(m, PeriodicHomogenizationSettings())
    assert d['domain_mode'] == PBC_DOMAIN_PUBLISHED and d['primitive_cell_count'] == 4
    # An oblique lattice with an irrational aspect has no rectangular periodic cell.
    a1 = np.array([10.0, 0.0]); a2 = np.array([math.sqrt(2.0), 7.3])
    nodes = np.array([[0.0, 0.0], [10.0, 0.0], [math.sqrt(2.0), 7.3]])
    obl = RVEModel(nodes, [Element(0, 1, 0, 1.0), Element(0, 2, 0, 1.0)], a1, a2, 20.0, name='oblique')
    with pytest.raises(ValueError):
        select_fe_rve(obl)
    dom, d = resolve_periodic_domain(obl, PeriodicHomogenizationSettings(supercell_repeat_x=1, supercell_repeat_y=1))
    assert d['fallback_from'] == PBC_DOMAIN_PUBLISHED and d['domain_mode'] == PBC_DOMAIN_LATTICE_SUPERCELL
