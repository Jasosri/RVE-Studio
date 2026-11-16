import numpy as np
import pytest

from rvestudio.examples_data import thesis_hexagonal_case
from rvestudio.fe_rve import build_published_fe_rve
from rvestudio.geometry import BUILTIN_GEOMETRIES
from rvestudio.models import RVEModel, Element
from rvestudio.abaqus_homogenization import (
    PeriodicHomogenizationSettings, PBC_DOMAIN_PUBLISHED, resolve_periodic_domain,
    generate_periodic_homogenization_input,
)
from rvestudio.shell_fe import element_stiffness, shell_section, run_dbc_shell_check
from rvestudio.homogenization import homogenize_rve
from rvestudio.validation import extract_six, PROPERTY_NAMES


def _family(gid, par, h=150.0):
    g = BUILTIN_GEOMETRIES[gid](par)
    return RVEModel(g.nodes, [Element(int(i), int(j), 0, 1.0) for i, j in g.connectivity], g.a1, g.a2, h)


@pytest.mark.parametrize('gid,par,cells', [
    ('HEX', {'L1': 40.0, 'L2': 60.0}, [2, 1]),
    ('SQUARE', {'Lx': 25.0, 'Ly': 30.0}, [2, 2]),
    ('TRI', {'Lh': 40.0, 'R': 1.0}, [2, 1]),
    ('MIXED', {'Lh': 40.0, 'R': 2.0}, [2, 2]),
])
def test_published_rule_holds_whole_cells_and_exact_material(gid, par, cells):
    dbc, d = build_published_fe_rve(_family(gid, par))
    assert d['cells_x1_x2'] == cells
    assert d['primitive_cells'] == 4
    assert d['wall_measure_relative_error'] < 1e-10
    np.testing.assert_allclose(d['rve_area'], 4 * d['primitive_area'], rtol=1e-12)


def test_hexagonal_block_scales_with_wall_lengths():
    _, d = build_published_fe_rve(_family('HEX', {'L1': 40.0, 'L2': 60.0}))
    np.testing.assert_allclose(d['p'], 2 * np.sqrt(3) * 60.0, rtol=1e-12)
    np.testing.assert_allclose(d['q'], 2 * 40.0 + 60.0, rtol=1e-12)
    assert sorted(d['face_wall_faces'].values()) == ['X1_0', 'X1_0', 'X1_P', 'X1_P']


def test_periodic_check_can_use_the_same_published_block():
    model, mats, lams = thesis_hexagonal_case(1)
    st = PeriodicHomogenizationSettings(mesh_size=10.0, domain_mode=PBC_DOMAIN_PUBLISHED)
    dm, dg = resolve_periodic_domain(model, st)
    assert dg['primitive_cell_count'] == 4 and np.isclose(dg['supercell_area'], 2 * 86.60254037844388 * 150)
    deck, meta = generate_periodic_homogenization_input(dm, mats, lams, st, np.array([1e-3, 0, 0, 0, 0, 0]),
                                                        resolved_domain=True, domain_diagnostics=dg)
    assert meta['pbc_audit']['passed'] and meta['pbc_audit']['unpaired_boundary_nodes'] == 0
    assert '_HALF' in deck  # both periodic images of a face wall carry half the moduli


def test_shell_element_passes_membrane_and_bending_patch_tests():
    model, mats, lams = thesis_hexagonal_case(12)  # unbalanced, unsymmetric: A16 and B non-zero
    A, B, D, As = shell_section(lams[0], mats)
    a, b = 7.0, 4.0
    K = element_stiffness(a, b, A, B, D, As, drill=0.0)
    nx, ny = 3, 2
    xy = np.array([[i * a, j * b] for j in range(ny + 1) for i in range(nx + 1)])
    nid = lambda i, j: j * (nx + 1) + i
    Kg = np.zeros((6 * len(xy), 6 * len(xy)))
    for j in range(ny):
        for i in range(nx):
            c = [nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)]
            d = np.array([[6 * n + k for k in range(6)] for n in c]).ravel()
            Kg[np.ix_(d, d)] += K
    eps, kap = np.array([1e-3, 2e-4, -3e-4]), np.array([1e-4, 2e-4, -1e-4])
    def field(x, y):
        bx = kap[0] * x + 0.5 * kap[2] * y; by = kap[1] * y + 0.5 * kap[2] * x
        w = -(0.5 * kap[0] * x * x + 0.5 * kap[1] * y * y + 0.5 * kap[2] * x * y)
        return [eps[0] * x + 0.5 * eps[2] * y, eps[1] * y + 0.5 * eps[2] * x, w, -by, bx, 0.0]
    U = np.concatenate([field(x, y) for x, y in xy])
    bnd = [k for k, (x, y) in enumerate(xy) if x in (0, nx * a) or y in (0, ny * b)]
    pres = sorted({6 * k + i for k in bnd for i in range(5)} | {6 * k + 5 for k in range(len(xy))})
    free = [i for i in range(len(U)) if i not in pres]
    u = U.copy(); u[free] = np.linalg.solve(Kg[np.ix_(free, free)], -Kg[np.ix_(free, pres)] @ U[pres])
    e = np.concatenate((eps, kap))
    W_exact = 0.5 * e @ np.block([[A, B], [B, D]]) @ e * (nx * a * ny * b)
    assert abs(0.5 * u @ Kg @ u / W_exact - 1) < 1e-10
    assert np.max(np.abs(u - U)) < 1e-12


def test_builtin_shell_check_reproduces_published_fe_for_a_cross_ply_core():
    # Config 2 ([0,0,90,0,0] E-glass): the published FE values (thesis Tables 3.5-3.10).
    model, mats, lams = thesis_hexagonal_case(2)
    r = run_dbc_shell_check(model, mats, lams, mesh_size=10.0)
    published = {'Q11': 591.66, 'Q22': 591.64, 'Q12': 584.02, 'Q44': 149.54, 'Q55': 149.44}
    for k, v in published.items():
        assert abs(r['Q'][k] / v - 1) < 0.01, (k, r['Q'][k], v)
    # Antiplane shear of a regular honeycomb under these DBCs is exactly affine, so
    # Q44 = Q55 = the analytical value - this is the check that face walls carry half.
    ana = dict(zip(PROPERTY_NAMES, extract_six(homogenize_rve(model, mats, lams))))
    assert abs(r['Q']['Q44'] / ana['Q44'] - 1) < 1e-6 and abs(r['Q']['Q55'] / ana['Q55'] - 1) < 1e-6
    # Full-thickness face walls break junction equilibrium of that exact field and
    # overshoot Q55 by about 10 %, which the published values rule out.
    full = run_dbc_shell_check(model, mats, lams, mesh_size=10.0, face_wall_section='full', cases=('Q55',))
    assert full['energies_Nmm']['Q55'] / r['energies_Nmm']['Q55'] > 1.08
