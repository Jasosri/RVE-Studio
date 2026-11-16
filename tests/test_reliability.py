from rvestudio.api import solve_example, summary
from rvestudio.reliability import LEVELS, assess


def test_cross_ply_symmetric_wall_is_reliable_with_shear_flexible_wall():
    r = assess(solve_example(2), wall_shear='timoshenko')
    assert r.in_plane == LEVELS[0] and r.out_of_plane == LEVELS[0]
    assert abs(r.walls[0].psi) < 1e-12 and r.walls[0].chi_B < 1e-12


def test_euler_bernoulli_wall_flags_shear_parameter():
    r = assess(solve_example(9), wall_shear='euler-bernoulli')
    assert r.walls[0].beta_max > 0.1
    assert r.in_plane != LEVELS[0]
    assert any('Timoshenko' in m for m in r.messages)


def test_strong_bend_twist_is_outside_verified_range_but_out_of_plane_is_reliable():
    r = assess(solve_example(7), wall_shear='timoshenko')
    assert abs(r.walls[0].psi - 0.715) < 1e-3
    assert r.in_plane == LEVELS[2] and r.out_of_plane == LEVELS[0]


def test_unsymmetric_wall_reports_membrane_bending_measure():
    r = assess(solve_example(13), wall_shear='timoshenko')
    assert abs(r.walls[0].chi_B - 0.0439) < 1e-3
    assert r.out_of_plane != LEVELS[0]


def test_summary_carries_reliability():
    s = summary(solve_example(1))
    assert s['reliability']['walls'][0]['psi'] > 0.3
    assert 'limits' in s['reliability']


def _solve(config, wall_type, mode, shear='timoshenko'):
    import copy
    from rvestudio.examples_data import paper_hexagonal_case
    from rvestudio.homogenization import homogenize_rve
    m, mats, lams = paper_hexagonal_case(config)
    m = copy.copy(m); m.wall_element_type = wall_type; m.strip_bending_mode = mode; m.wall_shear = shear
    return homogenize_rve(m, mats, lams)


def test_twist_corrected_option_reduces_bending_stiffness_by_calibrated_release():
    import math
    from rvestudio.laminate import twist_length, twist_release
    cyl = _solve(7, 'paper', 'cylindrical'); tc = _solve(7, 'anisotropicStrip', 'twistRelease')
    sec = cyl.sections[0]; D = sec.D; psi = D[0, 2] / math.sqrt(D[0, 0] * D[2, 2])
    eta = twist_release(psi, twist_length(sec, 50.0))
    assert abs(tc.element_data[0]['strip_meta']['Dstrip'] / D[0, 0] - (1 - eta * psi ** 2)) < 2e-3
    assert tc.props.G12 < cyl.props.G12
    r = assess(tc)
    assert r.option == 'twist-corrected' and r.in_plane == LEVELS[0]      # psi = 0.715 <= 0.85, beta = 0.055
    r1 = assess(_solve(1, 'anisotropicStrip', 'twistRelease'))
    assert r1.in_plane == LEVELS[0]                                       # psi = 0.325
    assert abs(twist_release(0.5, 0.03) - (2.0 + 2.15 * 0.25) * 0.03 ** 0.73) < 1e-12   # no reference configuration


def test_bound_options_are_not_rated_reliable():
    r = assess(_solve(2, 'anisotropicStrip', 'freeMoment'))
    assert r.in_plane == LEVELS[2] and r.option == 'free moment'


def test_bend_twist_index_is_eta_psi_squared():
    from rvestudio.laminate import twist_release
    from rvestudio.reliability import PI_GOOD
    w = assess(solve_example(1), wall_shear='timoshenko').walls[0]
    assert abs(w.twist_release - twist_release(w.psi, w.twist_length)) < 1e-12
    assert abs(w.bend_twist_index - w.twist_release * w.psi ** 2) < 1e-12
    assert w.bend_twist_index < PI_GOOD                                   # config 1: about 2 % at most
    w7 = assess(solve_example(7), wall_shear='timoshenko').walls[0]
    assert w7.bend_twist_index > 0.04 and w7.level == 'outside the verified range'


def test_user_twist_release_coefficients_change_stiffness_and_rating(tmp_path):
    import copy
    from rvestudio.examples_data import paper_hexagonal_case
    from rvestudio.homogenization import homogenize_rve
    from rvestudio.io import save_project, load_project
    m, mats, lams = paper_hexagonal_case(7)
    m = copy.copy(m); m.wall_element_type = 'anisotropicStrip'; m.strip_bending_mode = 'twistRelease'; m.wall_shear = 'timoshenko'
    base = homogenize_rve(m, mats, lams)
    m.twist_release = (4.0, 2.15, 0.73)
    more = homogenize_rve(m, mats, lams)
    assert more.props.G12 < base.props.G12                    # more twist released -> softer
    r = assess(more)
    assert r.in_plane == LEVELS[1] and any('user-defined' in n for n in r.messages)
    f = tmp_path / 'p.rvejson'; save_project(f, m, mats, lams)
    assert tuple(load_project(f)[0].twist_release) == (4.0, 2.15, 0.73)


def test_euler_bernoulli_cylindrical_wall_is_graded_on_pi_plus_beta():
    """Wall shear and bend-twist errors add for the Euler-Bernoulli cylindrical wall (Paper A, Table 19)."""
    import copy
    from rvestudio.examples_data import paper_hexagonal_case
    from rvestudio.homogenization import homogenize_rve
    from rvestudio.reliability import assess
    m, mats, lams = paper_hexagonal_case(1)      # Pi = 0.0178, beta = 0.0197: each alone below 0.02
    rel = assess(homogenize_rve(m, mats, lams))
    w = rel.walls[0]
    assert w.bend_twist_index < 0.018 and w.beta_max < 0.02
    assert rel.in_plane == "use with care" and any("beta + Pi" in n for n in w.notes)
    mt = copy.deepcopy(m); mt.wall_shear = "timoshenko"     # shear-flexible wall: Pi alone
    assert assess(homogenize_rve(mt, mats, lams)).in_plane == "reliable"
    m2, mats2, lams2 = paper_hexagonal_case(2)   # psi = 0, beta = 0.044 -> use with care
    assert assess(homogenize_rve(m2, mats2, lams2)).in_plane == "use with care"
