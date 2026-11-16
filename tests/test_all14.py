import numpy as np
import pytest
from rvestudio.validation import PAPER_PUBLISHED, PAPER_URL, run_all14_benchmark


def test_paper_reference_comparison_is_transparent():
    b = run_all14_benchmark(tolerance_pct=1.0)
    assert b['reference'] == PAPER_URL
    assert len(b['detail']) == 84
    assert b['total'] == 84
    # Euler-Bernoulli walls (the published model): all 84 values within 0.3 %.
    assert b['within_tolerance'] == 84 and b['passed']
    assert b['max_abs_paper_pct'] < 0.3
    assert list(b['summary']['Property']) == ['Q11','Q22','Q12','Q66','Q44','Q55']
    assert all(b['summary']['Within_tolerance'].to_numpy() == 14)
    assert np.allclose(b['paper_published'], PAPER_PUBLISHED)
    # Changing a reporting threshold must not alter any physical calculation.
    strict = run_all14_benchmark(tolerance_pct=0.1)
    assert not strict['passed'] and np.allclose(strict['detail']['Calculated_MPa'], b['detail']['Calculated_MPa'])


def test_benchmark_rejects_invalid_thresholds():
    for bad in [0.0,-3.0,float('nan'),float('inf')]:
        with pytest.raises(ValueError):
            run_all14_benchmark(tolerance_pct=bad)


def test_timoshenko_walls_change_only_q66_and_move_it_towards_the_published_fe():
    from rvestudio.examples_data import paper_hexagonal_case
    from rvestudio.homogenization import homogenize_rve
    from rvestudio.validation import extract_six
    from rvestudio.abaqus_homogenization import _THESIS_FE_Q
    err_eb, err_t = [], []
    for n in range(1, 15):
        m, mats, lams = paper_hexagonal_case(n)
        assert m.wall_shear == 'euler-bernoulli'
        eb = extract_six(homogenize_rve(m, mats, lams))
        m.wall_shear = 'timoshenko'
        t = extract_six(homogenize_rve(m, mats, lams))
        for k in (0, 1, 2, 4, 5):
            assert abs(t[k] / eb[k] - 1) < 0.003
        assert t[3] < eb[3]                                   # shear flexibility softens wall bending
        fe = _THESIS_FE_Q['Q66'][n - 1]
        err_eb.append(abs(eb[3] / fe - 1)); err_t.append(abs(t[3] / fe - 1))
    assert np.mean(err_t) < np.mean(err_eb)                   # 8.6 % against 11.4 %


def test_wall_beam_is_saved_with_the_project(tmp_path):
    from rvestudio.examples_data import paper_hexagonal_case
    from rvestudio.io import save_project, load_project
    m, mats, lams = paper_hexagonal_case(1)
    m.wall_shear = 'timoshenko'
    save_project(tmp_path / 'p.rvejson', m, mats, lams)
    assert load_project(tmp_path / 'p.rvejson')[0].wall_shear == 'timoshenko'


def test_benchmark_wall_beam_and_bending_options_without_published_fe():
    from rvestudio.validation import run_all14_benchmark
    eb = run_all14_benchmark(1.0)
    tim = run_all14_benchmark(1.0, wall_shear='timoshenko')
    tc = run_all14_benchmark(1.0, wall_shear='timoshenko', bending='twist-corrected')
    free = run_all14_benchmark(1.0, bending='free moment')
    assert eb['wall_shear'] == 'euler-bernoulli' and eb['bending'] == 'cylindrical' and eb['within_tolerance'] == 84
    assert tim['wall_shear'] == 'timoshenko' and tim['within_tolerance'] == 70
    assert tc['bending'] == 'twist-corrected' and free['bending'] == 'free moment'
    q66 = lambda out: out['calculated'][:, 3]
    assert (q66(tc) <= q66(tim) + 1e-12).all() and (q66(free) <= q66(eb) + 1e-12).all()
    assert not any('FE' in c for c in eb['detail'].columns) and 'error_vs_fe_pct' not in eb


def test_generate_db_writes_the_chosen_wall_beam(tmp_path):
    import pandas as pd
    from rvestudio.database import generate_parametric_database
    out = tmp_path / 'db.csv'
    generate_parametric_database(out, 4, 1, wall_shear='timoshenko')
    assert set(pd.read_csv(out)['WallShear']) == {'timoshenko'}
    generate_parametric_database(out, 4, 1)
    assert set(pd.read_csv(out)['WallShear']) == {'euler-bernoulli'}


def test_generate_db_writes_the_chosen_bending_option_and_twist_coefficients(tmp_path):
    import pandas as pd
    from rvestudio.database import generate_parametric_database, build_case_from_row
    out = tmp_path / 'db.csv'
    df = generate_parametric_database(out, 4, 1, wall_shear='timoshenko', bending='twist-corrected', twist_release=(3.0, 2.0, 0.7))
    d = pd.read_csv(out)
    assert set(d['StripMode']) == {'twistRelease'} and set(d['SolverType']) == {'anisotropicStrip'}
    model = build_case_from_row(d.iloc[0].to_dict())[0]
    assert model.twist_release == (3.0, 2.0, 0.7) and model.strip_bending_mode == 'twistRelease'
