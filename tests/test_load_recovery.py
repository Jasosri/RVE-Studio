import numpy as np

from rvestudio.examples_data import paper_hexagonal_case
from rvestudio.homogenization import homogenize_rve
from rvestudio.load_recovery import recover_wall_loads


def _result(config=1):
    model, mats, lams = paper_hexagonal_case(config)
    return model, homogenize_rve(model, mats, lams)


def test_wall_load_recovery_internal_identities():
    model, out = _result(1)
    E = np.array([-1.2e-3, -0.7e-3, -0.25e-3, 0.18e-3, -0.12e-3, 0.4e-3])
    rec = recover_wall_loads(model, out, macro_strain=E)
    assert len(rec.walls) == len(model.elements)
    assert rec.validation['equilibrium_relative_residual'] < 1e-9
    assert rec.validation['energy_relative_error'] < 1e-10
    assert rec.validation['macro_work_relative_error'] < 1e-10
    assert rec.validation['element_energy_sum_relative_error'] < 1e-10
    assert rec.validation['max_membrane_force_identity_error'] < 1e-10
    assert rec.validation['max_N33_conjugate_identity_error'] < 1e-10
    assert rec.validation['constraint_residual_norm'] < 1e-9


def test_stress_control_recovers_requested_macro_stress():
    model, out = _result(1)
    sigma = np.array([-2.0, -1.0, -0.2, 0.15, -0.10, 0.5])
    rec = recover_wall_loads(model, out, macro_stress=sigma)
    assert np.allclose(rec.macro_stress, sigma, rtol=1e-9, atol=1e-9)
    assert rec.validation['stress_control_relative_error'] < 1e-9
    assert any(abs(w['N11_N_per_mm']) > 1e-10 for w in rec.walls)
    assert any(abs(w['N33_N_per_mm']) > 1e-10 for w in rec.walls)
    assert any(abs(w['N13_N_per_mm']) > 1e-10 for w in rec.walls)


def test_wall_load_recovery_is_linear_in_macro_strain():
    model, out = _result(1)
    Ea = np.array([-1e-3, 0, 0, 0, 0, 0.25e-3])
    Eb = np.array([0, -0.5e-3, -0.2e-3, 0, 0.1e-3, 0])
    ra = recover_wall_loads(model, out, macro_strain=Ea)
    rb = recover_wall_loads(model, out, macro_strain=Eb)
    rab = recover_wall_loads(model, out, macro_strain=Ea+Eb)
    for wa, wb, wab in zip(ra.walls, rb.walls, rab.walls):
        for key in ('N11_N_per_mm','N33_N_per_mm','N13_N_per_mm'):
            assert np.isclose(wab[key], wa[key]+wb[key], rtol=1e-10, atol=1e-10)
