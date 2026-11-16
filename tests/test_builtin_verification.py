"""Built-in shell checks (no Abaqus), the one-design verification and its report."""
import numpy as np
import pytest

from rvestudio.examples_data import paper_hexagonal_case
from rvestudio.homogenization import homogenize_rve
from rvestudio.shell_fe import run_pbc_shell_check, run_dbc_shell_check
from rvestudio.verification import (verify_model, VerificationSettings, paper_configuration_of, signed,
                                    summary_table, solid_terms, SHELL_TERMS, SOLID_TERMS, pbc_values, q_matrix,
                                    full_c_terms, matrix_comparison, ENG_TERMS)
from rvestudio.abaqus_homogenization import (homogenization_cases, reconstruct_ch_from_energies,
                                             plane_stress_shell_reduction)
from rvestudio.fe_report import write_verification_report, write_verification_workbook
from rvestudio.cli import main as cli_main


def test_periodic_shell_check_gives_the_full_c_from_strain_energy():
    m, mats, lams = paper_hexagonal_case(1)
    A = homogenize_rve(m, mats, lams).CH
    r = run_pbc_shell_check(m, mats, lams, mesh_size=10.0)
    C = np.asarray(r['CH'], float)
    assert np.isfinite(C).all() and r['symmetry_error'] < 1e-10   # all six states, incl. E33
    for i, j in ((2, 2), (0, 2), (1, 2)):                         # through-height row: exact
        assert C[i, j] == pytest.approx(A[i, j], rel=1e-6)
    assert np.all(np.linalg.eigvalsh(C) > 0)
    for i in (3, 4):                                              # antiplane shear: exact (affine fields)
        assert C[i, i] == pytest.approx(A[i, i], rel=1e-6)
    assert C[0, 0] == pytest.approx(A[0, 0], rel=1e-3)
    assert C[0, 1] == pytest.approx(A[0, 1], rel=1e-3)
    assert C[5, 5] == pytest.approx(A[5, 5], rel=0.04)           # wall bending, 10 mm mesh
    assert set(r['C']) == set(SOLID_TERMS) and r['periodic_pairs'] > 0 and r['dofs_reduced'] > 0


def test_mutual_energies_equal_the_21_state_energy_identification():
    """Built-in C_ij = u_i K u_j / V is what Abaqus gets from the 21 ALLSE energies (W_ij - W_i - W_j)."""
    m, mats, lams = paper_hexagonal_case(7)
    C = np.asarray(run_pbc_shell_check(m, mats, lams, mesh_size=10.0)['CH'], float)
    a, V = 1e-3, 1000.0
    cases = [{'components': c['components'],
              'strain_energy': 0.5 * V * float(np.asarray(c['strain']) @ C @ np.asarray(c['strain']))}
             for c in homogenization_cases(a)]
    assert len(cases) == 21
    assert np.allclose(reconstruct_ch_from_energies(cases, V, a), C, rtol=1e-9, atol=1e-9)


def test_q_and_engineering_constants_follow_from_c():
    m, mats, lams = paper_hexagonal_case(7)
    CH = homogenize_rve(m, mats, lams).CH
    v = pbc_values(CH)
    Q5 = np.asarray(v['Qmatrix'])
    assert np.allclose(Q5, plane_stress_shell_reduction(CH)[0])
    assert v['Q15']['Q44'] == pytest.approx(Q5[3, 3]) and v['Q15']['Q55'] == pytest.approx(Q5[4, 4])   # 13, 23 shear
    assert v['Q15']['Q66'] == pytest.approx(Q5[2, 2]) and v['Q15']['Q16'] == pytest.approx(Q5[0, 2])
    S = np.linalg.inv(CH)
    assert v['eng']['E1'] == pytest.approx(1 / S[0, 0]) and v['eng']['G13'] == pytest.approx(1 / S[4, 4])
    assert v['eng']['nu13'] == pytest.approx(-S[2, 0] / S[0, 0]) and list(v['eng']) == list(ENG_TERMS)
    assert full_c_terms(CH)['C2313'] == pytest.approx(CH[3, 4])
    assert np.allclose(q_matrix(CH), Q5)


def test_dbc_shell_and_3d_terms_follow_the_top_bottom_conditions():
    """U3 free gives the reduced shell terms, U3 = 0 the eps33 = 0 terms (configuration 1: Q11 324, C1111 391)."""
    m, mats, lams = paper_hexagonal_case(1)
    res = homogenize_rve(m, mats, lams)
    C = solid_terms(res.CH)
    r = run_dbc_shell_check(m, mats, lams, mesh_size=10.0)
    assert r['Q']['Q11'] == pytest.approx(324.4, rel=2e-3) and r['C']['C1111'] == pytest.approx(C['C1111'], rel=2e-3)
    assert r['C']['C1122'] == pytest.approx(C['C1122'], rel=2e-3)
    assert r['C']['C1313'] == r['Q']['Q44'] and r['C']['C2323'] == r['Q']['Q55']
    assert r['C']['C1111'] / r['Q']['Q11'] > 1.15                   # Poisson contraction through the height


def test_verify_model_compares_every_source(tmp_path):
    m, mats, lams = paper_hexagonal_case(1)
    v = verify_model(m, mats, lams, VerificationSettings(mesh_size=10.0, return_fields=True))
    assert 'published_FE' not in v and not any('published' in k for r in v['rows'] for k in r)
    cols, rows, _notes = summary_table(v['rows'], v['project'], 'shell', v['analytical'])
    assert [r[0] for r in rows] == list(SHELL_TERMS) and 'built-in DBC' in cols and not any('published' in c for c in cols)
    cols3, rows3, _n = summary_table(v['rows'], v['project'], '3d', v['analytical'])
    assert [r[0] for r in rows3] == list(SOLID_TERMS) and 'built-in DBC (U3 = 0)' in cols3
    colsc, rowsc, notes = summary_table(v['rows'], v['project'], 'C', v['analytical'])
    assert 'built-in PBC' in colsc and {'C1111', 'C3333', 'C1133', 'C1212'} <= {r[0] for r in rowsc}
    assert 'C1123' not in {r[0] for r in rowsc} and '≈ 0' in notes                 # zero couplings not listed
    by = {(r['model'], r['group'], r['term']): r for r in v['rows']}
    assert abs(by[('PBC', 'C', 'C1111')]['analytical_vs_FE_pct']) < 0.5
    assert abs(by[('PBC', 'C', 'C1313')]['analytical_vs_FE_pct']) < 1e-6
    assert abs(by[('PBC', 'eng', 'G13')]['analytical_vs_FE_pct']) < 1e-6
    assert by[('PBC', 'eng', 'nu12')]['unit'] == '-' and by[('PBC', 'Q', 'Q11')]['unit'] == 'MPa'
    assert abs(by[('DBC', 'shell', 'Q11')]['analytical_vs_FE_pct']) < 1.0
    mc = matrix_comparison(v['rows'], v['project'], 'C', v['analytical'])
    assert mc['solver'] == 'built-in' and mc['FE'].shape == (6, 6) and np.isnan(mc['difference_pct'][0, 3])
    assert any('Full C' in f for f in v['findings']) and v['benchmark_configuration'] == 1
    assert len(v['fe']['fields']['DBC']['u']) == 10 and len(v['fe']['fields']['PBC']['u']) == 6
    html = write_verification_report(v, tmp_path / 'v.html').read_text(encoding='utf-8')
    assert 'RVE verification report' in html and 'Findings' in html and 'data:image/png;base64,' in html
    assert 'Full stiffness matrix C' in html and 'Analytical vs FE (%)' in html and 'U3 = 0 on both faces' in html
    assert 'published' not in html.lower()
    assert '>-0.00<' not in html  # round-off differences print as +0.00
    import openpyxl
    wb = openpyxl.load_workbook(write_verification_workbook(v, tmp_path / 'v.xlsx'))
    assert {'Summary', 'C matrix', 'Q matrix', 'Engineering constants', 'DBC shell terms', 'DBC 3D terms',
            'All results'} <= set(wb.sheetnames)


def test_benchmark_configurations_are_recognised_by_content():
    m, mats, lams = paper_hexagonal_case(9)
    m.name = 'RVE'
    assert paper_configuration_of(m, mats, lams) == 9
    m.height = 120.0
    assert paper_configuration_of(m, mats, lams) is None


def test_cli_fe_rve(tmp_path, capsys):
    assert cli_main(['fe', 'rve', '--paper-config', '1', '--png', str(tmp_path / 'rve.png')]) == 0
    assert (tmp_path / 'rve.png').exists()
    assert '173.2 x 150.0' in capsys.readouterr().out


def test_signed_never_prints_negative_zero():
    assert [signed(x) for x in (-2.6e-7, 0.0, -0.004, 0.006, -10.111)] == ['+0.00', '+0.00', '+0.00', '+0.01', '-10.11']
    assert signed(-0.04, 1) == '+0.0' and signed(float('nan')) == '–'
