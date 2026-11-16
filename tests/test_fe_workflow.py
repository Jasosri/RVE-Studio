import json
import math
import os
import stat
import sys
from pathlib import Path

import pytest

from rvestudio import fe_workflow as fw
from rvestudio.cli import main as cli_main

HERE = Path(__file__).parent


def _fake_abaqus(tmp_path: Path) -> str:
    """Executable wrapper around tests/fake_abaqus.py (POSIX only)."""
    exe = tmp_path / 'fake_abaqus'
    exe.write_text(f'#!/bin/sh\nexec "{sys.executable}" "{HERE / "fake_abaqus.py"}" "$@"\n', encoding='utf-8')
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    return str(exe)


def test_parse_config_list():
    assert fw.parse_config_list('1-3,7') == [1, 2, 3, 7]
    assert fw.parse_config_list([9, 1, 1]) == [1, 9]
    with pytest.raises(ValueError):
        fw.parse_config_list('0-2')


def test_richardson_recovers_second_order():
    hs = [10.0, 5.0, 2.5]
    fs = [100 * (1 + 0.04 * (h / 10) ** 2) for h in hs]
    r = fw.richardson(hs, fs)
    assert r['status'] == 'monotone'
    assert r['p'] == pytest.approx(2.0, rel=1e-8)
    assert r['extrapolated'] == pytest.approx(100.0, rel=1e-10)


def test_richardson_flags_oscillation():
    assert fw.richardson([4, 2, 1], [1.0, 1.2, 1.1])['status'] == 'oscillatory'


def test_prepare_is_incremental(tmp_path):
    root = tmp_path / 'study'
    s = fw.FEStudySettings(meshes=[10.0, 5.0], checks=['dbc'])
    st = fw.prepare_study(root, paper_configs='1', settings=s)
    assert set(st['projects']) == {'paper_config_01'}
    assert all(r['prepared'] for r in st['runs']) and len(st['runs']) == 2
    ana = st['projects']['paper_config_01']['analytical']
    assert ana['paper_config'] == 1 and 'published_FE' not in ana and 'paper_proposed' not in ana   # no published values
    assert ana['C']['C1111'] == pytest.approx(ana['CH'][0][0]) and ana['C']['C2323'] == pytest.approx(ana['CH'][3][3])
    assert ana['C21']['C3333'] == pytest.approx(ana['CH'][2][2]) and ana['C21']['C1133'] == pytest.approx(ana['CH'][0][2])
    assert len(ana['C21']) == 21 and len(ana['Q15']) == 15 and set(ana['eng']) >= {'E1', 'E3', 'G13', 'nu12'}
    # Extending adds a mesh and a project without touching existing runs.
    st2 = fw.prepare_study(root, paper_configs='2', settings=fw.FEStudySettings(meshes=[2.5], checks=['dbc']))
    assert len(st2['runs']) == 6
    assert st2['settings']['meshes'] == [10.0, 5.0, 2.5]
    jobs = fw.study_jobs(root)
    assert len(jobs) == 6 * 10 and all(j['state'] == fw.PENDING for j in jobs)   # ten DBC cases per mesh
    assert {j['case'] for j in jobs} == {'Q11', 'Q22', 'Q12', 'Q66', 'Q44', 'Q55', 'C1111', 'C2222', 'C1122', 'C1212'}
    assert (root / 'run_all_jobs.sh').exists() and (root / 'run_all_jobs.bat').exists()
    with pytest.raises(ValueError):
        fw.prepare_study(root, paper_configs='1', settings=fw.FEStudySettings(checks=['dbc'], dbc_displacement=2.0))


def test_u3_restrained_decks_and_21_state_energy_pbc(tmp_path):
    root = tmp_path / 'study'
    st = fw.prepare_study(root, paper_configs='1', settings=fw.FEStudySettings(meshes=[10.0]))
    dbc = next(r for r in st['runs'] if r['check'] == 'dbc')
    deck = (root / dbc['run_dir'] / '07_C1111' / 'DBC_C1111.inp').read_text()
    assert 'U3 = 0 on the top and bottom faces' in deck and 'RBM gauge' not in deck
    shell = (root / dbc['run_dir'] / '01_Q11' / 'DBC_Q11.inp').read_text()
    assert 'top/bottom U3 free' in shell and 'RBM gauge: node' in shell
    pbc = next(r for r in st['runs'] if r['check'] == 'pbc')
    cases = sorted(p.name for p in (root / pbc['run_dir']).iterdir() if p.is_dir())
    assert len(cases) == 21 and cases[:6] == ['01_E11', '02_E22', '03_E33', '04_G23', '05_G13', '06_G12']
    assert '07_E11_E22' in cases and '21_G13_G12' in cases                # the 15 paired states
    assert pbc['pbc_settings']['identification_method'] == 'energy21'


@pytest.mark.skipif(os.name == 'nt', reason='POSIX fake launcher')
def test_end_to_end_with_fake_abaqus(tmp_path):
    root = tmp_path / 'study'
    fw.prepare_study(root, paper_configs='1', settings=fw.FEStudySettings(meshes=[10.0, 5.0, 2.5]))
    exe = _fake_abaqus(tmp_path)
    seen = []
    out = fw.run_study(root, exe, progress=lambda i, n, j: seen.append((i, n)))
    st = fw.study_status(root)
    assert st['counts'] == {fw.DONE: 3 * (10 + 21)} and out['failed'] == 0
    rows = fw.read_rows(root / 'results' / 'abaqus_results.csv')
    assert {(r['model'], r['group']) for r in rows} == {('DBC', 'shell'), ('DBC', '3d'), ('PBC', 'C'), ('PBC', 'Q'),
                                                        ('PBC', 'eng')}
    assert {r['unit'] for r in rows if r['group'] == 'eng' and r['term'].startswith('nu')} == {'-'}
    conv = [c for c in fw.convergence(rows)]
    assert {(c['model'], c['group']) for c in conv} == {(m, g) for m, g in {(r['model'], r['group']) for r in rows}}
    for row in conv:
        # The fake solver scales C by (1 + 0.04 (h/10)^2): stiffness terms converge with p = 2 to the analytical
        # value; Poisson ratios are mesh-independent ("converged").
        assert row['status'] in ('monotone', 'converged'), row
        if row['status'] == 'monotone':
            assert row['p'] == pytest.approx(2.0, rel=1e-6)
        assert abs(row['analytical_vs_FE_pct']) < 1e-6, row
    assert (root / 'results' / 'paper_table.tex').read_text().count(r'\\') >= 3
    html = (root / 'FE_STUDY_REPORT.html').read_text(encoding='utf-8')
    assert 'FE verification study' in html and 'Abaqus DBC: deformed shapes' in html and 'FE RVE' in html
    assert 'Full stiffness matrix C' in html and 'Q matrix' in html and 'Engineering constants' in html
    assert 'DBC shell terms' in html and 'DBC 3D stiffness terms' in html and 'Mesh convergence' in html
    assert 'Abaqus PBC' in html and 'published' not in html.lower()
    assert (root / 'FE_STUDY_RESULTS.xlsx').exists()
    assert any(root.rglob('nodal_displacements.csv'))
    # A second run attempts nothing new.
    again = fw.run_study(root, exe)
    assert again['attempted'] == 0


@pytest.mark.skipif(os.name == 'nt', reason='POSIX fake launcher')
def test_failed_jobs_are_reported_and_the_study_carries_on(tmp_path, monkeypatch):
    root = tmp_path / 'study'
    fw.prepare_study(root, paper_configs='1', settings=fw.FEStudySettings(meshes=[10.0]))
    monkeypatch.setenv('FAKE_ABAQUS_PBC_FAIL', '1')
    out = fw.run_study(root, _fake_abaqus(tmp_path))
    st = fw.study_status(root)
    assert st['counts'] == {fw.DONE: 10, fw.FAILED: 21} and out['failed'] == 21
    html = (root / 'FE_STUDY_REPORT.html').read_text(encoding='utf-8')
    assert 'Problems' in html and 'fake PBC extractor failure' in html


def test_cli_fe_group(tmp_path, capsys):
    root = tmp_path / 's'
    assert cli_main(['fe', 'prepare', '--root', str(root), '--paper-configs', '1', '--models', 'dbc', '--meshes', '10']) == 0
    assert cli_main(['fe', 'status', '--root', str(root)]) == 0
    assert 'PENDING=10' in capsys.readouterr().out
    assert cli_main(['fe', 'check', '--root', str(root)]) == 0
    out = capsys.readouterr().out
    assert 'Q11' in out and 'C1111' in out and 'FE_STUDY_REPORT.html' in out


def test_doctor_without_abaqus():
    rep = fw.abaqus_doctor('definitely_not_abaqus_xyz')
    assert rep['found'] is False and rep['responds'] is False and rep['messages']


def test_builtin_check_writes_both_term_groups_and_the_report(tmp_path):
    root = tmp_path / 'study'
    fw.prepare_study(root, paper_configs='2', settings=fw.FEStudySettings(meshes=[10.0]))
    out = fw.check_study(root)
    assert out['rows'] == (21 + 15 + 12) + (6 + 6) and out['meshes'] == [10.0]
    rows = fw.read_rows(root / 'results' / 'builtin_results.csv')
    assert not any('published' in k for r in rows for k in r)
    by = {(r['model'], r['group'], r['term']): r for r in rows}
    for t in ('C1111', 'C1133', 'C3333', 'C1313'):
        assert abs(by[('PBC', 'C', t)]['analytical_vs_FE_pct']) < 0.5, t
    assert by[('PBC', 'C', 'C1123')]['negligible'] is True                 # 0/90 stack: no 11-23 coupling
    assert abs(by[('PBC', 'eng', 'E3')]['analytical_vs_FE_pct']) < 0.1
    assert abs(by[('DBC', 'shell', 'Q11')]['analytical_vs_FE_pct']) < 1.0
    assert (root / 'results' / 'builtin_fields' / 'paper_config_02.npz').exists()
    f = fw.load_fields(root / 'results' / 'builtin_fields' / 'paper_config_02.npz')
    assert len(f['DBC']['u']) == 10 and len(f['PBC']['u']) == 6
    html = (root / 'FE_STUDY_REPORT.html').read_text(encoding='utf-8')
    assert 'Built-in DBC: deformed shapes' in html and 'built-in PBC' in html and out.get('workbook')
    # Re-checking a mesh replaces its rows instead of duplicating them.
    assert fw.check_study(root)['rows'] == 60


def test_study_records_the_published_rve(tmp_path):
    root = tmp_path / 'study'
    st = fw.prepare_study(root, paper_configs='1', settings=fw.FEStudySettings(meshes=[10.0], checks=['dbc', 'pbc']))
    assert st['settings']['pbc_domain'] == 'published-rve'
    assert st['settings']['dbc_face_walls'] == 'half-thickness'
    dbc = next(r for r in st['runs'] if r['check'] == 'dbc')
    dom = json.loads((root / dbc['run_dir'] / 'thesis_dbc_snapshot.json').read_text())['domain']
    assert dom['primitive_cells'] == 4 and math.isclose(dom['p'], 173.20508075688772) and math.isclose(dom['q'], 150.0)
    pbc = next(r for r in st['runs'] if r['check'] == 'pbc')
    assert pbc['prepared'], pbc['error']


def test_stale_studies_are_refused_for_run_and_collect(tmp_path):
    """Study folders prepared before 1.11 (no U3-restrained cases) or 1.9.1 must not be run or collected."""
    root = tmp_path / 'study'
    fw.prepare_study(root, paper_configs='1', settings=fw.FEStudySettings(meshes=[10.0], checks=['dbc']))
    assert fw.stale_runs(root, fw.load_study(root)[1]) == []
    run_dir = fw.load_study(root)[1]['runs'][0]['run_dir']
    snap_path = root / run_dir / 'thesis_dbc_snapshot.json'
    snap = json.loads(snap_path.read_text())
    snap['cases'] = [c for c in snap['cases'] if not c['name'].startswith('C')]   # what a 1.10 snapshot looks like
    snap_path.write_text(json.dumps(snap))
    for call in (lambda: fw.run_study(root, 'abaqus'), lambda: fw.collect_study(root),
                 lambda: fw.prepare_study(root, paper_configs='1', settings=fw.FEStudySettings(meshes=[5.0], checks=['dbc']))):
        with pytest.raises(ValueError, match='older RVE Studio'):
            call()


def test_studies_from_before_1_13_are_refused(tmp_path):
    root = tmp_path / 'study'
    fw.add_projects(root, paper_configs='1', settings=fw.FEStudySettings(meshes=[10.0]))
    p = root / fw.STUDY_FILE
    st = json.loads(p.read_text())
    st['format'] = 2                                                     # what a 1.11-1.12 study file says
    p.write_text(json.dumps(st))
    with pytest.raises(ValueError, match='older RVE Studio'):
        fw.load_results(root)


def test_rigid_body_gauge_is_minimal_and_fast():
    """One U3 gauge for the free-top in-plane cases, none for the shear and U3-restrained cases."""
    import time
    from rvestudio.examples_data import paper_hexagonal_case
    from rvestudio.abaqus_thesis_dbc import generate_thesis_dbc_input, ThesisDBCSettings
    m, mats, lams = paper_hexagonal_case(1)
    t0 = time.time()
    for case, n_gauges in (('Q11', 1), ('Q66', 1), ('Q44', 0), ('Q55', 0), ('C1111', 0), ('C1212', 0)):
        _deck, meta = generate_thesis_dbc_input(m, mats, lams, ThesisDBCSettings(mesh_size=5.0), case)
        audit = meta['rigid_body_audit']
        assert audit['post_stabilization_rank'] == 6
        assert len(audit['added_gauge_constraints']) == n_gauges
        if n_gauges:
            assert audit['added_gauge_constraints'][0]['dof'] == 3
    assert time.time() - t0 < 30.0


def test_abaqus_wall_orientation_projects_about_the_wall_normal(tmp_path):
    from rvestudio.fe_workflow import prepare_study, _deck_orientation_stale
    prepare_study(tmp_path / 's', paper_configs='1', settings=None) if False else None
    import subprocess, sys
    subprocess.run([sys.executable, '-m', 'rvestudio.cli', 'fe', 'prepare', '--root', str(tmp_path / 's'),
                    '--paper-configs', '1', '--meshes', '10'], check=True, capture_output=True)
    decks = sorted((tmp_path / 's').rglob('*.inp'))
    assert decks
    for inp in decks:
        lines = [l.strip() for l in inp.read_text().splitlines() if not l.startswith('**')]
        idx = [i for i, l in enumerate(lines) if l.upper().startswith('*ORIENTATION')]
        assert idx and all(lines[i + 2].replace(' ', '') == '3,0.' for i in idx), inp
    assert not _deck_orientation_stale(decks[0].parent)
    old = decks[0].read_text().replace('\n3, 0.\n', '\n')
    decks[0].write_text(old)
    assert _deck_orientation_stale(decks[0].parent)


def test_mesh_verdict_from_three_meshes():
    from rvestudio.fe_workflow import richardson, mesh_verdict
    # second-order error, r = 2: in the asymptotic range, GCI small -> converged
    f = lambda h: 10.0 * (1 + 0.001 * h ** 2)
    r = richardson([10, 5, 2.5], [f(10), f(5), f(2.5)])
    assert abs(r['p'] - 2) < 1e-9 and abs(r['asymptotic_ratio'] - 1) < 0.05
    assert mesh_verdict(r) == 'converged'
    # same order, larger error -> refine (GCI)
    g = lambda h: 10.0 * (1 + 0.002 * h ** 2)
    assert mesh_verdict(richardson([10, 5, 2.5], [g(10), g(5), g(2.5)])) == 'refine (GCI)'
    # coarse mesh outside the asymptotic range (higher-order term dominates at 10 mm)
    q = lambda h: 10.0 * (1 + 0.001 * h ** 2 + 0.0004 * h ** 4)
    assert mesh_verdict(richardson([10, 5, 2.5], [q(10), q(5), q(2.5)])) == 'refine (not asymptotic)'
    # oscillating but within tolerance
    r = richardson([10, 5, 2.5], [10.02, 9.99, 10.0])
    assert r['status'] == 'oscillatory' and mesh_verdict(r, 10.0, 9.99) == 'converged (oscillating)'
