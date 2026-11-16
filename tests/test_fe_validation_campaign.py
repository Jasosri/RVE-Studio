import csv
import json
from pathlib import Path

import numpy as np
import pytest

from rvestudio.abaqus_homogenization import (
    PeriodicHomogenizationSettings, _build_shell_mesh,
    periodic_mesh_pairs, _fractional_coordinates, _mod1,
)
from rvestudio.fe_validation import prepare_campaign, collect_campaign, _target_metrics
from rvestudio.io import load_project
from rvestudio.homogenization import homogenize_rve

PROJECT = Path(__file__).resolve().parents[1] / 'examples' / 'fe_validation' / 'P00050.rvejson'
TARGET = PROJECT.with_name('C_target_MPa.csv')


def test_optimised_project_legacy_model_and_cost_metadata():
    m, materials, laminates = load_project(PROJECT)
    assert m.wall_element_type == 'paper'
    assert materials['Boron/Epoxy'].cost_per_kg == 3.0
    assert np.isclose(homogenize_rve(m, materials, laminates).CH[0, 0], 1469.7267779)


def test_fast_periodic_pairing_matches_brute_force_small_model():
    m, _, _ = load_project(PROJECT)
    nodes, _, _, _ = _build_shell_mesh(m, 15.0)
    fast = periodic_mesh_pairs(nodes, m)
    ids = sorted(nodes)
    xyz = np.asarray([nodes[i] for i in ids])
    frac = _fractional_coordinates(xyz, m)
    mod = _mod1(frac, 1e-8)
    assigned = np.zeros(len(ids), bool)
    brute = []
    for i in range(len(ids)):
        if assigned[i]:
            continue
        assigned[i] = True
        for j in range(i+1, len(ids)):
            if assigned[j]:
                continue
            delta = np.abs(mod[i]-mod[j]); delta = np.minimum(delta, 1.-delta)
            if np.max(delta) <= 1e-8 and np.max(np.abs(frac[j]-frac[i]-np.rint(frac[j]-frac[i]))) <= 1e-7:
                assigned[j] = True
                dv = xyz[j]-xyz[i]
                if np.linalg.norm(dv) > 1e-12:
                    brute.append((ids[i], ids[j]))
    assert {(a,b) for a,b,_ in fast} == set(brute)


def test_campaign_does_not_fabricate_fe_results(tmp_path):
    campaign = tmp_path / 'fe'
    m = prepare_campaign(campaign, [PROJECT], mesh_sizes=[20.], target_path=TARGET,
                         domain_mode='primitive')
    assert len(m['records']) == 1
    rec = m['records'][0]
    assert rec['pbc_pass']
    run = campaign / rec['run_dir']
    assert len(list(run.glob('*/*.inp'))) == 6
    assert not (run / 'homogenized_results.json').exists()
    rows = collect_campaign(campaign)
    assert rows[0]['Status'] == 'PENDING'
    assert not any('C11_FE_MPa' in r for r in rows)
    assert '0/1' in (campaign / 'FE_VALIDATION_REPORT.md').read_text()
    # A completed result must be explicitly provided by Abaqus, not the
    # native solver. Synthetic data here tests import logic only.
    C = homogenize_rve(*load_project(PROJECT)).CH
    (run / 'homogenized_results.json').write_text(json.dumps({
        'status': 'COMPLETED', 'CH': C.tolist(), 'identification_method': 'direct_reaction6',
        'direct_reaction_audit': {'max_energy_reaction_work_error': 0.0, 'macro_reaction_symmetry_error': 0.0},
        'pbc_audit': {'passed': True}, 'energy_vs_reaction_relative_error': 0.0
    }))
    rows = collect_campaign(campaign)
    assert rows[0]['Status'] == 'FE_COMPLETED'
    assert rows[0]['FE_vs_analytical_Frobenius_pct'] < 1e-9
    assert (campaign / 'figures' / 'P00050_FE_mesh_convergence.png').exists()


def test_campaign_rejects_invalid_targets():
    target = np.eye(6); target[-1,-1] = -1
    with pytest.raises(ValueError, match='positive definite'):
        _target_metrics(np.eye(6), target)


def test_reject_unverified_direct_reactions(tmp_path):
    root = tmp_path / 'fe'
    m = prepare_campaign(root, [PROJECT], mesh_sizes=[30.], domain_mode='primitive')
    run = root / m['records'][0]['run_dir']
    C = homogenize_rve(*load_project(PROJECT)).CH
    (run / 'homogenized_results.json').write_text(json.dumps({
        'status': 'COMPLETED', 'CH': C.tolist(),
        'identification_method': 'direct_reaction6',
        'direct_reaction_audit': {'max_energy_reaction_work_error': 0.15,
                                  'macro_reaction_symmetry_error': 0.001},
        'pbc_audit': {'passed': True}
    }))
    rows = collect_campaign(root)
    assert rows[0]['Status'] == 'FE_INVALID'
    assert 'work mismatch' in rows[0]['Error']
