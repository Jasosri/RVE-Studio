import json
import numpy as np

from rvestudio.examples_data import thesis_hexagonal_case
from rvestudio.abaqus_homogenization import (
    PeriodicHomogenizationSettings,
    homogenization_cases,
    generate_periodic_homogenization_input,
    prepare_periodic_homogenization,
    reconstruct_ch_from_energies,
    reconstruct_ch_from_macro_reactions,
    reconstruct_ch_direct_from_reactions, CH_METHOD_ENERGY, CH_METHOD_DIRECT,
    independent_terms,
    _macro_coefficients,
)


def test_21_energy_cases_are_complete():
    cases = homogenization_cases(0.001)
    assert len(cases) == 21
    assert sum(len(c['components']) == 1 for c in cases) == 6
    assert sum(len(c['components']) == 2 for c in cases) == 15
    pairs = {tuple(c['components']) for c in cases if len(c['components']) == 2}
    assert len(pairs) == 15


def test_engineering_shear_mapping_is_symmetric():
    # delta = [2,4,6], reference length 1. G12 produces [dy/2, dx/2, 0].
    B = _macro_coefficients(np.array([2.0, 4.0, 6.0]), 1.0)
    np.testing.assert_allclose(B[:, 5], [2.0, 1.0, 0.0])
    np.testing.assert_allclose(B[:, 4], [3.0, 0.0, 1.0])
    np.testing.assert_allclose(B[:, 3], [0.0, 3.0, 2.0])


def test_reconstruct_all_21_terms_from_synthetic_energy():
    # Symmetric positive-definite matrix with non-zero normal/shear couplings.
    A = np.array([
        [8.0,1.0,.8,.3,.2,.4],
        [0.0,7.0,.7,.2,.3,.5],
        [0.0,0.0,6.0,.4,.1,.2],
        [0.0,0.0,0.0,5.0,.5,.3],
        [0.0,0.0,0.0,0.0,4.0,.4],
        [0.0,0.0,0.0,0.0,0.0,3.0],
    ])
    C = A.T @ A
    amp = 0.002; V = 1234.5
    rows = []
    for case in homogenization_cases(amp):
        e = np.asarray(case['strain'])
        W = 0.5 * V * float(e @ C @ e)
        rows.append({'components': case['components'], 'strain_energy': W})
    got = reconstruct_ch_from_energies(rows, V, amp)
    np.testing.assert_allclose(got, C, rtol=1e-11, atol=1e-11)
    terms = independent_terms(got, 1e-12)
    assert len(terms) == 21
    assert any(t['term'] == 'C16' and t['significant'] for t in terms)


def test_macro_reaction_crosscheck_recovers_unsymmetric_columns():
    C = np.array([
        [10,2,1,.2,.3,.4], [2,9,1.5,.1,.2,.5], [1,1.5,8,.3,.4,.2],
        [.2,.1,.3,5,.4,.2], [.3,.2,.4,.4,4,.3], [.4,.5,.2,.2,.3,3]
    ], float)
    V=100.; amp=.001; Lref=1.0
    rows=[]
    # Simulate Abaqus returning the opposite sign for constraint reactions.
    for j in range(6):
        rf = -(C[:,j] * V * amp / Lref)
        rows.append({'components':[j], 'macro_reactions':rf.tolist()})
    got, err = reconstruct_ch_from_macro_reactions(rows,V,amp,Lref)
    np.testing.assert_allclose(got,C)
    assert err < 1e-12


def test_periodic_input_and_21_case_package(tmp_path):
    model,mats,lams = thesis_hexagonal_case(1)
    settings=PeriodicHomogenizationSettings(mesh_size=10.0,strain_amplitude=0.001,identification_method=CH_METHOD_ENERGY)
    case=homogenization_cases(settings.strain_amplitude)[0]
    text, meta = generate_periodic_homogenization_input(model,mats,lams,settings,case['strain'],job_name='H')
    assert '*ELEMENT, TYPE=S4' in text
    assert '*EQUATION' in text
    assert '*NSET, NSET=MACRO_ALL' in text
    assert '*ENERGY OUTPUT' in text and 'ALLSE' in text
    assert 'MACRO_1, 1, 1, 0.001' in text
    assert meta['periodic_pair_count'] > 0
    files = prepare_periodic_homogenization(tmp_path,model,mats,lams,settings)
    assert len(files.cases) == 21
    snap=json.loads(files.snapshot_file.read_text())
    assert snap['mode'] == 'general_anisotropic_periodic_energy_homogenization'
    assert len(snap['cases']) == 21


def test_pbc_anchor_does_not_clamp_shell_rotations():
    model,mats,lams = thesis_hexagonal_case(1)
    settings=PeriodicHomogenizationSettings(mesh_size=10.0,strain_amplitude=0.001)
    case=homogenization_cases(settings.strain_amplitude)[0]
    text,_ = generate_periodic_homogenization_input(model,mats,lams,settings,case['strain'],job_name='H')
    assert 'PBC_ANCHOR, 1, 3, 0.' in text
    assert 'PBC_ANCHOR, 1, 6, 0.' not in text


def test_chapter3_hex_rve_area_volume_definition():
    from rvestudio.abaqus_homogenization import rve_area_volume_diagnostics
    model,_,_ = thesis_hexagonal_case(1)
    d=rve_area_volume_diagnostics(model)
    L1=L2=50.0; theta=np.deg2rad(120.0)
    a=2*L2*np.sin(theta)
    b=L1+L2*np.cos(np.pi-theta)
    expected_area=a*b
    np.testing.assert_allclose(d['lattice_area_det'],expected_area,rtol=1e-13)
    np.testing.assert_allclose(d['chapter3_area_p_times_q'],expected_area,rtol=1e-13)
    np.testing.assert_allclose(d['rve_volume'],expected_area*150.0,rtol=1e-13)


def test_general_plane_stress_schur_reduction_and_thesis_terms():
    from rvestudio.abaqus_homogenization import plane_stress_shell_reduction
    C=np.array([
        [100.,20.,10.,2.,3.,4.],
        [20.,90.,8.,1.,2.,5.],
        [10.,8.,80.,6.,7.,9.],
        [2.,1.,6.,30.,3.,2.],
        [3.,2.,7.,3.,25.,1.],
        [4.,5.,9.,2.,1.,20.],
    ])
    Q,terms=plane_stress_shell_reduction(C)
    keep=[0,1,5,4,3]
    expected=C[np.ix_(keep,keep)]-np.outer(C[keep,2],C[2,keep])/C[2,2]
    np.testing.assert_allclose(Q,expected)
    assert np.isclose(terms['Q11'],expected[0,0])
    assert np.isclose(terms['Q66'],expected[2,2])
    assert np.isclose(terms['Q44'],expected[3,3])
    assert np.isclose(terms['Q55'],expected[4,4])


def test_thesis_fe_reference_config1_is_labelled_q_not_c():
    from rvestudio.abaqus_homogenization import thesis_fe_q_reference
    model,_,_ = thesis_hexagonal_case(1)
    ref=thesis_fe_q_reference(model)
    assert ref == {'Q11':322.63,'Q22':322.61,'Q12':319.03,'Q66':2.47,'Q44':335.65,'Q55':337.46}


def test_rectangular_supercell_config1_area_wall_measure_and_basis():
    from rvestudio.abaqus_homogenization import (
        build_rectangular_supercell, find_rectangular_supercell_basis,
    )
    model,_,_ = thesis_hexagonal_case(1)
    basis=find_rectangular_supercell_basis(model)
    assert basis['primitive_cells_in_minimum_rectangle'] == 2
    b1=np.asarray(basis['b1']); b2=np.asarray(basis['b2'])
    assert abs(float(np.dot(b1,b2))) < 1e-10
    super_model,d=build_rectangular_supercell(model,repeats_x=2,repeats_y=2)
    assert d['primitive_cell_count'] == 8
    np.testing.assert_allclose(d['supercell_area'],8*d['primitive_area'],rtol=1e-12)
    np.testing.assert_allclose(d['actual_supercell_wall_length'],d['expected_supercell_wall_length'],rtol=1e-12)
    assert d['area_consistency_relative_error'] < 1e-12
    assert d['wall_measure_relative_error'] < 1e-12
    assert abs(float(np.dot(super_model.a1,super_model.a2))) < 1e-10


def test_rectangular_supercell_deck_has_complete_strict_pbc_audit(tmp_path):
    from rvestudio.abaqus_homogenization import PBC_DOMAIN_RECTANGULAR_SUPERCELL
    model,mats,lams=thesis_hexagonal_case(1)
    settings=PeriodicHomogenizationSettings(
        mesh_size=100.0,strain_amplitude=0.001,
        domain_mode=PBC_DOMAIN_RECTANGULAR_SUPERCELL,
        supercell_repeat_x=2,supercell_repeat_y=2,
        strict_pbc_audit=True,
    )
    case=homogenization_cases(settings.strain_amplitude)[0]
    text,meta=generate_periodic_homogenization_input(model,mats,lams,settings,case['strain'],job_name='RECT')
    d=meta['periodic_domain']; audit=meta['pbc_audit']
    assert d['domain_mode'] == PBC_DOMAIN_RECTANGULAR_SUPERCELL
    assert d['primitive_cell_count'] >= 4
    assert audit['passed'] is True
    assert audit['periodic_translation_rank'] == 3
    assert audit['rigid_body_constraint_rank'] == 6
    assert audit['unpaired_boundary_nodes'] == 0
    assert all(x['complete'] for x in audit['boundary_coverage'].values())
    assert '** Periodic domain: rectangular_supercell' in text
    files=prepare_periodic_homogenization(tmp_path,model,mats,lams,settings)
    snap=json.loads(files.snapshot_file.read_text())
    assert snap['periodic_domain']['primitive_cell_count'] == d['primitive_cell_count']
    np.testing.assert_allclose(snap['rve_area'],d['supercell_area'],rtol=1e-12)


def test_rectangular_supercell_search_and_pbc_audit_for_builtin_families():
    from rvestudio.geometry import BUILTIN_GEOMETRIES
    from rvestudio.models import RVEModel, Element, Material, Laminate, Ply
    from rvestudio.abaqus_homogenization import PBC_DOMAIN_RECTANGULAR_SUPERCELL
    params={
        'HEX':{'L1':50.,'L2':50.}, 'SQUARE':{'Lx':50.,'Ly':40.},
        'MIXED':{'Lh':40.,'R':2.}, 'TRI':{'Lh':50.,'R':1.},
    }
    mats={'M':Material('M',10000.,5000.,3000.,0.25,3000.)}
    lams=[Laminate([Ply('M',0.,0.5)])]
    for gid,fn in BUILTIN_GEOMETRIES.items():
        g=fn(params[gid])
        model=RVEModel(g.nodes,[Element(int(i),int(j),0,float(w)) for (i,j),w in zip(g.connectivity,g.wall_weights)],g.a1,g.a2,100.,wall_groups=g.wall_groups,name=gid)
        settings=PeriodicHomogenizationSettings(mesh_size=25.,domain_mode=PBC_DOMAIN_RECTANGULAR_SUPERCELL,supercell_repeat_x=2,supercell_repeat_y=2)
        _,meta=generate_periodic_homogenization_input(model,mats,lams,settings,np.array([.001,0,0,0,0,0]),job_name='B')
        assert meta['pbc_audit']['passed'] is True
        assert meta['periodic_domain']['primitive_cell_count'] >= 4
        assert meta['periodic_domain']['area_consistency_relative_error'] < 1e-10
        assert meta['periodic_domain']['wall_measure_relative_error'] < 1e-7


def test_prepare_package_writes_exact_resolved_supercell_to_deck(tmp_path):
    from rvestudio.abaqus_homogenization import PBC_DOMAIN_RECTANGULAR_SUPERCELL
    model,mats,lams=thesis_hexagonal_case(1)
    settings=PeriodicHomogenizationSettings(
        mesh_size=100.0,strain_amplitude=0.001,
        domain_mode=PBC_DOMAIN_RECTANGULAR_SUPERCELL,
        supercell_repeat_x=2,supercell_repeat_y=2,strict_pbc_audit=True,
    )
    files=prepare_periodic_homogenization(tmp_path,model,mats,lams,settings,run_name='handoff')
    snap=json.loads(files.snapshot_file.read_text())
    analysis=snap['analysis_model']
    meta=snap['cases'][0]['deck_meta']
    assert snap['periodic_domain']['primitive_cell_count'] == 8
    assert len(analysis['elements']) > len(model.elements)
    assert meta['analysis_geometry_wall_count'] == len(analysis['elements'])
    assert meta['analysis_geometry_node_count'] == len(analysis['nodes'])
    np.testing.assert_allclose(meta['analysis_a1'],analysis['a1'])
    np.testing.assert_allclose(meta['analysis_a2'],analysis['a2'])
    # The actual first Abaqus input deck must contain the enlarged resolved domain.
    text=files.cases[0].input_file.read_text()
    node_lines=[]; reading=False
    for line in text.splitlines():
        if line.startswith('*NODE'):
            reading=True; continue
        if reading and line.startswith('*'):
            break
        if reading and line.strip() and not line.startswith('**'):
            parts=[x.strip() for x in line.split(',')]
            try:
                node_lines.append((float(parts[1]),float(parts[2]),float(parts[3])))
            except Exception:
                pass
    xs=[p[0] for p in node_lines]; ys=[p[1] for p in node_lines]
    assert max(xs)-min(xs) > 170.0
    assert max(ys)-min(ys) > 299.0
    assert text.count('*ELSET, ELSET=WALL_') == len(analysis['elements'])


def test_lattice_supercell_is_exact_translation_of_primitive_config1():
    from rvestudio.abaqus_homogenization import build_lattice_supercell
    model,_,_=thesis_hexagonal_case(1)
    sm,d=build_lattice_supercell(model,repeats_a1=2,repeats_a2=2)
    assert d['domain_mode']=='lattice_supercell'
    assert d['primitive_cell_count']==4
    np.testing.assert_allclose(sm.a1,2*model.a1)
    np.testing.assert_allclose(sm.a2,2*model.a2)
    np.testing.assert_allclose(d['supercell_area'],4*d['primitive_area'],rtol=1e-12)
    np.testing.assert_allclose(d['actual_supercell_wall_length'],4*d['primitive_wall_length'],rtol=1e-12)
    # Every supercell wall vector must be exactly one of the primitive wall vectors
    # (up to sign); clipping would create shorter boundary fragments and fail this.
    primitive_vectors=[]
    for e in model.elements:
        v=np.asarray(model.nodes[e.node_j])-np.asarray(model.nodes[e.node_i])
        primitive_vectors.append(v)
    for e in sm.elements:
        v=np.asarray(sm.nodes[e.node_j])-np.asarray(sm.nodes[e.node_i])
        assert any(np.allclose(v,pv,atol=1e-10) or np.allclose(v,-pv,atol=1e-10) for pv in primitive_vectors)


def test_lattice_supercell_strict_pbc_audit_and_exact_deck(tmp_path):
    from rvestudio.abaqus_homogenization import PBC_DOMAIN_LATTICE_SUPERCELL
    model,mats,lams=thesis_hexagonal_case(1)
    settings=PeriodicHomogenizationSettings(
        mesh_size=25.,strain_amplitude=.001,domain_mode=PBC_DOMAIN_LATTICE_SUPERCELL,
        supercell_repeat_x=2,supercell_repeat_y=2,strict_pbc_audit=True,
    )
    files=prepare_periodic_homogenization(tmp_path,model,mats,lams,settings,run_name='lattice')
    snap=json.loads(files.snapshot_file.read_text())
    dom=snap['periodic_domain']; meta=snap['cases'][0]['deck_meta']
    assert dom['primitive_cell_count']==4
    assert dom['construction'].startswith('exact primitive-RVE translations')
    assert meta['pbc_audit']['passed'] is True
    assert meta['pbc_audit']['periodic_translation_rank']==3
    assert meta['pbc_audit']['rigid_body_constraint_rank']==6
    assert meta['pbc_audit']['unpaired_boundary_nodes']==0
    analysis=snap['analysis_model']
    assert len(analysis['elements'])==4*len(model.elements)
    np.testing.assert_allclose(analysis['a1'],2*model.a1)
    np.testing.assert_allclose(analysis['a2'],2*model.a2)
    assert (files.run_dir/'periodic_domain_preview.svg').exists()


def test_lattice_supercell_pbc_audit_for_builtin_families():
    from rvestudio.geometry import BUILTIN_GEOMETRIES
    from rvestudio.models import RVEModel, Element, Material, Laminate, Ply
    from rvestudio.abaqus_homogenization import PBC_DOMAIN_LATTICE_SUPERCELL
    params={
        'HEX':{'L1':50.,'L2':50.}, 'SQUARE':{'Lx':50.,'Ly':40.},
        'MIXED':{'Lh':40.,'R':2.}, 'TRI':{'Lh':50.,'R':1.},
    }
    mats={'M':Material('M',10000.,5000.,3000.,0.25,3000.)}
    lams=[Laminate([Ply('M',0.,0.5)])]
    for gid,fn in BUILTIN_GEOMETRIES.items():
        g=fn(params[gid])
        model=RVEModel(g.nodes,[Element(int(i),int(j),0,float(w)) for (i,j),w in zip(g.connectivity,g.wall_weights)],g.a1,g.a2,100.,wall_groups=g.wall_groups,name=gid)
        settings=PeriodicHomogenizationSettings(mesh_size=25.,domain_mode=PBC_DOMAIN_LATTICE_SUPERCELL,supercell_repeat_x=2,supercell_repeat_y=2)
        _,meta=generate_periodic_homogenization_input(model,mats,lams,settings,np.array([.001,0,0,0,0,0]),job_name='LAT')
        assert meta['pbc_audit']['passed'] is True
        assert meta['periodic_domain']['primitive_cell_count'] == 4
        assert meta['periodic_domain']['wall_measure_relative_error'] < 1e-10


def test_direct_six_reaction_states_identify_fully_coupled_ch():
    rng = np.random.default_rng(1801)
    B = rng.normal(size=(6, 6))
    C = B.T @ B + 5.0*np.eye(6)
    V, amp, Lref = 1312.5, 0.001, 1.5
    rows = []
    for j in range(6):
        r = (-C[:, j] * V * amp / Lref).tolist() # globally reversed Abaqus RF convention
        W = 0.5 * V * amp**2 * C[j, j]
        rows.append({'components': [j], 'macro_reactions': r, 'strain_energy': W})
    out, audit = reconstruct_ch_direct_from_reactions(rows, V, amp, Lref)
    np.testing.assert_allclose(out, C, rtol=1e-12, atol=1e-12)
    assert audit['reaction_sign'] == -1
    assert audit['macro_reaction_symmetry_error'] < 1e-12
    assert audit['max_energy_reaction_work_error'] < 1e-12
    # Reject a column with opposite sign, rather than fixing it silently.
    rows[2]['macro_reactions'] = (-np.asarray(rows[2]['macro_reactions'])).tolist()
    import pytest
    with pytest.raises(ValueError, match='Inconsistent macro-reaction signs'):
        reconstruct_ch_direct_from_reactions(rows, V, amp, Lref)


def test_direct_ch_is_default_and_creates_only_six_decks(tmp_path):
    model, mats, lams = thesis_hexagonal_case(1)
    settings = PeriodicHomogenizationSettings(mesh_size=20.0)
    assert settings.identification_method == CH_METHOD_DIRECT
    files = prepare_periodic_homogenization(tmp_path, model, mats, lams, settings, run_name='Direct_CH')
    assert len(files.cases) == 6
    snap = json.loads(files.snapshot_file.read_text())
    assert snap['mode'] == 'general_anisotropic_periodic_direct_CH'
    assert 'macro reaction' in snap['identification']
    assert [c['name'] for c in snap['cases']] == ['E11','E22','E33','G23','G13','G12']


def test_full_direct_ch_runner_with_mocked_abaqus_odb(tmp_path, monkeypatch):
    """Exercise all six deck/result/assembly stages without pretending to run Abaqus."""
    from pathlib import Path
    from rvestudio import abaqus_homogenization as hm
    model, materials, laminates = thesis_hexagonal_case(1)
    settings = PeriodicHomogenizationSettings(mesh_size=25.0, strain_amplitude=0.001)
    files = prepare_periodic_homogenization(tmp_path, model, materials, laminates,
                                             settings, run_name='Mock_Direct')
    volume = json.loads(files.snapshot_file.read_text())['rve_volume']
    target = np.array([
        [24., 2., 1., 0.3, 0.1, 0.4],
        [2., 21., 1.1, .2, .3, 0.5],
        [1., 1.1, 18., .1, .2, .3],
        [.3, .2, .1, 12., .4, .1],
        [.1, .3, .2, .4, 11., .3],
        [.4, .5, .3, .1, .3, 10.]
    ])
    cases = {str(c.case_dir): (j, c) for j, c in enumerate(files.cases)}
    meta = json.loads(files.snapshot_file.read_text())['cases']

    def fake_abaqus(argv, cwd, log):
        j, c = cases[str(Path(cwd))]
        if len(argv) > 1 and str(argv[1]).startswith('job='):
            (c.case_dir / (c.job_name + '.odb')).write_text('test stub only')
        else:
            macro_ids = meta[j]['deck_meta']['macro_node_ids']
            force = target[:, j] * volume * settings.strain_amplitude / settings.reference_length
            (c.results_file).write_text(json.dumps({
                'status': 'COMPLETED',
                'strain_energy_ALLSE': 0.5*volume*settings.strain_amplitude**2*target[j,j],
                'macro_reaction_u1': {str(mid): float(v) for mid,v in zip(macro_ids, force)},
                'macro_displacement_u1': {},
                'odb': str(c.case_dir / (c.job_name+'.odb')),
            }))
        return 0
    monkeypatch.setattr(hm, '_run_shell_command', fake_abaqus)
    results = hm.run_periodic_homogenization(files, 'fake_abaqus', settings)
    assert results['status'] == 'COMPLETED'
    assert results['case_count'] == 6
    assert results['CH_source'] == 'six_macro_reaction_columns'
    np.testing.assert_allclose(results['CH'], target, atol=1e-12)
    assert 'CH_energy' not in results
    assert results['direct_reaction_audit']['max_energy_reaction_work_error'] < 1e-12
    assert (files.run_dir / 'CH_FE_MPa.csv').exists()
    assert not (files.run_dir / 'CH_energy_MPa.csv').exists()
    # This is a mocked output-pipeline test, NOT Abaqus validation.
