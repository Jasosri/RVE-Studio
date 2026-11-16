import json
import numpy as np

from rvestudio.examples_data import thesis_hexagonal_case
from rvestudio.abaqus_thesis_dbc import (
    ThesisDBCSettings, THESIS_DBC_CASES, build_thesis_dbc_rectangle,
    generate_thesis_dbc_input, reconstruct_thesis_q, prepare_thesis_dbc,
    _case_constraints,
)


def test_table_3_4_conditions_are_transcribed_exactly():
    d=1.0
    assert _case_constraints('Q11',d)==[('X1_0',1,0.),('X1_P',1,d),('X2_0',2,0.),('X2_Q',2,0.)]
    assert _case_constraints('Q22',d)==[('X1_0',1,0.),('X1_P',1,0.),('X2_0',2,0.),('X2_Q',2,d)]
    assert _case_constraints('Q12',d)==[('X1_0',1,0.),('X1_P',1,d),('X2_0',2,0.),('X2_Q',2,d)]
    assert _case_constraints('Q66',d)==[('X1_0',1,0.),('X1_0',2,0.),('X1_P',1,0.),('X1_P',2,d),('X2_0',1,0.),('X2_Q',1,0.)]
    assert _case_constraints('Q44',d)==[
        ('X1_0',1,0.),('X1_0',2,0.),('X1_0',3,0.),('X1_P',1,0.),('X1_P',2,0.),('X1_P',3,d),
        ('X2_0',1,0.),('X2_0',2,0.),('X2_Q',1,0.),('X2_Q',2,0.),('X3_0',1,0.),('X3_0',2,0.),('X3_H',1,0.),('X3_H',2,0.)]
    assert _case_constraints('Q55',d)==[
        ('X1_0',1,0.),('X1_0',2,0.),('X1_P',1,0.),('X1_P',2,0.),('X2_0',1,0.),('X2_0',2,0.),('X2_0',3,0.),
        ('X2_Q',1,0.),('X2_Q',2,0.),('X2_Q',3,d),('X3_0',1,0.),('X3_0',2,0.),('X3_H',1,0.),('X3_H',2,0.)]


def test_config1_dbc_rve_is_the_published_four_cell_block():
    # Thesis Fig. 3.6(b) / Table 3.3 and paper Fig. 8(b): p = 2|a1|, q = |a1 + 2 a2|.
    model,_,_=thesis_hexagonal_case(1)
    dbc,diag=build_thesis_dbc_rectangle(model)
    np.testing.assert_allclose(diag['p'],2*86.60254037844388,rtol=1e-12)
    np.testing.assert_allclose(diag['q'],150.0,rtol=1e-12)
    assert diag['primitive_cells']==4 and diag['cells_x1_x2']==[2,1]
    np.testing.assert_allclose(diag['rve_area'],4*diag['primitive_area'],rtol=1e-12)
    np.testing.assert_allclose(diag['rve_volume'],2*86.60254037844388*150.0*150.0,rtol=1e-12)
    # 8 inclined walls, 2 full and 2 half straight walls inside, 4 half straight walls in the x1 faces.
    assert len(dbc.elements)==16
    np.testing.assert_allclose(diag['actual_wall_length'],650.0,rtol=1e-12)
    np.testing.assert_allclose(diag['effective_wall_length'],600.0,rtol=1e-12)
    assert sorted(diag['face_wall_faces'].values())==['X1_0','X1_0','X1_P','X1_P']
    np.testing.assert_allclose(diag['face_wall_length'],100.0,rtol=1e-12)
    assert diag['wall_measure_relative_error'] < 1e-12


def test_config1_dbc_rve_faces_cut_where_the_thesis_figure_shows():
    model,_,_=thesis_hexagonal_case(1)
    dbc,diag=build_thesis_dbc_rectangle(model)
    p,q=diag['p'],diag['q']; X=np.asarray(dbc.nodes)
    segs=[(X[e.node_i],X[e.node_j]) for e in dbc.elements]
    # x2 faces: straight walls cut at mid-length (25 mm stubs) at x1 = 0, p/2, p.
    for face_y in (0.0,q):
        stubs=sorted(round(float(a[0]),6) for a,b in segs if abs(a[0]-b[0])<1e-9 and min(abs(a[1]-face_y),abs(b[1]-face_y))<1e-9)
        assert stubs==[0.0,round(p/2,6),round(p,6)]
    # One complete hexagon in the middle, centred in the block.
    np.testing.assert_allclose(diag['cell_centroid_rve'],[p/2,q/2],atol=1e-9)
    poly=np.asarray(diag['central_cell_polygon_rve'])
    assert len(poly)==6 and poly[:,0].min()>1 and poly[:,0].max()<p-1 and poly[:,1].min()>1 and poly[:,1].max()<q-1


def test_face_walls_get_half_plies_in_the_deck_and_full_option_is_available():
    model,mats,lams=thesis_hexagonal_case(1)
    tp=lams[0].plies[0].thickness
    text,meta=generate_thesis_dbc_input(model,mats,lams,ThesisDBCSettings(mesh_size=10.0),'Q11')
    assert meta['face_wall_section']=='half-thickness' and len(meta['face_walls'])==4
    assert f"{0.5*tp:.12g}, 3, " in text and f"{tp:.12g}, 3, " in text
    text_full,meta_full=generate_thesis_dbc_input(model,mats,lams,ThesisDBCSettings(mesh_size=10.0,face_wall_section='full'),'Q11')
    assert f"{0.5*tp:.12g}, 3, " not in text_full
    text_hs,_=generate_thesis_dbc_input(model,mats,lams,ThesisDBCSettings(mesh_size=10.0,face_wall_section='half-stiffness'),'Q11')
    assert '_HALF' in text_hs


def test_dbc_deck_leaves_shell_rotations_free_and_only_adds_u3_reference_for_inplane():
    model,mats,lams=thesis_hexagonal_case(1)
    dbc,diag=build_thesis_dbc_rectangle(model)
    settings=ThesisDBCSettings(mesh_size=10.0,displacement=1.0)
    text,meta=generate_thesis_dbc_input(dbc,mats,lams,settings,'Q11',resolved_domain=True,domain_diagnostics=diag)
    assert '*STEP, NAME=DBC' in text
    assert '1, 6, 0.' not in text
    assert len(meta['minimal_rigid_body_stabilization'])==1
    assert meta['minimal_rigid_body_stabilization'][0]['dof']==3
    assert meta['rigid_body_audit']['pre_stabilization_rank']==5
    assert meta['rigid_body_audit']['post_stabilization_rank']==6
    text44,meta44=generate_thesis_dbc_input(dbc,mats,lams,settings,'Q44',resolved_domain=True,domain_diagnostics=diag)
    assert meta44['minimal_rigid_body_stabilization'] is None
    assert meta44['rigid_body_audit']['pre_stabilization_rank']==6
    assert meta44['rigid_body_audit']['post_stabilization_rank']==6
    assert all(meta['boundary_sets'][k] > 0 for k in ['X1_0','X1_P','X2_0','X2_Q','X3_0','X3_H'])


def test_thesis_table_3_3_config1_energies_recover_published_fe_q_values():
    # Table 3.3 energies for configuration 1 (N mm). The Q12 OCR in the PDF is 96749.70.
    model,_,_=thesis_hexagonal_case(1)
    _,diag=build_thesis_dbc_rectangle(model)
    energies={'Q11':20951.1,'Q22':27936.0,'Q12':96749.70,'Q66':160.6,'Q44':21801.7,'Q55':29224.4}
    q=reconstruct_thesis_q(energies,diag['p'],diag['q'],diag['height'],1.0)
    ref={'Q11':322.63,'Q22':322.61,'Q12':319.03,'Q66':2.47,'Q44':335.65,'Q55':337.46}
    for k in ref:
        assert abs(q[k]-ref[k])/ref[k] < 0.002


def test_prepare_thesis_dbc_writes_ten_case_package(tmp_path):
    """Six Table 3.4 shell-term cases (top/bottom U3 free) and four U3-restrained 3D-term cases."""
    from rvestudio.abaqus_thesis_dbc import DBC_CASES, DBC_C_CASES
    model,mats,lams=thesis_hexagonal_case(1)
    files=prepare_thesis_dbc(tmp_path,model,mats,lams,ThesisDBCSettings(mesh_size=10.0))
    assert len(files.cases)==10
    assert [c.name for c in files.cases]==list(DBC_CASES)==list(THESIS_DBC_CASES)+list(DBC_C_CASES)
    snap=json.loads(files.snapshot_file.read_text())
    assert snap['mode']=='thesis_table_3_4_prescribed_displacement_validation'
    assert np.isclose(snap['domain']['rve_area'],2*86.60254037844388*150.0)
    assert (files.run_dir/'thesis_dbc_rve.svg').exists()
    assert files.cases[0].input_file.exists()


def test_all_thesis_dbc_cases_have_full_rigid_body_rank_after_minimal_gauge():
    for cfg in range(1,15):
        model,mats,lams=thesis_hexagonal_case(cfg)
        dbc,diag=build_thesis_dbc_rectangle(model)
        for case in THESIS_DBC_CASES+('C1111','C1212'):
            _text,meta=generate_thesis_dbc_input(dbc,mats,lams,ThesisDBCSettings(mesh_size=10.0),case,resolved_domain=True,domain_diagnostics=diag)
            audit=meta['rigid_body_audit']
            assert audit['connected_shell_components']==1
            assert audit['post_stabilization_rank']==6
            if case in {'Q11','Q22','Q12','Q66'}:
                assert audit['pre_stabilization_rank']==5
                assert len(audit['added_gauge_constraints'])==1
                assert audit['added_gauge_constraints'][0]['dof']==3
                # The only unconstrained rigid motion is pure translation in x3.
                n=np.asarray(audit['pre_nullspace'][0],float)
                assert abs(abs(n[2])-1.0)<1e-10
                assert np.linalg.norm(np.delete(n,2))<1e-10
            else:
                assert audit['pre_stabilization_rank']==6
                assert audit['added_gauge_constraints']==[]
