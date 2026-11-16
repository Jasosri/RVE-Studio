import numpy as np
from rvestudio.examples_data import thesis_hexagonal_case
from rvestudio.abaqus_thesis_dbc import _case_constraints, build_thesis_dbc_rectangle, generate_thesis_dbc_input, ThesisDBCSettings
from rvestudio.abaqus_homogenization import generate_periodic_homogenization_input, PeriodicHomogenizationSettings

def test_physical_top_bottom_are_never_periodically_tied_in_publication_dbc():
    model,mats,lams=thesis_hexagonal_case(1)
    dbc,diag=build_thesis_dbc_rectangle(model)
    assert diag['top_bottom_periodic'] is False
    for case in ('Q11','Q22','Q12','Q66'):
        assert not any(t[0].startswith('X3_') for t in _case_constraints(case,1.0))
    for case in ('Q44','Q55'):
        tops=[x for x in _case_constraints(case,1.0) if x[0] in ('X3_0','X3_H')]
        assert {(t[0],t[1]) for t in tops} == {(x,i) for x in ('X3_0','X3_H') for i in (1,2)}
        deck,meta=generate_thesis_dbc_input(dbc,mats,lams,ThesisDBCSettings(mesh_size=10),case,resolved_domain=True,domain_diagnostics=diag)
        assert meta['top_bottom_periodic'] is False
        assert '*EQUATION' not in deck.upper()

def test_triply_periodic_case_is_explicitly_labelled_as_different_model():
    model,mats,lams=thesis_hexagonal_case(1)
    deck,meta=generate_periodic_homogenization_input(model,mats,lams,PeriodicHomogenizationSettings(mesh_size=10.0),np.array([1e-3,0,0,0,0,0]))
    assert meta['top_bottom_periodic'] is True
    assert meta['finite_height_core_validation'] is False
    assert '3-D PERIODIC MATERIAL SURROGATE' in deck
