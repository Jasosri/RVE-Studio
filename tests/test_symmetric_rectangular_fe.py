"""Integration regression for cylindrical ply optimisation + rectangular FE geometry."""
from pathlib import Path
import json
import numpy as np
import pytest
from rvestudio.io import load_project
from rvestudio.rectangular_selection import inspect_rectangular_rve
from rvestudio.abaqus_homogenization import (
    PeriodicHomogenizationSettings, PBC_DOMAIN_SYMMETRIC_RECTANGULAR,
    PBC_DOMAIN_RECTANGULAR_SUPERCELL, resolve_periodic_domain,
    generate_periodic_homogenization_input,
)

EXAMPLES=Path(__file__).resolve().parents[1]/'examples'/'symmetric_rectangular_fe'

@pytest.mark.parametrize('name', ['SQUARE_symmetric_0_90_0','MIXED_symmetric_0_90_0'])
def test_rectangular_selection_reflection_material_and_measure(name):
    model,materials,laminates=load_project(EXAMPLES/f'{name}.rvejson')
    rect,diag=inspect_rectangular_rve(model,laminates)
    assert diag['mirror_geometry_symmetric']
    assert diag['same_stack_on_mirror_walls']
    assert diag['area_consistency_relative_error']<1e-10
    assert diag['wall_measure_relative_error']<1e-8
    assert abs(float(np.dot(rect.a1,rect.a2))) < 1e-8
    assert diag['primitive_cell_count']>=1


def test_symmetry_selected_fe_generates_real_abaqus_pbc_deck():
    model,materials,laminates=load_project(EXAMPLES/'SQUARE_symmetric_0_90_0.rvejson')
    settings=PeriodicHomogenizationSettings(domain_mode=PBC_DOMAIN_SYMMETRIC_RECTANGULAR,
                                           supercell_repeat_x=1,supercell_repeat_y=1,mesh_size=20.)
    deck,meta=generate_periodic_homogenization_input(model,materials,laminates,settings,
                                                     np.array([0.001,0,0,0,0,0.]),job_name='SYM_RECT_TEST')
    assert '*ELEMENT, TYPE=S4' in deck.upper()
    assert '*EQUATION' in deck.upper()
    assert meta['pbc_audit']['passed']
    assert meta['periodic_domain']['mirror_geometry_symmetric']
    assert meta['periodic_domain']['mirror_stack_symmetric']


def test_off_axis_optimum_still_allows_unrestricted_rectangle_but_not_symmetric_mode():
    path=Path(__file__).resolve().parents[1]/'examples'/'fe_validation'/'P00050.rvejson'
    model,materials,laminates=load_project(path)
    rect,_=resolve_periodic_domain(model,PeriodicHomogenizationSettings(domain_mode=PBC_DOMAIN_RECTANGULAR_SUPERCELL,supercell_repeat_x=1,supercell_repeat_y=1))
    assert len(rect.elements)>0
    with pytest.raises(ValueError,match='symmetric rectangular'):
        resolve_periodic_domain(model,PeriodicHomogenizationSettings(domain_mode=PBC_DOMAIN_SYMMETRIC_RECTANGULAR,supercell_repeat_x=1,supercell_repeat_y=1))


def test_invalid_off_axis_layup_is_blocked_even_for_symmetric_wall_graph():
    model,materials,laminates=load_project(EXAMPLES/'SQUARE_symmetric_0_90_0.rvejson')
    laminates[0].plies[1].angle_deg=45
    settings=PeriodicHomogenizationSettings(domain_mode=PBC_DOMAIN_SYMMETRIC_RECTANGULAR,
                                           supercell_repeat_x=1,supercell_repeat_y=1,mesh_size=20.)
    with pytest.raises(ValueError,match='0/90-degree'):
        generate_periodic_homogenization_input(model,materials,laminates,settings,
                                                 np.array([0.001,0,0,0,0,0.]))
