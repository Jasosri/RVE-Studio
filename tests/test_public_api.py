import json
import numpy as np
import pytest
from rvestudio.api import solve, solve_example, solve_project, summary, benchmark
from rvestudio.examples_data import paper_hexagonal_case
from rvestudio.io import save_project
from rvestudio.properties import effective_properties
from rvestudio.validation import extract_six


def test_paper_q44_q55_order_with_distinct_shears():
    # Internal Voigt C[3,3]=C2323, C[4,4]=C1313.
    c = np.diag([50.,60.,70., 5., 11., 9.])
    props = effective_properties(c)
    assert props.G23 == pytest.approx(5)
    assert props.G13 == pytest.approx(11)
    assert props.Qshear.shape == (2,2)
    assert props.Qshear[0,0] == pytest.approx(11)  # Paper S44
    assert props.Qshear[1,1] == pytest.approx(5)   # Paper S55
    assert props.Qshear[0,1] == pytest.approx(0)


def test_public_api_solve_and_json_safe_summary(tmp_path):
    model, mats, stacks = paper_hexagonal_case(2)
    out = solve(model,mats,stacks)
    item = summary(out, include_ch=True)
    assert item['software'] == 'RVE Studio'
    assert item['transverse_shear_order'] == ['13','23']
    assert item['voigt_order_3d'] == ['11','22','33','23','13','12']
    assert item['stiffness']['Q44'] == pytest.approx(out.props.Qshear[0,0])
    assert item['stiffness']['Q55'] == pytest.approx(out.props.Qshear[1,1])
    json.dumps(item, allow_nan=False)
    assert np.allclose(extract_six(out), extract_six(solve_example(2)))
    path = tmp_path / 'model.rvejson'
    save_project(path, model, mats, stacks)
    assert np.allclose(solve_project(path).CH, out.CH)


def test_legacy_model_identifier_compatibility():
    model, mats, lams = paper_hexagonal_case(1)
    ref = solve(model,mats,lams).CH
    model.wall_element_type = 'thesis'
    assert np.allclose(solve(model,mats,lams).CH, ref)


def test_solve_validates_api_inputs():
    model, mats, lams = paper_hexagonal_case(1)
    with pytest.raises(TypeError): solve(None,mats,lams)
    with pytest.raises(ValueError): solve(model,{},lams)
    with pytest.raises(ValueError): solve(model,mats,[])
    with pytest.raises(ValueError): paper_hexagonal_case(0)
