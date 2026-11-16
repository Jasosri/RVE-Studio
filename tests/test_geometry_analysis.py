import math
import numpy as np
import pytest

from rvestudio.examples_data import paper_hexagonal_case
from rvestudio.geometry_analysis import analyse_geometry, characteristic_pqr
from rvestudio.models import Element


def test_paper_hex_geometry_parameters():
    model, _, _ = paper_hexagonal_case(1)
    s = analyse_geometry(model.nodes, model.elements, model.a1, model.a2, periodic_tol=model.periodic_tol)
    assert s.is_valid
    assert not s.issues
    assert s.n_nodes == 5
    assert s.n_walls == 3
    assert np.allclose([w.length for w in s.walls], [50.0, 50.0, 50.0])
    assert math.isclose(s.p, 86.6025403784, rel_tol=1e-10)
    assert math.isclose(s.q, 75.0, rel_tol=1e-10)
    assert math.isclose(s.r, 43.3012701892, rel_tol=1e-10)
    assert math.isclose(s.cell_area, s.p * s.q, rel_tol=1e-10)
    assert math.isclose(s.lattice_angle_deg, 120.0, abs_tol=1e-10)
    assert len(s.periodic_pairs) >= 1


def test_characteristic_pqr_is_rotation_invariant():
    theta = math.radians(37)
    R = np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])
    a1 = R @ np.array([10.0, 0.0])
    a2 = R @ np.array([-3.0, 7.0])
    p, q, rs, r = characteristic_pqr(a1, a2)
    assert np.allclose([p, q, rs, r], [10.0, 7.0, -3.0, 3.0])


def test_geometry_validation_catches_bad_walls():
    nodes = np.array([[0.0, 0.0], [10.0, 0.0], [10.0, 0.0]])
    elements = [Element(0, 1), Element(1, 2), Element(0, 1)]
    s = analyse_geometry(nodes, elements, [10, 0], [0, 10])
    messages = " ".join(i.message for i in s.issues)
    assert "coincident" in messages
    assert "zero" in messages.lower()
    assert "duplicates" in messages.lower()
    assert not s.is_valid


def test_builtin_geometry_detection_and_parameter_extraction():
    from rvestudio.geometry import hexagonal, square_family, mixed_rhombus_triangular, triangular_family
    from rvestudio.models import Element

    cases=[
        (hexagonal,{'L1':50.0,'L2':60.0},'HEX',{'L1':50.0,'L2':60.0}),
        (square_family,{'Lx':40.0,'Ly':55.0},'SQUARE',{'Lx':40.0,'Ly':55.0}),
        (mixed_rhombus_triangular,{'Lh':40.0,'R':2.0},'MIXED',{'Lh':40.0,'R':2.0}),
        (triangular_family,{'Lh':45.0,'R':1.0},'TRI',{'Lh':45.0,'R':1.0}),
    ]
    for fn,params,gid,expected in cases:
        g=fn(params)
        elems=[Element(int(i),int(j),0,float(w)) for (i,j),w in zip(g.connectivity,g.wall_weights)]
        out=analyse_geometry(g.nodes,elems,g.a1,g.a2)
        assert out.detected_geometry_id==gid
        assert 0.55 <= out.detection_confidence <= 1.0
        assert set(out.parameter_sources) == set(expected)
        assert all(out.parameter_sources[k] for k in expected)
        for key,value in expected.items():
            assert out.detected_parameters[key]==pytest.approx(value)


def test_project_wall_groups_round_trip(tmp_path):
    from rvestudio.examples_data import paper_hexagonal_case
    from rvestudio.io import save_project, load_project
    model, mats, lams = paper_hexagonal_case(1)
    model.wall_groups = ['L1','L2','L2']
    path = tmp_path/'grouped.rvejson'
    save_project(path, model, mats, lams)
    loaded, _, _ = load_project(path)
    assert loaded.wall_groups == ['L1','L2','L2']
