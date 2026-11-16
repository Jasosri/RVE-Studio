from pathlib import Path
import numpy as np

from rvestudio.parametric_geometry import scale_xy_reference_geometry, make_scaled_xy_plugin_source
from rvestudio.plugins import discover_geometry_plugins


def _reference():
    nodes=np.array([[0.,0.],[2.,0.],[0.,3.],[2.,3.]])
    conn=np.array([[0,1],[0,2],[1,3]])
    a1=np.array([2.,0.]); a2=np.array([0.,3.])
    groups=['h','v1','v2']; weights=np.ones(3)
    return nodes,conn,a1,a2,groups,weights


def test_scale_xy_reference_geometry_preserves_topology_and_target_extents():
    nodes,conn,a1,a2,groups,weights=_reference()
    g=scale_xy_reference_geometry(nodes,conn,a1,a2,groups,weights,geometry_id='X',name='X',width=4,height=6)
    assert np.allclose(np.ptp(g.nodes,axis=0),[4.,6.])
    assert np.array_equal(g.connectivity,conn)
    assert np.allclose(g.a1,[4.,0.])
    assert np.allclose(g.a2,[0.,6.])
    assert g.wall_groups==groups


def test_wizard_plugin_source_is_discoverable_and_parametric(tmp_path: Path):
    nodes,conn,a1,a2,groups,weights=_reference()
    source=make_scaled_xy_plugin_source(
        geometry_id='WIZTEST',name='Wizard test',nodes=nodes,connectivity=conn,a1=a1,a2=a2,
        wall_groups=groups,wall_weights=weights,width_parameter='W',height_parameter='H')
    (tmp_path/'wiztest.py').write_text(source,encoding='utf-8')
    registry=discover_geometry_plugins(tmp_path)
    assert 'WIZTEST' in registry
    g=registry['WIZTEST']({'W':5.0,'H':9.0})
    assert np.allclose(np.ptp(g.nodes,axis=0),[5.0,9.0])
    assert np.allclose(g.a1,[5.0,0.0])
    assert np.allclose(g.a2,[0.0,9.0])
