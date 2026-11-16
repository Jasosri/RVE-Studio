"""Parameter diagrams of the Geometry Library (v1.12)."""
import json
from pathlib import Path

import numpy as np
import pytest
from matplotlib.figure import Figure

from rvestudio.geometry import DEFAULT_GEOMETRY_LIBRARY
from rvestudio.geometry_diagram import (entry_diagram, builtin_guide_figure, caption_html, reference_values,
                                        library_geometry, parameter_walls, save_library_diagrams)
from rvestudio.parametric_geometry import scale_xy_reference_geometry
from rvestudio.plugins import discover_geometry_plugins

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('entry', DEFAULT_GEOMETRY_LIBRARY, ids=lambda e: e['id'])
def test_builtin_families_draw_every_parameter_on_its_walls(entry):
    geom, vals = library_geometry(entry)
    kind, walls, fam = parameter_walls(geom, vals)
    assert kind == 'walls' and fam == entry['id']
    assert set(walls) == {p['name'] for p in entry['parameters']}
    assert all(walls.values())
    # Every wall named for a length parameter really has that length.
    for name in ('L1', 'L2', 'Lx', 'Ly', 'Lh'):
        for w in walls.get(name, []):
            i, j = geom.connectivity[w]
            assert np.linalg.norm(geom.nodes[j] - geom.nodes[i]) == pytest.approx(vals[name])
    fig, info = entry_diagram(entry)
    assert [n for n, _t in info['parameters']] == [p['name'] for p in entry['parameters']]
    html = caption_html(info)
    assert all(f"<b>{p['name']}</b>" in html for p in entry['parameters']) and 'range' in html


def test_reference_values_use_reference_or_mid_range():
    assert reference_values({'parameters': [{'name': 'A', 'min': 1, 'max': 3}, {'name': 'B', 'min': 0, 'max': 9,
                                                                                 'reference': 2}]}) == {'A': 2.0, 'B': 2.0}


def test_wizard_width_height_family_and_plugins(tmp_path):
    reg = discover_geometry_plugins(ROOT / 'plugins')
    g = reg['NODESELECT']({'NodeCount': 6})
    reg['WH1'] = lambda p: scale_xy_reference_geometry(g.nodes, g.connectivity, g.a1, g.a2, g.wall_groups,
                                                       g.wall_weights, geometry_id='WH1', name='Wizard',
                                                       width=p['Width'], height=p['Height'])
    _fig, info = entry_diagram({'id': 'WH1', 'name': 'Wizard', 'parameters': [
        {'name': 'Width', 'min': 80, 'max': 160}, {'name': 'Height', 'min': 60, 'max': 140}]}, reg)
    assert info['kind'] == 'box' and [n for n, _ in info['parameters']] == ['Width', 'Height']
    lib = json.loads((ROOT / 'data' / 'geometries.json').read_text())
    paths = save_library_diagrams(lib, tmp_path, reg)
    assert len(paths) == len(lib) and all(p.stat().st_size > 5000 for p in paths)


def test_unknown_generator_is_a_clear_error():
    with pytest.raises(ValueError, match='No generator'):
        library_geometry({'id': 'NOPE', 'parameters': []})


def test_guide_figure_has_four_panels():
    fig = builtin_guide_figure(Figure(figsize=(8, 6)))
    assert len([a for a in fig.axes if a.get_title()]) == 4


def test_plugin_walls_are_coloured_by_the_parameter_named_in_their_group():
    from pathlib import Path
    from rvestudio.plugins import discover_geometry_plugins
    from rvestudio.geometry_diagram import parameter_walls
    reg = discover_geometry_plugins(Path(__file__).resolve().parents[1] / "plugins")
    kind, walls, fam = parameter_walls(reg["KAGOME"]({"L1": 50, "L2": 50}), {"L1": 50, "L2": 50})
    assert kind == "walls" and walls == {"L1": [0, 1], "L2": [2, 3, 4, 5]} and fam == "KAGOME"
