"""User-defined FE RVE (v1.12): a block of n1 x n2 cells, or a rectangle drawn over the replicated design."""
import copy
import json

import numpy as np
import pytest

from rvestudio.geometry import BUILTIN_GEOMETRIES
from rvestudio.models import RVEModel, Element
from rvestudio.examples_data import paper_hexagonal_case
from rvestudio.io import save_project, load_project
from rvestudio.fe_rve import (select_fe_rve, fe_rve_layout, window_from_drag, centred_window, spec_rectangle,
                              normalise_fe_rve_spec, describe_fe_rve_spec, snap_to_cell, fe_rve_summary, MAX_RVE_CELLS)
from rvestudio.abaqus import _divisions
from rvestudio.abaqus_homogenization import (PeriodicHomogenizationSettings, prepare_periodic_homogenization,
                                             build_published_periodic_domain, PBC_DOMAIN_PUBLISHED, CH_METHOD_ENERGY)
from rvestudio.abaqus_thesis_dbc import prepare_thesis_dbc, ThesisDBCSettings
from rvestudio.verification import builtin_check
from rvestudio import fe_workflow as fw


def _model(gid, par, h=150.0):
    g = BUILTIN_GEOMETRIES[gid](par)
    els = [Element(int(i), int(j), 0, float(w)) for (i, j), w in zip(g.connectivity, g.wall_weights)]
    return RVEModel(g.nodes, els, g.a1, g.a2, h, name=gid, wall_groups=g.wall_groups)


FAMILIES = [('HEX', {'L1': 50, 'L2': 50}), ('SQUARE', {'Lx': 40, 'Ly': 60}), ('TRI', {'Lh': 50, 'R': 1.0}),
            ('MIXED', {'Lh': 40, 'R': 2.0})]


# ------------------------------------------------------------------ specification
def test_spec_validation_and_labels():
    assert normalise_fe_rve_spec(None) is None
    assert normalise_fe_rve_spec({'mode': 'auto'}) is None
    assert normalise_fe_rve_spec({'mode': 'cells', 'cells': [3, 2]}) == {'mode': 'cells', 'cells': [3, 2]}
    w = normalise_fe_rve_spec({'mode': 'window', 'cells': (2.0, 1), 'centre': (1, 2)})
    assert w == {'mode': 'window', 'cells': [2, 1], 'centre': [1.0, 2.0]}
    for bad in ({'mode': 'x'}, {'mode': 'cells', 'cells': [0, 1]}, {'mode': 'cells', 'cells': [MAX_RVE_CELLS + 1, 1]},
                {'mode': 'window', 'cells': [1, 1]}, {'mode': 'window', 'cells': [1, 1], 'centre': [float('nan'), 0]}):
        with pytest.raises(ValueError):
            normalise_fe_rve_spec(bad)
    assert describe_fe_rve_spec(None) == 'automatic'
    assert describe_fe_rve_spec(w).startswith('drawn 2 x 1')


def test_spec_is_saved_with_the_project(tmp_path):
    m, mats, lams = paper_hexagonal_case(1)
    m.fe_rve = {'mode': 'window', 'cells': [3, 1], 'centre': [86.6, 100.0]}
    save_project(tmp_path / 'p.rvejson', m, mats, lams)
    m2, _, _ = load_project(tmp_path / 'p.rvejson')
    assert m2.fe_rve == m.fe_rve
    m.fe_rve = None
    save_project(tmp_path / 'q.rvejson', m, mats, lams)
    assert 'fe_rve' not in json.loads((tmp_path / 'q.rvejson').read_text())
    assert load_project(tmp_path / 'q.rvejson')[0].fe_rve is None


# ------------------------------------------------------------------ selection
@pytest.mark.parametrize('gid,par', FAMILIES)
def test_drawn_rectangle_is_whole_cells_centred_on_the_geometry(gid, par):
    m = _model(gid, par)
    lay = fe_rve_layout(m)
    rw, rh = lay['rect_width'], lay['rect_height']
    c = lay['auto_centre']
    # A sloppy drag of about 2.7 x 1.2 cells, off any cell centre.
    drag = window_from_drag(lay, c[0] - 1.1 * rw, c[1] - 0.4 * rh, c[0] + 1.6 * rw, c[1] + 0.8 * rh, snap_centre=False)
    assert drag['cells'] == [3, 1]
    spec, note = centred_window(m, lay, drag['cells'], drag['centre'])
    assert 'mirror-symmetric' in note                                            # centred on the geometry
    assert spec['centre'] in [[x, y] for x, y, _k in lay['centre_candidates']]
    m.fe_rve = spec
    _, d = select_fe_rve(m)
    assert d['rve_mode'] == 'window' and d['selection']['user_defined']
    assert d['p'] == pytest.approx(3 * rw) and d['q'] == pytest.approx(rh)
    assert d['primitive_cells'] == 3 * lay['rectangular_cell']['primitive_cells']
    assert d['wall_measure_relative_error'] < 1e-7                               # same material as 3 x 1 cells
    assert d['mirror_symmetry']['about_x1_mid'] and d['mirror_symmetry']['about_x2_mid']
    assert d['central_cell_polygon_rve']                                          # a complete cell at the centre
    assert not d['selection']['warnings']
    x0, y0, w, h = spec_rectangle(lay, spec)
    assert [x0, y0] == pytest.approx(d['window_origin_local']) and (w, h) == pytest.approx((d['p'], d['q']))
    assert fe_rve_summary(d).startswith('User-drawn FE RVE')


def test_free_centre_is_kept_and_warned_about():
    m, _, _ = paper_hexagonal_case(1)
    lay = fe_rve_layout(m)
    rw, rh = lay['rect_width'], lay['rect_height']
    spec, note = centred_window(m, lay, [2, 1], [13 + rw, 7 + rh / 2], snap=False)
    assert spec['centre'] == [13 + rw, 7 + rh / 2] and 'kept' in note
    m.fe_rve = spec
    _, d = select_fe_rve(m)
    assert any('mirror-symmetric' in w for w in d['selection']['warnings'])
    snapped, note = centred_window(m, lay, [2, 1], [13 + rw, 7 + rh / 2])
    assert 'Centred on the geometry' in note
    m.fe_rve = snapped
    assert not select_fe_rve(m)[1]['selection']['warnings']


def test_cells_mode_and_model_spec_reach_every_entry_point():
    m, _, _ = paper_hexagonal_case(1)
    m.fe_rve = {'mode': 'cells', 'cells': [3, 2]}
    _, d = select_fe_rve(m)
    assert d['cells_x1_x2'] == [3, 2] and d['rve_mode'] == 'cells'
    _, auto = select_fe_rve(m, use_model_spec=False)
    assert auto['cells_x1_x2'] == [2, 1] and auto['rve_mode'] == 'auto'
    _, forced = select_fe_rve(m, cells=(1, 1))                    # explicit arguments win over the design's choice
    assert forced['cells_x1_x2'] == [1, 1]
    dom, dd = build_published_periodic_domain(m)                  # the PBC domain follows the design's choice
    assert dd['repeat_x'] == 3 and dd['repeat_y'] == 2
    lay = fe_rve_layout(m)
    x0, y0, w, h = spec_rectangle(lay, m.fe_rve)
    assert [x0, y0] == pytest.approx(d['window_origin_local'])


def test_window_needs_a_size():
    m, _, _ = paper_hexagonal_case(1)
    with pytest.raises(ValueError):
        select_fe_rve(m, window=(0.0, 0.0))


# ------------------------------------------------------------------ FE: periodic results do not depend on the RVE
def test_pbc_tensor_is_independent_of_the_drawn_rve_and_dbc_matches_when_centred():
    m, mats, lams = paper_hexagonal_case(1)
    lay = fe_rve_layout(m)
    rw, rh = lay['rect_width'], lay['rect_height']
    ref = builtin_check(m, mats, lams, 10.0)
    for spec in (centred_window(m, lay, [2, 2], [0, 0])[0],                            # centred, 2 x 2
                 centred_window(m, lay, [2, 1], [13 + rw, 7 + rh / 2], snap=False)[0]):  # off-centre
        mm = copy.deepcopy(m)
        mm.fe_rve = spec
        out = builtin_check(mm, mats, lams, 10.0, models=('PBC',))
        for t, v in ref['PBC']['C'].items():
            assert out['PBC']['C'][t] == pytest.approx(v, rel=2e-3, abs=0.05), t
    mm = copy.deepcopy(m)
    mm.fe_rve = centred_window(m, lay, [2, 1], [rw, rh / 2])[0]   # centred 2 x 1 = the automatic RVE
    out = builtin_check(mm, mats, lams, 10.0, models=('DBC',))
    for t, v in ref['DBC']['Q'].items():
        assert out['DBC']['Q'][t] == pytest.approx(v, rel=1e-6, abs=1e-6), t


def test_mesh_divisions_ignore_round_off():
    assert _divisions(50.0, 10.0) == 5
    assert _divisions(50.00000000000003, 10.0) == 5
    assert _divisions(50.001, 10.0) == 6
    assert _divisions(1e-3, 10.0) == 1


@pytest.mark.parametrize('cells', [(2, 1), (3, 1), (3, 2)])
def test_abaqus_decks_build_for_drawn_rves(tmp_path, cells):
    m, mats, lams = paper_hexagonal_case(1)
    lay = fe_rve_layout(m)
    m.fe_rve = centred_window(m, lay, list(cells), [0.0, 0.0])[0]
    ps = PeriodicHomogenizationSettings(mesh_size=10, identification_method=CH_METHOD_ENERGY,
                                        domain_mode=PBC_DOMAIN_PUBLISHED, strict_pbc_audit=True)
    prepare_periodic_homogenization(tmp_path / 'pbc', m, mats, lams, ps, run_name='t', resume=True)  # audit passes
    files = prepare_thesis_dbc(tmp_path / 'dbc', m, mats, lams, ThesisDBCSettings(mesh_size=10), run_name='t')
    dom = json.loads(files.snapshot_file.read_text())['domain']
    assert dom['rve_mode'] == 'window' and dom['cells_x1_x2'] == list(cells)


def test_study_keeps_designs_with_different_rves_apart(tmp_path):
    m, mats, lams = paper_hexagonal_case(1)
    save_project(tmp_path / 'a.rvejson', m, mats, lams)
    m.fe_rve = {'mode': 'cells', 'cells': [3, 2]}
    (tmp_path / 'b').mkdir()
    save_project(tmp_path / 'b' / 'a.rvejson', m, mats, lams)
    study = fw.add_projects(tmp_path / 's', [tmp_path / 'a.rvejson', tmp_path / 'b' / 'a.rvejson'],
                            settings=fw.FEStudySettings(meshes=[10.0], checks=['pbc']))
    assert len(study['projects']) == 2
    specs = [load_project(tmp_path / 's' / v['project'])[0].fe_rve for v in study['projects'].values()]
    assert None in specs and {'mode': 'cells', 'cells': [3, 2]} in specs


def test_cli_sets_and_saves_the_rve(tmp_path, capsys):
    from rvestudio.cli import main
    out = tmp_path / 'c.rvejson'
    assert main(['fe', 'rve', '--paper-config', '1', '--cells', '3', '1', '--centre', '80', '95', '--save', str(out)]) == 0
    spec = load_project(out)[0].fe_rve
    assert spec['mode'] == 'window' and spec['cells'] == [3, 1]
    assert 'nearest cell centre' in capsys.readouterr().out
    assert main(['fe', 'rve', '--project', str(out), '--auto', '--save', str(out)]) == 0
    assert load_project(out)[0].fe_rve is None
