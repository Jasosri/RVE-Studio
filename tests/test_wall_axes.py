"""Wall direction: Node I -> Node J defines the positive local 1bar axis (v1.12.1)."""
import copy
import math

import numpy as np
import pytest
from matplotlib.figure import Figure

from rvestudio.examples_data import paper_hexagonal_case
from rvestudio.homogenization import homogenize_rve
from rvestudio.models import Element
from rvestudio.validation import PAPER_PUBLISHED, extract_six
from rvestudio.wall_axes import (default_honeycomb_orientation, default_orientation_text, draw_wall_axes,
                                 html_note, long_note, reverse_walls, wall_frames)


def test_default_honeycomb_walls_point_into_the_junction():
    assert default_honeycomb_orientation() == [(1, 5), (2, 5), (4, 5)]
    t = default_orientation_text()
    assert 'W1 = N1 → N5' in t and 'junction node 5' in t
    assert 'Node I → Node J' in html_note() and 'overline' in html_note()
    assert '0.3 %' in long_note()


def test_wall_frames_follow_node_i_to_node_j():
    nodes = [[0, 0], [10, 0], [0, 10]]
    f = wall_frames(nodes, [(1, 2), (2, 1), (1, 3), (1, 1), (1, 9)])
    assert f[0]['angle_deg'] == pytest.approx(0) and np.allclose(f[0]['e1'], [1, 0]) and np.allclose(f[0]['e2'], [0, 1])
    assert f[1]['angle_deg'] == pytest.approx(180) and np.allclose(f[1]['e2'], [0, -1])     # reversed: 1bar and 2bar flip
    assert f[2]['angle_deg'] == pytest.approx(90) and f[2]['length'] == pytest.approx(10)
    assert not f[3]['ok'] and not f[4]['ok']                                                  # zero length, missing node


def test_frames_match_the_solver_rotation():
    """The drawn 1bar angle is the angle the wall element is rotated by (atan2 of J - I)."""
    m, _mats, _lams = paper_hexagonal_case(1)
    f = wall_frames(m.nodes, [(e.node_i + 1, e.node_j + 1) for e in m.elements])
    for fr, e in zip(f, m.elements):
        d = m.nodes[e.node_j] - m.nodes[e.node_i]
        assert math.radians(fr['angle_deg']) == pytest.approx(math.atan2(d[1], d[0]))


def test_reverse_walls():
    rows = [[1, 5, 1.0, 1, 'a'], [2, 5, 1.0, 1, ''], [4, 5, 1.0, 2, '']]
    out = reverse_walls(rows, [0, 2])
    assert out == [[5, 1, 1.0, 1, 'a'], [2, 5, 1.0, 1, ''], [5, 4, 1.0, 2, '']]
    assert rows[0][:2] == [1, 5]                                                              # input untouched


def _orient(m, pairs):
    m = copy.deepcopy(m)
    m.elements = [Element(i - 1, j - 1, e.laminate, e.weight) for (i, j), e in zip(pairs, m.elements)]
    return m


@pytest.mark.parametrize('config,changes', [(2, False), (3, False), (7, True), (12, True)])
def test_reversing_a_wall_matters_only_for_off_axis_stacks(config, changes):
    m, mats, lams = paper_hexagonal_case(config)
    a = homogenize_rve(m, mats, lams).CH
    b = homogenize_rve(_orient(m, [(5, 1), (2, 5), (4, 5)]), mats, lams).CH
    rel = np.abs(a - b).max() / np.abs(a).max()
    assert (rel > 1e-2) if changes else (rel < 1e-10)


def test_default_orientation_reproduces_the_paper_and_reversal_does_not():
    worst = {}
    for label, pairs in (('default', [(1, 5), (2, 5), (4, 5)]), ('reversed', [(5, 1), (2, 5), (4, 5)])):
        w = 0.0
        for c in range(1, 15):
            m, mats, lams = paper_hexagonal_case(c)
            e = 100 * (extract_six(homogenize_rve(_orient(m, pairs), mats, lams)) - PAPER_PUBLISHED[c - 1]) / PAPER_PUBLISHED[c - 1]
            w = max(w, float(np.abs(e).max()))
        worst[label] = w
    assert worst['default'] < 0.3 and worst['reversed'] > 5.0


def test_drawing_labels_every_wall(tmp_path):
    m, _mats, _lams = paper_hexagonal_case(1)
    fig = Figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    frames = draw_wall_axes(ax, m.nodes, [(1, 5), (2, 5), (4, 5), (0, 0)], m.a1, m.a2, highlight=[1])
    texts = [t.get_text() for t in ax.texts]
    assert {'W1: N1→N5', 'W2: N2→N5', 'W3: N4→N5'} <= set(texts) and not frames[3]['ok']
    fig.savefig(tmp_path / 'w.png')
