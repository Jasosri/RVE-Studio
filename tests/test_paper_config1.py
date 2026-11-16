from rvestudio.examples_data import paper_hexagonal_case
from rvestudio.homogenization import homogenize_rve
from rvestudio.validation import extract_six, PAPER_PUBLISHED
import numpy as np


def test_paper_configuration_1():
    model, mats, lams = paper_hexagonal_case(1)
    assert model.wall_element_type == 'paper'
    out = homogenize_rve(model,mats,lams)
    vals = extract_six(out)
    # Eq (5) order [13,23]; this symmetric honeycomb happens to give equal shear.
    assert np.max(np.abs(100*(vals - PAPER_PUBLISHED[0])/PAPER_PUBLISHED[0])) < 2.0
    assert abs(out.relative_density - 0.0722) < 1e-12
