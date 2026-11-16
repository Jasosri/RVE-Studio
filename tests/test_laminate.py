import numpy as np
from rvestudio.examples_data import default_materials
from rvestudio.models import Ply,Laminate
from rvestudio.laminate import analyze_laminate


def test_laminate_tools_reference():
    mats=default_materials(); tp=0.625270341532304
    lam=Laminate([Ply('E-glass/Epoxy',a,tp) for a in [45,-45,0,-45,45]])
    out=analyze_laminate(lam,mats)
    assert np.linalg.norm(out['B']) < 1e-8
    assert abs(out['D'][0,0]-44045.585835) < 1e-3
    assert abs(out['strip_free'].Dstrip-31289.252167) < 1e-3
