from rvestudio.examples_data import default_materials
from rvestudio.models import Ply,Laminate
from rvestudio.laminate import laminate_section,reduced_bending_stiffness


def test_free_moment_reduces_bending_for_benchmark_laminate():
    mats=default_materials(); tp=0.625270341532304
    sec=laminate_section(Laminate([Ply('E-glass/Epoxy',a,tp) for a in [45,-45,0,-45,45]]),mats)
    cyl=reduced_bending_stiffness(sec,'cylindrical'); free=reduced_bending_stiffness(sec,'freeMoment')
    assert free.Dstrip < cyl.Dstrip
    assert 0.70 < free.Dstrip/cyl.Dstrip < 0.72
