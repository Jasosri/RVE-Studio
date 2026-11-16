import math
import numpy as np

from rvestudio.models import LaminateSection
from rvestudio.buckling import laminated_plate_buckling, evaluate_wall_buckling
from rvestudio.examples_data import paper_hexagonal_case
from rvestudio.homogenization import homogenize_rve
from rvestudio.load_recovery import recover_wall_loads


def _isotropic_section(D0=1000.0, nu=0.3):
    D=np.array([[D0,nu*D0,0.0],[nu*D0,D0,0.0],[0.0,0.0,(1.0-nu)*D0/2.0]])
    return LaminateSection(A=np.eye(3)*1e5,B=np.zeros((3,3)),Bneutral=np.zeros((3,3)),
        Doriginal=D,D=D,Seq=1e4,A44=1e4,kappa_shear=1.0,d_neutral=0.0,
        thickness=1.0,z=np.array([-.5,.5]),Qbar=np.zeros((0,3,3)),symmetry_indicator=0.0)


def test_square_isotropic_uniaxial_buckling_matches_closed_form():
    sec=_isotropic_section(); a=b=100.0
    out=laminated_plate_buckling(sec,a,b,-1.0,0.0,0.0,m_max=4,n_max=4)
    expected=4.0*math.pi**2*1000.0/a**2
    assert np.isclose(out.lambda_cr,expected,rtol=1e-11,atol=1e-12)
    assert (out.dominant_m,out.dominant_n)==(1,1)


def test_buckling_factor_scales_inverse_with_load():
    sec=_isotropic_section(); a=b=100.0
    one=laminated_plate_buckling(sec,a,b,-1.0,-0.4,0.2,m_max=6,n_max=6)
    two=laminated_plate_buckling(sec,a,b,-2.0,-0.8,0.4,m_max=6,n_max=6)
    assert np.isfinite(one.lambda_cr) and np.isfinite(two.lambda_cr)
    assert np.isclose(two.lambda_cr,0.5*one.lambda_cr,rtol=1e-10,atol=1e-12)


def test_wall_buckling_uses_recovered_combined_resultants():
    model,mats,lams=paper_hexagonal_case(1)
    result=homogenize_rve(model,mats,lams)
    # Includes global normal and transverse-shear components, producing local
    # biaxial + shear resultants in the inclined wall plates.
    rec=recover_wall_loads(model,result,macro_stress=[-1.0,0.0,-1.0,0.0,0.25,0.0])
    assert any(abs(w['N13_N_per_mm'])>1e-6 for w in rec.walls)
    buck=evaluate_wall_buckling(model,result,rec,m_max=6,n_max=6,safety_factor=1.0)
    assert len(buck.walls)==len(model.elements)
    assert np.isfinite(buck.minimum_lambda) and buck.minimum_lambda>0
    assert buck.critical_wall in range(1,len(model.elements)+1)
    assert all('mode_convergence_relative_error' in w for w in buck.walls)


def test_square_isotropic_pure_shear_converges_to_classical_ssss_value():
    sec=_isotropic_section(); a=b=100.0
    out=laminated_plate_buckling(sec,a,b,0.0,0.0,1.0,m_max=10,n_max=10)
    # Classical simply-supported square-plate coefficient k_s ~= 9.34.
    reference=9.34*math.pi**2*1000.0/b**2
    assert np.isclose(out.lambda_cr,reference,rtol=3e-3)
