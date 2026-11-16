from __future__ import annotations
import numpy as np
from .models import EffectiveProperties


def effective_properties(CH:np.ndarray)->EffectiveProperties:
    S=np.linalg.solve(CH,np.eye(6)); S=0.5*(S+S.T)
    ip=[0,1,5]; Qp=np.linalg.solve(S[np.ix_(ip,ip)],np.eye(3)); Qp=0.5*(Qp+Qp.T)
    # Paper (Sriharan et al., 2024), Eq. (5): S44 = C1313 and S55 = C2323.
    # Internal CH uses engineering Voigt [11,22,33,23,13,12], hence [4,3].
    ish=[4,3]; Qs=np.linalg.solve(S[np.ix_(ish,ish)],np.eye(2)); Qs=0.5*(Qs+Qs.T)
    mask=np.ones((6,6),bool); mask[0:3,0:3]=False; mask[3:6,3:6]=False
    norm=float(np.linalg.norm(CH[mask])/max(np.linalg.norm(CH),np.finfo(float).eps))
    return EffectiveProperties(E1=1/S[0,0],E2=1/S[1,1],E3=1/S[2,2],G23=1/S[3,3],G13=1/S[4,4],G12=1/S[5,5],
        nu12=-S[1,0]/S[0,0],nu21=-S[0,1]/S[1,1],nu13=-S[2,0]/S[0,0],nu31=-S[0,2]/S[2,2],
        nu23=-S[2,1]/S[1,1],nu32=-S[1,2]/S[2,2],Qplane=Qp,Qshear=Qs,normal_shear_coupling_norm=norm)
