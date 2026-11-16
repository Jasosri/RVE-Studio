from __future__ import annotations
import warnings
import numpy as np
from scipy.sparse import csc_matrix, bmat
from scipy.sparse.linalg import splu
from typing import Mapping
from .models import Material,Laminate,RVEModel,HomogenizationResult
from .laminate import laminate_section
from .walls import wall_element_paper,wall_element_anisotropic_strip,wall_transform,is_timoshenko
from .periodic import build_periodic_constraints
from .properties import effective_properties
from .physical import compute_physical_properties


DENSE_KKT_LIMIT = 600


def _dense_kkt_solve(Z,RHS):
    sol=np.linalg.solve(Z,RHS)
    # np.linalg.solve accepts a numerically singular matrix silently; reject it as the sparse LU would.
    if not np.all(np.isfinite(sol)) or np.linalg.norm(Z@sol-RHS)>1e-8*max(1.0,np.linalg.norm(RHS))*max(1.0,np.abs(Z).max()):
        raise np.linalg.LinAlgError("matrix is exactly singular")
    return sol


def homogenize_rve(model:RVEModel,materials:Mapping[str,Material],laminates:list[Laminate])->HomogenizationResult:
    nodes=model.nodes; n_node=nodes.shape[0]; ndof=4*n_node; area=abs(np.linalg.det(np.column_stack((model.a1,model.a2)))); volume=area*model.height
    sections=[]
    for lam in laminates:
        sections.append(lam.section_override if lam.section_override is not None else laminate_section(lam,materials))
    K=np.zeros((ndof,ndof)); H=np.zeros((ndof,6)); R=np.zeros((6,6)); element_data=[]
    for ei,e in enumerate(model.elements):
        if e.node_i<0 or e.node_i>=n_node or e.node_j<0 or e.node_j>=n_node or e.node_i==e.node_j: raise ValueError(f"Invalid connectivity in element {ei+1}")
        if e.laminate<0 or e.laminate>=len(laminates): raise ValueError(f"Invalid laminate ID in element {ei+1}")
        dxy=nodes[e.node_j]-nodes[e.node_i]; L=float(np.hypot(*dxy)); alpha=float(np.arctan2(dxy[1],dxy[0])); sec=sections[e.laminate]
        typ=model.wall_element_type.lower().replace('_',''); shear=is_timoshenko(getattr(model,'wall_shear','euler-bernoulli'))
        if typ in ('anisotropicstrip','anisotropic','fullanisotropic'):
            Kloc,g,r,meta=wall_element_anisotropic_strip(sec,L,model.height,model.strip_bending_mode,shear,getattr(model,'twist_release',None))
        elif typ in ('paper','thesis','original','cylindrical'):
            Kloc,g,r,meta=wall_element_paper(sec.A,sec.D,sec.Seq,L,model.height,shear)
        else: raise ValueError(f"Unknown wall_element_type: {model.wall_element_type}")
        T=wall_transform(alpha); Kg=T@Kloc@T.T; gg=T@g
        edofs=np.r_[np.arange(4*e.node_i,4*e.node_i+4),np.arange(4*e.node_j,4*e.node_j+4)]
        K[np.ix_(edofs,edofs)]+=e.weight*Kg; H[edofs,2]+=e.weight*gg; R[2,2]+=e.weight*r
        element_data.append({"L":L,"alpha":alpha,"edofs":edofs,"laminate":e.laminate,"weight":e.weight,
            "Kglobal":Kg,"Klocal":Kloc,"glocal":g,"rlocal":r,"T":T,"strip_meta":meta})
    C,Dmacro,cinfo=build_periodic_constraints(model); ncon=C.shape[0]
    RHS=np.vstack((-H,Dmacro))
    Cd=C.toarray() if hasattr(C,'toarray') else np.asarray(C,float)
    try:
        if ndof+ncon<=DENSE_KKT_LIMIT:
            # Small RVEs (every practical cell): a dense LU is several times faster than building a sparse system.
            Z=np.zeros((ndof+ncon,ndof+ncon)); Z[:ndof,:ndof]=K; Z[:ndof,ndof:]=Cd.T; Z[ndof:,:ndof]=Cd
            SOL=_dense_kkt_solve(Z,RHS)
        else:
            Z=bmat([[csc_matrix(K),csc_matrix(Cd.T)],[csc_matrix(Cd),csc_matrix((ncon,ncon))]],format='csc')
            lu=splu(Z); SOL=lu.solve(RHS)
    except Exception as exc:
        raise np.linalg.LinAlgError(f"RVE KKT solve failed; check connectivity/periodicity: {exc}") from exc
    U=SOL[:ndof]; Lam=SOL[ndof:]
    KH=U.T@K@U+U.T@H+H.T@U+R; CH=KH/volume; CH=0.5*(CH+CH.T)
    S=np.linalg.solve(CH,np.eye(6)); S=0.5*(S+S.T); props=effective_properties(CH); phys=compute_physical_properties(model,materials,laminates)
    return HomogenizationResult(CH=CH,S=S,props=props,Umacro=U,Lambda=Lam,K=K,H=H,R=R,constraints_C=C,constraints_Dmacro=Dmacro,
        constraint_info=cinfo,sections=sections,element_data=element_data,area=area,volume=volume,physical=phys)
