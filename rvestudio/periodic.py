from __future__ import annotations
import numpy as np
from scipy.linalg import qr
from .models import RVEModel


def _independent_rows(C:np.ndarray,D:np.ndarray):
    if C.shape[0]==0: raise ValueError("Constraint matrix is empty")
    _,R,piv=qr(C.T,pivoting=True,mode='economic')
    dr=np.abs(np.diag(R))
    if dr.size==0: raise ValueError("Constraint matrix is empty")
    tol=max(C.T.shape)*np.finfo(float).eps*max(dr.max(),1.0)
    rank=int(np.sum(dr>tol)); keep=np.sort(piv[:rank])
    return C[keep,:],D[keep,:],keep


def build_periodic_constraints(model:RVEModel):
    nodes=model.nodes; n_node=nodes.shape[0]; ndof=4*n_node; Acell=np.column_stack((model.a1,model.a2)); pairs=[]
    pp=model.periodic_pairs
    if pp is not None and np.asarray(pp).size:
        pp=np.asarray(pp,float)
        if pp.ndim!=2 or pp.shape[1] not in (2,4): raise ValueError("periodic_pairs must have 2 or 4 columns")
        if pp.shape[1]==2:
            for row in pp:
                m,s=map(int,row[:2]); dx=nodes[s]-nodes[m]; shift=np.linalg.solve(Acell,dx)
                pairs.append([m,s,*shift])
        else:
            pairs=pp.tolist()
    else:
        xi=np.linalg.solve(Acell,nodes.T).T; assigned=np.zeros(n_node,bool)
        for i in range(n_node):
            if assigned[i]: continue
            group=[i]; assigned[i]=True
            for j in range(i+1,n_node):
                if assigned[j]: continue
                dd=xi[j]-xi[i]
                if np.linalg.norm(dd-np.rint(dd),ord=np.inf)<=model.periodic_tol:
                    group.append(j); assigned[j]=True
            if len(group)>1:
                master=group[0]
                for slave in group[1:]:
                    shift=np.rint(xi[slave]-xi[master]).astype(int)
                    pairs.append([master,slave,*shift])
    if not pairs: raise ValueError("No periodic node pairs found; provide duplicate boundary nodes or explicit periodic_pairs")
    Crows=[]; Drows=[]
    for row in pairs:
        m,s=int(row[0]),int(row[1]); nshift=np.array(row[2:4],float); dx=Acell@nshift
        if pp is not None and np.asarray(pp).shape[1]==2: dx=nodes[s]-nodes[m]
        r=np.zeros(ndof); r[4*s]=1; r[4*m]=-1; Crows.append(r); Drows.append([dx[0],0,0,0,0,dx[1]/2])
        r=np.zeros(ndof); r[4*s+1]=1; r[4*m+1]=-1; Crows.append(r); Drows.append([0,dx[1],0,0,0,dx[0]/2])
        r=np.zeros(ndof); r[4*s+2]=1; r[4*m+2]=-1; Crows.append(r); Drows.append([0,0,0,dx[1],dx[0],0])
        r=np.zeros(ndof); r[4*s+3]=1; r[4*m+3]=-1; Crows.append(r); Drows.append([0,0,0,0,0,0])
    anchor=int(model.anchor_node)
    for d in range(3):
        r=np.zeros(ndof); r[4*anchor+d]=1; Crows.append(r); Drows.append([0]*6)
    if model.fix_anchor_rotation:
        r=np.zeros(ndof); r[4*anchor+3]=1; Crows.append(r); Drows.append([0]*6)
    C=np.vstack(Crows); D=np.asarray(Drows,float); C,D,keep=_independent_rows(C,D)
    info={"pairs":np.asarray(pairs,float),"anchor_node":anchor,"n_constraints":C.shape[0],"rank_C":C.shape[0],"kept_constraint_rows":keep}
    return C,D,info
