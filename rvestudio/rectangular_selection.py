"""Checks geometric mirror symmetry of an exact, commensurate rectangular FE cut.

A rectangular outer boundary and periodicity do NOT imply reflection symmetry.
This module reports the geometry check separately from the laminate symmetry check.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .models import RVEModel, Laminate
from .abaqus_homogenization import build_rectangular_supercell


def _signature(lam: Laminate) -> tuple:
    return tuple((p.material,round(float(p.angle_deg),9),round(float(p.thickness),9)) for p in lam.plies)


def _reflection_error(model: RVEModel, laminate_signatures: list[tuple] | None = None,
                      *, samples_per_wall: int = 6) -> dict:
    """Numerical wall-set Hausdorff check under two axis reflections, periodic modulo B.

    Segments can be clipped at the rectangular boundary and do not have to share
    matching endpoint indices. The error is measured in physical mm.
    """
    B=np.column_stack((np.asarray(model.a1,float),np.asarray(model.a2,float)))
    if abs(float(B[:,0] @ B[:,1])) > 1e-7 * np.linalg.norm(B[:,0])*np.linalg.norm(B[:,1]):
        raise ValueError('Symmetry check requires an orthogonal FE cell')
    uv=np.linalg.solve(B,np.asarray(model.nodes,float).T).T
    start=np.array([uv[e.node_i] for e in model.elements]); end=np.array([uv[e.node_j] for e in model.elements]);
    # Include periodic copies; match real segment distance, never a closest-node shortcut.
    shifts=np.array([(i,j) for i in (-1,0,1) for j in (-1,0,1)],float)
    dst0=(start[None,:,:]+shifts[:,None,:]).reshape(-1,2)
    dst1=(end[None,:,:]+shifts[:,None,:]).reshape(-1,2)
    labels=np.asarray([e.laminate for e in model.elements],int)
    label_targets=np.tile(labels,len(shifts))
    sigs=laminate_signatures
    physical_scale=max(np.linalg.norm(B[:,0]),np.linalg.norm(B[:,1]))
    dvec=(dst1-dst0)@B.T
    denom=np.einsum('ij,ij->i',dvec,dvec)
    errors={}
    for ax in ('x','y'):
        max_geometry=0.; max_mat=0.
        for k,(a,b) in enumerate(zip(start,end)):
            for t in np.linspace(.07,.93,samples_per_wall):
                u=(1-t)*a+t*b
                u[0 if ax=='x' else 1]*=-1
                q=u@B.T
                dv=q-dst0@B.T
                s=np.clip(np.einsum('ij,ij->i',dv,dvec)/np.maximum(denom,1e-30),0,1)
                dist=np.linalg.norm(dv-s[:,None]*dvec,axis=1)
                max_geometry=max(max_geometry,float(dist.min()))
                if sigs is not None:
                    compatible=np.array([sigs[j]==sigs[int(labels[k])] for j in label_targets],bool)
                    max_mat=max(max_mat,float(dist[compatible].min()) if compatible.any() else float('inf'))
        errors[ax+'_geometry_max_mm']=max_geometry
        errors[ax+'_geometry_relative']=max_geometry/physical_scale
        if sigs is not None:
            errors[ax+'_same_stack_max_mm']=max_mat
            errors[ax+'_same_stack_relative']=max_mat/physical_scale
    return errors


def inspect_rectangular_rve(model: RVEModel, laminates: list[Laminate] | None = None,
                            *, repeats_x: int = 1,repeats_y: int = 1,
                            search_limit: int = 6,
                            prefer_mirror_symmetry: bool = True,
                            symmetry_tolerance: float = 1e-6) -> tuple[RVEModel,dict]:
    rectangular,diag=build_rectangular_supercell(model,repeats_x=repeats_x,repeats_y=repeats_y,search_limit=search_limit,
                                                prefer_mirror_symmetry=prefer_mirror_symmetry)
    sigs=[_signature(l) for l in laminates] if laminates is not None else None
    err=_reflection_error(rectangular,sigs)
    g=all(err[ax+'_geometry_relative']<=symmetry_tolerance for ax in ('x','y'))
    s=all(err.get(ax+'_same_stack_relative',float('inf'))<=symmetry_tolerance for ax in ('x','y')) if sigs is not None else None
    # Mirror behaviour of inclined fibres has NOT been verified by this stack-identity test.
    diag={**diag,**err,'mirror_symmetry_tolerance':symmetry_tolerance,
          'mirror_geometry_symmetric':bool(g),'same_stack_on_mirror_walls':s,
          'laminate_orientation_mirror_verified':False,
          'symmetry_scope':'geometry and same numerical stacking sequences only; ply-angle handedness under reflection is not certified'}
    return rectangular,diag
