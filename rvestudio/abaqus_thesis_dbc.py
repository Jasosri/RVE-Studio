from __future__ import annotations

from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
import json
import math
import re
from typing import Callable

import numpy as np

from .models import RVEModel, Material, Laminate, Element
from .fe_rve import build_published_fe_rve, FACE_WALL_SECTIONS
from .abaqus import _safe_name, _ids_lines, _abaqus_material_blocks, _build_shell_mesh, _run_shell_command

# Shell terms (sigma33 = 0): thesis Table 3.4 = paper Table 3, U3 free on the top and bottom faces.
THESIS_DBC_CASES = ("Q11", "Q22", "Q12", "Q66", "Q44", "Q55")
# 3D stiffness terms (eps33 = 0): the same side-face displacements with U3 = 0 on the top and bottom faces.
# C1313 and C2323 come from the Q44 and Q55 cases (U3 is prescribed on the side faces there).
DBC_C_CASES = ("C1111", "C2222", "C1122", "C1212")
DBC_CASES = THESIS_DBC_CASES + DBC_C_CASES
SHELL_TERMS = THESIS_DBC_CASES
SOLID_TERMS = ("C1111", "C2222", "C1122", "C1212", "C1313", "C2323")
_C_SIDE_FACES = {"C1111": "Q11", "C2222": "Q22", "C1122": "Q12", "C1212": "Q66"}


@dataclass(slots=True)
class ThesisDBCSettings:
    element_type: str = "S4"
    mesh_size: float = 10.0
    displacement: float = 1.0
    cpus: int = 1
    save_odb: bool = True
    minimal_rigid_body_stabilization: bool = True
    face_wall_section: str = "half-thickness"
    rve_cells: tuple[int, int] | None = None

    def validate(self) -> None:
        if self.face_wall_section not in FACE_WALL_SECTIONS:
            raise ValueError(f"face_wall_section must be one of {FACE_WALL_SECTIONS}.")
        if self.rve_cells is not None and (len(self.rve_cells) != 2 or min(int(c) for c in self.rve_cells) < 1):
            raise ValueError("rve_cells must be two positive integers (cells along x1, x2) or None for the rule of the source study.")
        if self.element_type.upper() != "S4":
            raise ValueError("Thesis DBC validation currently supports S4 only.")
        if not math.isfinite(self.mesh_size) or self.mesh_size <= 0:
            raise ValueError("Mesh size must be positive and finite.")
        if not math.isfinite(self.displacement) or self.displacement == 0:
            raise ValueError("DBC displacement must be finite and non-zero.")
        if int(self.cpus) < 1:
            raise ValueError("CPU count must be at least 1.")


@dataclass(slots=True)
class ThesisDBCCaseFiles:
    index: int
    name: str
    case_dir: Path
    job_name: str
    input_file: Path
    extractor_file: Path
    results_file: Path
    log_file: Path


@dataclass(slots=True)
class ThesisDBCFiles:
    run_dir: Path
    snapshot_file: Path
    summary_file: Path
    log_file: Path
    cases: list[ThesisDBCCaseFiles]


def _wall_measure(model: RVEModel) -> float:
    return float(sum(np.linalg.norm(np.asarray(model.nodes[e.node_j], float) - np.asarray(model.nodes[e.node_i], float)) for e in model.elements))


def build_thesis_dbc_rectangle(model: RVEModel, *, cells: tuple[int, int] | None = None,
                               face_wall_section: str = "half-thickness") -> tuple[RVEModel, dict]:
    """Construct the source-study finite-height FE RVE used with thesis Table 3.4 / paper Table 3.

    For the honeycomb this is the p x q x h block of thesis Fig. 3.6(b) and paper
    Fig. 8(b): p = 2|a1|, q = |a1 + 2 a2| (173.2 x 150 mm for 50 mm walls), four
    primitive cells, one complete hexagon in the middle, the x1 faces through
    junction columns and the x2 faces at mid-length of the straight walls.  Walls
    lying in a face are shared with the neighbouring block and receive half of the
    laminate (``face_wall_section``).  See :mod:`rvestudio.fe_rve` for the general rule.
    """
    return build_published_fe_rve(model, cells=cells, face_wall_section=face_wall_section)


def _boundary_node_sets(nodes: dict[int,tuple[float,float,float]], p:float,q:float,h:float) -> dict[str,list[int]]:
    scale=max(abs(p),abs(q),abs(h),1.0); tol=max(1e-7*scale,1e-7)
    sets={"X1_0":[],"X1_P":[],"X2_0":[],"X2_Q":[],"X3_0":[],"X3_H":[]}
    for nid,(x,y,z) in nodes.items():
        if abs(x)<=tol: sets["X1_0"].append(nid)
        if abs(x-p)<=tol: sets["X1_P"].append(nid)
        if abs(y)<=tol: sets["X2_0"].append(nid)
        if abs(y-q)<=tol: sets["X2_Q"].append(nid)
        if abs(z)<=tol: sets["X3_0"].append(nid)
        if abs(z-h)<=tol: sets["X3_H"].append(nid)
    missing=[k for k,v in sets.items() if not v]
    if missing:
        raise ValueError("Thesis DBC mesh has empty boundary set(s): "+", ".join(missing))
    return sets


def _case_constraints(case:str, d:float) -> list[tuple[str,int,float]]:
    # Shell terms: exact translational conditions of thesis Table 3.4 (top/bottom U3 free).
    # 3D terms: the same side faces plus U3 = 0 on the top and bottom faces (eps33 = 0).
    if case in _C_SIDE_FACES:
        return _case_constraints(_C_SIDE_FACES[case],d)+[("X3_0",3,0.),("X3_H",3,0.)]
    if case=="Q11":
        return [("X1_0",1,0.),("X1_P",1,d),("X2_0",2,0.),("X2_Q",2,0.)]
    if case=="Q22":
        return [("X1_0",1,0.),("X1_P",1,0.),("X2_0",2,0.),("X2_Q",2,d)]
    if case=="Q12":
        return [("X1_0",1,0.),("X1_P",1,d),("X2_0",2,0.),("X2_Q",2,d)]
    if case=="Q66":
        return [("X1_0",1,0.),("X1_0",2,0.),("X1_P",1,0.),("X1_P",2,d),("X2_0",1,0.),("X2_Q",1,0.)]
    if case=="Q44":
        return [
            ("X1_0",1,0.),("X1_0",2,0.),("X1_0",3,0.),
            ("X1_P",1,0.),("X1_P",2,0.),("X1_P",3,d),
            ("X2_0",1,0.),("X2_0",2,0.),("X2_Q",1,0.),("X2_Q",2,0.),
            ("X3_0",1,0.),("X3_0",2,0.),("X3_H",1,0.),("X3_H",2,0.),
        ]
    if case=="Q55":
        return [
            ("X1_0",1,0.),("X1_0",2,0.),("X1_P",1,0.),("X1_P",2,0.),
            ("X2_0",1,0.),("X2_0",2,0.),("X2_0",3,0.),
            ("X2_Q",1,0.),("X2_Q",2,0.),("X2_Q",3,d),
            ("X3_0",1,0.),("X3_0",2,0.),("X3_H",1,0.),("X3_H",2,0.),
        ]
    raise ValueError(f"Unknown thesis DBC case: {case}")




def _rigid_body_row(xyz: tuple[float,float,float], dof: int) -> np.ndarray:
    """Map [Tx,Ty,Tz,Rx,Ry,Rz] to one constrained translational DOF."""
    x,y,z=(float(xyz[0]),float(xyz[1]),float(xyz[2]))
    if int(dof)==1:
        return np.array([1.,0.,0.,0.,z,-y],float)
    if int(dof)==2:
        return np.array([0.,1.,0.,-z,0.,x],float)
    if int(dof)==3:
        return np.array([0.,0.,1.,y,-x,0.],float)
    raise ValueError("Rigid-body audit only supports translational DOFs 1, 2 and 3.")

def _rank_and_null(A: np.ndarray) -> tuple[int,list[float],list[list[float]]]:
    """Rank and rigid-body null space of the (m x 6) constraint map A.

    Uses the economy SVD: only the 6 x 6 right factor is needed, so the cost is
    O(36 m) instead of O(m^3) for the full m x m left factor.
    """
    _u,s,vh=np.linalg.svd(A,full_matrices=A.shape[0]<6)
    tol=max(A.shape)*np.finfo(float).eps*(float(s[0]) if len(s) else 1.0)*100.0
    rank=int(np.sum(s>tol))
    null=vh[rank:,:] if rank<6 else np.zeros((0,6),float)
    return rank,[float(x) for x in s],[[float(x) for x in row] for row in null]

def _rigid_body_rank(nodes: dict[int,tuple[float,float,float]], bc_map: dict[tuple[int,int],float]) -> tuple[int,list[float],list[list[float]]]:
    rows=[]
    for (nid,dof) in sorted(bc_map):
        if int(dof) in (1,2,3):
            rows.append(_rigid_body_row(nodes[int(nid)],int(dof)))
    if not rows:
        return 0,[],np.eye(6).tolist()
    return _rank_and_null(np.vstack(rows))

def _mesh_connected_components(nodes: dict[int,tuple[float,float,float]], wall_elements: dict[int,list[tuple[int,tuple[int,int,int,int]]]]) -> list[list[int]]:
    adj={int(n):set() for n in nodes}
    for elems in wall_elements.values():
        for _eid,conn in elems:
            cc=[int(x) for x in conn]
            for i,a in enumerate(cc):
                for b in cc[i+1:]:
                    adj[a].add(b);adj[b].add(a)
    unseen=set(adj); comps=[]
    while unseen:
        start=min(unseen); unseen.remove(start); stack=[start]; comp=[]
        while stack:
            u=stack.pop();comp.append(u)
            for v in adj[u]:
                if v in unseen:
                    unseen.remove(v);stack.append(v)
        comps.append(sorted(comp))
    return comps

def _add_minimal_rbm_gauge(nodes: dict[int,tuple[float,float,float]], bc_map: dict[tuple[int,int],float], preferred_nodes: list[int]|None=None) -> dict:
    """Remove only rigid-body null modes by adding zero translational gauges.

    Existing prescribed-displacement conditions are left untouched. A candidate is
    accepted only when it increases the six-parameter rigid-body constraint rank.
    """
    pre_rank,svals,null=_rigid_body_rank(nodes,bc_map)
    added=[]
    rank=pre_rank
    if rank<6:
        order=[]
        for nid in list(preferred_nodes or [])+sorted(nodes):
            if int(nid) not in order: order.append(int(nid))
        rows=[_rigid_body_row(nodes[int(n)],int(d)) for (n,d) in sorted(bc_map) if int(d) in (1,2,3)]
        A=np.vstack(rows) if rows else np.zeros((0,6),float)
        N=np.asarray(null,float).reshape(-1,6)
        for nid in order:
            for dof in (1,2,3):
                key=(nid,dof)
                if key in bc_map: continue
                r=_rigid_body_row(nodes[int(nid)],dof)
                # A row raises the rank only if it is not orthogonal to the current null space.
                if N.shape[0]==0 or float(np.linalg.norm(N@r))<=1e-9*max(1.0,float(np.linalg.norm(r))):
                    continue
                A_try=np.vstack((A,r))
                new_rank,_,new_null=_rank_and_null(A_try)
                if new_rank>rank:
                    bc_map[key]=0.0
                    added.append({"node":nid,"dof":dof,"value":0.0,"rank_before":rank,"rank_after":new_rank})
                    rank=new_rank; A=A_try; N=np.asarray(new_null,float).reshape(-1,6)
                    if rank==6: break
            if rank==6: break
    post_rank,post_s,post_null=_rigid_body_rank(nodes,bc_map)
    return {
        "pre_stabilization_rank":int(pre_rank),
        "post_stabilization_rank":int(post_rank),
        "required_rank":6,
        "pre_nullspace":null,
        "post_nullspace":post_null,
        "added_gauge_constraints":added,
        "singular_values_before":svals,
        "singular_values_after":post_s,
    }

def generate_thesis_dbc_input(model:RVEModel, materials:dict[str,Material], laminates:list[Laminate], settings:ThesisDBCSettings, case:str, *, job_name:str="RVE_DBC", resolved_domain:bool=False, domain_diagnostics:dict|None=None) -> tuple[str,dict]:
    settings.validate()
    if case not in DBC_CASES:
        raise ValueError(f"Unknown thesis DBC case: {case}")
    if not resolved_domain:
        model,domain_diag=build_thesis_dbc_rectangle(model,cells=settings.rve_cells,face_wall_section=settings.face_wall_section)
    else:
        domain_diag=dict(domain_diagnostics or {})
    p=float(domain_diag.get("p",model.a1[0])); q=float(domain_diag.get("q",model.a2[1])); h=float(model.height)
    nodes,wall_elements,_,_=_build_shell_mesh(model,settings.mesh_size)
    bsets=_boundary_node_sets(nodes,p,q,h)
    factors=[float(x) for x in (domain_diag.get("wall_section_factors") or [1.0]*len(model.elements))]
    if len(factors)!=len(model.elements):
        raise ValueError("Wall section factors do not match the DBC wall list.")
    face_mode=str(domain_diag.get("face_wall_section","full"))
    face_walls=[int(i) for i in domain_diag.get("face_walls",[])]
    mats_all=dict(materials)
    if face_mode=="half-stiffness" and any(f<1.0 for f in factors):
        for key,mat in materials.items():
            dens=float(mat.density_kg_m3)
            mats_all[key+"__HALF"]=Material(name=f"{mat.name}_HALF",E1=0.5*mat.E1,E2=0.5*mat.E2,G12=0.5*mat.G12,nu12=mat.nu12,
                                           Gth=0.5*float(mat.Gth if mat.Gth is not None else mat.G12),
                                           density_kg_m3=0.5*dens if math.isfinite(dens) else math.nan)
    material_text,mat_names=_abaqus_material_blocks(mats_all)

    out=[
        "** RVE Studio Python - thesis Table 3.4 / paper Table 3 prescribed-displacement validation",
        f"** Case: {case}",
        f"** RVE: source-study finite-height block p x q x h = {p:.6g} x {q:.6g} x {h:.6g} mm "
        f"({domain_diag.get('primitive_cells','?')} primitive cells; thesis Fig. 3.6b, paper Fig. 8b)",
        f"** Walls lying in a face (shared with the next block): {face_walls if face_walls else 'none'}; section treatment: {face_mode}",
        ("** 3D stiffness term (eps33 = 0): side faces as thesis Table 3.4 case "+_C_SIDE_FACES[case]+", U3 = 0 on the top and bottom faces; shell rotations free."
         if case in _C_SIDE_FACES else
         "** Shell term (sigma33 = 0): translational DBCs of thesis Table 3.4, top/bottom U3 free; shell rotations free."),
        "** Any zero-translation gauge below only removes a rigid-body mode; it imposes no strain.",
        "*HEADING",f"RVE Studio thesis DBC - {model.name} - {case}","*NODE",
    ]
    for nid,xyz in sorted(nodes.items()): out.append(f"{nid}, {xyz[0]:.12g}, {xyz[1]:.12g}, {xyz[2]:.12g}")
    out.append("*ELEMENT, TYPE=S4")
    for wi in sorted(wall_elements):
        for eid,conn in wall_elements[wi]: out.append(f"{eid}, {conn[0]}, {conn[1]}, {conn[2]}, {conn[3]}")
    for wi,elems in wall_elements.items(): out.extend([f"*ELSET, ELSET=WALL_{wi}",_ids_lines(eid for eid,_ in elems)])
    if face_walls:
        out.extend(["*ELSET, ELSET=FACE_WALLS",_ids_lines(eid for wi in face_walls for eid,_ in wall_elements[wi])])
    for name,ids in bsets.items(): out.extend([f"*NSET, NSET={name}",_ids_lines(ids)])
    out.append(material_text)
    for wi,e in enumerate(model.elements,1):
        lam=laminates[e.laminate]
        f=factors[wi-1]
        p0=np.asarray(model.nodes[e.node_i],float);p1=np.asarray(model.nodes[e.node_j],float);dv=p1-p0;dv/=np.linalg.norm(dv)
        if f<1.0:
            out.append(f"** Wall W{wi} lies in face {domain_diag.get('face_wall_faces',{}).get(str(wi),'?')}: {face_mode} section (factor {f:g})")
        out.extend([f"*ORIENTATION, NAME=ORI_W{wi}, SYSTEM=RECTANGULAR",f"{dv[0]:.12g}, {dv[1]:.12g}, 0., 0., 0., 1.","3, 0.",f"*SHELL SECTION, ELSET=WALL_{wi}, COMPOSITE, ORIENTATION=ORI_W{wi}"])
        for pi,ply in enumerate(lam.plies,1):
            if ply.material not in mat_names: raise ValueError(f"Wall W{wi} ply {pi} references missing material {ply.material!r}.")
            thick=float(ply.thickness)*(f if face_mode=="half-thickness" else 1.0)
            mname=mat_names[ply.material+"__HALF"] if (face_mode=="half-stiffness" and f<1.0) else mat_names[ply.material]
            out.append(f"{thick:.12g}, 3, {mname}, {float(ply.angle_deg):.12g}, PLY{pi}")

    # Expand face-set conditions to unique node/DOF values so shared corners never receive conflicting duplicate constraints.
    bc_map:dict[tuple[int,int],float]={}
    for set_name,dof,val in _case_constraints(case,float(settings.displacement)):
        for nid in bsets[set_name]:
            key=(int(nid),int(dof))
            if key in bc_map and not math.isclose(bc_map[key],float(val),rel_tol=0,abs_tol=1e-12):
                raise ValueError(f"Conflicting Table 3.4 DBC at node {nid}, DOF {dof}: {bc_map[key]} vs {val}")
            bc_map[key]=float(val)
    # Audit the actual meshed shell model rather than assuming the face DBCs
    # suppress rigid-body motion. Disconnected shell components are unsafe because
    # each component could carry its own rigid-body modes.
    components=_mesh_connected_components(nodes,wall_elements)
    if len(components)!=1:
        raise ValueError(f"Thesis DBC mesh contains {len(components)} disconnected shell components; refusing to run a potentially singular model.")

    rbm_pre_rank,_,rbm_pre_null=_rigid_body_rank(nodes,bc_map)
    if settings.minimal_rigid_body_stabilization:
        # Prefer nodes on x3=0 so the gauge is easy to identify in Abaqus. The
        # greedy rank test adds a zero translation only if it removes a genuine
        # rigid-body null mode. No shell rotation is ever clamped here.
        rbm_audit=_add_minimal_rbm_gauge(nodes,bc_map,bsets.get("X3_0",[]))
    else:
        rbm_audit={
            "pre_stabilization_rank":rbm_pre_rank,"post_stabilization_rank":rbm_pre_rank,"required_rank":6,
            "pre_nullspace":rbm_pre_null,"post_nullspace":rbm_pre_null,"added_gauge_constraints":[],
        }
    rbm_audit["connected_shell_components"]=len(components)
    rbm_audit["component_node_counts"]=[len(c) for c in components]
    if int(rbm_audit["post_stabilization_rank"])<6:
        raise ValueError(
            f"Thesis DBC case {case} remains rigid-body under-constrained after minimal stabilization "
            f"(rank {rbm_audit['post_stabilization_rank']}/6)."
        )
    stabilization=(rbm_audit["added_gauge_constraints"] or None)

    out.append(f"** Rigid-body audit: rank {rbm_audit['pre_stabilization_rank']}/6 before gauge; {rbm_audit['post_stabilization_rank']}/6 after gauge")
    if stabilization:
        for g in stabilization:
            out.append(f"** RBM gauge: node {g['node']} U{g['dof']}=0; rank {g['rank_before']} -> {g['rank_after']}")
    out.extend(["*STEP, NAME=DBC, NLGEOM=NO","*STATIC","0.1, 1.0, 1e-08, 0.1","*BOUNDARY"])
    for (nid,dof),val in sorted(bc_map.items()): out.append(f"{nid}, {dof}, {dof}, {val:.16g}")
    out.extend(["*OUTPUT, FIELD, VARIABLE=PRESELECT","*ELEMENT OUTPUT","S, E","*OUTPUT, HISTORY","*ENERGY OUTPUT","ALLSE","*END STEP"])
    meta={
        "case":case,"p":p,"q":q,"height":h,"rve_area":p*q,"rve_volume":p*q*h,
        "face_walls":face_walls,"face_wall_section":face_mode,"primitive_cells":domain_diag.get("primitive_cells"),
        "boundary_sets":{k:len(v) for k,v in bsets.items()},"specified_translation_constraints":len(bc_map),
        "minimal_rigid_body_stabilization":stabilization,"rigid_body_audit":rbm_audit,"domain":domain_diag,
        "mesh_node_count":len(nodes),"mesh_element_count":sum(len(v) for v in wall_elements.values()),
        "table_3_4_statement":"all unspecified degrees of freedom remain free; shell rotations are not prescribed",
        "boundary_physics":"finite_height_core_published_DBC",
        "top_bottom_periodic":False,
    }
    return "\n".join(out)+"\n",meta


def _extractor_source() -> str:
    return r'''from __future__ import print_function
import json, os, sys, traceback
from odbAccess import openOdb
odb=None
out_path='case_results.json'
odb_path=sys.argv[-2] if len(sys.argv)>=3 else sys.argv[-1]
if len(sys.argv)>=3: out_path=sys.argv[-1]
def repo_get(repo,name):
    try:return repo[name]
    except Exception:pass
    try:
        target=str(name).upper()
        for key in repo.keys():
            if str(key).upper()==target:return repo[key]
    except Exception:pass
    return None
def write(data):
    with open(out_path,'w') as f: json.dump(data,f,indent=2)
def find_allse(step):
    vals=[]
    for rk in step.historyRegions.keys():
        reg=step.historyRegions[rk]
        for hk in reg.historyOutputs.keys():
            if str(hk).upper()=='ALLSE':
                data=reg.historyOutputs[hk].data
                if data:vals.append((str(rk),float(data[-1][1])))
    if not vals: raise RuntimeError('ALLSE history output not found')
    vals.sort(key=lambda x:abs(x[1]),reverse=True)
    return vals[0][1],vals
def write_nodal(odb,step,path='nodal_displacements.csv'):
    """Last-frame nodal translations for RVE Studio's deformed-shape plots."""
    frame=step.frames[-1]
    fo=repo_get(frame.fieldOutputs,'U')
    if fo is None: return 0
    coords={}
    for iname in odb.rootAssembly.instances.keys():
        for n in odb.rootAssembly.instances[iname].nodes:
            coords[(str(iname),int(n.label))]=n.coordinates
    count=0
    with open(path,'w') as f:
        f.write('label,x,y,z,u1,u2,u3\n')
        for v in fo.values:
            iname=str(v.instance.name) if v.instance is not None else ''
            c=coords.get((iname,int(v.nodeLabel)))
            if c is None: continue
            d=list(v.data)+[0.0,0.0,0.0]
            f.write('%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n'%(int(v.nodeLabel),c[0],c[1],c[2],d[0],d[1],d[2]))
            count+=1
    return count
try:
    odb=openOdb(path=odb_path,readOnly=True)
    step=repo_get(odb.steps,'DBC')
    if step is None: raise RuntimeError('DBC step not found')
    energy,candidates=find_allse(step)
    result={'status':'COMPLETED','odb':odb_path,'step':'DBC','strain_energy_ALLSE':float(energy),'allse_candidates':candidates,'frames':len(step.frames)}
    try:
        result['nodal_displacements']=write_nodal(odb,step)
    except Exception as nexc:
        result['nodal_displacements_error']=str(nexc)
    write(result); print(json.dumps(result,indent=2))
except Exception as exc:
    err={'status':'ERROR','odb':locals().get('odb_path',''),'error':str(exc),'traceback':traceback.format_exc()}
    try:write(err)
    except Exception:pass
    print(json.dumps(err,indent=2)); sys.exit(2)
finally:
    if odb is not None:
        try:odb.close()
        except Exception:pass
'''


def reconstruct_dbc_c(energies:dict[str,float], p:float,q:float,h:float, displacement:float=1.0) -> dict[str,float]:
    """3D stiffness terms (eps33 = 0) from the U3-restrained cases and the Q44/Q55 cases."""
    V=float(p)*float(q)*float(h); d=float(displacement)
    e11=d/float(p); e22=d/float(q); g12=d/float(p); g13=d/float(p); g23=d/float(q)
    W={k:float(energies[k]) for k in DBC_C_CASES+("Q44","Q55")}
    return {
        "C1111":2*W["C1111"]/(V*e11*e11),
        "C2222":2*W["C2222"]/(V*e22*e22),
        "C1122":(W["C1122"]-W["C1111"]-W["C2222"])/(V*e11*e22),
        "C1212":2*W["C1212"]/(V*g12*g12),
        "C1313":2*W["Q44"]/(V*g13*g13),
        "C2323":2*W["Q55"]/(V*g23*g23),
    }


def reconstruct_dbc_terms(energies:dict[str,float], p:float,q:float,h:float, displacement:float=1.0) -> dict[str,dict]:
    """{'Q': shell terms, 'C': 3D terms}; a group is present only when all of its cases are."""
    out={}
    if all(k in energies for k in THESIS_DBC_CASES):
        out["Q"]=reconstruct_thesis_q(energies,p,q,h,displacement)
    if all(k in energies for k in DBC_C_CASES+("Q44","Q55")):
        out["C"]=reconstruct_dbc_c(energies,p,q,h,displacement)
    return out


def reconstruct_thesis_q(energies:dict[str,float], p:float,q:float,h:float, displacement:float=1.0) -> dict[str,float]:
    V=float(p)*float(q)*float(h); d=float(displacement)
    e11=d/float(p); e22=d/float(q); g12=d/float(p); g13=d/float(p); g23=d/float(q)
    W={k:float(energies[k]) for k in THESIS_DBC_CASES}
    return {
        "Q11":2*W["Q11"]/(V*e11*e11),
        "Q22":2*W["Q22"]/(V*e22*e22),
        "Q12":(W["Q12"]-W["Q11"]-W["Q22"])/(V*e11*e22),
        "Q66":2*W["Q66"]/(V*g12*g12),
        "Q44":2*W["Q44"]/(V*g13*g13),
        "Q55":2*W["Q55"]/(V*g23*g23),
    }


def write_dbc_rve_svg(path:str|Path, dbc_model:RVEModel, diag:dict) -> None:
    """Plan view of the DBC RVE: interior walls, half-section face walls and face labels."""
    p=float(diag['p']); q=float(diag['q']); pad=0.12*max(p,q); s=480.0/max(p+2*pad,q+2*pad)
    W=(p+2*pad)*s; H=(q+2*pad)*s
    X=lambda x:(x+pad)*s; Y=lambda y:H-(y+pad)*s
    face=set(int(i) for i in diag.get('face_walls',[]))
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.0f}" height="{H:.0f}" viewBox="0 0 {W:.1f} {H:.1f}" font-family="sans-serif" font-size="12">',
         '<rect width="100%" height="100%" fill="white"/>',
         f'<rect x="{X(0):.2f}" y="{Y(q):.2f}" width="{p*s:.2f}" height="{q*s:.2f}" fill="#f3f6fb" stroke="#555" stroke-dasharray="5,4"/>']
    poly=diag.get('central_cell_polygon_rve')
    if poly:
        pts=' '.join(f'{X(x):.2f},{Y(y):.2f}' for x,y in poly)
        out.append(f'<polygon points="{pts}" fill="#dce8f7" stroke="none"/>')
    for i,e in enumerate(dbc_model.elements,1):
        a=dbc_model.nodes[e.node_i]; b=dbc_model.nodes[e.node_j]
        col,wd=('#d9480f',3.0) if i in face else ('#1c4e80',4.5)
        out.append(f'<line x1="{X(a[0]):.2f}" y1="{Y(a[1]):.2f}" x2="{X(b[0]):.2f}" y2="{Y(b[1]):.2f}" stroke="{col}" stroke-width="{wd}" stroke-linecap="round"/>')
    out.append(f'<text x="{X(p/2):.1f}" y="{Y(-pad*0.45):.1f}" text-anchor="middle">x1 = 0 ... p = {p:.2f} mm</text>')
    out.append(f'<text x="{X(-pad*0.55):.1f}" y="{Y(q/2):.1f}" text-anchor="middle" transform="rotate(-90 {X(-pad*0.55):.1f} {Y(q/2):.1f})">x2 = 0 ... q = {q:.2f} mm</text>')
    out.append(f'<text x="{X(p/2):.1f}" y="{Y(q+pad*0.55):.1f}" text-anchor="middle">{diag.get("primitive_cells","?")} primitive cells; orange = wall in a face ({diag.get("face_wall_section","")})</text>')
    out.append('</svg>')
    Path(path).write_text('\n'.join(out),encoding='utf-8')


def prepare_thesis_dbc(run_root:str|Path, model:RVEModel, materials:dict[str,Material], laminates:list[Laminate], settings:ThesisDBCSettings, *, run_name:str|None=None) -> ThesisDBCFiles:
    settings.validate(); dbc_model,diag=build_thesis_dbc_rectangle(model,cells=settings.rve_cells,face_wall_section=settings.face_wall_section)
    root=Path(run_root).expanduser().resolve();root.mkdir(parents=True,exist_ok=True)
    stamp=datetime.now().strftime('%Y%m%d_%H%M%S'); name=_safe_name(run_name or f"{model.name}_Thesis_DBC_{stamp}","DBC")
    run_dir=root/name;run_dir.mkdir(parents=False,exist_ok=False)
    cases=[]; metadata=[]; extractor=_extractor_source()
    for i,case in enumerate(DBC_CASES,1):
        cdir=run_dir/f"{i:02d}_{case}";cdir.mkdir()
        job=_safe_name(f"DBC_{case}","DBC")
        inp=cdir/f"{job}.inp"; ex=cdir/'extract_dbc_odb.py'; res=cdir/'case_results.json'; log=cdir/'abaqus_case.log'
        text,meta=generate_thesis_dbc_input(dbc_model,materials,laminates,settings,case,job_name=job,resolved_domain=True,domain_diagnostics=diag)
        inp.write_text(text,encoding='utf-8');ex.write_text(extractor,encoding='utf-8')
        cases.append(ThesisDBCCaseFiles(i,case,cdir,job,inp,ex,res,log));metadata.append({"index":i,"name":case,"deck_meta":meta,"directory":cdir.name})
    snapshot={
        "mode":"thesis_table_3_4_prescribed_displacement_validation",
        "source":"Chapter 3, Table 3.4",
        "model":{"name":model.name,"nodes":model.nodes.tolist(),"elements":[asdict(e) for e in model.elements],"a1":model.a1.tolist(),"a2":model.a2.tolist(),"height":model.height},
        "dbc_model":{"nodes":dbc_model.nodes.tolist(),"elements":[asdict(e) for e in dbc_model.elements],"a1":dbc_model.a1.tolist(),"a2":dbc_model.a2.tolist(),"height":dbc_model.height},
        "domain":diag,"settings":asdict(settings),"cases":metadata,
    }
    snapshot_file=run_dir/'thesis_dbc_snapshot.json';snapshot_file.write_text(json.dumps(snapshot,indent=2),encoding='utf-8')
    rbm_summary={c['name']:c['deck_meta'].get('rigid_body_audit',{}) for c in metadata}
    (run_dir/'dbc_rigid_body_audit.json').write_text(json.dumps(rbm_summary,indent=2),encoding='utf-8')
    # Simple preview as CSV, intentionally easy to inspect/import.
    with (run_dir/'thesis_dbc_geometry_nodes.csv').open('w',encoding='utf-8') as f:
        f.write('node,x1_mm,x2_mm\n')
        for i,pnt in enumerate(dbc_model.nodes,1):f.write(f"{i},{pnt[0]:.16g},{pnt[1]:.16g}\n")
    faces=diag.get('face_wall_faces',{}); fac=diag.get('wall_section_factors',[1.0]*len(dbc_model.elements))
    with (run_dir/'thesis_dbc_geometry_walls.csv').open('w',encoding='utf-8') as f:
        f.write('wall,node_i,node_j,laminate,length_mm,face,section_factor\n')
        for i,e in enumerate(dbc_model.elements,1):
            ln=float(np.linalg.norm(dbc_model.nodes[e.node_j]-dbc_model.nodes[e.node_i]))
            f.write(f"{i},{e.node_i+1},{e.node_j+1},{e.laminate},{ln:.12g},{faces.get(str(i),'')},{fac[i-1]:g}\n")
    try:
        write_dbc_rve_svg(run_dir/'thesis_dbc_rve.svg',dbc_model,diag)
    except Exception:
        pass
    return ThesisDBCFiles(run_dir,snapshot_file,run_dir/'thesis_dbc_results.json',run_dir/'thesis_dbc_run.log',cases)


def run_thesis_dbc(files:ThesisDBCFiles, abaqus_command:str, settings:ThesisDBCSettings, progress:Callable[[int,int,str],None]|None=None) -> dict:
    settings.validate();exe=str(abaqus_command).strip()
    if not exe:raise ValueError("Abaqus executable/command is empty.")
    snap=json.loads(files.snapshot_file.read_text(encoding='utf-8'));diag=snap['domain'];energies={};completed=[]
    with files.log_file.open('w',encoding='utf-8',errors='replace') as master:
        master.write('RVE Studio thesis Table 3.4 DBC validation\n')
        for i,case in enumerate(files.cases,1):
            if progress:progress(i,len(files.cases),case.name)
            master.write(f"[{i}/{len(files.cases)}] {case.name}\n");master.flush()
            with case.log_file.open('w',encoding='utf-8',errors='replace') as log:
                rc=_run_shell_command([exe,f"job={case.job_name}",f"input={case.input_file.name}",f"cpus={int(settings.cpus)}","interactive"],case.case_dir,log)
                if rc!=0:raise RuntimeError(f"Abaqus thesis DBC case {case.name} returned exit code {rc}. See {case.log_file}")
                odb=case.case_dir/f"{case.job_name}.odb"
                if not odb.exists():raise RuntimeError(f"Abaqus thesis DBC case {case.name} produced no ODB. See {case.log_file}")
                rc2=_run_shell_command([exe,'python',case.extractor_file.name,odb.name,case.results_file.name],case.case_dir,log)
                if not case.results_file.exists():raise RuntimeError(f"DBC extractor produced no result for {case.name} (wrapper {rc2}).")
            raw=json.loads(case.results_file.read_text(encoding='utf-8'))
            if str(raw.get('status','')).upper()!='COMPLETED':raise RuntimeError(raw.get('traceback') or raw.get('error') or f"DBC extraction failed for {case.name}")
            energy=float(raw['strain_energy_ALLSE']);energies[case.name]=energy;completed.append({'index':i,'name':case.name,'strain_energy':energy,'odb':raw.get('odb','')})
            if not settings.save_odb:
                try:odb.unlink()
                except Exception:pass
    terms=reconstruct_dbc_terms(energies,float(diag['p']),float(diag['q']),float(diag['height']),float(settings.displacement))
    Q=terms['Q']
    result={
        'status':'COMPLETED','mode':'Thesis Table 3.4 prescribed-displacement FE validation','Q_DBC':Q,'C_DBC':terms.get('C'),'case_energies':completed,
        'p':float(diag['p']),'q':float(diag['q']),'height':float(diag['height']),'rve_area':float(diag['rve_area']),'rve_volume':float(diag['rve_volume']),
        'domain':diag,'boundary_condition_note':'Shell terms: prescribed translations of thesis Table 3.4, top/bottom U3 free. 3D terms: the same side faces with U3 = 0 on the top and bottom faces. Unspecified DOFs and all shell rotations remain free; only the minimum zero translational gauge needed to reach rigid-body rank 6/6 is added.',
        'run_dir':str(files.run_dir),
    }
    files.summary_file.write_text(json.dumps(result,indent=2),encoding='utf-8')
    with (files.run_dir/'Q_thesis_DBC_MPa.csv').open('w',encoding='utf-8') as f:
        f.write('term,value_MPa\n')
        for k in THESIS_DBC_CASES:
            f.write(f"{k},{Q[k]:.16g}\n")
    if terms.get('C'):
        with (files.run_dir/'C_DBC_MPa.csv').open('w',encoding='utf-8') as f:
            f.write('term,value_MPa\n')
            for k in SOLID_TERMS:f.write(f"{k},{terms['C'][k]:.16g}\n")
    with (files.run_dir/'case_energies.csv').open('w',encoding='utf-8') as f:
        f.write('case,strain_energy_Nmm\n')
        for k in DBC_CASES:f.write(f"{k},{energies[k]:.16g}\n")
    return result
