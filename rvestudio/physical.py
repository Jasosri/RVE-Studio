from __future__ import annotations
import math
import numpy as np
from typing import Mapping
from .models import Material,Laminate,RVEModel,PhysicalProperties


def _volume_to_m3(unit:str)->float:
    u=unit.lower().strip()
    if u=='mm': return 1e-9
    if u=='cm': return 1e-6
    if u=='m': return 1.0
    raise ValueError(f"Unsupported length unit {unit!r}; use mm, cm or m")


def compute_physical_properties(model:RVEModel,materials:Mapping[str,Material],laminates:list[Laminate])->PhysicalProperties:
    area=abs(np.linalg.det(np.column_stack((model.a1,model.a2)))); VRVE=area*model.height; conv=_volume_to_m3(model.length_unit)
    mat_vol={n:0.0 for n in materials}; mat_mass={n:0.0 for n in materials}; mat_cost={n:0.0 for n in materials}; used=set(); wall_rows=[]
    complete=True; complete_cost=True
    for eidx,e in enumerate(model.elements,1):
        L=float(np.linalg.norm(model.nodes[e.node_j]-model.nodes[e.node_i])); lam=laminates[e.laminate]; Vwall=0.; Mwall=0.; Cwall=0.; known=True; cost_known=True
        for ply in lam.plies:
            V=e.weight*L*model.height*ply.thickness; Vwall+=V; used.add(ply.material); mat_vol[ply.material]+=V
            mat=materials[ply.material]; rho=mat.density_kg_m3
            if math.isfinite(rho) and rho>0:
                m=V*conv*rho; Mwall+=m; mat_mass[ply.material]+=m
                cpk=float(mat.cost_per_kg)
                if math.isfinite(cpk) and cpk>=0:
                    c=m*cpk; Cwall+=c; mat_cost[ply.material]+=c
                else:
                    cost_known=False; complete_cost=False
            else:
                known=False; complete=False; cost_known=False; complete_cost=False
        wall_rows.append({"wall_id":eidx,"length":L,"thickness":lam.thickness,"solid_volume":Vwall,
                          "average_density_kg_m3":Mwall/(Vwall*conv) if known and Vwall>0 else math.nan,
                          "mass_kg":Mwall if known else math.nan,"weight_n":Mwall*model.gravity if known else math.nan,
                          "material_cost":Cwall if cost_known else math.nan})
    solid=sum(r["solid_volume"] for r in wall_rows); void=VRVE-solid; rho_rel=solid/VRVE
    missing=[n for n in used if not (math.isfinite(materials[n].density_kg_m3) and materials[n].density_kg_m3>0)]
    missing_cost=[n for n in used if not (math.isfinite(float(materials[n].cost_per_kg)) and float(materials[n].cost_per_kg)>=0)]
    if complete:
        total=sum(mat_mass[n] for n in used); weight=total*model.gravity; avg=total/(VRVE*conv); eff=total/(solid*conv); relm=avg/eff
    else:
        total=weight=avg=eff=relm=math.nan
    if complete and complete_cost:
        total_cost=sum(mat_cost[n] for n in used); effective_cost=total_cost/total if total>0 else math.nan; cost_per_rve_m3=total_cost/(VRVE*conv) if VRVE>0 else math.nan
    else:
        total_cost=effective_cost=cost_per_rve_m3=math.nan
    breakdown=[]
    for n in sorted(used):
        rho=materials[n].density_kg_m3; mass=mat_mass[n] if math.isfinite(rho) and rho>0 else math.nan; cpk=float(materials[n].cost_per_kg); cost=mat_cost[n] if math.isfinite(mass) and math.isfinite(cpk) and cpk>=0 else math.nan
        breakdown.append({"material":n,"density_kg_m3":rho,"cost_per_kg":cpk,"solid_volume":mat_vol[n],"mass_kg":mass,"material_cost":cost,
                          "volume_fraction_of_solid":mat_vol[n]/solid if solid>0 else 0.0,
                          "mass_fraction":mass/total if complete and total>0 else math.nan})
    return PhysicalProperties(model.length_unit,area,VRVE,solid,void,rho_rel,total,weight,avg,eff,relm,complete,missing,breakdown,wall_rows,
                              total_cost,effective_cost,cost_per_rve_m3,complete and complete_cost,missing_cost)
