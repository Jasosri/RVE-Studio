"""Cylindrical-only inverse design across built-in geometries, ply counts, stacking sequences and materials.
Reference script; machine-generated candidates are finite search only, not proven global optimum.
"""
import os,sys,time,json,math
from pathlib import Path
import numpy as np,pandas as pd
from scipy.linalg import eigh
P=Path(__file__).resolve().parent
for candidate in (P/'RVE_Studio_Node_Explorer',P.parent.parent,P.parent.parent.parent):
 if (candidate/'rvestudio').is_dir(): sys.path.insert(0,str(candidate)); break
from rvestudio.models import Material,Ply,Laminate,Element,RVEModel
from rvestudio.geometry import BUILTIN_GEOMETRIES
from rvestudio.library import load_material_library,material_mapping
from rvestudio.homogenization import homogenize_rve
from rvestudio.load_recovery import recover_wall_loads
from rvestudio.buckling import evaluate_wall_buckling
from rvestudio.io import save_project
from rvestudio.api import solve_example

OUT=P/'results'; OUT.mkdir(exist_ok=True)
TARGET=np.loadtxt(P/'target/results/C_target_6x6_MPa.csv',delimiter=',')
assert TARGET.shape==(6,6) and np.linalg.eigvalsh(TARGET)[0]>0
h=150.; RHO_MAX=.15; TMIN=.35; PLY_MIN=.08; MAX_RATIO=1.9; SAFETY=1.5; STIFFNESS_FACTOR=1.05; COST_CAP=50. # illustrative units per m2
mats=material_mapping(load_material_library(P/'source_libraries/materials.json' if (P/'source_libraries/materials.json').exists() else None))
MATERIALS={'Glass':['E-glass/Epoxy'], 'Boron':['Boron/Epoxy'], 'Hybrid_GBG':['E-glass/Epoxy','Boron/Epoxy','E-glass/Epoxy'], 'Hybrid_BGB':['Boron/Epoxy','E-glass/Epoxy','Boron/Epoxy']}
COLORS=['HEX','SQUARE','MIXED','TRI']
GEO_BOUNDS={'HEX':((35,70),(35,70)),'SQUARE':((35,70),(35,70)),'MIXED':((30,60),(1.5,2.5)),'TRI':((35,70),(.8,1.2))}
PARAM_NAMES={'HEX':('L1','L2'),'SQUARE':('Lx','Ly'),'MIXED':('Lh','R'),'TRI':('Lh','R')}
LOADS={'uniaxial':[-1.,0.,0.,0.,0.,0.], 'biaxial':[-1.,-.5,0.,0.,0.,0.], 'combined':[-1.,-.5,0.,0.,0.,.25]}

# Symmetric stacking families based on 0, +45, -45 and 90 deg; distinct integer ply counts.
def stacks(n):
    half=(n+1)//2
    base=[(0,), (90,), (45,), (0,90),(90,0), (45,-45), (0,45),(45,0),(0,90,45),(90,45,0),(-45,0)]
    patterns=[]
    for it in range(len(base)):
        if it<3: a=[base[it][0]]*half
        else:
            src=list(base[it]);a=[src[j%len(src)] for j in range(half)]
        pattern=tuple(a+a[-2::-1]);
        if pattern not in patterns:patterns.append(pattern)
    return patterns

rng=np.random.default_rng(20261030)
# low-dimensional geometry sweep, fixed seed, includes central and 4 LHS-ish randomized points
GEOPARAMS={}
for gid in COLORS:
    b=GEO_BOUNDS[gid]; points=[[(x+y)/2 for x,y in b]]
    points.extend([[float(rng.uniform(*b[0])),float(rng.uniform(*b[1]))] for _ in range(5)])
    GEOPARAMS[gid]=points

def material_names(kind,n):
    scheme=MATERIALS[kind]
    if len(scheme)==1:return [scheme[0]]*n
    if kind=='Hybrid_GBG':return [scheme[0] if i%2==0 else scheme[1] for i in range(n)]
    return [scheme[0] if i%2==0 else scheme[1] for i in range(n)]

def geometry_meta(gid,vals):
    return BUILTIN_GEOMETRIES[gid](dict(zip(PARAM_NAMES[gid],vals)))

def construct(meta,rho):
    gid=meta['Geometry']; p=[meta['Param1'],meta['Param2']];geo=geometry_meta(gid,p)
    lengths=np.linalg.norm(geo.nodes[geo.connectivity[:,1]]-geo.nodes[geo.connectivity[:,0]],axis=1)
    weights=geo.wall_weights;area=abs(np.linalg.det(np.column_stack((geo.a1,geo.a2))))
    factors=np.array(meta['WallFactors'],float);ts=rho*area*factors/np.sum(lengths*weights*factors)
    n=int(meta['PlyCount']);angles=meta['Angles']; plym=material_names(meta['MaterialScheme'],n)
    laminates=[Laminate([Ply(plym[i],float(angles[i]),float(t/n)) for i in range(n)],'Wall %d'%(j+1)) for j,t in enumerate(ts)]
    elems=[Element(int(a),int(b),k,float(weights[k])) for k,(a,b) in enumerate(geo.connectivity)]
    model=RVEModel(geo.nodes,elems,geo.a1,geo.a2,h,wall_groups=list(geo.wall_groups),wall_element_type='paper',strip_bending_mode='cylindrical',name=meta['ID'])
    return model,laminates,ts

def calc(meta,rho,with_buckling=False,mmax=4):
    model,lams,ts=construct(meta,rho)
    if min(ts)<TMIN-1e-8 or min(ts)/len(lams[0].plies)<PLY_MIN-1e-8 or max(ts)/min(ts)>MAX_RATIO+1e-8:
        return None
    out=homogenize_rve(model,mats,lams)
    adequate=float(eigh(out.CH,TARGET,eigvals_only=True,check_finite=False)[0])
    cost=float(out.physical.cost_per_rve_m3*.15) if out.physical.has_complete_cost else math.inf
    rec=dict(Adequacy=adequate, Density=float(out.physical.average_density_kg_m3),ArealMass=float(out.physical.average_density_kg_m3*.15),
             RelDensity=float(out.physical.relative_density),TotalCostPerM2=cost,
             MatrixGap=float(np.linalg.norm(out.CH-TARGET,'fro')/np.linalg.norm(TARGET,'fro')),MinThickness=float(ts.min()),MinPlyThickness=float(ts.min()/meta['PlyCount']),ThicknessRatio=float(ts.max()/ts.min()),
             Q11=float(out.props.Qplane[0,0]),Q22=float(out.props.Qplane[1,1]),Q66=float(out.props.Qplane[2,2]),Q44=float(out.props.Qshear[0,0]),Q55=float(out.props.Qshear[1,1]))
    if with_buckling:
        for name,stress in LOADS.items():
            wrec=recover_wall_loads(model,out,macro_stress=stress)
            buck=evaluate_wall_buckling(model,out,wrec,m_max=mmax,n_max=mmax,safety_factor=SAFETY,convergence_check=False)
            rec['Buck_'+name]=float(buck.minimum_lambda)
            rec['Wall_'+name]=int(buck.critical_wall) if buck.critical_wall is not None else -1
            row=next((r for r in buck.walls if r['wall']==buck.critical_wall),{})
            rec['Mode_'+name]=f"({row.get('dominant_m','?')},{row.get('dominant_n','?')})"
    return rec,out,model,lams

# create candidate list, include individual wall thickness factor patterns
allmeta=[]; count=0
for gid in COLORS:
 for k,p in enumerate(GEOPARAMS[gid]):
  for mat in MATERIALS:
   for n in [3,5,7,9]:
    for angles in stacks(n):
     count+=1
     geo=geometry_meta(gid,p)
     nw=len(geo.connectivity)
     # standard uniform allocation; wall factor search follows top solutions
     meta={'ID':f'P{count:05d}','Geometry':gid,'Param1':p[0],'Param2':p[1], 'ParamNames':PARAM_NAMES[gid],
           'MaterialScheme':mat,'PlyCount':n,'Angles':list(angles),'WallFactors':[1.]*nw,'WallCount':nw,'NodeCount':len(geo.nodes)}
     allmeta.append(meta)

(OUT/'search_spec.json').write_text(json.dumps({'seed':20261030,'target_matrix':TARGET.tolist(),'geometry_ranges':GEO_BOUNDS,
     'materials':MATERIALS,'ply_counts':[3,5,7,9],'angle_options':[0,45,-45,90],
     'minimum_wall_mm':TMIN,'minimum_ply_mm':PLY_MIN,'max_thickness_ratio':MAX_RATIO,'max_relative_density':RHO_MAX,
     'cost_cap_arbitrary_per_m2':COST_CAP,'minimum_stiffness_factor':STIFFNESS_FACTOR,'buckling_safety_factor':SAFETY,'loads_MPa':LOADS,'solver':'paper cylindrical ONLY'},indent=2))
(OUT/'search_profiles.json').write_text(json.dumps(allmeta),encoding='utf8')

if __name__=='__main__':
 phase=sys.argv[1] if len(sys.argv)>1 else 'screen'
 if phase=='screen':
  t0=time.time();rows=[]; failures=[]
  for idx,meta in enumerate(allmeta):
   try:
    val=calc(meta,RHO_MAX)
    if val is None:rows.append({'ID':meta['ID'],'Status':'THICKNESS_INFEASIBLE'});continue
    r,_,_,_=val
    rows.append({**{k:meta[k] for k in ['ID','Geometry','MaterialScheme','PlyCount','NodeCount','WallCount','Param1','Param2']},'Angles':str(meta['Angles']),'Status':'OK',**r})
   except Exception as e:
    rows.append({'ID':meta['ID'],'Status':'FAILED','Error':str(e)})
   if (idx+1)%500==0:print('SCREEN',idx+1,'/',len(allmeta),'t',round(time.time()-t0,1),flush=True)
  df=pd.DataFrame(rows);df.to_csv(OUT/'initial_screen.csv',index=False)
  print('DONE',len(df),df.Status.value_counts().to_dict(),'max rho adequate',sum(df.Adequacy>=1.05),'time',round(time.time()-t0,1),flush=True)
 elif phase=='refine':
  screen=pd.read_csv(OUT/'initial_screen.csv');ok=screen[(screen.Status=='OK') & (screen.Adequacy>=STIFFNESS_FACTOR)].copy()
  ok['ApproxDensity']=ok['Density']*np.minimum(1.0,STIFFNESS_FACTOR/ok['Adequacy'])
  ids=set(ok.nsmallest(700,'ApproxDensity').ID)
  for gid in COLORS:
   for mm in MATERIALS:
    ids.update(ok[(ok.Geometry==gid)&(ok.MaterialScheme==mm)].nsmallest(35,'ApproxDensity').ID)
  for nn in (3,5,7,9):ids.update(ok[ok.PlyCount==nn].nsmallest(60,'ApproxDensity').ID)
  ok=ok[ok.ID.isin(ids)].copy()
  print('REFINE CANDIDATES',len(ok),flush=True)
  meta_by_id={m['ID']:m for m in allmeta}
  # First bound density needed for stiffness, with feasibility and thickness constraints
  t0=time.time();results=[]
  sl=int(os.environ.get('RVE_START','0'));sr=int(os.environ.get('RVE_STOP',str(len(ok))))
  for idx,z in enumerate(ok.to_dict('records')[sl:sr]):
   meta=meta_by_id[z['ID']];rmax=calc(meta,RHO_MAX)
   if rmax is None:continue
   rho0=RHO_MAX; tsmin=rmax[0]['MinThickness'];n=meta['PlyCount']
   lo=max(.003, rho0*TMIN/tsmin, rho0*(n*PLY_MIN)/tsmin)
   if lo>rho0:continue
   def adequate(rho):
    rr=calc(meta,rho)
    if rr is None:return False,None
    return rr[0]['Adequacy']>=STIFFNESS_FACTOR,rr[0]
   success,lowrec=adequate(lo)
   if success:rho=lo
   else:
    hi=rho0
    for _ in range(7):
     mid=(lo+hi)/2;passed,_=adequate(mid)
     if passed:hi=mid
     else:lo=mid
    rho=min(rho0,hi*1.0002)
   v=calc(meta,rho)
   if v is None:continue
   rr=v[0];results.append({**{k:meta[k] for k in ['ID','Geometry','MaterialScheme','PlyCount','Param1','Param2','NodeCount','WallCount']},'Angles':str(meta['Angles']),'Status':'OK' if rr['Adequacy']>=STIFFNESS_FACTOR else 'INFEASIBLE',**rr})
   if (idx+1)%200==0:print('REFINE',idx+1,'/',len(ok),'t',round(time.time()-t0,1),flush=True)
  df=pd.DataFrame(results).sort_values('Density');df.to_csv(OUT/f'stiffness_refined_{sl}_{sr}.csv',index=False)
  print('REFINE DONE',len(df),'best',df.head(10)[['ID','Geometry','MaterialScheme','PlyCount','Density','Adequacy','TotalCostPerM2']].to_string(index=False),flush=True)
 elif phase=='buckling':
  refined=pd.read_csv(OUT/'stiffness_refined.csv');refined=refined[(refined.Status=='OK') & (refined.TotalCostPerM2<=COST_CAP)]
  meta_by_id={m['ID']:m for m in allmeta}
  # Selection: strongest candidate mix low mass, per topology and material and ply count
  ids=list(refined.head(85).ID)
  for gid in COLORS:
   for mat in MATERIALS:
    ids += list(refined[(refined.Geometry==gid)&(refined.MaterialScheme==mat)].head(5).ID)
  for n in [3,5,7,9]:ids+=list(refined[refined.PlyCount==n].head(6).ID)
  ids=list(dict.fromkeys(ids))
  t0=time.time();records=[]
  sl=int(os.environ.get('RVE_START','0'));sr=int(os.environ.get('RVE_STOP',str(len(ids))))
  print('BUCKLING SHORTLIST',len(ids),'batch',sl,sr,flush=True)
  for ii,cid in enumerate(ids[sl:sr]):
   meta=meta_by_id[cid]
   rho_stiff=float(refined[refined.ID==cid].iloc[0].RelDensity)
   rmaxall=calc(meta,RHO_MAX,with_buckling=False)
   if rmaxall is None:continue
   caprho=min(RHO_MAX,RHO_MAX*COST_CAP/max(rmaxall[0]['TotalCostPerM2'],1e-15)*0.999999)
   if caprho<rho_stiff:
    records.append(dict(ID=cid,Geometry=meta['Geometry'],MaterialScheme=meta['MaterialScheme'],PlyCount=meta['PlyCount'],Status='COST_INFEASIBLE'));continue
   val=calc(meta,caprho,with_buckling=True,mmax=4)
   if val is None:continue
   rmax=val[0]
   def grade(rec):
    return min(rec['Adequacy']/STIFFNESS_FACTOR, COST_CAP/rec['TotalCostPerM2'],min(rec[f'Buck_{n}']/SAFETY for n in LOADS))
   if grade(rmax)<1-1e-7:
    records.append(dict(ID=cid,Geometry=meta['Geometry'],MaterialScheme=meta['MaterialScheme'],PlyCount=meta['PlyCount'],Status='NO_FEASIBLE_AT_ALLOWED_MAX',UpperGrade=grade(rmax)));continue
   lowval=calc(meta,rho_stiff,with_buckling=True,mmax=4)
   if lowval is not None and grade(lowval[0])>=1:
    rho=rho_stiff
   else:
    lo=rho_stiff;hi=caprho
    for _ in range(9):
     mid=(lo+hi)/2
     value=calc(meta,mid,with_buckling=True,mmax=4)
     if value is not None and grade(value[0])>=1:hi=mid
     else:lo=mid
    rho=min(caprho,hi*1.00005)
   value=calc(meta,rho,with_buckling=True,mmax=5)
   if value is None:continue
   rr=value[0];status='PASS' if grade(rr)>=1 else 'RECHECK_FAIL'
   records.append({**{k:meta[k] for k in ['ID','Geometry','MaterialScheme','PlyCount','Param1','Param2','NodeCount','WallCount']},'Angles':str(meta['Angles']), 'Status':status,**rr})
   if (ii+1)%20==0:print('BUCKLING',ii+1,'/',len(ids),'t',round(time.time()-t0,1),flush=True)
  df=pd.DataFrame(records);df.sort_values('Density',na_position='last').to_csv(OUT/f'buckling_refined_{sl}_{sr}.csv',index=False)
  winners=df[df.Status=='PASS'];print('BUCKLING COMPLETE',len(df),'PASS',len(winners),'best\n',winners.head(12)[['ID','Geometry','MaterialScheme','PlyCount','Density','Adequacy','Buck_uniaxial','Buck_biaxial','Buck_combined']].to_string(index=False),flush=True)
 elif phase=='verify':
  bk=pd.read_csv(OUT/'buckling_refined.csv'); winners=bk[bk.Status=='PASS'].nsmallest(8,'Density'); profiles={m['ID']:m for m in allmeta}
  for _,r in winners.iterrows():
   meta=profiles[r.ID];value=calc(meta,float(r.RelDensity),True,6)
   rec,out,model,lams=value
   filename=OUT/f'project_{r.ID}.rvejson';save_project(filename,model,mats,lams)
   np.savetxt(OUT/f'C_achieved_{r.ID}_MPa.csv',out.CH,delimiter=',',fmt='%.9f')
   print(r.ID,r.Geometry,r.MaterialScheme,r.PlyCount,'rho_eff',rec['Density'],'stiff',rec['Adequacy'],'bucklings',[round(rec[f'Buck_{n}'],3) for n in LOADS], 'cost',rec['TotalCostPerM2'],flush=True)
  np.savetxt(OUT/'C_target_MPa.csv',TARGET,delimiter=',',fmt='%.9f')
