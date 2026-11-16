"""Verify saved numerical results using the real RVE Studio cylindrical-only mechanics API."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.linalg import eigh
P=Path(__file__).resolve().parent
for candidate in (P/'RVE_Studio_Node_Explorer',P.parent.parent,P.parent.parent.parent):
 if (candidate/'rvestudio').is_dir(): sys.path.insert(0,str(candidate)); break
from rvestudio.io import load_project
from rvestudio.homogenization import homogenize_rve
from rvestudio.load_recovery import recover_wall_loads
from rvestudio.buckling import evaluate_wall_buckling
R=P/'results'
T=np.loadtxt(R/'C_target_MPa.csv',delimiter=',')
assert T.shape==(6,6) and np.linalg.eigvalsh(T)[0]>0
assert np.allclose(T,T.T,atol=1e-8)
screen=pd.read_csv(R/'initial_screen.csv');refined=pd.read_csv(R/'stiffness_refined.csv');bk=pd.read_csv(R/'buckling_refined.csv')
assert len(screen)==4128 and screen.Status.eq('OK').all()
assert len(refined)==1177 and len(bk)==136 and (bk.Status=='PASS').sum()==110
passed=bk[bk.Status=='PASS'].sort_values('Density')
assert passed.iloc[0].ID=='P00050'
loads={'uniaxial':[-1,0,0,0,0,0],'biaxial':[-1,-.5,0,0,0,0],'combined':[-1,-.5,0,0,0,.25]}
checks=[]
for src in sorted(R.glob('project_P*.rvejson')):
 model,mats,lams=load_project(src)
 assert model.wall_element_type.lower()=='paper' and model.strip_bending_mode.lower()=='cylindrical'
 rs=homogenize_rve(model,mats,lams)
 cid=src.stem.removeprefix('project_');f=R/f'C_achieved_{cid}_MPa.csv'
 assert f.exists(),f
 mat=np.loadtxt(f,delimiter=',')
 assert np.allclose(rs.CH,mat,rtol=1e-7,atol=1e-6)
 eig=float(eigh(rs.CH,T,eigvals_only=True)[0]);assert eig >= 1.05-1e-6
 assert min(l.thickness for l in lams)>=.35-1e-7
 assert min(min(p.thickness for p in l.plies) for l in lams)>=.08-1e-7
 cost=float(rs.physical.cost_per_rve_m3*.15);assert cost<=50+1e-6
 bd={}
 for name,values in loads.items():
  rec=recover_wall_loads(model,rs,macro_stress=values)
  buck=evaluate_wall_buckling(model,rs,rec,m_max=6,n_max=6,safety_factor=1.5,convergence_check=False)
  bd[name]=buck.minimum_lambda
  assert buck.minimum_lambda>=1.5-1e-5,(cid,name,buck.minimum_lambda)
 checks.append(dict(ID=cid,RhoKgM3=rs.physical.average_density_kg_m3,ArealMassKgM2=rs.physical.average_density_kg_m3*.15,StiffnessFactor=eig,CostIndex=cost,**{'Buck_'+k:v for k,v in bd.items()}))
pd.DataFrame(checks).to_csv(R/'verified_saved_projects.csv',index=False)
print('VERIFIED',len(checks),'projects; 4128 screens, 1177 refined, 136 buckling screens, 110 pass.')
for c in checks:print(c['ID'],round(c['RhoKgM3'],3),'lambda=',round(c['StiffnessFactor'],6),'buck_min=',round(min(c['Buck_'+k] for k in loads),4))
