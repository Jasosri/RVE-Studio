"""Use the Node-count geometry selector with the standard RVE Studio solver."""
import numpy as np
from pathlib import Path
from rvestudio.plugins import discover_geometry_plugins
from rvestudio.models import Material,Ply,Laminate,Element,RVEModel
from rvestudio.library import PACKAGE_ROOT,load_material_library,material_mapping
from rvestudio.homogenization import homogenize_rve
from scipy.linalg import eigh

# The NODESELECT library row is opt-in. This direct API example works regardless.
selector=discover_geometry_plugins(PACKAGE_ROOT/'plugins')['NODESELECT']
target=np.loadtxt(Path(__file__).resolve().parents[1]/'data'/'C_target_6x6_MPa.csv',delimiter=',')
mats=material_mapping(load_material_library())
for n in range(3,11):
    geo=selector({'NodeCount':n,'StretchX':1,'StretchY':1,'Shear':0})
    lengths=np.array([np.linalg.norm(geo.nodes[j]-geo.nodes[i]) for i,j in geo.connectivity])
    rho=.08; thickness=rho*abs(np.linalg.det(np.column_stack([geo.a1,geo.a2])))/np.dot(lengths,geo.wall_weights)
    laminates=[Laminate([Ply('Boron/Epoxy',a,thickness/5) for a in [0,45,90,45,0]]) for _ in geo.connectivity]
    model=RVEModel(geo.nodes,[Element(int(i),int(j),k,float(geo.wall_weights[k])) for k,(i,j) in enumerate(geo.connectivity)],geo.a1,geo.a2,150)
    scores=[]
    for mode in ['paper','freeMoment']:
        model.wall_element_type='paper' if mode=='paper' else 'anisotropicStrip'
        model.strip_bending_mode='cylindrical' if mode=='paper' else 'freeMoment'
        result=homogenize_rve(model,mats,laminates)
        scores.append(float(eigh(result.CH,target,eigvals_only=True)[0]))
    print(n,'nodes',len(geo.connectivity),'walls','adequacy',min(scores))
