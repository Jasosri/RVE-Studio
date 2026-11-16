from __future__ import annotations
import json
from dataclasses import asdict
from pathlib import Path
import numpy as np
from .models import Material,Ply,Laminate,Element,RVEModel


def save_project(path,model:RVEModel,materials:dict[str,Material],laminates:list[Laminate]):
    data={'model':{'name':model.name,'nodes':model.nodes.tolist(),'a1':model.a1.tolist(),'a2':model.a2.tolist(),'height':model.height,'periodic_tol':model.periodic_tol,
                   'anchor_node':model.anchor_node,'fix_anchor_rotation':model.fix_anchor_rotation,'periodic_pairs':model.periodic_pairs.tolist() if model.periodic_pairs is not None else None,'wall_element_type':model.wall_element_type,'strip_bending_mode':model.strip_bending_mode,'wall_shear':model.wall_shear,'twist_release':list(model.twist_release) if model.twist_release is not None else None,'length_unit':model.length_unit,'gravity':model.gravity,
                   'elements':[{'node_i':e.node_i,'node_j':e.node_j,'laminate':e.laminate,'weight':e.weight} for e in model.elements],
                   'wall_groups':list(model.wall_groups)},
          **({'fe_rve':model.fe_rve} if model.fe_rve else {}),
          'materials':[asdict(m) for m in materials.values()],
          'laminates':[{'name':l.name,'plies':[{'material':p.material,'angle_deg':p.angle_deg,'thickness':p.thickness} for p in l.plies]} for l in laminates]}
    Path(path).write_text(json.dumps(data,indent=2),encoding='utf-8')


def load_project(path):
    d=json.loads(Path(path).read_text(encoding='utf-8')); materials={m['name']:Material(**m) for m in d['materials']}; laminates=[Laminate([Ply(**p) for p in l['plies']],l.get('name','Laminate')) for l in d['laminates']]; md=d['model']
    model=RVEModel(np.array(md['nodes']),[Element(**e) for e in md['elements']],np.array(md['a1']),np.array(md['a2']),md['height'],periodic_pairs=np.array(md['periodic_pairs'],float) if md.get('periodic_pairs') is not None else None,periodic_tol=md.get('periodic_tol',1e-8),anchor_node=md.get('anchor_node',0),fix_anchor_rotation=md.get('fix_anchor_rotation',False),wall_element_type=md.get('wall_element_type','paper'),strip_bending_mode=md.get('strip_bending_mode','cylindrical'),length_unit=md.get('length_unit','mm'),gravity=md.get('gravity',9.80665),name=md.get('name','RVE'),wall_groups=list(md.get('wall_groups',[])),wall_shear=md.get('wall_shear','timoshenko'),twist_release=tuple(md['twist_release']) if md.get('twist_release') is not None else None)  # files without the key predate 1.11 and used Timoshenko walls
    if d.get('fe_rve'):
        from .fe_rve import normalise_fe_rve_spec
        model.fe_rve=normalise_fe_rve_spec(d['fe_rve'])
    return model,materials,laminates
