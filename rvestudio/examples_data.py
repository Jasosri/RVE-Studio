from __future__ import annotations
import math
import numpy as np
from .models import Material,Ply,Laminate,Element,RVEModel
from .geometry import hexagonal


def default_materials():
    return {
        'Boron/Epoxy': Material('Boron/Epoxy',204000,18500,5590,0.23,5590,2080),
        'E-glass/Epoxy': Material('E-glass/Epoxy',38600,8270,4140,0.26,4140,1800),
        'Carbon/Epoxy example': Material('Carbon/Epoxy example',131000,10300,6900,0.22,6900,math.nan),
    }

LAYUPS={
  1:[45,-45,0,-45,45],2:[0,0,90,0,0],3:[90,0,0,0,90],4:[0,45,90,45,0],5:[45,0,90,0,45],6:[45,-45,45,-45,45]
}


def paper_hexagonal_case(config:int=1):
    mats=default_materials(); g=hexagonal({'L1':50.,'L2':50.}); h=150.; rho=0.0722; area=abs(np.linalg.det(np.column_stack((g.a1,g.a2)))); twall=rho*area/(50+2*50); tp=twall/5
    G='E-glass/Epoxy'; B='Boron/Epoxy'
    matseq={1:[G]*5,2:[G]*5,3:[G]*5,4:[G]*5,5:[G]*5,6:[G]*5,7:[B,G,G,G,B],8:[G,B,G,B,G],9:[B,G,G,G,B],10:[G,B,G,B,G],11:[G,G,B,G,G],12:[G,G,G,G,B],13:[B,G,G,B,G],14:[G,G,G,G,B]}
    ang={1:[45,-45,0,-45,45],2:[0,0,90,0,0],3:[90,0,0,0,90],4:[0,45,90,45,0],5:[45,0,90,0,45],6:[45,-45,45,-45,45],7:[45,-45,0,-45,45],8:[0,0,90,0,0],9:[0,0,90,0,0],10:[90,0,0,0,90],11:[0,45,90,45,0],12:[45,0,90,0,45],13:[0,0,90,0,0],14:[45,-45,0,-45,45]}
    if config not in matseq: raise ValueError('config must be 1..14')
    lam=Laminate([Ply(m,a,tp) for m,a in zip(matseq[config],ang[config])],name=f'Paper configuration {config}')
    elems=[Element(int(i),int(j),0,float(w)) for (i,j),w in zip(g.connectivity,g.wall_weights)]
    model=RVEModel(g.nodes,elems,g.a1,g.a2,h,periodic_tol=1e-9,anchor_node=0,wall_element_type='paper',strip_bending_mode='cylindrical',length_unit='mm',name=f'Hexagonal config {config}')
    return model,mats,[lam]

# Compatibility with earlier Python project integrations.
thesis_hexagonal_case = paper_hexagonal_case
