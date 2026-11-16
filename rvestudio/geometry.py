from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, Mapping
import math
import numpy as np

@dataclass(slots=True)
class Geometry:
    id: str
    name: str
    nodes: np.ndarray
    connectivity: np.ndarray
    a1: np.ndarray
    a2: np.ndarray
    wall_groups: list[str]
    wall_weights: np.ndarray
    note: str=""
    length_a: float=math.nan
    length_b: float=math.nan


def hexagonal(p:Mapping[str,float])->Geometry:
    L1=float(p['L1']); L2=float(p['L2']); th=math.radians(120); a=2*L2*math.sin(th); c=L2*math.sin(th); b=L1+L2*math.cos(math.pi-th)
    nodes=np.array([[0,0],[-c,b],[a,0],[a-c,b],[0,L1]],float); conn=np.array([[0,4],[1,4],[3,4]],int)
    return Geometry('HEX','Hexagonal',nodes,conn,np.array([a,0.]),np.array([-c,b]),['vertical-L1','inclined-L2-A','inclined-L2-B'],np.ones(3),'Unequal-wall hexagonal; theta=120 deg',L1,L2)


def square_family(p:Mapping[str,float])->Geometry:
    Lx=float(p['Lx']); Ly=float(p['Ly']); nodes=np.array([[0,0],[-2*Lx,2*Ly],[2*Lx,0],[0,2*Ly],[-Lx,Ly],[Lx,Ly],[0,Ly]],float)
    conn=np.array([[0,6],[6,3],[4,6],[6,5]],int)
    return Geometry('SQUARE','Square-family',nodes,conn,np.array([2*Lx,0.]),np.array([-2*Lx,2*Ly]),['vertical-bottom','vertical-top','horizontal-left','horizontal-right'],np.ones(4),'Square/rectangular family',Lx,Ly)


def mixed_rhombus_triangular(p:Mapping[str,float])->Geometry:
    Lh=float(p['Lh']); R=float(p['R']); Ld=R*Lh
    if R<=1: raise ValueError('Mixed rhombus-triangular requires R=Ld/Lh > 1')
    q=2*math.sqrt(Ld*Ld-Lh*Lh); nodes=np.array([[0,0],[0,q],[2*Lh,0],[2*Lh,q],[0,q/2],[2*Lh,q/2],[Lh,q/2]],float)
    conn=np.array([[1,6],[3,6],[0,6],[2,6],[4,6],[6,5]],int)
    return Geometry('MIXED','Mixed rhombus-triangular',nodes,conn,np.array([2*Lh,0.]),np.array([0.,q]),['diag-top-left','diag-top-right','diag-bottom-left','diag-bottom-right','horizontal-left','horizontal-right'],np.ones(6),'Mixed rhombus-triangular',Lh,Ld)


def triangular_family(p:Mapping[str,float])->Geometry:
    Lh=float(p['Lh']); R=float(p['R']); Ld=R*Lh
    if R<=0.5: raise ValueError('Triangular family requires R=Ld/Lh > 0.5')
    q=2*math.sqrt(Ld*Ld-(0.5*Lh)**2); nodes=np.array([[0,0],[-Lh,q],[Lh,0],[0,q],[2*Lh,0],[Lh,q],[-0.5*Lh,q/2],[1.5*Lh,q/2],[0.5*Lh,q/2]],float)
    conn=np.array([[3,8],[5,8],[0,8],[2,8],[6,8],[8,7]],int)
    return Geometry('TRI','Triangular-family',nodes,conn,np.array([2*Lh,0.]),np.array([-Lh,q]),['diag-top-left','diag-top-right','diag-bottom-left','diag-bottom-right','horizontal-left','horizontal-right'],np.ones(6),'Triangular family',Lh,Ld)

BUILTIN_GEOMETRIES={'HEX':hexagonal,'SQUARE':square_family,'MIXED':mixed_rhombus_triangular,'TRI':triangular_family}
DEFAULT_GEOMETRY_LIBRARY=[
    {'active':True,'id':'HEX','name':'Hexagonal','function':'HEX','parameters':[{'name':'L1','min':35.,'max':70.},{'name':'L2','min':35.,'max':70.}],'weight':1.0,'solver':'paper','strip_mode':'cylindrical','description':'Built-in hexagonal family. L1 = wall length normal to a1; L2 = inclined-wall length.'},
    {'active':True,'id':'SQUARE','name':'Square-family','function':'SQUARE','parameters':[{'name':'Lx','min':35.,'max':70.},{'name':'Ly','min':35.,'max':70.}],'weight':1.0,'solver':'paper','strip_mode':'cylindrical','description':'Built-in square/rectangular family. Lx = wall length parallel to a1; Ly = wall length normal to a1.'},
    {'active':True,'id':'MIXED','name':'Mixed rhombus-triangular','function':'MIXED','parameters':[{'name':'Lh','min':30.,'max':60.},{'name':'R','min':1.5,'max':2.5}],'weight':1.0,'solver':'paper','strip_mode':'cylindrical','description':'Built-in mixed rhombus-triangular family. Lh = reference wall length parallel to a1; R = Ld/Lh.'},
    {'active':True,'id':'TRI','name':'Triangular-family','function':'TRI','parameters':[{'name':'Lh','min':35.,'max':70.},{'name':'R','min':0.8,'max':1.2}],'weight':1.0,'solver':'paper','strip_mode':'cylindrical','description':'Built-in triangular family. Lh = reference wall length parallel to a1; R = Ld/Lh.'},
]
