"""Copy this file to create a custom periodic RVE geometry plugin."""
from __future__ import annotations
import numpy as np
from rvestudio.geometry import Geometry

GEOMETRY_ID='CUSTOM'

def generate(p):
    # Example parameters. Add matching entries to data/geometries.json.
    Lx=float(p.get('Lx',50.0)); Ly=float(p.get('Ly',50.0))
    nodes=np.array([[0,0],[-2*Lx,2*Ly],[2*Lx,0],[0,2*Ly],[-Lx,Ly],[Lx,Ly],[0,Ly]],float)
    connectivity=np.array([[0,6],[6,3],[4,6],[6,5]],int)  # ZERO-BASED node IDs
    return Geometry(
        id=GEOMETRY_ID,
        name='Custom periodic RVE',
        nodes=nodes,
        connectivity=connectivity,
        a1=np.array([2*Lx,0.0]),
        a2=np.array([-2*Lx,2*Ly]),
        wall_groups=['wall-1','wall-2','wall-3','wall-4'],
        wall_weights=np.ones(4),
        note='Replace this template with your topology.',
        length_a=Lx,
        length_b=Ly,
    )
