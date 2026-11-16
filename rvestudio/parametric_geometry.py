from __future__ import annotations

from typing import Sequence
import numpy as np

from .geometry import Geometry


def scale_xy_reference_geometry(
    nodes: Sequence[Sequence[float]],
    connectivity: Sequence[Sequence[int]],
    a1: Sequence[float],
    a2: Sequence[float],
    wall_groups: Sequence[str],
    wall_weights: Sequence[float],
    *,
    geometry_id: str,
    name: str,
    width: float,
    height: float,
    reference_width: float | None = None,
    reference_height: float | None = None,
) -> Geometry:
    """Return an x/y-scaled copy of a reference periodic geometry.

    This is the pure-Python operation used by the v1.3 Geometry Definition Wizard
    for arbitrary custom topologies.  Scaling acts on coordinate differences from
    the reference bounding-box minimum and on the corresponding lattice-vector
    components.  Connectivity and wall weights are unchanged.
    """
    pts=np.asarray(nodes,float)
    if pts.ndim != 2 or pts.shape[1] != 2 or len(pts)==0:
        raise ValueError("nodes must be a non-empty N x 2 array")
    xmin,ymin=np.min(pts,axis=0); xmax,ymax=np.max(pts,axis=0)
    w0=float(reference_width if reference_width is not None else xmax-xmin)
    h0=float(reference_height if reference_height is not None else ymax-ymin)
    width=float(width); height=float(height)
    if w0<=0 or h0<=0 or width<=0 or height<=0:
        raise ValueError("Reference and target width/height must be positive")
    sx=width/w0; sy=height/h0
    out=pts.copy(); out[:,0]=xmin+(out[:,0]-xmin)*sx; out[:,1]=ymin+(out[:,1]-ymin)*sy
    u=np.asarray(a1,float); v=np.asarray(a2,float)
    u2=np.array([u[0]*sx,u[1]*sy]); v2=np.array([v[0]*sx,v[1]*sy])
    return Geometry(
        str(geometry_id),str(name),out,np.asarray(connectivity,int),u2,v2,
        list(wall_groups),np.asarray(wall_weights,float),
        "Wizard-defined W/H parametric family",width,height,
    )


def make_scaled_xy_plugin_source(
    *,
    geometry_id: str,
    name: str,
    nodes: Sequence[Sequence[float]],
    connectivity: Sequence[Sequence[int]],
    a1: Sequence[float],
    a2: Sequence[float],
    wall_groups: Sequence[str],
    wall_weights: Sequence[float],
    width_parameter: str,
    height_parameter: str,
) -> str:
    """Create the importable plugin source used by the Geometry Definition Wizard."""
    pts=np.asarray(nodes,float)
    if pts.ndim != 2 or pts.shape[1] != 2 or len(pts)==0:
        raise ValueError("nodes must be a non-empty N x 2 array")
    xmin,ymin=np.min(pts,axis=0); xmax,ymax=np.max(pts,axis=0)
    width=float(xmax-xmin); height=float(ymax-ymin)
    if width<=0 or height<=0:
        raise ValueError("Custom width/height parameterisation requires non-zero width and height")
    return (
        '"""Parametric geometry created by the RVE Studio Geometry Definition Wizard."""\n'
        'from __future__ import annotations\n'
        'import numpy as np\n'
        'from rvestudio.parametric_geometry import scale_xy_reference_geometry\n\n'
        f'GEOMETRY_ID={str(geometry_id)!r}\n'
        f'_REF_NODES={pts.tolist()!r}\n'
        f'_CONNECTIVITY={np.asarray(connectivity,int).tolist()!r}\n'
        f'_A1={np.asarray(a1,float).tolist()!r}\n'
        f'_A2={np.asarray(a2,float).tolist()!r}\n'
        f'_GROUPS={list(wall_groups)!r}\n'
        f'_WEIGHTS={np.asarray(wall_weights,float).tolist()!r}\n'
        f'_W0={width!r}\n'
        f'_H0={height!r}\n\n'
        'def generate(p):\n'
        f'    W=float(p.get({width_parameter!r},_W0))\n'
        f'    H=float(p.get({height_parameter!r},_H0))\n'
        '    return scale_xy_reference_geometry(\n'
        '        _REF_NODES,_CONNECTIVITY,_A1,_A2,_GROUPS,_WEIGHTS,\n'
        f'        geometry_id=GEOMETRY_ID,name={str(name)!r},width=W,height=H,\n'
        '        reference_width=_W0,reference_height=_H0,\n'
        '    )\n'
    )
