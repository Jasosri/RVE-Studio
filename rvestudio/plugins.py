"""Geometry plug-ins: built-in families plus every *.py in the plug-in folder that defines GEOMETRY_ID and generate().

The registry is cached per folder and rebuilt only when a plug-in file is added, removed or changed, so large
batches do not re-import the plug-ins for every case."""
from __future__ import annotations
import importlib.util
from pathlib import Path
from typing import Callable,Mapping
from .geometry import Geometry,BUILTIN_GEOMETRIES

_CACHE: dict = {}


def _signature(plugin_dir:Path):
    return tuple(sorted((f.name, f.stat().st_mtime_ns, f.stat().st_size) for f in plugin_dir.glob('*.py')))


def discover_geometry_plugins(plugin_dir:Path) -> dict[str,Callable[[Mapping[str,float]],Geometry]]:
    plugin_dir=Path(plugin_dir)
    if not plugin_dir.exists(): return dict(BUILTIN_GEOMETRIES)
    key=str(plugin_dir.resolve()); sig=_signature(plugin_dir)
    hit=_CACHE.get(key)
    if hit is not None and hit[0]==sig: return dict(hit[1])
    registry=dict(BUILTIN_GEOMETRIES)
    for file in sorted(plugin_dir.glob('*.py')):
        if file.name.startswith('_') or file.stem=='custom_geometry_template': continue
        spec=importlib.util.spec_from_file_location(f'rvestudio_user_{file.stem}',file)
        if spec is None or spec.loader is None: continue
        module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        geom_id=getattr(module,'GEOMETRY_ID',None); generate=getattr(module,'generate',None)
        if geom_id and callable(generate): registry[str(geom_id).upper()]=generate
    _CACHE[key]=(sig,registry)
    return dict(registry)
