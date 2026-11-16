"""Editable per-project or per-user material/geometry/layup libraries.

Editable source checkouts keep their original `data/` and `plugins/` folders.
Installed wheels initialise editable libraries under ~/.rve_studio_v1 on first
library use, so a non-admin installation never attempts to write site-packages.
"""
from __future__ import annotations
import json
import shutil
from pathlib import Path
from .models import Material

import sys

_SOURCE_ROOT = Path(__file__).resolve().parents[1]
_DEFAULTS_ROOT = Path(__file__).resolve().parent / 'defaults'
USER_ROOT = Path.home() / '.rve_studio_v1'


def library_root(source_root: Path = _SOURCE_ROOT, frozen: bool | None = None) -> Path:
    """Folder holding the editable data/ and plugins/ libraries.

    A source checkout edits its own data/ folder. An installed package, and the Windows application built with
    PyInstaller (whose files may sit in a read-only Program Files folder), use ~/.rve_studio_v1."""
    frozen = bool(getattr(sys, 'frozen', False)) if frozen is None else frozen
    if not frozen and (source_root / 'data' / 'materials.json').is_file():
        return source_root
    return USER_ROOT


PACKAGE_ROOT = library_root()
DATA_DIR = PACKAGE_ROOT / 'data'


def ensure_library_dirs() -> None:
    """Copy packaged defaults to a writable location without overwriting edits."""
    if PACKAGE_ROOT == _SOURCE_ROOT:
        return
    for subdir in ('data', 'plugins'):
        dest = PACKAGE_ROOT / subdir
        dest.mkdir(parents=True, exist_ok=True)
        source = _DEFAULTS_ROOT / subdir
        if source.exists():
            for file in source.iterdir():
                if file.is_file() and not (dest / file.name).exists():
                    shutil.copy2(file, dest / file.name)


def load_json(path: Path):
    ensure_library_dirs()
    return json.loads(Path(path).read_text(encoding='utf-8'))


def save_json(path: Path, data):
    ensure_library_dirs()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding='utf-8')


def load_material_library(path: Path | None = None) -> list[dict]:
    return load_json(path or DATA_DIR / 'materials.json')


def load_geometry_library(path: Path | None = None) -> list[dict]:
    return load_json(path or DATA_DIR / 'geometries.json')


def load_layup_library(path: Path | None = None) -> list[dict]:
    return load_json(path or DATA_DIR / 'layups.json')


def material_mapping(rows: list[dict] | None = None) -> dict[str, Material]:
    if rows is None:
        rows = load_material_library()
    out = {}
    for r in rows:
        out[r['name']] = Material(
            name=r['name'], E1=float(r['E1']), E2=float(r['E2']), G12=float(r['G12']), nu12=float(r['nu12']),
            Gth=float(r.get('Gth', r['G12'])),
            density_kg_m3=float(r['density_kg_m3']) if r.get('density_kg_m3') is not None else float('nan'),
            source=r.get('source', ''), source_url=r.get('source_url', ''),
            cost_per_kg=float(r['cost_per_kg']) if r.get('cost_per_kg') not in (None, '') else float('nan'))
    return out
