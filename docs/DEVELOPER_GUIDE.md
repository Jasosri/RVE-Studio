# RVE Studio — Developer Guide

## Package layout

```text
rvestudio/
  models.py              data classes (Material, Ply, Laminate, Element, RVEModel, results)
  laminate.py            Qbar, ABD, neutral plane, S_eq (exact Gauss), lamination parameters, invariants,
                         bending options, twist length and twist release (TWIST_RELEASE_DEFAULT)
  walls.py               strain-energy wall element (Euler-Bernoulli / Timoshenko), anisotropic strip
  periodic.py            periodic-pair detection and constraint matrices
  homogenization.py      assembly, KKT solve (dense LU <= 600 unknowns, sparse above), condensation
  properties.py          compliance, Q matrix, engineering constants
  physical.py            volume, mass, density, cost
  reliability.py         bend-twist index, chi_B, beta; ratings and notes (Theory, Section 9)
  validation.py          14 benchmark designs with any wall model (BENDING_OPTIONS)
  database.py            LHS generation, batch solving (cached libraries, chunked process pool)
  postprocess.py         derived metrics, Pareto, fits, sensitivity
  load_recovery.py       wall membrane resultants
  buckling.py            CLPT Ritz/Navier wall buckling
  shell_fe.py            built-in periodic shell FE (PBC strain energy, DBC)
  fe_rve.py              FE RVE selection
  fe_workflow.py         FE study folders, Abaqus orchestration, Richardson / GCI / mesh_verdict
  fe_report.py           HTML report and workbook
  abaqus*.py             Abaqus deck writers and extractors
  geometry*.py, plugins.py, parametric_geometry.py   geometry families, analysis, plug-ins
  io.py                  .rvejson projects
  api.py, cli.py         public Python API and command line
  gui/                   PySide6 interface (main_window, fe_page, sketch, wall_editor, ui_helpers, ...)
tests/                   pytest suite (GUI tests skip without PySide6)
docs/                    user guide, theory, API, FE validation, developer guide
```

## Set-up

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pytest                      # full suite; set QT_QPA_PLATFORM=offscreen for the GUI tests on a headless machine
ruff check --select F rvestudio tests
```

## Regression discipline

Any change to laminate, wall, periodic or homogenisation code is accepted only after

```bash
pytest
rvestudio benchmark                      # 84/84 within 0.3 % with the default wall model
rvestudio benchmark --wall-beam timoshenko --bending twist-corrected
```

`tests/test_all14.py` asserts the 0.3 % reproduction of the original model and the effect of every wall model;
`tests/test_reliability.py` the ratings and the twist-release calibration; `tests/test_fe_workflow.py` the
Richardson extrapolation and mesh verdicts.

## Conventions that must not change silently

* Voigt order `[11, 22, 33, 23, 13, 12]`; shell `Q44 = C1313`, `Q55 = C2323` (paper Eq. 5).
* Wall local axis 1-bar from Node I to Node J; ply angles from 1-bar.
* `wall_element_type="paper"` is the cylindrical wall; `"thesis"` and `"original"` are accepted aliases.
* `strip_bending_mode`: `cylindrical`, `freeMoment`, `twistRelease` (UI: Anisotropic strip, Free moment,
  Twist-corrected).
* Twist-release coefficients: `laminate.TWIST_RELEASE_DEFAULT = (2.0, 2.15, 0.73)`; `RVEModel.twist_release`
  is `None` for the default.

## Geometry plug-in contract

A file in `plugins/` defines

```python
GEOMETRY_ID = "MYCELL"

def generate(p):            # p: dict of parameter values
    return Geometry(...)    # nodes, zero-based connectivity, a1, a2, wall_groups, wall_weights
```

and `data/geometries.json` lists the parameter names and sampling bounds. The registry is cached and rebuilt when
a plug-in file changes.

## GUI architecture

`RVESketchView` edits only coordinates, connectivity, weights, stack IDs and lattice vectors; it never calls the
solver. `geometry_analysis.py` is the GUI-independent extraction and validation layer, and `RVEModel` is the
boundary between GUI and numerical core. `ui_helpers.WallModelOptions` is the shared wall-model selector of
Benchmark Designs and Generate DB. Global styling is `APP_QSS` in `main_window.py`; page accent colours are in
`PAGE_COLOURS`.

## Building the Windows application and installer

`build_windows_installer.bat` runs PyInstaller with `RVE_Studio.spec` (one-folder build `dist/RVE_Studio/RVE_Studio.exe`;
bundles `rvestudio/defaults`, `rvestudio/gui/assets` and the PDF guides in `rvestudio/docs`; unused Qt modules are
left out) and then Inno Setup 6 with `installer/RVE_Studio.iss` (`dist/RVE_Studio_Setup_<version>.exe`: per-user or
all-users install, Start-menu shortcuts, optional desktop icon and `.rvejson` association, uninstaller). The frozen
application keeps its editable libraries in `~/.rve_studio_v1` (`library.library_root`), because the install folder
may be read-only, and opens a `.rvejson` passed on the command line. Update `installer/version_info.txt` with the
version. The *Release* workflow builds, smoke-tests and attaches the installer on `windows-latest`.

## Releasing

1. Update `rvestudio/__init__.py` (`__version__`), `pyproject.toml`, `CITATION.cff`, `.zenodo.json` and `CHANGELOG.md`.
2. `pytest`, `rvestudio benchmark`.
3. Commit, tag `vX.Y.Z`, push the tag. The *Release* workflow builds the wheel and sdist and attaches them to a
   GitHub release; with the Zenodo–GitHub integration enabled, Zenodo archives the release and mints a DOI.
