# RVE Studio — Composite cellular core design

[![CI](https://github.com/Jasosri/RVE-Studio/actions/workflows/ci.yml/badge.svg)](https://github.com/Jasosri/RVE-Studio/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/ZENODO_DOI.svg)](https://doi.org/ZENODO_DOI)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue)

**RVE Studio** computes the effective stiffness of laminated-composite cellular cores — honeycombs and any other
periodic wall network — from a periodic representative volume element, tells you whether the wall model can be
trusted for your design, validates it against periodic shell finite elements with a mesh-convergence verdict, and
generates and screens large design databases. Desktop application (PySide6), command line and Python API.

![RVE Studio](rvestudio/gui/assets/rve3d_banner.png)

## What it does

| | |
|---|---|
| **Periodic RVE designer** | draw any wall network or load a parametric family; full 6×6 Cᴴ, shell Q matrix and engineering constants in about 1 ms |
| **Wall models** | Euler–Bernoulli or Timoshenko walls; cylindrical, **twist-corrected** (calibrated on 320 walls, editable coefficients), anisotropic-strip and free-moment bending |
| **Reliability rating** | every solve states whether the wall model is inside its verified range (bend–twist index Π, χ_B, β) |
| **FE design & validation** | periodic shell FE of the same design — built-in solver in seconds or Abaqus — full C from strain energy, DBC shell and 3D terms, Richardson extrapolation, GCI and a mesh verdict for every term |
| **Benchmark designs** | the 14 benchmark designs of the source paper with any wall model (84/84 reference values within 0.3 %) |
| **Design databases** | Latin-hypercube databases with exact relative density, parallel batch solving (≈ 20 000 designs per minute per core), optional wall-buckling screening |
| **Results explorer** | Ashby plots and power-law scaling, anisotropy, E(θ), sensitivity, Pareto and robust multi-objective screening, manufacturing rules |
| **Laminate tools** | ABD, neutral plane, S_eq, lamination parameters, invariants, bend–twist ratio |

## Install

```bash
pip install .                 # from this folder; Python >= 3.10
rvestudio-gui                 # desktop application
rvestudio --help              # command line
```

**Windows installer (no Python needed):** download `RVE_Studio_Setup_<version>.exe` from the
[releases page](https://github.com/Jasosri/RVE-Studio/releases), run it, and start RVE Studio from the Start menu. It installs
for the current user without administrator rights (or for all users), adds Start-menu shortcuts to the application and the
two guides, and can open `.rvejson` project files by double-click. Build it yourself with `build_windows_installer.bat`.

From source on Windows: run `install_windows.bat` once, then `START_RVE_STUDIO.bat`.

## Quick start

```bash
rvestudio benchmark                                           # reproduce the original model: 84/84 within 0.3 %
rvestudio benchmark --wall-beam timoshenko --bending twist-corrected
rvestudio single --config 7 --json                            # one design with its reliability rating
rvestudio generate-db cases.csv --cases 5000 --wall-beam timoshenko --bending twist-corrected
rvestudio batch cases.csv results.csv --workers 8
rvestudio fe check --root Study --paper-configs 1 --meshes 10 5 2.5   # shell FE with mesh verdicts
```

```python
from rvestudio.api import solve_example, summary, benchmark
res = solve_example(config=1)
print(summary(res)["engineering_constants"], summary(res)["reliability"]["level"])
print(benchmark(wall_shear="timoshenko", bending="twist-corrected")["within_tolerance"])
```

## Documentation

* [User Guide](docs/USER_GUIDE.md) — every workspace step by step, command line, API, troubleshooting
* [Theory and Modelling Guidance](docs/THEORY.md) — equations, wall models, verified limits, choosing a wall
  model, mesh convergence
* [API reference](docs/API_REFERENCE.md) · [FE design & validation](docs/FE_VERIFICATION.md) ·
  [Developer guide](docs/DEVELOPER_GUIDE.md) · [Changelog](CHANGELOG.md)

PDF versions of the user and theory guides are attached to every release (`docs/build_pdfs.sh`).

## Choosing a wall model

| Goal | Wall beam | Bending option |
|---|---|---|
| Reproduce the source paper | Euler–Bernoulli | Cylindrical |
| Most accurate for any stacking sequence (within 2 % of converged shell FE for \|ψ\| ≤ 0.85, β ≤ 0.10) | Timoshenko | Twist-corrected |
| Bounds on wall bending | any | Cylindrical (upper) / Free moment (lower) |

## Tests

```bash
pip install -e ".[dev]"
pytest                         # GUI tests run when PySide6 is available (QT_QPA_PLATFORM=offscreen headless)
```

Continuous integration runs the suite on Linux and Windows with Python 3.10 and 3.12.

## Citing

RVE Studio builds on earlier papers. If you use it, please cite the RVE Studio paper together with them, and the
software release through its Zenodo DOI (see [CITATION.cff](CITATION.cff)):

1. J. Sriharan et al. RVE Studio: open-source strain-energy homogenisation of laminated-composite cellular cores —
   verification of the full stiffness tensor, a bend–twist reliability limit and design tools. *In preparation.*
2. J. Sriharan, M. Dias, S. Adhikari, D. Fernando (2025). A semi-analytical model for predicting the shear buckling
   of laminated composite honeycomb cores in sandwich panels. *Composite Structures* 351, 118629.
   doi:[10.1016/j.compstruct.2024.118629](https://doi.org/10.1016/j.compstruct.2024.118629)
3. J. Sriharan, M. Dias, D. Fernando, S. Adhikari (2024). A unified approach for the prediction of the effective
   properties of laminated composite cellular core. *Engineering Structures* 299, 117106.
   doi:[10.1016/j.engstruct.2023.117106](https://doi.org/10.1016/j.engstruct.2023.117106)

## Licence

MIT — see [LICENSE](LICENSE).
