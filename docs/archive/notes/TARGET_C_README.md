# RVE Studio v1 — Cylindrical target-C optimisation release

The core RVE Studio API and GUI are preserved. A reproducible cylindrical-only target-C inverse-design study is in `research/cylindrical_target_C/`. It selects from four built-in topologies, E-glass/Boron/hybrid ply stacks, 3/5/7/9 plies and symmetric fibre orientations; it ranks minimum areal mass subject to the full 6x6 stiffness target, explicit thickness limits, illustrative cost bounds and three load-dependent buckling screens.

To verify: first install dependencies (`pip install -e .`), then run `python research/cylindrical_target_C/verify_cylindrical_results.py`. This has been tested by re-solving nine archived projects. To regenerate fifteen figures: `python research/cylindrical_target_C/make_figures.py`. The core GUI may still offer other wall models, but **every reported research case study intentionally selects `paper` + `cylindrical`**. This software release has not been visually smoke-tested in PySide6.
