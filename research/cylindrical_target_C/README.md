# RVE Studio v1 — Cylindrical-only five case studies / minimum-mass full-C inverse design

## Purpose and scope

This package supports the **18-page manuscript** `RVE_Studio_Cylindrical_Only_Five_Case_Studies_Manuscript.docx`. **All analysed stiffness and wall-recovery results are paper-style cylindrical bending only** (`wall_element_type='paper'`, `strip_bending_mode='cylindrical'`). The installed RVE Studio program may support other wall formulations, but no other formulation enters these results.

The target 6x6 stiffness matrix (MPa, internal engineering Voigt `[11,22,33,23,13,12]`) is in `results/C_target_MPa.csv`. The target originates from 80% of the journal configuration-2 example. The search requires the **smallest generalised eigenvalue** of (`C_achieved`,`C_target`) >= **1.05**, i.e. the achieved tensor has at least 5% more strain energy for *every macroscopic strain state*. **It does not require matching individual stiffness entries**. Matrix Frobenius mismatch is reported as a separate diagnostic.

## Decision variables

- Built-in geometry functions: HEX (5 nodes, 3 walls), SQUARE (7 nodes, 4 walls), MIXED (7 nodes, 6 walls), TRI (9 nodes, 6 walls).
- Geometry parameters: L1/L2, Lx/Ly, Lh/R according to source `data/geometries.json` and `source_libraries/geometries.json`.
- Materials: all E-glass/epoxy; all boron/epoxy; alternating glass/boron; alternating boron/glass.
- Individual equal-thickness ply count 3, 5, 7 or 9, with symmetric angle sequences from the discrete palette 0,+45,-45,90 degrees.
- Relative density scaled by a common positive wall-thickness factor; there is no arbitrary free optimisation of every wall thickness in this release.

## Constraints

- Full-matrix minimum generalised stiffness factor >= **1.05**.
- Relative density <= **0.15**.
- Total wall thickness >= **0.35 mm** and every ply thickness >= **0.08 mm**.
- Maximum wall-thickness ratio <= **1.9**.
- Illustrative core-material cost index <= **50 units/m2** (glass=1 per kg, boron=3 per kg; **not monetary pricing**).
- Linear reduced-order buckling multiplier >= **1.5** for each of three macroscopic stress vectors (MPa):
  - uniaxial `[-1,0,0,0,0,0]`
  - biaxial `[-1,-0.5,0,0,0,0]`
  - combined `[-1,-0.5,0,0,0,0.25]`
- Core height **150 mm**; *core-only mass*, excluding facesheets, adhesive and joining material.

## Computed results

- **4,128** solved initial combinations at rho*=0.15, **3,297** above the 1.05 stiffness margin at that density.
- **1,177** thickness-refined profiles: 958 prioritised profiles plus 219 additional exploratory profiles archived from the earlier broader refinement, with overlapping records de-duplicated.
- **136** candidates underwent three-load buckling refinement; **110** passed every declared screen. The shortlist is not a global search over the complete design space.
- **Minimum sampled mass among the 110 passing designs: P00050**, built-in HEX, all boron/epoxy, 3-ply [0/45/0], L1=L2=52.5 mm, relative density 0.0399463, density 83.0883 kg/m3, areal mass 12.4632 kg/m2, minimum stiffness factor 1.050524, cost index 37.3897 units/m2 and 6x6-basis buckling factors 5.3635, 4.6040 and 3.6860. The relative matrix Frobenius gap is **157.7%**: this design meets a **minimum stiffness inequality**, not the desired matrix coefficients exactly.
- **Near-mass hybrid alternative P00148**: HEX, boron/glass alternating, 5-ply [0/90/45/90/0], 83.5457 kg/m3, cost 28.4259 units/m2, minimum stiffness factor 1.05720, minimum buckling factor 3.6022, matrix gap 87.8%.

## Verification and reanalysis

1. Obtain the previously delivered software source ZIP (`RVE_Studio_Node_Explorer_3_to_10.zip`) or install `RVE Studio` as an editable pip package. The **bundled cylindrical optimisation release** also includes this software.
2. Install dependencies: `pip install -e /path/to/RVE_Studio_...` or `pip install /path/to/RVE Studio-1.0.0-py3-none-any.whl` plus required numerical libraries.
3. Run `python verify_cylindrical_results.py` from this research root: it reloads the nine selected `.rvejson` projects and recomputes each complete stiffness matrix and all three 6x6-basis buckling factors. It checks the 4,128/1,177/136/110 dataset counts and reports the winner.
4. Run `python make_figures.py` to regenerate all fifteen included figures from **cylindrical-only CSV data** (`prior_cyl` and `results`).
5. Run `python build_report.py` to generate a new Word document within this folder (requires `python-docx`).
6. The full search driver `run_cylindrical_study.py` has `screen`, `refine`, `buckling`, and `verify` stages. Repeating the long stages overwrites data. For reproducible *batch* refinement the original selected 958 profiles are processed in ordered index ranges using `RVE_START`, `RVE_STOP` environmental variables and combined with `drop_duplicates('ID')`. The supplemental auxiliary 219 records are retained in `results/stiffness_refined_0_350.csv` and the supplied merged full table. The screened shortlist is based on the supplied full merged dataset.

## File inventory

- `run_cylindrical_study.py` — actual solver-calling search driver and design palette.
- `verify_cylindrical_results.py` — independent direct project re-solves and tensor checks.
- `make_figures.py`, `build_report.py` — figure and manuscript source.
- `source_libraries/` — frozen constituent properties and geometry metadata.
- `target/results/C_target_6x6_MPa.csv`, `results/C_target_MPa.csv` — fixed target, duplicated for convenient use.
- `results/search_profiles.json`, `results/search_spec.json` — seeded design inputs and all constraint definitions.
- `results/initial_screen.csv`, `stiffness_refined.csv`, `buckling_refined.csv`, `verified_saved_projects.csv` — numerical tables.
- `results/project_*.rvejson` and `C_achieved_*.csv` — saved designs and achieved 6x6 tensors.
- `prior_cyl/` — cylindrical-only subsets of the previous 14-configuration/atlas/robust/buckling/custom-geometry datasets.
- `figures/F01_*.png` to `F15_*.png` — report figures.

## Limitations and research reference

The minimisation is a **finite-sample optimisation** and 136-profile stability shortlist, **not** a mathematical proof of global optimality. The wall buckling calculation uses simply-supported laminated CLPT/Ritz plates and is **not** validated whole-core/panel instability. Manufacturing and cost limits are declared assumptions. The paper benchmark reproduces **70 of 84** analytical stiffness entries within 1%, with unresolved Q66 discrepancy up to **14.54%**.

Sriharan, J., Dias, M., Fernando, D., Adhikari, S. (2024). *A unified approach for the prediction of the effective properties of laminated composite cellular core*. **Engineering Structures**, 299, 117106. DOI: https://doi.org/10.1016/j.engstruct.2023.117106.
