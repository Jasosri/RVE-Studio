> **v1.8.1 update:** The default campaign uses six PBC macro-reaction load cases to
> calculate the full C^H directly, not 21 cases. Set `--identification energy21`
> to use the original 21-case energy method. See `DIRECT_CH_PBC.md`.

# FE validation: cylindrical native RVE versus independent Abaqus S4 shell model

## Status and scope

RVE Studio Python **1.7.5** adds a reproducible FE validation *campaign* around the existing
v1.7.4 Abaqus/Standard periodic CH21 backend. **Preparing a campaign does not
run Abaqus.** An Abaqus executable, licence and job environment must be available
before independent FE matrices can be claimed. The shipped campaign has `PENDING`
FE results. No manufactured FE values, artificial convergence curves or solved
buckling mode shapes are included.

The FE path uses conforming **S4 layered shell** vertical wall midsurfaces, with
ply-specific material and orientation. Each primitive cell or exact tiled
parallelogram supercell uses a strict periodic-mesh-pair and rigid-body audit.
Six pure macroscopic engineering-strain cases and fifteen paired cases identify
all 21 independent components of the symmetric effective 6x6 matrix from Abaqus
`ALLSE` energy. The six-case macro-reaction reconstruction provides an independent
energy/work consistency check.

This is **shell FE**, not 3D solid FE. It is not, by itself, a non-linear buckling
study. The native cylindrical/Timoshenko wall calculation uses different height
and junction assumptions than the periodic Abaqus model. A disagreement may arise
from mechanical model-form differences as well as mesh error; it is not automatically
an implementation bug.

## Windows setup and execution

1. Install Python and dependencies, with Abaqus/Standard accessible (e.g. `abaqus` or `abq2021`).
2. In the main software directory run `python -m pip install -e .`.
3. The three optimised projects are in `examples/fe_validation/`.
4. Prepare the full independent FE mesh campaign:

```powershell
python -m rvestudio.cli fe-campaign prepare --root FE_Validation_Runs `
  --project examples/fe_validation/P00050.rvejson `
  --project examples/fe_validation/P00148.rvejson `
  --project examples/fe_validation/P00004.rvejson `
  --target examples/fe_validation/C_target_MPa.csv `
  --meshes 10 5 2.5 --pbc-domain lattice-supercell `
  --supercell-x 2 --supercell-y 2 --strain 0.001 --cpus 4
```

Or run on the pre-generated package supplied with this release (skip the prepare step):

```powershell
python -m rvestudio.cli fe-campaign run --root FE_Validation_Runs --abaqus-command abaqus
python -m rvestudio.cli fe-campaign collect --root FE_Validation_Runs
```

Each mesh level involves 21 separate linear elastic Abaqus jobs. The packaged
campaign comprises **3 project designs × 3 mesh sizes × 21 strain states = 189
Abaqus jobs** and a potentially significant computational workload. Run on a
licensed machine with adequate memory and disk space. Completed strain-state
results are restartable; `run` reuses valid saved checkpoints.

## Outputs and verification

- `projects/`: exact input `.rvejson` snapshots; `C_target_MPa.csv`.
- `analytical/<design>/`: native cylindrical `CH_cylindrical_MPa.csv` and physical properties.
- `abaqus_runs/<design>_<mesh>/`: 21 ready-to-run `.inp` decks, ODB extraction
  scripts, a strict PBC audit and a lattice-domain preview SVG.
- `FE_VALIDATION_REPORT.md`: provenance-aware status table and completed comparisons.
- `fe_convergence_summary.csv`: 6x6 diagonal components, relative FE/native
  matrix errors, and smallest target-generalised eigenvalue (only when solved).
- `fe_matrices/<design>/`: the full Abaqus `C_FE` 6x6 matrix **only when solved**.
- `figures/`: actual RVE geometry and native-vs-target figures. FE convergence
  figures appear only when genuine completed FE matrices are present.

The FE strain order is `[11,22,33,23,13,12]` with engineering shear.
Internal `C44` means **2–3 shear**, `C55` means **1–3 shear**. The
publication-oriented equivalent shell convention has `Q44=1–3`, `Q55=2–3`.
Compare *like against like* (raw 6x6 C against raw 6x6 C, not C against Q).

Recommended acceptance checks (define tolerances before analysing results):
(1) all PBC audit tests pass; (2) `ALLSE > 0` for every nonzero applied strain;
(3) the identified symmetric C is positive definite; (4) energy- and
reaction-derived C agree to an acceptably small numerical residual; (5) each
critical stiffness component converges from 10→5→2.5 mm; (6) the FE-native
matrix discrepancy is presented as a **model comparison**, not as a mesh error.

An external 3D-solid FE and eigenbuckling study is still needed before asserting
full model validation of the optimiser's buckling screening.
