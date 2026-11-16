> **v1.10.0 — one automatic FE RVE for every core.** Every DBC and PBC model (single run, study, campaign,
> built-in check) now picks the FE RVE with the published honeycomb rule, for square, triangular, mixed and
> user cores too. **Abaqus FE → Preview FE RVE** shows it; **Verify without Abaqus** checks any design against
> the built-in DBC and 3D-periodic shell FE in seconds. Studies write `FE_STUDY_REPORT.html` and
> `FE_STUDY_RESULTS.xlsx`. See docs/VERIFICATION.md.
>
> **v1.9.1 — finite-height FE verification on the published RVE.** For the Paper A study (14 configurations,
> DBC + PBC, three meshes, built-in shell check) use `fe-study` (docs/FE_STUDY.md) or
> `RUN_FE_STUDY_PAPER_A_WINDOWS.bat`. Both FE problems now use the source study's RVE: 173.2 × 150 × 150 mm
> for the 50 mm honeycomb, four primitive cells, walls lying in the x1 faces at half plies
> (docs/FINITE_HEIGHT_RVE_BOUNDARY_CONDITIONS.md). Study folders prepared with 1.9.0 are refused.
> The `fe-campaign` workflows below are unchanged.

# v1.8.1 direct PBC C^H update

The **default** periodic FE analysis now uses **six** independent macroscopic strain
states and extracts all six macro reactions from each case. The full homogenised
stiffness `C^H` is calculated directly from periodic FE reactions; elastic moduli
are optional derived quantities. Select **Direct C^H from PBC macro reactions**
in the GUI, or use `--identification direct_reaction6` in CLI. The original
21-case energy reconstruction remains available via `--identification energy21`.
For new runs use [docs/DIRECT_CH_PBC.md](docs/DIRECT_CH_PBC.md).

The v1.8.0 instructions below refer to the former 21-case default, which is
now an optional mode. Abaqus simulations remain pending until executed.

---

# RVE Studio v1.8.0 — Cylindrical Ply Optimisation + Abaqus FE

This release merges the `RVE_Studio_Cylindrical_Ply_Optimisation` numerical GUI/API/research examples with the Abaqus S4 periodic FE module from RVE Studio v1.7.5.

## What's integrated

* Existing 14-configuration paper benchmark, cylindrical analytical solver, CLT layups, geometry and materials libraries, multi-objective Results Explorer, costs, constraints, wall loads and reduced-order buckling are retained.
* Abaqus FE workspace: S4 composite shells, explicit per-wall plies, 21-state complete 6×6 stiffness identification, strict periodic PBC audits, mesh convergence campaigns and post-solver comparisons.
* **NEW symmetric rectangular RVE choice**: finds a commensurate rectangular cut from an existing periodic lattice and searches phase shifts for both horizontal and vertical mirror symmetry. Checks identical stacks on mirrored walls and conservatively permits only 0°/90° ply axes. Unsupported geometries/layups are **rejected**, not silently assumed symmetric.
* **Unrestricted rectangular periodic cut** remains available for published/optimised off-axis layups, including P00050 (0°/45°/0°), even though they are not mirror-certified.
* The primitive and the exact parallelogram-lattice supercell remain available for arbitrary lattice vectors.

## GUI steps

1. On Windows run `install_windows.bat`, then `START_RVE_STUDIO.bat`.
2. Open **Abaqus FE**. Choose **Select symmetric in-built RVE...** and select `SQUARE` or `MIXED`, each preloaded with a symmetric glass/epoxy `[0/90/0]` stacking sequence and cylindrical analytical settings. The active Single RVE Designer is updated.
3. Select **Periodic homogenisation — full 21-term Cᴴ** (the new FE tab defaults to it), and **Symmetric rectangular RVE** for PBC domain.
4. Choose `1×1` rectangular repeats and an S4 mesh size, then **Generate 21-case package**. The periodic boundary audit must say PASS.
5. For mesh convergence use **Prepare FE mesh-convergence campaign**. Run on a licensed Abaqus workstation; use the matching `fe-campaign run` and `fe-campaign collect` commands below.
6. To analyse the original minimum-mass **P00050** candidate: load `examples/fe_validation/P00050.rvejson` and use **Rectangular periodic cut — clips boundary cells**, NOT *Symmetric rectangular*. Its off-axis plies and geometry are not mirror-certified.

## Reproducible command line

```bash
python -m pip install -e .
python -m pytest -q

# Create 21 Abaqus jobs for the symmetric SQUARE example (no Abaqus needed for preparation):
python -m rvestudio.cli fe-campaign prepare \
  --project examples/symmetric_rectangular_fe/SQUARE_symmetric_0_90_0.rvejson \
  --root FE_Runs_SymmetricSquare --meshes 15 7.5 \
  --pbc-domain symmetric-rectangular --supercell-x 1 --supercell-y 1

# This candidate is not mirror-symmetric; generate an unrestricted rectangular FE domain:
python -m rvestudio.cli fe-campaign prepare \
  --project examples/fe_validation/P00050.rvejson \
  --root FE_Runs_OptimalP00050 --meshes 15 7.5 \
  --pbc-domain rectangular-supercell --supercell-x 1 --supercell-y 1

# On a licensed Abaqus computer:
python -m rvestudio.cli fe-campaign run --root FE_Runs_SymmetricSquare --abaqus-command abaqus
python -m rvestudio.cli fe-campaign collect --root FE_Runs_SymmetricSquare
```

The self-contained FE deck examples and geometry sketch are under `examples/symmetric_rectangular_fe/`.

## Numerical and physical distinctions

- **Cylindrical analytical bending** is the only wall reduction used as the FE-comparison reference. FE shell bending remains an independent structural formulation; selecting a symmetric geometry does not impose a cylindrical constraint on Abaqus.
- Rectangular FE means an **orthogonal periodic analysis window**, not just the rectangular bounding box of a skew primitive cell. The primitive lattice must have a commensurate orthogonal integer supercell. Wall segments intersected by the boundary are clipped, preserving total material measure and periodic correspondence.
- The new check is an **in-plane geometry/material mirror-screen**; a 0°/90° stack does not guarantee every anisotropic response decouples, and no symmetry reduction of the Abaqus degrees of freedom is applied. Full periodic PBCs are still applied on all rectangular faces and the vertical direction.
- **Full 3D Cᴴ internal Voigt order:** `[11,22,33,23,13,12]`. The **journal shell convention** is `Q44 = 13 shear` and `Q55 = 23 shear`. Do not interchange the internal C44 and C55 entries.
- The FE backend is **Abaqus/Standard S4 shell midsurfaces**, not a 3D solid continuum or experimental validation.
- All FE runs shipped in the demo campaign are **PREPARED/PENDING**; Abaqus is not available in this environment, and no FE agreement claims are made.

## Evidence shipped

`examples/symmetric_rectangular_fe/rectangular_fe_geometries.png` illustrates the rectangular wall networks (not Abaqus results). The accompanying optional campaign ZIP contains SQUARE/MIXED symmetry-audited input decks and a P00050 unrestricted rectangle, plus project JSON, PBC audit, macro-strain decks and analytical baseline files.

## Limitations

Reflective geometrical symmetry alone does not establish **mechanical** mirror symmetry when off-axis composite plies, spatially varying stacks or wall-junction orientations are involved. The strict selector conservatively accepts only the supported 0°/90°-ply layouts with identical ply stacks on mirrored segments; all other layouts can still be evaluated with unrestricted periodic PBCs.

Run genuine Abaqus FE results and a mesh-convergence check before using results for journal validation. Buckling predictions in the optimisation module are still preliminary reduced-order screens.
