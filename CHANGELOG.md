# Changelog

## 1.0.0 — RVE Studio, first public release

First public release of RVE Studio, following the development builds 1.0–1.15.1 (history below). Projects,
databases and scripts of the development builds open and run unchanged.

**Reorganised**
- Package `rvestudio`; commands `rvestudio` (alias `rve-studio`) and `rvestudio-gui`; window title *RVE Studio*;
  libraries in `~/.rve_studio_v1` as before.
- *Validation & Benchmarks* is now **Benchmark Designs**; *FE Verification* is now **FE Design & Validation**
  (steps *Design the FE model*, *Validate now*, *Validate with Abaqus*).

**Wall models everywhere**
- **Editable twist-corrected wall:** the coefficients of η = (a + b ψ²) λⁿ (calibrated a = 2.0, b = 2.15,
  n = 0.73) can be edited in the Solver tab, on Benchmark Designs and in Generate DB, saved with the project
  (`twist_release`), passed in the API (`RVEModel(twist_release=(a, b, n))`, `benchmark(twist_release=...)`) and on
  the command line (`--twist A B N`). Designs solved with edited coefficients are rated *use with care*.
- **Benchmark Designs and Generate DB** take any wall beam (Euler–Bernoulli, Timoshenko) and any bending option
  (Cylindrical, Twist-corrected, Anisotropic strip, Free moment): GUI selector, `rvestudio benchmark|generate-db
  --wall-beam ... --bending ... --twist ...`, `api.benchmark(wall_shear=..., bending=...)`. Database rows carry
  `SolverType`, `StripMode`, `WallShear`, `TwistA/B/N`.
- Benchmark Designs no longer shows published FE values; it compares the chosen wall model with the reference
  values of the cylindrical Euler–Bernoulli model. Shell FE validation is under FE Design & Validation.

**Reliability rating**
- The Euler–Bernoulli cylindrical wall is graded on Π + β: its wall-shear error (≈ 100β %) and bend–twist error
  (≈ 100Π %) act on the same terms and add. Rated against the converged shell FE of the 320 calibration walls,
  every design rated *reliable* is within 2.0 % and every design rated *use with care* within 5.2 %, for every
  wall beam and bending option. (Grading the two measures separately had let errors of up to 3.4 % through as
  *reliable*.) Benchmark design 1 with the default wall is now rated *use with care* (3.9 % against FE).

**Windows installer**
- `RVE_Studio_Setup_<version>.exe` (PyInstaller + Inno Setup): no Python needed, per-user install without
  administrator rights, Start-menu shortcuts to the application and both guides, optional desktop icon and `.rvejson`
  file association, uninstaller. Built by `build_windows_installer.bat` and by the release workflow; new application
  icon. The installed application keeps its libraries in `~/.rve_studio_v1` and opens a project given on the command line.

**Help and citation**
- Help (page, menu, F1) now opens the User Guide and the Theory Guide (PDF, shipped with the package).
- The About page, README, guides and `CITATION.cff` ask users to cite the RVE Studio paper together with the papers
  it builds on (Engineering Structures 2024; Composite Structures 2025).
- The cylindrical Euler–Bernoulli wall is labelled "default" (and "original model" in the documentation) instead of
  "published model".

**Geometry**
- The example plug-in is now a Kagome (trihexagonal) cell, `plugins/kagome_geometry.py` (parameters L1, L2),
  replacing the re-entrant hexagon example.

**Interface**
- Dashboard fits the window without a scroll bar (pictures scale with the window; library counts in the banner).
- Dashboard tiles always open their own tab: *Open RVE designer* after *Open wall buckling* opens the design tab.
- New look for every workspace: coloured page headers and navigation badges, styled group boxes, tables with
  alternating rows, rounded tabs, focus highlights, slim scroll bars, coloured result banners.

**Faster**
- One RVE solve 1.6 → 1.0 ms (dense KKT solve for small RVEs; exact Gauss integration of S_eq instead of adaptive
  quadrature, identical to round-off).
- Database batches 4.5 → 2.9 ms per case (libraries and geometry plug-ins cached per process; chunked process
  pool, results in input order).
- Built-in shell FE 4–7 × faster (symmetric-mode sparse LU): PBC at 2.5 mm 29 s → 4.5 s, DBC at 5 mm 6.5 s → 2.8 s;
  results identical to round-off.

**Documentation and release**
- New detailed [User Guide](docs/USER_GUIDE.md) and [Theory and Modelling Guidance](docs/THEORY.md) (also as PDF),
  rewritten README, API reference and developer guide.
- Open-source licence, `CITATION.cff`, `.zenodo.json`, continuous integration (tests on Linux and Windows, Python
  3.10 and 3.12, lint) and a release workflow that builds the wheel and sdist for each tagged version and is
  archived by Zenodo with a DOI.

---

# Development history

## 1.15.1 — twist-release factor and bend–twist limit without a reference configuration

- The twist release is written in absolute form, **η = (2.0 + 2.15 ψ²) λ^0.73** with the twist length
  λ = √(D66/A44)/L (same fit to the 320 walls; no reference wall). `laminate.TWIST_LENGTH_REF` is gone;
  `laminate.TWIST_LENGTH_RANGE` = (0.010, 0.100) is the calibrated range.
- The cylindrical-wall limit uses the **bend–twist index Π = η ψ²**, the fraction by which the cylindrical wall
  over-estimates the wall bending stiffness (the in-plane properties are over-predicted by about 100 Π %):
  reliable for Π ≤ 0.018 (≈2 %), use with care for Π ≤ 0.045 (≈5 %). Replaces ψ_e. Wall reports carry
  `twist_length`, `twist_release` and `bend_twist_index`; database rows gain `BendTwistIndexMax`.

## 1.15.0 — twist release calibrated on 320 walls; mesh-convergence verdict

- **Twist-corrected wall recalibrated** on 320 random symmetric walls (4–8 plies; glass, carbon and boron;
  12 ply angles; t/L 0.031–0.124) with mesh-converged periodic shell FE (2.5/1.25/0.625 mm):
  η = (0.150 + 0.162 ψ²)(λ/λ_ref)^0.73 with the **twist length** λ = √(D66/A44)/L (λ_ref = 0.0291, the
  configuration-1 walls). New helpers `laminate.twist_length`, `laminate.twist_release`. Within 2 % of the FE in
  every in-plane constant up to |ψ| = 0.85 for walls with β ≤ 0.10 (five-fold cross-validated; 0.7 % on 48
  independent walls).
- **Reliability limits updated:** cylindrical wall ψ_e = |ψ|(λ/λ_ref)^0.37 ≤ 0.30 (≈2 %) and ≤ 0.45 (≈5 %);
  twist-corrected wall reliable for |ψ| ≤ 0.85 within 0.010 ≤ λ ≤ 0.100; shear-flexible walls with β > 0.10 are
  rated *use with care* (shear model ≈5 %), β > 0.5 outside the tested range. Wall reports carry the twist length.
- **Mesh-convergence verdict** for every FE term (`fe_workflow.mesh_verdict`): converged / refine (GCI) / refine
  (not asymptotic) / converged or refine (oscillating) / refine (diverging), from the observed order, GCI_fine,
  Roache's asymptotic ratio GCI32/(r^p GCI21) and a 1 % tolerance; changes below tol/100 count as converged.
  `richardson` also returns `gci_coarse_pct`, `asymptotic_ratio` and `max_change_pct`. Shown in the HTML report
  (columns *asymptotic ratio*, *mesh verdict* and an explanation) and summarised by `fe check|report` with the
  terms that need a finer mesh.
- Tests: verdict cases (asymptotic, GCI too large, non-asymptotic, oscillating), twist length, updated ratings.

## 1.14.1 — interface: table editor window, illustrated dashboard, no clipped text

- **Geometry tables in their own window.** The Node table, Wall table and Explicit periodic pairs no longer share
  the space under the sketch. A bar below the drawing shows "5 nodes · 3 walls · …" with **Edit nodes…**,
  **Edit walls…**, **Edit periodic pairs…** and **Wall editor (drawing)…**. The editor is a resizable, non-modal
  window (it can stay open beside the sketch; edits update the drawing immediately), with *Refresh sketch from
  tables* and *Close*. The sketch now uses the full height of the page.
- **Dashboard redesigned.** The "Paper benchmark – 14 configurations" card is gone (the benchmark stays under
  *Validation & Benchmarks*). New banner with quick-start buttons, coloured statistic cards and six illustrated
  tiles that open the matching workspace: periodic 3D RVE, laminate ply stack, FE deformed shape, design
  database, Pareto front (multi-objective optimisation) and wall buckling mode. The pictures are generated from
  the software itself by `tools/make_dashboard_art.py` and shipped in `rvestudio/gui/assets`.
- **Text is never cut off.** Every button is at least as wide as its text, tab titles are not elided (scroll
  arrows appear instead), and every workspace scrolls in a small window instead of squeezing its controls.
- Checkable buttons (sketch modes) show their state clearly; window title shows v1.14.
- Tests: dashboard tiles and assets, no benchmark card, button widths, table editor tabs and summary.

## 1.14.0 — reliability check, twist-corrected wall, load-recovery fix (Paper A, Sections 5–7)

- **Reliability check for every design** (`rvestudio.reliability`). Per wall laminate: bend–twist ratio
  ψ = D16/√(D11 D66), effective ratio ψ_e = |ψ|·((t/L)/0.0625)^¼ (shortest wall), membrane–bending measure
  χ_B = B11²/(A11 D11) and wall shear parameter β = 12 D_b/(L² S_eq). Ratings *reliable* / *use with care* /
  *outside the verified range* for the in-plane (E1, E2, G12, ν12, ν13) and out-of-plane (E3, G13, G23, ν31)
  properties, with the reason. Limits from a 58-wall periodic shell-FE study: cylindrical wall ψ_e ≤ 0.30 (≈2 %),
  ≤ 0.40 (≈5 %); twist-corrected wall |ψ| ≤ 0.50 (≈2 %), ≤ 0.62 (≈4 %) for 0.031 ≤ t/L ≤ 0.125; χ_B ≤ 0.02 / 0.05;
  Euler–Bernoulli walls β ≤ 0.02 / 0.05.
  Shown after every solve in the GUI (status line and tooltip), in `summary(result)['reliability']`, printed by
  `rvestudio.cli single`, and written to every database row (PsiMax, ChiBMax, BetaMax, ReliabilityInPlane,
  ReliabilityOutOfPlane).
- **Twist-corrected wall bending option** (GUI *Bending option → Twist-corrected*; model
  `wall_element_type='anisotropicStrip'`, `strip_bending_mode='twistRelease'`): D_b = D11 (1 − η ψ²) with the
  calibrated twist release η = (0.14 + 0.38 ψ²)·√((t/L)/0.0625), capped at 1. Also reported by the FE
  verification as the analytical option "twist-corrected". The default remains the published cylindrical wall.
- **Fix:** wall load recovery — the ε33-conjugate identity of a wall carrying (almost) no N33 compared two
  round-off numbers and could report spurious 1e-2 errors; it is now scaled by the wall's largest resultant.
- Plug-in example `plugins/reentrant_geometry.py` (hexagonal cell with free cell angle; re-entrant for θ < 0).
- Tests: reliability ratings, twist-corrected stiffness, slenderness scaling, database columns (179 tests).

## 1.13.0 — FE verification by strain energy: full C, Q matrix and engineering constants; no published FE

- **PBC = strain-energy homogenisation of the full stiffness.** The top and bottom faces are now periodic through the
  height (U3⁺ − U3⁻ = ε33 h) and all six unit strains are applied, so the PBC model gives the **full 6×6 C**
  (21 terms) instead of six 3D terms. Abaqus: 21 load cases per design and mesh (six single strains, 15 pairs);
  Cᵢᵢ = 2Wᵢ/(V e²), Cᵢⱼ = (Wᵢⱼ − Wᵢ − Wⱼ)/(V e²) from ALLSE; macro reactions only as a cross-check. Built-in
  check: Cᵢⱼ = uᵢᵀ K uⱼ / V, the same quantity (tested to 1e-9 against the 21-state formula), under a second per
  design.
- From C: the **Q matrix** (σ33 = 0, Schur complement; 15 terms in the paper's notation Q11, Q22, Q12, Q66,
  Q44 = 13, Q55 = 23 and the couplings) and the **engineering constants** E1, E2, E3, G12, G13, G23, ν12, ν13,
  ν23, ν21, ν31, ν32 (from S = C⁻¹). The analytical model gives the same three from its Cᴴ.
- **No published FE in FE verification.** Published FE values, columns, findings and the "vs published" charts are
  gone from the page, the report, the workbook, the CLI and the Abaqus summaries (DBC and PBC). Every FE value is
  compared with the analytical model only. The DBC model stays (shell terms with U3 free, 3D terms with U3 = 0). The
  Validation & Benchmarks page (paper comparison) is unchanged.
- GUI: Results → *Stiffness matrices* (C 6×6 and Q 5×5 as analytical | FE | difference), *C, Q and engineering
  constants*, *DBC terms*; deformed shapes of the six PBC unit strains. Report: the matrices side by side, term
  tables, overview heat maps for Q and the engineering constants; workbook sheets *C matrix*, *Q matrix*,
  *Engineering constants*, *DBC shell terms*, *DBC 3D terms*. CLI: `fe check|report` print the matrices.
- Result rows carry a **unit** (MPa, or - for Poisson ratios) and a **negligible** flag: matrix terms zero to
  round-off in both models (e.g. C1123 of a 0/90 stack) get no % difference and are listed under the tables.
  C terms are named by tensor indices (C1133, C2313, …).
- Richardson convergence treats changes below 1e-9 (round-off of the 21-state subtraction) as converged.
- Study folders from before 1.13 are refused (new PBC decks, new result format); prepare a new folder.
- Tests: full C incl. the E33 row, mutual energies = 21-state identification, Q and engineering constants from C,
  21-case Abaqus decks, end-to-end with the fake Abaqus (C, Q and engineering constants converge to the analytical
  values), no published values in rows, report or workbook, refusal of old studies.

## 1.12.1 — wall direction (Node I → Node J, local 1̄ axis); wall editor

- **The local 1̄ axis is now explained where walls are defined.** Every wall runs from Node I (start) to Node J
  (end); that order sets its positive local 1̄ axis (3̄ = core height, 2̄ = wall normal), and ply angles are
  measured from 1̄. The Wall table shows this above the table with the default honeycomb orientation
  (W1 = N1 → N5, W2 = N2 → N5, W3 = N4 → N5, every wall into the junction) and has header tooltips; columns are
  now *Node I (start)*, *Node J (end)* and *Angle of 1̄*.
- Why it matters (checked): reversing a wall leaves 0/90 symmetric stacks unchanged but changes off-axis ones by
  up to 40 % of the largest C^H entry (configurations 4–7, 11, 12, 14). The default orientation reproduces all 84
  published values within 0.3 %; with W1 reversed (N5 → N1) the largest difference is 14.6 % (6 of 84 entries
  beyond 1 %).
- **Sketch:** an arrowhead on every wall points from Node I to Node J (toggle *Wall direction I → J*); hovering a
  wall shows its nodes. **Reverse I ↔ J** swaps the nodes of the selected walls.
- **Wall editor…** (Wall table, sketch tools): a large window with the walls in a full-size table beside a live
  drawing of every wall's I → J arrow and normal 2̄, neighbouring cells and the unit cell; add, delete, reverse;
  invalid rows are marked and block *Apply to Wall table*.
- **More room for the tables:** the drawing and the Node/Wall/periodic tables are now split by a draggable
  divider (the tables were limited to 245 px).
- New module `rvestudio.wall_axes` (convention text, wall frames, `reverse_walls`, `draw_wall_axes`); Help tab
  *Wall direction*; Theory §5 states the wall-local frame.
- Tests: wall frames, default orientation, reversal, the effect of reversal on the benchmark, drawing; GUI smoke
  tests of the editor and the Wall table actions (with Qt).

## 1.12.0 — user-defined FE RVE; parameter diagrams in the Geometry Library

- **Define FE RVE…** (FE Verification page). The selected design is shown replicated in the FE frame (x1 along
  a1) with the rectangular periodic cells, the cell centres and the automatic RVE. Choose **Automatic**
  (published rule, unchanged default), a **Block** of n1 × n2 rectangular periodic cells centred on a complete
  cell, or **drag a rectangle** over the replicated design. A drawn rectangle is rounded to whole rectangular
  periodic cells (always exactly periodic, a whole number of primitive cells, same material per area) and is
  **centred on the geometry wherever possible**: the nearest cell centre, junction or wall mid-point that gives a
  mirror-symmetric block without slivers (triangular networks centre on a junction). Untick *Centre on the
  geometry* to keep the drawn centre; a warning then says when the faces cut the walls at unmatched points.
- The choice is stored on the design (`RVEModel.fe_rve`, saved in `.rvejson` as `"fe_rve"`), so **every** DBC and
  PBC model of that design — built-in check, Abaqus decks, report, CLI — uses it. Each listed design keeps its own
  choice; designs that differ only in their FE RVE get separate entries in a study.
- Built-in check, configuration 1, 10 mm: the PBC C^H is the same (within 0.01 MPa) for the automatic RVE, a
  3 × 2 block and drawn 2 × 2 and off-centre 2 × 1 rectangles, as periodicity requires. DBC terms depend on where
  the faces cut the walls: a geometry-centred drawn RVE reproduces the automatic one exactly, an off-centre one
  gave Q11 = 13 MPa instead of 324 MPa — hence centring by default.
- API: `fe_rve.select_fe_rve(model, window=(x1, x2), cells=(n1, n2), use_model_spec=...)`,
  `normalise_fe_rve_spec`, `fe_rve_layout`, `centred_window`, `window_from_drag`, `spec_rectangle`;
  `fe_plots.plot_rve_layout`. CLI: `rvestudio fe rve --cells N1 N2 [--centre X1 X2 [--free-centre]] [--auto]
  --save design.rvejson`.
- **Fix (all meshers).** Wall and height divisions used `ceil(length / mesh)`, so a 50 mm wall that came out of a
  cut as 50.00000000000003 mm got 6 elements instead of 5 and its periodic copy on the opposite face had
  unmatched nodes (PBC audit failure). Divisions are now robust to round-off (`abaqus._divisions`).
- **Geometry Library: parameter diagrams.** Selecting a family draws it at its reference values (middle of each
  range): each parameter's walls in its own colour with a dimension line, wall numbers W1, W2… (as in the
  Geometry Inspector), neighbouring cells, the unit cell centred on the geometry and a1, a2. Built-in families
  use the Inspector's parameter→wall evidence (R drawn as Ld = R·Lh); wizard W/H families show W and H as the
  drawing width and height; other plug-ins show the geometry with their parameter values. **Parameter guide**
  now opens these diagrams for all built-in families and for the current sketch. Module `geometry_diagram`.
- Tests: user-defined RVE (spec, save/load, every family, centring, PBC invariance, Abaqus decks, study, CLI),
  mesh divisions, parameter diagrams; GUI smoke tests of the dialog, the page and the library diagrams (with Qt).

## 1.11.1 — Abaqus wall orientation fix

- **Fix (Abaqus decks).** Each wall's `*ORIENTATION` now has the data line `3, 0.`, so Abaqus projects the
  material axes about the wall normal (local 3). Without it Abaqus projected about local 1, which put the ply
  1-axis along the core height: the wall's axial stiffness became A22 instead of A11, and for configuration 1
  ([45/-45/0/-45/45], A22/A11 = 0.714) Abaqus gave Q11, Q22, Q12 and C1111, C2222, C1122 about 28% low while Q66,
  Q44 and Q55 were almost unchanged. The built-in solver and the analytical model were not affected.
- `fe run` and `fe report` refuse study folders whose decks lack the line; prepare a new study folder.
- Wall-in-face comments in PBC decks now sit before `*ORIENTATION` instead of inside it.
- GUI: *Anisotropic strip* added to the Bending option list.

## 1.11.0 — FE verification in one place; top/bottom conditions per term group

- **Boundary conditions.** The DBC model keeps the six thesis Table 3.4 cases (top/bottom U3 free) for the
  shell terms Q11, Q22, Q12, Q66, Q44, Q55 (σ33 = 0) and adds four cases with the same side faces and
  **U3 = 0 on the top and bottom faces** for the 3D terms C1111, C2222, C1122, C1212 (ε33 = 0); C1313 and
  C2323 come from the Q44/Q55 cases. The PBC model drops the E33 state, so **U3⁺ − U3⁻ = 0** between top and
  bottom in every state, and gives C1111, C2222, C1122, C1212, C1313, C2323 directly (no σ33 condensation).
  Configuration 1, built-in, 10 mm: Q11 = 324.4 MPa (published FE 322.6), C1111 = 390.6 MPa (analytical 391.0).
- New Abaqus PBC method `direct_reaction5` (`reconstruct_c_direct5`); DBC: `DBC_CASES`, `DBC_C_CASES`,
  `reconstruct_dbc_c`, `reconstruct_dbc_terms`. The built-in DBC check reuses factorisations across cases.
- **One FE page** in the GUI (*FE Verification*) replaces *Abaqus FE (single run)* and *FE Verification
  Study*: designs to add (current design, files, benchmark configurations), meshes, models, study folder;
  **Check now** (built-in) and **Run with Abaqus** (prepare, run, collect; resumable); tabs Results (shell-term
  and 3D-term tables), Charts, FE RVE, Deformed shapes (built-in DBC/PBC, Abaqus DBC), Abaqus jobs, Log.
- **One command group**: `rvestudio fe check | prepare | run | status | report | scripts | rve | doctor`. Removed
  from the CLI: `fe-dbc`, `fe-campaign`, `fe-study`, `fe-doctor`, `verify`, `fe-rve`. Removed from the GUI:
  other periodic domains, 21-case energy recovery, structural load cases, symmetric examples, the campaign.
  Their Python functions remain.
- **One result format**: long rows (project, solver, model, mesh, group, term, FE, analytical, published FE,
  differences) for the built-in check and Abaqus alike; `fe_workflow.load_results`, `verification.summary_table`
  and `findings` feed the GUI, the report and the CLI. One report (`FE_STUDY_REPORT.html`: boundary conditions,
  overview heat maps, every design, convergence, PBC audit, problems) and one workbook; the Markdown report and
  separate PNG figures are gone. Built-in displacement fields are stored per design for the deformed shapes.
- **Wall beam option.** Walls bend as **Euler–Bernoulli** beams by default, as in the published model
  (`RVEModel.wall_shear = "euler-bernoulli"`): the benchmark now reproduces all 84 published values within
  0.3% (was 70/84 at ±1%). **Timoshenko** (`"timoshenko"`, the behaviour before 1.11) adds wall shear
  flexibility; it changes Q66 by 1.4–14.7% and the other terms by < 0.2%, and is closer to shell FE (mean Q66
  error against the published FE 8.6% vs 11.4%). GUI: *Wall beam* in the *Wall model* box of the Single RVE Designer's Solver tab, on
  *Validation & Benchmarks* (with published-FE columns) and on *Generate DB*; CLI: `benchmark --wall-beam`,
  `generate-db --wall-beam`. Saved in project files (files without it load as Timoshenko); `WallShear`
  column in parametric databases. FE verification
  tables and findings show both beams.
- Study folders prepared before 1.11 are refused by `run` and `report` (no U3-restrained cases; six-state PBC).
- Docs: one FE guide, `docs/FE_VERIFICATION.md`, a new FE Verification section in `docs/USER_GUIDE.md` and an FE tab in the in-app Help; older FE documents in `docs/archive/`. One Windows script,
  `RUN_FE_VERIFICATION_WINDOWS.bat`, replaces five.
- Tests: ten-case DBC, U3-restrained decks without gauge, five-state PBC, end-to-end Abaqus study with a working
  fake PBC extractor, failure handling, the `fe` CLI, GUI smoke tests of the new page (with Qt).

## 1.10.0 — automatic FE RVE for every core, verification without Abaqus, reports

- **Automatic FE RVE everywhere.** `rvestudio.fe_rve.select_fe_rve` finds every cell type of a wall network,
  centres the block on the largest complete cell, tries 1×1 … 3×3 rectangular cells, skips blocks that leave
  slivers (< 5 % of the median wall) and reports cell types, mirror symmetry and the shortest wall piece. It
  reproduces the published honeycomb block exactly and gives 4-cell blocks for the square, triangular and mixed
  families (checked on random designs of every family). Overrides: `cells=(n1, n2)`, `centre_cell=k`.
- The periodic default (`PeriodicHomogenizationSettings.domain_mode`, `fe-campaign prepare`, the single-run
  GUI) is now this RVE (`published_rve`); lattices without a rectangular cell fall back to exact lattice copies
  with a stated reason. Older window modes remain selectable.
- **Built-in 3D-periodic check** `shell_fe.run_pbc_shell_check`: the full 6×6 C^H on the same block and mesh as
  the Abaqus PBC decks (master–slave elimination, one factorisation, six strain states). For configuration 1 it
  reproduces the analytical C^H to 0.01 % (Q66 0.5 % at 5 mm).
- **One-call verification** `rvestudio.verification.verify_model` and `rvestudio verify`: analytical (three wall
  options) vs built-in periodic C^H vs built-in finite-height DBC vs published FE, with findings in plain
  language. Benchmark configurations are recognised by content, so published FE values appear even for an
  unnamed model loaded in the designer.
- **Reports** (`rvestudio.fe_report`): self-contained HTML (figures embedded) and Excel workbooks for a single
  verification and for FE studies (overview heat map, colour-coded tables, parity/convergence, every FE RVE,
  built-in checks, Abaqus deformed shapes). `fe-study collect` and `fe-study check` write both.
- The DBC ODB extractor exports `nodal_displacements.csv` (never fails the energy extraction).
- `fe-study check` adds the periodic C^H check (`--no-pbc` to skip); `rvestudio fe-rve` prints/plots the FE RVE.
- **GUI.** Abaqus FE page: *Preview FE RVE*, *Verify without Abaqus*, new *FE RVE*, *Verification* (colour-coded)
  and *Deformed shapes* tabs, HTML/Excel export, automatic RVE as the periodic default. FE study page: *FE RVEs*
  tab, *Open HTML report*, *Open Excel*, colour-coded results, built-in check with the periodic tensor.
- `rvestudio.fe_plots`: shared drawings (RVE plan, boundary conditions, deformed shapes, comparisons).
- `rvestudio verify` prints an aligned table; every table, chart and report shows round-off differences as
  +0.00 (never −0.00).
- Tests: automatic RVE for every family, oblique-lattice fallback, periodic check vs analytical C^H, verification,
  reports, CLI, HTML/Excel in the fake-Abaqus study, GUI smoke test of the new page actions (with Qt).

## 1.9.1 — published FE RVE for the finite-height DBC check

- **Fixed: the DBC decks used the wrong RVE.** Up to 1.9.0 the finite-height check clipped a one-primitive-cell
  window (86.6 × 75 mm for the 50 mm honeycomb) at an arbitrary phase. The thesis (Fig. 3.6b, Tables 3.3–3.4) and
  the paper (Fig. 8b, Table 3) use a 173.2 × 150 × 150 mm block: four primitive cells, one complete hexagon in the
  middle, x1 faces through junction columns and x2 faces at mid-length of the straight walls. New module
  `rvestudio.fe_rve` builds exactly this block (p = 2|a1|, q = |a1 + 2a2|) and, for other wall networks, the
  smallest 1×1/2×1/1×2/2×2 block of rectangular periodic cells centred on a cell that holds that cell whole.
- Walls lying in a face are shared with the next block and get half of the laminate (`face_wall_section`:
  `half-thickness` default, `half-stiffness`, `full`). The published values require it: with half face walls the
  cross-ply cores give Q44 = Q55 exactly as published; full face walls overshoot Q55 by about 10 %.
- The Table 3.4 boundary sets, prescribed displacements and energy-to-Q formulae are unchanged (they matched the
  thesis and paper); V = p·q·h now uses the published block.
- Periodic check: new default window `published-rve` — the same block, with both periodic images of a face wall
  at half moduli so that they act as one wall. Both FE checks now run on the paper's RVE.
- New `rvestudio.shell_fe`: a built-in four-node laminated shell (incompatible-mode membrane, MITC4 bending and
  shear, full CLT A/B/D) that solves the six DBC cases without Abaqus. `fe-study check` and the
  **Built-in check (no Abaqus)** button write `results/shell_check.csv`, which the report includes.
- Studies prepared with 1.9.0 are refused for extension (their DBC decks used the old window); prepare a new one.
- DBC run folders now include `thesis_dbc_rve.svg` (plan view with face walls marked) and a wall table with
  lengths, faces and section factors; decks state the RVE, its cell count and the face-wall treatment.
- `fe-study run` and `collect` (and the Windows study script) now refuse a study whose DBC decks were prepared
  with 1.9.0, so stale one-cell decks cannot be solved or reported by mistake.
- The rigid-body audit of each DBC deck uses an economy SVD and a null-space test for the gauge search:
  identical decks and gauges, about 7x faster (504 decks for the 14 configurations in about a minute).

## 1.9.0 — unified FE verification study and interface

- New `rvestudio.fe_workflow` and `fe-study prepare|status|run|collect|scripts` CLI: any number of projects
  (plus published configurations 1–14) × meshes × {finite-height DBC, periodic PBC} in one study folder.
- Resumable job-level execution: one failed job no longer stops a study; finished jobs are skipped; SOLVED jobs
  are only re-extracted. Job state is read from disk, so studies can be run elsewhere and collected later.
- `run_all_jobs.bat` / `.sh` solve and extract every job without Python (clusters, licensed workstations).
- DBC results are now compared with the analytical model (previously only with thesis FE values), with
  Richardson extrapolation and GCI, a paper-ready CSV/LaTeX table, parity and convergence figures, and a report.
- `fe-doctor` checks that the Abaqus command exists and answers before a long run.
- GUI: new **FE Verification Study** workspace (four guided steps, live job table, stop/resume, results tab).
  The previous Abaqus page remains as **Abaqus FE (single run)**.
- GUI navigation now uses page names; fixes the Single RVE Designer Help button opening the Abaqus page.
- Tests: fake-Abaqus end-to-end test of prepare → run → collect, Richardson checks, CLI and GUI smoke tests.
- No Abaqus solutions are included or claimed in this release.

## 1.8.2 — finite-height top/bottom boundary-condition clarification

- Distinguished physical finite-height DBC validation from fully periodic 3D material-surrogate CH.
- Added the `fe-dbc prepare` and `fe-dbc execute` CLI workflows for published Table 3.4 shell-Q calculations.
- GUI starts in the source-matched finite-height DBC option; periodic mode carries a prominent warning.
- Added explicit height-face BC flags to FE metadata, documentation, and regression tests.
- No Abaqus FE solves are claimed for this release.

## 1.8.1 — Direct 6-state PBC stiffness identification
- Six pure engineering-strain Abaqus cases directly recover full 6x6 CH from macro-control reactions; no intermediate moduli.
- Global force-sign and energy-work checks; raw reciprocity audit; complete coupling terms retained.
- Full 21-case energy method retained as selectable independent cross-check.
- GUI selection, CLI `--identification`, validation campaign, CSV/JSON and tests updated.
- No Abaqus FE solution performed by this release process.

# RVE Studio — Release notes

## 1.0.0 (28 October 2026)

- Renamed the distributable and application to **RVE Studio**.
- Updated the primary research reference to Sriharan et al. (2024), *Engineering Structures* 299, 117106, DOI 10.1016/j.engstruct.2023.117106.
- Aligned the reduced shell transverse shear output with paper Eq. (5): **Q44 = C1313**, **Q55 = C2323**. The internal 3D engineering Voigt index order is unchanged.
- Rebased the 14-configuration regression on **the published proposed-model values** instead of another comparison dataset; users can adjust the acceptance tolerance and inspect or export all 84 differences.
- Added a GUI paper-reference link, clear shear labels, tolerance control, outlier filtering and CSV export.
- Added the stable `rvestudio.api` interface (`solve`, `solve_example`, `solve_project`, `summary`, `benchmark`) and CLI JSON/CSV export.
- Retained the existing geometry editor, material and laminate workflows, batch solver, postprocessing, load recovery, and screening-level wall-buckling functionality.
- Existing saved projects with the previous solver identifier continue to load. They are mapped to the same cylindrical strain-energy wall element.

### Benchmark status

At the default ±1% paper-comparison tolerance, 70/84 entries agree. The 14 Q66 entries differ by up to approximately 14.5363%; this is displayed as **REVIEW**, not hidden or relabeled as a successful paper reproduction. This comparison does not establish FE-level accuracy.

### Scope

The revised publication citation and transverse-shear mapping **do not** imply a new independent derivation or FE validation of the existing numerical kernel. Local-wall buckling remains a linear, simply-supported plate screening model.
- Made wheel installs self-contained: default data and plugin templates are bundled, with editable library files initialised under the user's home directory.
