# FE verification study (v1.10.0)

One folder, one workflow, every Abaqus check needed to verify the analytical model.

| Check | Boundary conditions | Output | Use it for |
|---|---|---|---|
| `dbc` | Finite-height core, thesis Table 3.4. Physical top/bottom faces, not paired. | Q11, Q22, Q12, Q66, Q44 (13), Q55 (23) | Like-for-like verification of the analytical model and the published values |
| `pbc` | Fully periodic 3D surrogate. Top/bottom faces paired through a3 = (0, 0, h). | Full 6×6 C^H, work and reciprocity audits | Full tensor, coupling terms; a different boundary-value problem |

## Command line

```bash
# 1. Prepare (no Abaqus needed). Extend later with more projects or meshes using the same command.
python -m rvestudio.cli fe-study prepare --root FE_Study_PaperA --paper-configs 1-14 \
    --checks dbc pbc --meshes 10 5 2.5 --cpus 4
python -m rvestudio.cli fe-study prepare --root FE_Study_PaperA --project examples/fe_validation/P00050.rvejson

# 2. Check Abaqus, then run. Interrupt at any time; re-running skips finished jobs.
python -m rvestudio.cli fe-doctor --abaqus-command abaqus
python -m rvestudio.cli fe-study run --root FE_Study_PaperA --abaqus-command abaqus --only dbc
python -m rvestudio.cli fe-study run --root FE_Study_PaperA --abaqus-command abaqus
python -m rvestudio.cli fe-study status --root FE_Study_PaperA

# 3. Collect.
python -m rvestudio.cli fe-study collect --root FE_Study_PaperA
```

Every DBC and PBC deck uses the automatic FE RVE (docs/VERIFICATION.md). Before Abaqus,
`fe-study check --root FE_Study_PaperA --meshes 10` runs the built-in shell checks (finite-height
DBC and 3D-periodic C^H) for every project; `collect` and `check` both refresh
`FE_STUDY_REPORT.html` and `FE_STUDY_RESULTS.xlsx`.

On a machine without RVE Studio (for example a cluster), copy the folder and run
`run_all_jobs.sh` or `run_all_jobs.bat`; then copy it back and run `fe-study collect`.
`fe-study scripts --abaqus-command <cmd>` rewrites both scripts for another launcher name.

## GUI

**FE Verification Study** in the navigation bar follows the same four steps: study folder,
what to verify (projects, published configurations, checks, meshes), run (with an Abaqus
check, stop/resume and script export), and a live job table plus collected results. The
job table re-reads the folder every five seconds, so runs started elsewhere show up too.
Double-click a job to open its folder.

## Job states

| State | Meaning |
|---|---|
| PENDING | Deck written, not solved |
| RUNNING | Abaqus lock file present |
| SOLVED | ODB present, results not yet extracted (the next run only extracts) |
| DONE | `case_results.json` written by the ODB extractor with status COMPLETED |
| FAILED | Abaqus or the extractor reported an error; see `abaqus_case.log`, `.sta`, `.msg` |
| NOT PREPARED | The deck could not be generated for this geometry (reason in the note) |

## Results (`results/`)

| File | Content |
|---|---|
| `dbc_results.csv` | Every project × mesh × term: FE, analytical, difference, thesis FE and published proposed-model values where available |
| `dbc_convergence.csv` | Richardson extrapolation (Celik et al. 2008): observed order p, extrapolated value, fine-grid GCI |
| `paper_table_dbc.csv`, `.tex` | One row per RVE: FE reference, analytical value and difference for the six terms |
| `pbc_results.csv`, `pbc_matrices/` | Full FE C^H per project and mesh, Frobenius difference to the analytical tensor, audits |
| `figures/dbc_parity.png`, `dbc_convergence.png` | Parity and mesh-convergence plots |
| `job_status.csv` | Every job and its state |
| `shell_check.csv`, `shell_check_pbc.csv` | Built-in shell checks: DBC shell terms; periodic C^H and terms vs the analytical model |
| `../FE_STUDY_REPORT.html` | Self-contained report: overview heat map, tables, plots, FE RVEs, Abaqus deformed shapes |
| `../FE_STUDY_RESULTS.xlsx` | Every table above as an Excel sheet |

The FE reference used for comparisons is the extrapolated value when three meshes converge
monotonically, otherwise the finest mesh. Nothing is estimated for jobs that have not finished.


## 1.9.1: published RVE and built-in check

Both checks now run on the published RVE (thesis Fig. 3.6b, paper Fig. 8b): for the honeycomb a 173.2 x 150 x
150 mm block of four primitive cells. DBC face walls default to half ply thickness (`--face-walls`); the periodic
window defaults to the same block (`--pbc-domain published-rve`). Studies prepared with 1.9.0 used a one-cell
window for the DBC decks and cannot be extended: prepare a new folder.

`fe-study check` (GUI: **Built-in check (no Abaqus)**) solves the six DBC cases with a built-in four-node
laminated shell on the same block and writes `results/shell_check.csv`; `collect` adds it to the report. Meshes
finer than 2.5 mm are skipped. It is a cross-check of the decks and the energy formulae, not the FE result:
Q66 in particular is mesh- and element-sensitive.
