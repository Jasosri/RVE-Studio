# Abaqus database batch analysis

## Purpose

The Abaqus batch workflow uses the same parametric input database as the analytical batch solver. Each active row is reconstructed into an `RVEModel`, material mapping and wall laminate list before being passed to the selected Abaqus backend.

## Backends

- **Abaqus — Structural FE**: one S4 verification job per RVE for the selected XZ shear, YZ shear or Z-compression load case.
- **Abaqus — Periodic CH21**: 21 periodic macro-strain jobs per RVE, identifying the complete symmetric 6×6 effective stiffness matrix by strain energy.

## Resume logic

Each batch creates a deterministic directory based on the input filename and backend, and a deterministic `Case_<CaseID>` directory inside it. A model/material/laminate/settings signature is stored with each run. Resume is refused if the underlying inputs changed.

For CH21, each of the 21 strain states has its own `case_results.json`. A state marked `COMPLETED` is reused on restart. If the final `homogenized_results.json` is complete, the entire RVE is reused immediately. The master results file is checkpointed after every RVE.

## Input columns

Existing columns such as `GeometryID`, geometry parameter values, `Height`, `Material`, `LayupID`, `T1_mm`, `T2_mm`, ... remain supported. Optional per-wall overrides are:

- `W1Material`, `W2Material`, ...
- `W1LayupID`, `W2LayupID`, ...
- `W1LayupAngles`, `W2LayupAngles`, ... using strings such as `45/-45/0`.

If an override is blank, the global row value is used.

## CH21 master results

The compact master table contains `C11` through the 21 unique upper-triangular `Cij` terms, engineering constants, conditioning/positive-definiteness information, macro-reaction cross-check metrics, runtime and the detailed run-folder path. When analytical comparison is enabled, `Analytical_Cij` and `DifferencePct_Cij` are added.

## Concurrency and licences

`Concurrent RVEs` controls how many independent RVE cases are submitted at the same time. It defaults to 1. `CPUs / job` controls Abaqus CPU allocation within each job. Both settings can affect Abaqus token consumption.

## CLI

Example periodic batch:

```bat
rvestudio abaqus-batch data\RVE_5000_Python_Input.csv abaqus_results.csv --backend ch21 --abaqus-command abq2021 --run-root Abaqus_Batch_Runs --mesh 10 --strain 0.001 --cpus 4 --concurrent 1 --limit 10
```

Run selected CaseIDs by repeating `--case-id`:

```bat
rvestudio abaqus-batch input.csv selected_fe.csv --case-id C000021 --case-id C000114
```

## Rectangular supercell batches (v1.7)

CH21 batches can use the same PBC-domain selection as a single Abaqus run. **Rectangular repeated supercell** is useful when you want clean opposite rectangular boundaries and a direct supercell-size convergence check. The X/Y repeat counts are stored in the batch manifest and in each case signature, so changing them prevents an incompatible resume from reusing stale results.
