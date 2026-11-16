# RVE Studio v1.8.1 — Direct FE stiffness from PBC: implementation QA

## Modification
- Default FE identification: `direct_reaction6` instead of the previous `energy21`.
- Six macroscopically independent engineering-strain states, with six conjugate Abaqus generalized reactions per state.
- Full raw 6×6 tensor recovered directly from PBC reactions, without using or fitting engineering moduli.
- Raw reciprocity defect and diagonal energy-work balance retained; one global reaction sign rather than per-column sign flipping.
- Major-symmetric tensor exported as `CH`, `CH_macro_reaction` and `CH_FE_MPa.csv`.
- Compliance and engineering constants derived *after* CH recovery.
- Optional `energy21` reconstructs CH independently from six pure and fifteen coupled strain-energy cases.
- FE campaign collector rejects completed direct runs with >1% work mismatch or >2% raw reciprocity error.

## Local tests (no Abaqus executable)
- 75 pure-Python tests passed; 7 GUI-dependent PySide6 tests skipped.
- End-to-end FE deck/result/CH collection exercised using a synthetic 6×6 positive-definite matrix and **mock ODB outputs** (not FE evidence).
- Real FE decks prepared for a symmetry-checked rectangular SQUARE cell and unrestricted rectangular P00050.
- Both periodic boundary-condition audits passed.
- Each example has 6 `.inp` decks rather than 21.
- No Abaqus solutions have been run, and **no FE CH numerical values have been validated**.

## Required external validation
1. Solve the prepared six-case model using a licensed Abaqus/Standard installation.
2. Extract six reactions and ALLSE for each load case; check work and symmetry audits.
3. Repeat at successively refined meshes and at multiple linear strain amplitudes.
4. Compare the full FE CH against the cylindrical analytical CH and the prescribed target CH.
5. Independently run 21-state energy mode as a high-confidence cross-check on representative designs.

Internal 3D Voigt ordering `[11,22,33,23,13,12]`; shell output convention `Q44=13`, `Q55=23`.
