# FE design & validation (built-in shell FE and Abaqus)

RVE Studio checks the analytical model against shell finite elements on one FE RVE per design, with
one of two solvers:

* **built-in** — a four-node laminated shell (Wilson incompatible-mode membrane, MITC4 bending and
  transverse shear, full A/B/D). Seconds per design, no licence. A cross-check, not a replacement for Abaqus.
* **Abaqus** — S4 composite shells, the same RVE, node sets, boundary conditions and energy formulae.

Both write into one **study folder** and one report.

## The FE RVE

Every FE model uses the block the source study used for the honeycomb (thesis Fig. 3.6b, paper Fig. 8b),
with the same rule applied to any wall network: x1 along a1; the smallest rectangular block of whole
periodic cells (1×1 … 3×3) that holds the largest cell complete in its centre; no wall piece shorter than
5 % of the median wall. Walls lying in a face are shared with the next block and carry half the laminate.
For the 50 mm honeycomb: 173.2 × 150 × 150 mm, four primitive cells. `rvestudio fe rve --paper-config 1`
prints it; the GUI shows it in the **FE RVE** tab.

### Choosing or drawing the FE RVE (v1.12)

The automatic RVE is the default. For any design you can instead choose a **block** of n1 × n2 rectangular
periodic cells (centred on a complete cell) or **draw a rectangle** over the design replicated in the FE frame:
*FE Design & Validation → select a design → Define FE RVE…*, or `rvestudio fe rve --cells N1 N2 [--centre X1 X2]
--save design.rvejson`. The size of a drawn rectangle is rounded to whole rectangular periodic cells, so the RVE
is always exactly periodic and holds a whole number of primitive cells. Its centre is moved to the **centre of the
geometry** wherever possible — the nearest cell centre, junction or wall mid-point that gives a mirror-symmetric
block without slivers — unless you untick *Centre on the geometry* (`--free-centre`).

The choice is saved with the design (`"fe_rve"` in the `.rvejson`) and used by every DBC and PBC model, built-in
and Abaqus. **PBC** terms do not depend on the RVE (configuration 1: identical within 0.01 MPa for the automatic,
3 × 2, drawn 2 × 2 and off-centre 2 × 1 RVEs). **DBC** terms do: they depend on where the faces cut the walls
(an off-centre RVE gave Q11 = 13 MPa against 324 MPa), which is why the centre is placed on the geometry by
default.

## Boundary-value problems

| Model | Side faces | Top and bottom faces | Terms |
|---|---|---|---|
| PBC, strain energy | periodic | periodic through the height: U3⁺ − U3⁻ = ε33 h; U1, U2 and rotations periodic | full 6×6 C (21 terms), Q matrix (σ33 = 0), engineering constants |
| DBC, shell terms | prescribed, thesis Table 3.4 = paper Table 3 | U3 free (in-plane cases; one U3 gauge node); U1 = U2 = 0 (Q44, Q55) | Q11, Q22, Q12, Q66, Q44, Q55 (σ33 = 0) |
| DBC, 3D terms | the same | U3 = 0 on both faces | C1111, C2222, C1122, C1212 (ε33 = 0) |

### PBC: full C, Q and engineering constants from strain energy

Six unit macro strains ε = e·eₖ (Voigt order 11, 22, 33, 23, 13, 12) are applied to the periodic RVE of
volume V = p·q·h. The stiffness follows from the strain energy W:

* **Abaqus** — 21 load cases per mesh: the six single strains and the 15 pairs eᵢ + eⱼ.
  Cᵢᵢ = 2 Wᵢ / (V e²) and Cᵢⱼ = (Wᵢⱼ − Wᵢ − Wⱼ) / (V e²) from the ALLSE of each case, so all 21 independent
  terms come from energies; the macro reactions are kept only as a cross-check (*PBC audit* in the report).
* **Built-in check** — the same six states with one factorisation; Cᵢⱼ = uᵢᵀ K uⱼ / V is the mutual strain
  energy, identical to the 21-state formula (a test checks this to 1e-9).

From C follow the **Q matrix** — C condensed to σ33 = 0 (Schur complement), 5×5 in the order 11, 22, 12, 13, 23
with the paper's names Q11, Q22, Q12, Q66, Q44 (1–3), Q55 (2–3) and the couplings Q16, Q26, Q45, … — and the
**engineering constants** from S = C⁻¹: E1 = 1/S11, E2, E3, G23 = 1/S44, G13 = 1/S55, G12 = 1/S66 and
νij = −Sji/Sii. The analytical model gives the same three from its Cᴴ, so FE and analytical are compared like for
like. C terms are named by tensor indices (C1133, C2313, …) to avoid confusing Voigt C44 (= C2323) with the
paper's Q44 (= C1313). Terms that are zero to round-off in both (e.g. C1123 for a 0/90 stack) are shown as ≈ 0
and get no % difference. The PBC results do not depend on the FE RVE.

### DBC: shell and 3D terms

C1313 and C2323 of the DBC model come from the Q44 and Q55 cases: U3 is prescribed on the side faces
there, so it cannot be restrained on the top and bottom. The DBC model has ten cases per mesh.

The stiffness follows from the strain energy W, V = p·q·h, e11 = γ12 = γ13 = d/p, e22 = γ23 = d/q:

| | U3 free | U3 = 0 |
|---|---|---|
| normal | Q11 = 2W/(V e11²), Q22 = 2W/(V e22²) | C1111, C2222 likewise |
| coupling | Q12 = (W12 − W11 − W22)/(V e11 e22) | C1122 likewise |
| in-plane shear | Q66 = 2W/(V γ12²) | C1212 likewise |
| transverse shear | Q44 = 2W/(V γ13²) (1–3), Q55 = 2W/(V γ23²) (2–3) | C1313 = Q44, C2323 = Q55 |

What is compared — every FE value with the same quantity of the analytical model, nothing else:

* **PBC** — the full C matrix, the Q matrix and the engineering constants (shown as side-by-side matrices
  analytical | FE | difference, and as term tables).
* **DBC** — the shell terms with the analytical Q, the 3D terms with the analytical Cᴴ entries.
* **Both wall beams** — the analytical columns show the design's wall beam and the other one: Euler–Bernoulli
  and Timoshenko (with wall shear flexibility). The shell elements are shear-flexible, so Timoshenko is the
  like-for-like comparison with FE; the findings name the analytical wall option closest to the PBC result.

Restraining U3 matters in the DBC model: for configuration 1 the built-in DBC gives Q11 = 324.4 MPa with U3 free
and C1111 = 390.6 MPa with U3 = 0 (analytical C1111 = 391.0). Q66, Q44 and Q55 hardly change.

## GUI: FE Design & Validation

1. **Design the FE model** — add the current design and/or project files (select one and press
   **Define FE RVE…** to choose or draw its FE RVE; **Automatic** resets it), type benchmark configurations
   (e.g. `1-14`), mesh sizes (`10`, or `10, 5, 2.5` for a Richardson extrapolation), models and a study folder.
2. **Validate now** — built-in shell FE of every design and mesh (meshes finer than 2.5 mm are skipped).
3. **Validate with Abaqus** — writes the decks if needed, solves every pending job one at a time and collects.
   Stop and press it again to resume; finished jobs are kept. *Scripts for another machine* writes
   `run_all_jobs.bat` / `.sh` that need no Python; copy the folder back and press *Reload results*.

The right side shows the selected design: the two term tables (colour-coded differences), charts, the
FE RVE, deformed shapes (built-in DBC and PBC, and Abaqus DBC from its nodal export), the Abaqus jobs and
a log. **Open report** opens `FE_STUDY_REPORT.html`; **Open Excel** `FE_STUDY_RESULTS.xlsx`.

## Command line

```bash
python -m rvestudio.cli fe check   --root MyStudy --project my_core.rvejson --meshes 10
python -m rvestudio.cli fe check   --root PaperA --paper-configs 1-14 --meshes 10 5
python -m rvestudio.cli fe run     --root PaperA --paper-configs 1-14 --meshes 10 5 2.5 --abaqus-command abaqus --cpus 4
python -m rvestudio.cli fe status  --root PaperA
python -m rvestudio.cli fe report  --root PaperA        # collect jobs finished elsewhere, rewrite the report
python -m rvestudio.cli fe prepare --root PaperA ...    # decks only (no Abaqus)
python -m rvestudio.cli fe scripts --root PaperA --abaqus-command "C:\SIMULIA\Commands\abaqus.bat"
python -m rvestudio.cli fe rve     --project my_core.rvejson --png rve.png
python -m rvestudio.cli fe rve     --project my_core.rvejson --cells 3 2 --centre 80 95 --save my_core_3x2.rvejson
python -m rvestudio.cli fe doctor  --abaqus-command abaqus
```

On Windows, `RUN_FE_VERIFICATION_WINDOWS.bat` runs the Paper A study (14 configurations, both models,
10/5/2.5 mm): built-in check, then Abaqus, then the report. Run it from the Abaqus command prompt, or pass
the path to `abaqus.bat` as its argument.

## Study folder

```
fe_study.json                 projects, analytical references, settings, runs
projects/                     the .rvejson of every design
dbc/<design>_h<mesh>/NN_<case>/   Abaqus deck, extractor, case_results.json, nodal_displacements.csv
pbc/<design>_h<mesh>/NN_<state>/  Abaqus deck, extractor, case_results.json
results/builtin_results.csv   built-in rows: project, solver, model, mesh, group, term, FE, analytical, differences
results/abaqus_results.csv    the same for Abaqus
results/convergence.csv       Richardson order, extrapolated value, GCI, asymptotic ratio and mesh verdict per term (three meshes)
results/paper_table.csv/.tex  Abaqus DBC shell terms vs analytical (extrapolated where monotone)
results/pbc_audit.csv         energy-work and reciprocity errors of the Abaqus PBC states
FE_STUDY_REPORT.html          one self-contained file: boundary conditions, overview, every design, convergence
FE_STUDY_RESULTS.xlsx         the same tables, one sheet each
run_all_jobs.bat / .sh        solve and extract every pending job without Python
```

Nothing is estimated: an Abaqus value counts only when its extractor wrote `status=COMPLETED`. Study
folders prepared before 1.11 (no U3-restrained DBC cases, PBC with an E33 state) are refused by `run` and
`report`; prepare a new folder.

## Python

```python
from rvestudio import fe_workflow as fw
from rvestudio.verification import verify_model, VerificationSettings

v = verify_model(model, materials, laminates, VerificationSettings(mesh_size=5.0))   # one design, no files
fw.prepare_study('PaperA', paper_configs='1-14', settings=fw.FEStudySettings(meshes=[10, 5, 2.5]))
fw.check_study('PaperA')                  # built-in
fw.run_study('PaperA', 'abaqus')          # Abaqus, then collect
R = fw.load_results('PaperA')             # rows, convergence, jobs
```

Older FE front ends (single-run Abaqus page, `fe-dbc`, `fe-campaign`, `fe-study`, other periodic domains,
21-case energy recovery, structural load cases) were removed from the GUI and the command line in 1.11;
their Python functions remain (`abaqus_thesis_dbc`, `abaqus_homogenization`, `fe_validation`) and their
documents are in `docs/archive/`.


## Is the mesh converged?

Give three meshes halving each time. Every term then gets the observed order p, the extrapolated value, GCI_fine,
the asymptotic ratio R_as = GCI_32/(r^p GCI_fine) and a verdict: *converged* (GCI_fine <= 1 %, 0.5 <= p <= 4,
|R_as - 1| <= 0.15, or the meshes agree to 0.01 %), *refine (GCI)*, *refine (not asymptotic)*,
*converged/refine (oscillating)* or *refine (diverging)*. Compare the analytical model only with converged terms,
using the extrapolated FE value. PBC terms do not depend on the RVE height, so a thin RVE (e.g. 10 mm) makes
2.5/1.25/0.625 mm meshes affordable for bend-twist coupled walls. Full explanation: [THEORY.md, Section 13](THEORY.md#13-mesh-convergence-is-my-fe-result-converged).
