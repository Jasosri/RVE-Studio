# Direct C^H from periodic Abaqus FE — v1.8.1

## Why this mode

The RVE Studio FE calculation does **not** infer homogenised stiffness by first
estimating E1, E2, E3, G23, G13, G12 and Poisson ratios. It obtains **the complete
6×6 C^H directly** from the force conjugates to the six imposed macroscopic
engineering-strain DOFs. The pre-existing 21-case energy-based full-C route is
still independently available.

## Kinematics and extraction

The internal vector ordering is
`epsilon_bar=[epsilon11,epsilon22,epsilon33,gamma23,gamma13,gamma12]`.
For corresponding periodic shell mesh nodes x+/x-:

`u(x+) - u(x-) = eps_tensor_bar (x+ - x-)`;
`rotation(x+) - rotation(x-) = 0`.

Engineering shear is converted to tensor strain using `epsilon_ij=gamma_ij/2`.
The six macro control nodes use scalar U1 values `U_i = Lref * epsilon_i`.
The 3D periodic lattice is `a1, a2, [0,0,h]`. Each of six FE cases imposes
`epsilon_j = a` and the other five macro strains zero.

Abaqus extracts the six conjugate reactions R_i^(j) from the macro control
nodes. The raw full 6×6 matrix is, in MPa when N and mm are used:

`Craw_ij = sign * Lref * R_i^(j) / (V_RVE * a)`.

`V_RVE = |det([a1,a2])|*h` for the resolved periodic FE domain.
The global reaction-sign convention is selected ONCE from the work balance
`sign * R_j^(j) * Lref*a = 2*ALLSE_j` for all six cases. Independent
sign changes between columns are rejected. The raw matrix is retained and
its reciprocity defect `||Craw-Craw.T||_F/||Craw||_F` reported.
The major-symmetric output is `C^H=(Craw+Craw.T)/2`; this is reported only
after recording the unsymmetrised matrix and its reciprocity error.

There is **no orthotropic projection** and no deletion of off-diagonal
normal–shear couplings. Output order is `[11,22,33,23,13,12]`; the published
shell names `Q44=13 shear, Q55=23 shear` are mapped separately in the shell
plane-stress output, not interchanged in the internal C matrix.

## Diagnostics and verification

A prepared FE deck must pass the exact-lattice/rectangular domain audit.
After Abaqus extraction, the direct route reports:

- `CH` and `CH_macro_reaction`: primary symmetric full stiffness;
- `CH_macro_reaction_raw`: unsymmetrised direct reaction matrix;
- `direct_reaction_audit`: maximum relative diagonal energy/work error,
  per-case work errors, global sign, and reciprocity error;
- `positive_definite`, condition number, all 21 independent coefficients;
- `compliance` and `engineering_constants`: **derived only after CH**.

The FE campaign collector marks the FE matrix **INVALID**, not verified, if the
maximum energy-work error exceeds 1% or the raw reciprocity error exceeds 2%.
These are initial numerical screening tolerances, not proof of physical accuracy.
Mesh convergence, strain amplitude and independent validation are still required.

## Command-line workflow

```bash
python -m pip install -e .
python -m rvestudio.cli fe-campaign prepare \
  --project examples/fe_validation/P00050.rvejson \
  --root FE_Direct_CH --meshes 10 5 2.5 \
  --pbc-domain rectangular-supercell --supercell-x 1 --supercell-y 1 \
  --identification direct_reaction6
python -m rvestudio.cli fe-campaign run --root FE_Direct_CH --abaqus-command abaqus
python -m rvestudio.cli fe-campaign collect --root FE_Direct_CH
```

For the independent full-energy cross-check run a **different output root** with
`--identification energy21`. For symmetric certified built-in SQUARE/MIXED
0/90 stacks you may use `--pbc-domain symmetric-rectangular`. Off-axis HEX
optimised layouts such as P00050 should use unrestricted `rectangular-supercell`
(or an exact parallelogram lattice supercell) instead.

Results are written only after actual Abaqus ODB extraction. Preparing input
decks never creates a numerical FE stiffness matrix.
