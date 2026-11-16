> **Current v1.8.1:** This document describes the original 21-case energy
> identification, available as an option. The **default direct 6-case C^H** method
> is described in [`DIRECT_CH_PBC.md`](DIRECT_CH_PBC.md).

# Abaqus periodic homogenisation — full 21-term C^H

## Purpose

RVE Studio v1.5 adds a general-anisotropic Abaqus/Standard homogenisation path. It is an extension of the Chapter-3 strain-energy-equivalence idea: the microscopic shell RVE and its effective continuum are required to store the same strain energy for prescribed macroscopic strain states.

The thesis introduces the general 3-D effective elasticity tensor through the second derivative of the homogenised strain-energy density and notes that this gives 21 independent elastic coefficients before an orthotropic reduction is imposed. RVE Studio therefore identifies the complete symmetric 6 x 6 stiffness first and does not prescribe any orthotropic zero pattern.

## Strain convention

The Abaqus periodic module uses the standard engineering strain vector

`E = [eps11, eps22, eps33, gamma23, gamma13, gamma12]^T`.

For a periodic mesh-node pair separated by `dx = x+ - x-`, the affine part of the displacement jump is

`u+ - u- = eps_tensor dx`,

where engineering shear is converted to tensor shear by `eps23=gamma23/2`, `eps13=gamma13/2`, and `eps12=gamma12/2`.

This symmetric 3-D strain convention is required for a conventional 6 x 6 continuum stiffness matrix. It generalises the Chapter-3 energy method; it should not be confused with the reduced FSDT transverse-shear kinematics used by the native analytical wall model.

## Periodicity

The cell uses

- `a1 = [a1x,a1y,0]`,
- `a2 = [a2x,a2y,0]`,
- `a3 = [0,0,h]`.

Mesh nodes whose fractional coordinates differ by integer lattice translations are paired automatically. For each pair, three displacement equations impose the affine macroscopic jump and three equations impose periodic shell rotations. One physical shell node fixes translations U1-U3 only. Anchor shell rotations remain free; the periodic translation equations remove the remaining global rigid rotations when the periodic translation vectors span 3-D space.

Six scalar control nodes represent the six engineering-strain components. Their `U1` values are prescribed to `strain * reference_length`. This also permits recovery of the generalized macro reactions for an independent stress-based check.

## 21 energy cases

For a symmetric elastic 6 x 6 stiffness there are 21 independent entries. RVE Studio runs:

1. six pure states: `E11`, `E22`, `E33`, `G23`, `G13`, `G12`;
2. all fifteen pair combinations of those six states.

For a pure state `i` with amplitude `d`,

`Cii = 2 Wi / (V d^2)`.

For a paired state `(i,j)`,

`Cij = (Wij - Wi - Wj) / (V d^2)`.

`Wi` is the Abaqus whole-model elastic strain energy `ALLSE` and `V = |a1 x a2| h` is the envelope RVE volume.

All 21 values are stored. The GUI's coupling threshold is only a significance flag based on

`eta_ij = |Cij| / sqrt(|Cii Cjj|)`.

It never modifies the raw stiffness matrix.

## Six-case macro-reaction cross-check

For each of the six pure cases, Abaqus reactions at all six macro control nodes provide a second identification of a complete stiffness column. The raw 6 x 6 reaction matrix is stored, along with:

- reciprocity/symmetry error `||C-C^T||/||C||`;
- relative difference between the symmetric reaction matrix and the 21-case energy matrix.

The energy method remains the primary result.

## Output

A successful run folder contains:

- `homogenization_snapshot.json` — exact geometry, materials, layups, settings, case definitions and mesh metadata;
- 21 numbered case folders containing the `.inp`, ODB extractor, Abaqus job files and case JSON;
- `homogenized_results.json` — complete result record;
- `CH_energy_MPa.csv` — full 6 x 6 effective stiffness;
- `compliance_1_per_MPa.csv` — inverse/pseudo-inverse compliance;
- `independent_21_terms.csv` — all independent terms and significance flags;
- `case_energies.csv` — the 21 strain energies;
- `homogenization_run.log`.

If `C^H` is singular or nearly singular, the full stiffness is still retained. The software uses a pseudo-inverse for display of compliance and reports the conditioning/positive-definiteness diagnostics rather than discarding the stiffness result.

## FE modelling assumptions

- Abaqus/Standard S4 shell elements.
- Each RVE wall is a vertical shell midsurface spanning the core height.
- Each wall uses its assigned composite Stack ID.
- Wall-local material 0 degrees follows the wall direction; local 2 follows the core-height direction.
- Weighted/multiplicity analytical walls are rejected by this FE path unless first expanded into explicit physical walls.
- The periodic cell is treated as a 3-D prismatic periodic medium with `a3=(0,0,h)` for the purpose of full 3-D constitutive identification.

## Validation status

The pure-Python tests verify the 21-case enumeration, engineering-shear mapping, periodic deck generation, run packaging, and exact reconstruction of a synthetic fully coupled symmetric stiffness matrix. The original 14-case Chapter-3 analytical benchmark remains unchanged.

A licensed Abaqus executable is not available in the package-building environment. Therefore each new Abaqus release should be validated on the target Abaqus version with at least:

1. a homogeneous material cell with a known stiffness matrix;
2. a mesh-refinement study;
3. strain-amplitude linearity checks;
4. energy-vs-macro-reaction agreement;
5. comparison against the thesis configurations for the stiffness terms for which the modelling assumptions coincide.

## v1.6.1 boundary-condition and area audit

The periodic equations retain equal shell rotations on corresponding periodic nodes, but the single anchor fixes **translations U1-U3 only**. Anchor rotations are not fixed. This matches the intent of the Chapter-3 periodic kinematics more closely and avoids artificially suppressing wall/junction rotation.

Energy is normalised by the periodic envelope volume

`V_RVE = |det([a1,a2])| * h`.

In a local frame aligned with `a1`, `p=|a1|`, `q=|det([a1,a2])|/p` and hence `A_RVE=p*q`. For the regular Chapter-3 hexagon this is exactly the appendix formula `V=a*b*h`.

The raw 3-D `C^H` and the equivalent-shell properties are distinct outputs. v1.6.1 adds the general anisotropic plane-stress Schur complement (eliminating epsilon_33 under sigma_33=0) and exports the resulting 5x5 shell matrix plus `Q11,Q22,Q12,Q66,Q44,Q55`. The FE values in thesis Tables 3.5-3.10 used prescribed-displacement boundary conditions (Table 3.4), not periodic PBC; their comparison is therefore labelled as a cross-formulation diagnostic rather than a same-BC regression.


## v1.7 rectangular repeated-supercell mode

The CH21 solver can use either the primitive periodic cell or a repeated rectangular supercell. In rectangular mode RVE Studio searches for integer coefficient vectors `(m,n)` and `(p,q)` such that

`b1 = m a1 + n a2`,  `b2 = p a1 + q a2`

and `b1 . b2 = 0` within the configured numerical tolerance. Because the coefficients are integers, the rectangle is commensurate with the primitive lattice. If `M=[[m,p],[n,q]]`, the minimum rectangle contains `|det(M)|` primitive cells. Additional X/Y repeats multiply this cell count.

The generated RVE area is therefore checked independently as

`A_super = |b1 x b2| = N |a1 x a2|`.

The repeated wall network is phase-shifted, clipped to the outer rectangle, and merged. The phase shift prevents a complete zero-thickness shell wall from being represented simultaneously on both opposite boundaries and therefore double-counted in the strain energy. As a second consistency check, the total explicit wall length in the supercell must equal `N` times the primitive-cell wall length.

Before any `.inp` is written, a strict PBC audit checks:

- complete lower/upper node pairing for the two in-plane supercell directions and the height direction;
- periodic translation rank = 3;
- rigid-body constraint rank = 6;
- no duplicate master/slave pair equations;
- connectivity after periodic ties;
- affine-field residual for the implemented engineering-strain mapping.

A failed audit aborts generation. `periodic_domain.json` and `pbc_audit.json` are saved in every periodic run folder. Rectangular runs additionally save `supercell_nodes.csv` and `supercell_walls.csv` for independent inspection.

## v1.7.4 exact parallelogram lattice supercell

When the objective is to repeat the original RVE without changing its geometry, select **Parallelogram lattice supercell — exact copies**. The primitive wall network is translated by

`x_ij = x + i*a1 + j*a2`,  `i=0..Na-1`, `j=0..Nb-1`.

The analysis vectors are `B1=Na*a1` and `B2=Nb*a2`. Every copied wall therefore has exactly the same vector, length, orientation, laminate assignment, and topology as the corresponding primitive wall. Shared nodes/walls are merged only for FE conformity. No geometric clipping is performed.

The expected checks are:

- primitive-cell count `N = Na*Nb`;
- `A_super = N*A_primitive`;
- total wall measure `L_super = N*L_primitive` for the explicit periodic wall representation;
- complete PBC pairing on `B1`, `B2`, and `a3=(0,0,h)`;
- periodic-translation rank 3/3 and rigid-body rank 6/6;
- no duplicate periodic pair equations and negligible affine-field residual.

The **Rectangular periodic cut** is a different model construction. It forms an orthogonal commensurate window through the repeated infinite lattice and can split wall segments where the rectangle intersects cells. It remains available for studies that specifically require an orthogonal FE window, but it should not be confused with direct repetition of complete primitive parallelogram RVEs.
