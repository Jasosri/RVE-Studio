# Abaqus FE integration

RVE Studio Python 1.7.4 provides three Abaqus/Standard S4 workflows.

## 1. Structural verification

The original v1.4 structural checks remain available: XZ shear, YZ shear and Z compression. They generate a single displacement-controlled S4 job and store ODB-derived reactions/displacements.

## 2. Periodic homogenisation — full 21-term C^H

The v1.5 workflow identifies the complete symmetric 6 x 6 effective stiffness matrix without imposing orthotropic zeros. It creates six pure and fifteen paired macro-strain jobs, extracts whole-model elastic strain energy `ALLSE`, and reconstructs all 21 independent coefficients.

Periodic equations are based on the RVE lattice vectors `a1`, `a2` and prismatic height vector `a3=(0,0,h)`. Corresponding S4 mesh nodes satisfy the affine displacement jump and periodic shell rotations. Six scalar macro-control nodes encode `[eps11,eps22,eps33,gamma23,gamma13,gamma12]`.

The same wall-specific composite shell sections and local material orientations used by structural verification are retained. A 0-degree ply follows the local wall direction.

See `ABAQUS_HOMOGENIZATION.md` for the equations, 21-case identification, reaction cross-check, output files and validation requirements.

## 3. Thesis Chapter-3 prescribed-displacement validation

Choose **Thesis DBC validation — Table 3.4** to reproduce the six Abaqus FE boundary-condition cases used in Chapter 3 for the effective shell stiffnesses `Q11`, `Q22`, `Q12`, `Q66`, `Q44`, and `Q55`. This is intentionally separate from the fully periodic CH21 formulation.

The geometry is rotated into the Chapter-3 local frame (x1 along a1) and clipped to the published `p x q x h` RVE of thesis Fig. 3.6(b) / paper Fig. 8(b): for the honeycomb p = 2|a1|, q = |a1 + 2a2| (173.2 x 150 mm for 50 mm walls), four primitive cells with one complete hexagon in the middle. Walls lying in a face get half of the laminate. The code verifies that the area is a whole number of primitive cells and that the shared-wall-weighted wall length matches before writing the decks (see FINITE_HEIGHT_RVE_BOUNDARY_CONDITIONS.md). The thesis FE setup used Abaqus 2021, S4 shells and a 10 mm mesh; the applied displacement in Table 3.4 is 1 mm.

The six translational DBCs are:

| Property | x1=0 | x1=p | x2=0 | x2=q | x3=0 | x3=h |
|---|---|---|---|---|---|---|
| Q11 | U1=0 | U1=1 | U2=0 | U2=0 | free | free |
| Q22 | U1=0 | U1=0 | U2=0 | U2=1 | free | free |
| Q12 | U1=0 | U1=1 | U2=0 | U2=1 | free | free |
| Q66 | U1=U2=0 | U1=0, U2=1 | U1=0 | U1=0 | free | free |
| Q44 | U1=U2=U3=0 | U1=U2=0, U3=1 | U1=U2=0 | U1=U2=0 | U1=U2=0 | U1=U2=0 |
| Q55 | U1=U2=0 | U1=U2=0 | U1=U2=U3=0 | U1=U2=0, U3=1 | U1=U2=0 | U1=U2=0 |

Every unlisted translational DOF and every shell rotational DOF remains free. For `Q11/Q22/Q12/Q66`, one physical node additionally receives `U3=0` solely to remove the otherwise free rigid-body translation in the third direction.

The ODB extractor reads whole-model elastic strain energy `ALLSE`, and the six shell-Q terms are reconstructed using the Chapter-3 characteristic dimensions and prescribed displacement.

## Running Abaqus

Open **Abaqus FE**, choose the analysis type, enter/detect the Abaqus command (`abaqus`, `abq2021`, etc.), set mesh size and CPUs, then generate or run. Generation does not require Abaqus; solution and ODB extraction require a licensed local Abaqus installation.
