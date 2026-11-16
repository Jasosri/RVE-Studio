# Finite-height cellular-core FE: RVE and boundary conditions (v1.9.1)

## The RVE used for the finite-height (DBC) check

The thesis (Sriharan 2024, Fig. 3.6 and Tables 3.3–3.4) and the paper (Sriharan et al., Eng. Struct. 299 (2024)
117106, Fig. 8 and Table 3) use, for the honeycomb, a rectangular block of the core:

| Item | Value (l1 = l2 = l3 = l4 = 50 mm, h = 150 mm) |
|---|---|
| p (x1, zig-zag direction, along a1) | 2·\|a1\| = 173.21 mm |
| q (x2, along the straight walls) | \|a1 + 2·a2\| = 150.00 mm |
| h (x3, core height) | 150 mm |
| Primitive cells in the block | 4 (one complete hexagon in the middle) |
| x1 = 0 and x1 = p faces | pass through junction columns; four 25 mm straight-wall halves lie in these faces |
| x2 = 0 and x2 = q faces | cut the straight walls at mid-length (25 mm stubs at x1 = 0, p/2, p) |
| Walls | 8 inclined (50 mm), 2 straight (50 mm), 2 straight halves inside, 4 straight halves in the faces |
| Wall length | 650 mm drawn; 600 mm = 4 × 150 mm once face walls count half |

RVE Studio up to 1.9.0 used a one-cell window (86.6 × 75 mm) at an arbitrary phase. That window has the same
aspect ratio, so the energy-to-Q formulae gave the same numbers for the same energies, but the boundary cuts the
walls in different places and the computed energies are not those of the published model. **Use 1.9.1 or later.**

`rvestudio.fe_rve.build_published_fe_rve` builds the block for any periodic wall network: the smallest of
1×1, 2×1, 1×2 or 2×2 rectangular periodic cells, centred on the centroid of a cell, that holds that cell whole.
For the built-in hexagonal family this is exactly the published block (2×1); square and mixed cores give 2×2,
triangular cores 2×1, all with four primitive cells. Only the honeycomb block has been checked against the
thesis figures.

### Walls lying in a face

A wall that lies in a face belongs half to this block and half to the next. It is given half of the laminate:

- `half-thickness` (default): every ply at half thickness — the natural Abaqus/CAE section for a shared wall;
- `half-stiffness`: full plies, all moduli halved;
- `full`: as drawn (overcounts the material by 8 % for the honeycomb; not recommended).

The published values settle this. Under the Q44 and Q55 boundary conditions the antiplane-shear field of a regular
honeycomb is exactly affine when the face walls carry half, so Q44 = Q55 and both equal the analytical value. The
thesis FE gives Q55/Q44 = 1.000 for the cross-ply cores (configurations 2, 3, 8, 9, 10, 13). With full face walls
Q55 comes out about 10 % higher.

## Boundary conditions (thesis Table 3.4 = paper Table 3)

| Property | x1 = 0 | x1 = p | x2 = 0 | x2 = q | x3 = 0 | x3 = h |
|---|---|---|---|---|---|---|
| Q11 | U1 = 0 | U1 = 1 | U2 = 0 | U2 = 0 | free | free |
| Q22 | U1 = 0 | U1 = 0 | U2 = 0 | U2 = 1 | free | free |
| Q12 | U1 = 0 | U1 = 1 | U2 = 0 | U2 = 1 | free | free |
| Q66 | U1 = U2 = 0 | U1 = 0, U2 = 1 | U1 = 0 | U1 = 0 | free | free |
| Q44 | U1 = U2 = U3 = 0 | U1 = U2 = 0, U3 = 1 | U1 = U2 = 0 | U1 = U2 = 0 | U1 = U2 = 0 | U1 = U2 = 0 |
| Q55 | U1 = U2 = 0 | U1 = U2 = 0 | U1 = U2 = U3 = 0 | U1 = U2 = 0, U3 = 1 | U1 = U2 = 0 | U1 = U2 = 0 |

Displacements in mm; every unlisted translation and every shell rotation is free. For the four in-plane cases
one node also gets U3 = 0, only to remove the rigid translation in x3. The stiffness follows from the total
strain energy W with V = p·q·h:

Q11 = 2W11/(V e11²), Q22 = 2W22/(V e22²), Q12 = (W12 − W11 − W22)/(V e11 e22), Q66 = 2W66/(V γ12²),
Q44 = 2W44/(V γ13²), Q55 = 2W55/(V γ23²), with e11 = γ12 = γ13 = 1/p and e22 = γ23 = 1/q.

Published Q44 is the 1–3 shear and Q55 the 2–3 shear; the internal 3D Voigt order is [11, 22, 33, 23, 13, 12].

## Built-in check before running Abaqus

`rvestudio.shell_fe.run_dbc_shell_check` solves the six cases on the same block with a four-node laminated
shell (Wilson incompatible-mode membrane, MITC4 bending and transverse shear, full CLT A/B/D, transverse shear
5/6·ΣGt). It passes the membrane, bending and coupled patch tests exactly. At the thesis mesh (10 mm) it gives,
against the thesis FE values for the 14 configurations: Q11, Q22, Q12 within about 1 % (2.2 % at most);
Q44 and Q55 within 0.3 % for the cross-ply cores; larger differences for the unbalanced or unsymmetric
±45° cores (configurations 7, 12, 14); and Q66 6–15 % higher. Q66 is bending dominated and converges slowly
(config 1: 2.655, 2.589, 2.558 MPa at 10, 5, 2.5 mm). Treat the check as a cross-check, not as the FE result.

```bash
python -m rvestudio.cli fe-study prepare --root FE_Study --paper-configs 1-14 --meshes 10 5 2.5
python -m rvestudio.cli fe-study check   --root FE_Study --meshes 10 5      # no Abaqus needed
python -m rvestudio.cli fe-study run     --root FE_Study --abaqus-command abaqus
python -m rvestudio.cli fe-study collect --root FE_Study
```

## Periodic surrogate

The periodic check (`--pbc-domain published-rve`, default) uses the same block with x1–x2 periodicity and the
top and bottom faces paired by a3 = (0, 0, h). Both periodic images of a face wall are kept, each with half the
moduli, so together they act as one wall and every face node has its mate. It is a different boundary-value
problem from the DBC check and gives the full 6×6 C^H.
