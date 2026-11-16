# Automatic FE RVE and verification without Abaqus (v1.10)

## One FE RVE for every case

Every FE model in RVE Studio — the finite-height DBC check, the 3D-periodic PBC check, single
runs, FE studies, FE campaigns and the built-in shell checks — now uses the same automatically
selected FE RVE. It follows the rule of the source study's honeycomb RVE (thesis Fig. 3.6b,
paper Fig. 8b) and applies it to any periodic wall network:

1. x1 is aligned with a1 and the smallest rectangular periodic cell is found.
2. The distinct open cells of the network are found and ranked by inscribed circle; the block is
   centred on the largest complete cell.
3. The block is the smallest of 1×1, 2×1, 1×2, 2×2, 3×1, 1×3, 3×2, 2×3, 3×3 rectangular cells
   that holds that cell strictly inside.
4. A block that would leave a wall piece shorter than 5 % of the median wall is skipped.
5. The wall measure is checked: wall length (face walls counted half) = cells × primitive wall
   length, and area = cells × primitive area.

| Core | Block | Cells | Centre | Walls in faces |
|---|---|---|---|---|
| Honeycomb (50 mm walls) | 173.2 × 150 mm (2 × 1) | 4 | complete hexagon | 4 half straight walls in the x1 faces (published RVE) |
| Square | 2 × 2 rectangular cells | 4 | complete square | none |
| Triangular | 2 × 1 | 4 | complete triangle | none |
| Mixed rhombus–triangular | 2 × 2 | 4 | complete rhombus | none |

Walls lying in a face are shared with the next block and carry half the laminate (half ply
thickness in the DBC model; half moduli on both periodic images in the PBC model). Lattices with no
rectangular periodic cell (oblique lattices) cannot have a DBC block; the PBC model then falls
back to exact lattice copies and says so.

```bash
python -m rvestudio.cli fe-rve --project my_core.rvejson --png my_core_rve.png
python -m rvestudio.cli fe-rve --paper-config 1            # prints the block and its checks
python -m rvestudio.cli fe-rve --project my_core.rvejson --cells 3 2   # force a block size
```

In the GUI: **Abaqus FE (single run) → Preview FE RVE**, and the **FE RVEs** tab of the FE
Verification Study.

## Verify without Abaqus

`verify` compares, for any project, in seconds:

| Source | Boundary-value problem | Terms |
|---|---|---|
| Analytical model | periodic wall network, three wall-bending options | full C^H, six shell terms |
| Built-in 3D-periodic shell FE | same block, every boundary node tied to its image (top/bottom too) | full C^H, six shell terms |
| Built-in finite-height DBC shell FE | same block, thesis Table 3.4 = paper Table 3 | six shell terms |
| Published FE (benchmark configurations 1–14) | Abaqus S4, same RVE and DBCs | six shell terms |

The periodic check solves the same boundary-value problem as the analytical model, so it verifies
the implementation term by term (couplings included). The DBC check is the like-for-like
reproduction of the source study's FE. The built-in element is a flat four-node laminated shell
(Wilson incompatible-mode membrane, MITC4 bending and shear, full A/B/D). It is a check, not a
replacement for Abaqus.

```bash
python -m rvestudio.cli verify --project my_core.rvejson --mesh 5 --html report.html --xlsx results.xlsx
python -m rvestudio.cli verify --paper-config 12 --mesh 10
```

In the GUI: **Abaqus FE (single run) → Verify without Abaqus** fills the **Verification** tab
(colour-coded table, findings, difference charts, C^H map), the **Deformed shapes** tab (all six
DBC cases or the six periodic unit strains) and saves an HTML report or Excel workbook.

What the 14 benchmark configurations show (5 mm mesh): for the ten balanced symmetric stacks the
analytical C^H equals the 3D-periodic shell result to 0.05 % in Q11, Q22, Q12, Q44 and Q55 and
within 3.2 % in Q66. The other four stacks differ through wall kinematics (bend–twist coupling in
configuration 7, extension–bending coupling in 12, 13 and 14); for the unsymmetric stacks the
anisotropic-strip option brings Q66 to within 4.3 %.

## Reports

`fe-study collect` and `fe-study check` now also write `FE_STUDY_REPORT.html` (one self-contained
file: overview heat map, tables with green/amber/red differences, parity and convergence plots,
every project's FE RVE, the built-in checks and the Abaqus deformed shapes) and
`FE_STUDY_RESULTS.xlsx` (one sheet per result table). The DBC ODB extractor writes
`nodal_displacements.csv` for each case so that the deformed shapes come from Abaqus itself.
