# RVE Studio — Public Python and command-line API

Research reference: Sriharan et al., *Engineering Structures*, 299 (2024) 117106; [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0141029623015213).

## Imports and calls

```python
from rvestudio.api import solve, solve_example, solve_project, summary, benchmark
from rvestudio.models import Material, Ply, Laminate, Element, RVEModel

result = solve_example(config=1)
public = summary(result, include_ch=True)

# Equivalent explicit solve:
# model, materials, laminates = paper_hexagonal_case(1)
# result = solve(model, materials, laminates)

# A saved project:
# result = solve_project('my_project.rvejson')

report = benchmark(tolerance_pct=1.0)                                   # original model
report = benchmark(wall_shear='timoshenko', bending='twist-corrected')   # any wall model
report = benchmark(bending='twist-corrected', twist_release=(2.0, 2.15, 0.73))
report['detail'].to_csv('benchmark.csv', index=False)
```

`solve(model, materials, laminates)` validates the model, nonempty material map and nonempty laminate sequence; returns `HomogenizationResult` with `CH`, `S`, effective properties, wall data, and physical properties. `solve_example(config)` accepts integer `1..14` for the paper's example wall stacks, and returns the same result type. `solve_project(path)` uses the existing `.rvejson` schema.

`summary(result, include_ch=False)` is a JSON-compatible dictionary with keys `software`, `reference`, `stress_unit`, `voigt_order_3d`, `transverse_shear_order`, `stiffness`, `engineering_constants`, `rve` and, when `include_ch` is true, the matrices `CH`, `Qplane` and `Qshear`.

`benchmark(*, tolerance_pct=1.0, wall_shear='euler-bernoulli', bending='cylindrical', twist_release=None)` solves the 14 benchmark designs with the chosen wall model and returns `detail` (84 rows: Config, Property, Calculated_MPa, Paper_MPa, Error_pct, Within_tolerance) and `summary` as pandas tables, with `wall_shear`, `bending`, `passed`, `total`, `within_tolerance`, `max_abs_paper_pct`, `calculated` and `paper_published`. `bending` is one of `cylindrical`, `twist-corrected`, `anisotropic strip`, `free moment`. With the default wall model all 84 values are within 0.3 %.

`summary(result)['reliability']` gives the rating (`level`, `in_plane`, `out_of_plane`), the wall model, every wall's `psi`, `twist_length`, `twist_release`, `bend_twist_index`, `chi_B`, `beta_max` and notes, and the limits used.

### Wall model of an `RVEModel`

| UI | `wall_element_type` | `strip_bending_mode` |
|---|---|---|
| Cylindrical | `paper` | — |
| Anisotropic strip | `anisotropicStrip` | `cylindrical` |
| Free moment | `anisotropicStrip` | `freeMoment` |
| Twist-corrected | `anisotropicStrip` | `twistRelease` |

`wall_shear='euler-bernoulli' | 'timoshenko'`; `twist_release=(a, b, n)` overrides the calibrated coefficients
`rvestudio.laminate.TWIST_RELEASE_DEFAULT = (2.0, 2.15, 0.73)`.

### Databases

`rvestudio.database.generate_parametric_database(output_file, num_cases=5000, seed=20260906, relative_density_range=(0.03, 0.15), ..., wall_shear='euler-bernoulli', bending=None, twist_release=None)` and `run_batch(input, output, workers=1, case_limit=None, progress=None, buckling_settings=None)`.

### FE validation

`rvestudio.shell_fe.run_pbc_shell_check(model, materials, laminates, mesh_size=10.0)` and `run_dbc_shell_check(...)`; `rvestudio.fe_workflow.richardson(hs, fs)` and `mesh_verdict(rich, f_fine, f_medium, tol_pct=1.0)`.

## Tensor and units

- Full stiffness `CH` and compliance `S`: engineering Voigt `[11,22,33,23,13,12]`.
- In-plane plane-stress `Qplane`: `[11,22,12]`.
- Shell `Qshear`: `[13,23]`; Q44 corresponds to C1313, Q55 to C2323, as in paper Eq. (5).
- Modulus/stress units follow material input (MPa for example library); length follows the model input (mm for examples). Density is kg/m³.
- `Q44` and `Q55` in **legacy CSV output** may have the opposite convention. Rerun the model or explicitly remap the old columns before merging data.

## Command-line interface

```bash
python -m rvestudio.cli --help
python -m rvestudio.cli --version
python -m rvestudio.cli single --config 2 --json --output case2.json
python -m rvestudio.cli solve-project example.rvejson --full --output output.json
python -m rvestudio.cli benchmark --tolerance 1 --csv comparison.csv
python -m rvestudio.cli benchmark --wall-beam timoshenko --bending twist-corrected --twist 2.0 2.15 0.73
python -m rvestudio.cli benchmark --json
python -m rvestudio.cli generate-db inputs.csv --cases 1000 --wall-beam timoshenko --bending twist-corrected
python -m rvestudio.cli batch inputs.csv outputs.csv --workers 2 --limit 50
```

### Model/solver compatibility

`RVEModel.wall_element_type='paper'` selects the cylindrical strain-energy wall implementation. Previously saved `'thesis'` and `'original'` identifiers are accepted for backward compatibility; new exports use the preferred `'paper'` identifier. The separate `'anisotropicStrip'` solver remains supported.

Buckling is a reduced-order wall screening calculation, not a direct reproduction of the paper's model or experimental results.
