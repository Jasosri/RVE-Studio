#!/usr/bin/env python3
"""Stand-in for the Abaqus launcher, used only by the test-suite.

``fake_abaqus job=J input=J.inp cpus=N interactive`` writes J.odb.
``fake_abaqus python extract_dbc_odb.py J.odb case_results.json`` writes an ALLSE
energy built from the study's analytical Q (shell cases) or C (U3-restrained cases)
times (1 + 0.04 (h/10)^2), so a three-mesh Richardson extrapolation must recover the
analytical value with p = 2.  ``extract_hom_odb.py`` (PBC) writes the energy and the
six macro reactions of the analytical C^H with the same mesh error.  With
FAKE_ABAQUS_PBC_FAIL=1 the PBC extraction fails instead, to exercise failure handling.
"""
import json
import os
import sys
from pathlib import Path

args = sys.argv[1:]
cwd = Path.cwd()
if args and args[0].startswith('job='):
    job = args[0].split('=', 1)[1]
    (cwd / f'{job}.odb').write_text('fake odb', encoding='utf-8')
    sys.exit(0)


def mesh_factor(h):
    return 1 + 0.04 * (h / 10.0) ** 2


if args and args[0] == 'python' and args[1] == 'extract_hom_odb.py':
    out = args[-1]
    if os.environ.get('FAKE_ABAQUS_PBC_FAIL') == '1':
        Path(out).write_text(json.dumps({'status': 'ERROR', 'error': 'fake PBC extractor failure'}), encoding='utf-8')
        sys.exit(2)
    run_dir = cwd.parent
    root = run_dir.parent.parent
    snap = json.loads((run_dir / 'homogenization_snapshot.json').read_text(encoding='utf-8'))
    study = json.loads((root / 'fe_study.json').read_text(encoding='utf-8'))
    pid = run_dir.name.rsplit('_h', 1)[0]
    h = float(snap['settings']['mesh_size'])
    L = float(snap['settings']['reference_length'])
    V = float(snap['rve_volume'])
    C = [[c * mesh_factor(h) for c in row] for row in study['projects'][pid]['analytical']['CH']]
    case = next(c for c in snap['cases'] if c['directory'] == cwd.name)
    e = [float(x) for x in case['strain']]
    sig = [sum(C[i][j] * e[j] for j in range(6)) for i in range(6)]
    W = 0.5 * V * sum(sig[i] * e[i] for i in range(6))
    ids = case['deck_meta']['macro_node_ids']
    Path(out).write_text(json.dumps({'status': 'COMPLETED', 'strain_energy_ALLSE': W,
                                     'macro_reaction_u1': {str(n): V * sig[k] / L for k, n in enumerate(ids)},
                                     'macro_displacement_u1': {str(n): L * e[k] for k, n in enumerate(ids)}}),
                         encoding='utf-8')
    sys.exit(0)

if args and args[0] == 'python':
    extractor, out = args[1], args[-1]
    run_dir = cwd.parent
    root = run_dir.parent.parent
    snap = json.loads((run_dir / 'thesis_dbc_snapshot.json').read_text(encoding='utf-8'))
    study = json.loads((root / 'fe_study.json').read_text(encoding='utf-8'))
    pid = run_dir.name.rsplit('_h', 1)[0]
    h = float(snap['settings']['mesh_size'])
    d = float(snap['settings']['displacement'])
    ana = study['projects'][pid]['analytical']
    Q = {k: v * mesh_factor(h) for k, v in ana['Q'].items()}
    C = {k: v * mesh_factor(h) for k, v in ana['C'].items()}
    dom = snap['domain']
    p, q, H = float(dom['p']), float(dom['q']), float(dom['height'])
    V = p * q * H
    e11, e22 = d / p, d / q
    W = {'Q11': 0.5 * V * Q['Q11'] * e11 ** 2, 'Q22': 0.5 * V * Q['Q22'] * e22 ** 2,
         'Q66': 0.5 * V * Q['Q66'] * e11 ** 2, 'Q44': 0.5 * V * Q['Q44'] * e11 ** 2,
         'Q55': 0.5 * V * Q['Q55'] * e22 ** 2,
         'C1111': 0.5 * V * C['C1111'] * e11 ** 2, 'C2222': 0.5 * V * C['C2222'] * e22 ** 2,
         'C1212': 0.5 * V * C['C1212'] * e11 ** 2}
    W['Q12'] = W['Q11'] + W['Q22'] + Q['Q12'] * V * e11 * e22
    W['C1122'] = W['C1111'] + W['C2222'] + C['C1122'] * V * e11 * e22
    case = cwd.name.split('_', 1)[1]
    side = {'C1111': 'Q11', 'C2222': 'Q22', 'C1122': 'Q12', 'C1212': 'Q66'}.get(case, case)
    # Nodal export, as the real extractor writes it: an affine field of the case.
    lines, mode = ['label,x,y,z,u1,u2,u3'], None
    for line in next(cwd.glob('*.inp')).read_text(encoding='utf-8').splitlines():
        if line.startswith('*'):
            mode = 'n' if line.upper().startswith('*NODE') and 'OUTPUT' not in line.upper() else None
            continue
        if mode == 'n' and line.strip() and not line.startswith('**'):
            n, x, y, z = (float(v) for v in line.split(',')[:4])
            ux = d / p * x if side in ('Q11', 'Q12') else 0.0
            uy = d / q * y if side in ('Q22', 'Q12') else (d / p * x if side == 'Q66' else 0.0)
            uz = (d / p * x if side == 'Q44' else d / q * y) if side in ('Q44', 'Q55') else 0.0
            lines.append(f'{int(n)},{x},{y},{z},{ux},{uy},{uz}')
    (cwd / 'nodal_displacements.csv').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    Path(out).write_text(json.dumps({'status': 'COMPLETED', 'strain_energy_ALLSE': W[case],
                                     'nodal_displacements': len(lines) - 1}), encoding='utf-8')
    sys.exit(0)
sys.exit(1)
