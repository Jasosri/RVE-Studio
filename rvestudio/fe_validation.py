"""Reproducible FE validation campaigns using the Abaqus S4 direct-CH/PBC backend.

The preparation and report steps do *not* require Abaqus; only ``run`` executes
external FE jobs. Missing FE results are explicitly marked PENDING (never
substituted with analytical predictions or synthetic energies).
"""
from __future__ import annotations

import csv
import json
import math
import re
import shutil
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.linalg import eigh

from .io import load_project
from .homogenization import homogenize_rve
from .abaqus_homogenization import (
    PeriodicHomogenizationSettings, PBC_DOMAIN_LATTICE_SUPERCELL, PBC_DOMAIN_PUBLISHED, PBC_DOMAIN_SYMMETRIC_RECTANGULAR,
    prepare_periodic_homogenization, run_periodic_homogenization, CH_METHOD_DIRECT, CH_METHOD_ENERGY,
)

VOIGT = ('11', '22', '33', '44', '55', '66')  # matrix *indices*, not shear-plane names
MANIFEST = 'campaign_manifest.json'


def _matrix(A, *, name='matrix') -> np.ndarray:
    C = np.asarray(A, dtype=float)
    if C.shape != (6, 6) or not np.all(np.isfinite(C)):
        raise ValueError(f'{name} must be a finite 6x6 matrix')
    if not np.allclose(C, C.T, atol=1e-8, rtol=1e-8):
        raise ValueError(f'{name} must be symmetric')
    return (C + C.T) / 2


def _safe_id(s: str) -> str:
    value = re.sub(r'[^A-Za-z0-9_.-]+', '_', s).strip('._')
    if not value or value in ('.', '..'):
        raise ValueError(f'Invalid case ID: {s!r}')
    return value


def _target_metrics(C: np.ndarray, target: np.ndarray) -> dict:
    C = _matrix(C, name='stiffness')
    target = _matrix(target, name='target')
    eig_target = np.linalg.eigvalsh(target)
    if eig_target.min() <= 0:
        raise ValueError('target stiffness is not positive definite')
    vals = eigh(C, target, eigvals_only=True)
    # Congruence-normalised mismatch gives every strain mode comparable weight.
    w, V = np.linalg.eigh(target)
    invsqrt = (V / np.sqrt(w)) @ V.T
    gap = invsqrt @ (C - target) @ invsqrt
    return {
        'min_generalized_eigenvalue': float(vals.min()),
        'max_generalized_eigenvalue': float(vals.max()),
        'relative_frobenius_difference': float(np.linalg.norm(C - target) / np.linalg.norm(target)),
        'target_scaled_difference': float(np.linalg.norm(gap) / np.sqrt(6)),
    }


def _fe_comparison(C: np.ndarray, analytical: np.ndarray, target: np.ndarray | None) -> dict:
    C = _matrix(C, name='FE stiffness')
    analytical = _matrix(analytical, name='analytical stiffness')
    diff = C - analytical
    with np.errstate(divide='ignore', invalid='ignore'):
        rel = 100 * np.diag(diff) / np.diag(analytical)
    d = {
        'fe_vs_analytical_relative_frobenius': float(np.linalg.norm(diff) / np.linalg.norm(analytical)),
        'fe_vs_analytical_max_abs_coefficient_MPa': float(np.max(np.abs(diff))),
        'fe_vs_analytical_diagonal_error_percent': rel.tolist(),
    }
    if target is not None:
        d.update({f'fe_target_{k}': v for k, v in _target_metrics(C, target).items()})
    return d


def prepare_campaign(root: str | Path, projects: Sequence[str | Path], *,
                     mesh_sizes: Sequence[float] = (10., 5., 2.5),
                     target_path: str | Path | None = None,
                     domain_mode: str = PBC_DOMAIN_PUBLISHED,
                     repeats: tuple[int, int] = (2, 2),
                     strain_amplitude: float = 1e-3,
                     cpus: int = 1,
                     identification_method: str = CH_METHOD_DIRECT,
                     resume: bool = False) -> dict:
    """Create direct-CH6 (or optional energy-CH21) decks and analytical baselines.

    Notes: ``prepare`` does not run FE, and never writes a completed FE matrix.
    """
    root = Path(root).expanduser().resolve()
    if not projects:
        raise ValueError('At least one RVE project is required')
    meshes = [float(h) for h in mesh_sizes]
    if not meshes or any(not math.isfinite(h) or h <= 0 for h in meshes) or len(set(meshes)) != len(meshes):
        raise ValueError('Mesh sizes must be distinct positive finite values')
    if root.exists() and (root / MANIFEST).exists():
        if not resume:
            raise FileExistsError(f'Campaign already exists: {root}; use --resume')
        existing = json.loads((root / MANIFEST).read_text(encoding='utf-8'))
        if {Path(x).stem for x in projects} != {Path(x).stem for x in existing['project_paths']}:
            raise ValueError('Cannot resume with a different project selection')
        if meshes != existing['mesh_sizes'] or domain_mode != existing['domain_mode'] or list(repeats) != existing['repeats'] or identification_method != existing.get('identification_method', CH_METHOD_ENERGY):
            raise ValueError('Cannot resume with changed mesh/domain settings')
    root.mkdir(parents=True, exist_ok=True)
    target = np.loadtxt(target_path, delimiter=',') if target_path else None
    if target is not None:
        target = _matrix(target, name='target')
        if np.linalg.eigvalsh(target).min() <= 0:
            raise ValueError('Target stiffness must be positive definite')
        np.savetxt(root / 'C_target_MPa.csv', target, delimiter=',', fmt='%.14g')
    records = []
    names = set()
    for p in projects:
        p = Path(p).expanduser().resolve()
        model, materials, laminates = load_project(p)
        case_id = _safe_id(p.stem)
        if case_id in names:
            raise ValueError(f'Duplicate project ID: {case_id}')
        names.add(case_id)
        project_copy = root / 'projects' / f'{case_id}.rvejson'
        project_copy.parent.mkdir(parents=True, exist_ok=True)
        if p != project_copy:
            shutil.copy2(p, project_copy)
        # The user's prescribed comparison is strictly against cylindrical bending.
        if model.wall_element_type.lower() not in ('thesis','original','paper','cylindrical'):
            raise ValueError(f'{p}: non-cylindrical wall model {model.wall_element_type!r} is not allowed')
        if model.strip_bending_mode.lower() != 'cylindrical':
            raise ValueError(f'{p}: non-cylindrical strip mode is not allowed')
        a = homogenize_rve(model, materials, laminates)
        analytic = _matrix(a.CH, name='analytical stiffness')
        ref_dir = root / 'analytical' / case_id
        ref_dir.mkdir(parents=True, exist_ok=True)
        np.savetxt(ref_dir / 'CH_cylindrical_MPa.csv', analytic, delimiter=',', fmt='%.14g')
        ref = {
            'project': str(Path('projects') / f'{case_id}.rvejson'), 'case_id': case_id, 'source_model': model.name,
            'wall_model': 'cylindrical', 'CH': analytic.tolist(),
            'relative_density': float(a.relative_density),
            'density_kg_m3': float(a.physical.average_density_kg_m3),
            'mass_kg': float(a.mass_kg),
            'area_mm2': float(a.area), 'height_mm': float(model.height),
            'areal_mass_kg_m2': float(a.physical.average_density_kg_m3 * model.height * 1e-3),
            'target_metrics': (_target_metrics(analytic, target) if target is not None else None),
        }
        (ref_dir / 'analytical_reference.json').write_text(json.dumps(ref, indent=2), encoding='utf-8')
        for h in meshes:
            s = PeriodicHomogenizationSettings(
                mesh_size=h, strain_amplitude=strain_amplitude, cpus=cpus, identification_method=identification_method,
                domain_mode=domain_mode, supercell_repeat_x=repeats[0],
                supercell_repeat_y=repeats[1], strict_pbc_audit=True,
            )
            # Fixed directory is deterministic across restarts and makes ODB
            # checkpoints reusable even if the campaign is interrupted.
            label = f'{case_id}_h{format(h, ".8g").replace(".", "p")}'
            files = prepare_periodic_homogenization(
                root / 'abaqus_runs', model, materials, laminates, s,
                run_name=label, resume=resume,
            )
            snap = json.loads(files.snapshot_file.read_text(encoding='utf-8'))
            meta = snap['cases'][0]['deck_meta']
            if not meta['pbc_audit']['passed']:
                raise ValueError(f'PBC audit failed: {label}')
            records.append({
                'case_id': case_id, 'project': str(Path('projects') / f'{case_id}.rvejson'), 'mesh_mm': h,
                'run_dir': str(files.run_dir.relative_to(root)),
                'mesh_nodes': int(meta['mesh_node_count']),
                'mesh_elements': int(meta['mesh_element_count']),
                'periodic_pairs': int(meta['periodic_pair_count']),
                'periodic_cell_count': int(snap['periodic_domain'].get('primitive_cell_count', 1)),
                'pbc_pass': True, 'status': 'PREPARED', 'error': '',
                'settings': asdict(s),
            })
    manifest = {
        'schema': 'rve_fe_validation_campaign_v1',
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'execution_note': 'PREPARED only; no FE result until Abaqus ODB jobs are solved and extracted',
        'identification_method': identification_method,
        'solver': ('Abaqus/Standard layered S4; six macro-strain reactions → direct CH' if identification_method == CH_METHOD_DIRECT else 'Abaqus/Standard layered S4; 21 strain-energy states → CH') + '; periodic affine PBC',
        'analytical_reference': 'Native RVE Studio cylindrical wall model only',
        'strain_order': ['11','22','33','23','13','12'],
        'Q_output_convention': 'Q44=13 shear, Q55=23 shear; internal C44=23, C55=13',
        'project_paths': [str(Path('projects') / (Path(x).stem + '.rvejson')) for x in projects],
        'mesh_sizes': meshes, 'domain_mode': domain_mode,
        'repeats': list(repeats), 'target': (target.tolist() if target is not None else None),
        'records': records,
    }
    (root / MANIFEST).write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    collect_campaign(root)
    return manifest


def _load_manifest(root):
    root = Path(root).expanduser().resolve()
    f = root / MANIFEST
    if not f.exists():
        raise FileNotFoundError(f'Missing FE campaign manifest: {f}')
    m = json.loads(f.read_text(encoding='utf-8'))
    if m.get('schema') != 'rve_fe_validation_campaign_v1':
        raise ValueError('Unsupported campaign manifest')
    return root, m


def run_campaign(root: str | Path, *, abaqus_command: str = 'abaqus',
                 resume: bool = True) -> dict:
    root, m = _load_manifest(root)
    if not abaqus_command.strip():
        raise ValueError('Specify an installed Abaqus command')
    for rec in m['records']:
        from .abaqus_homogenization import _periodic_files_from_directory
        rundir = root / rec['run_dir']
        snapshot = json.loads((rundir / 'homogenization_snapshot.json').read_text(encoding='utf-8'))
        files = _periodic_files_from_directory(rundir, snapshot)
        settings = PeriodicHomogenizationSettings(**rec['settings'])
        try:
            output = run_periodic_homogenization(files, abaqus_command, settings, resume=resume)
            rec['status'] = 'COMPLETED' if output.get('status') == 'COMPLETED' else 'FAILED'
            rec['error'] = ''
        except Exception as exc:
            rec['status'] = 'FAILED'
            rec['error'] = str(exc)
        (root / MANIFEST).write_text(json.dumps(m, indent=2), encoding='utf-8')
    collect_campaign(root)
    return m


def collect_campaign(root: str | Path) -> list[dict]:
    """Create CSV, Markdown and plots from completed *real* Abaqus result JSONs.

    A matrix is accepted only from ``homogenized_results.json`` with
    ``status=COMPLETED``; never substitute an analytical result for missing FE.
    """
    root, manifest = _load_manifest(root)
    target = _matrix(manifest['target'], name='target') if manifest['target'] is not None else None
    rows = []
    for rec in manifest['records']:
        base = json.loads((root / 'analytical' / rec['case_id'] / 'analytical_reference.json').read_text(encoding='utf-8'))
        A = _matrix(base['CH'], name='analytical reference')
        result_path = root / rec['run_dir'] / 'homogenized_results.json'
        row = {
            'CaseID': rec['case_id'], 'Mesh_mm': rec['mesh_mm'],
            'MeshNodes': rec['mesh_nodes'], 'ElementsS4': rec['mesh_elements'],
            'CellCount': rec['periodic_cell_count'], 'PBCAuditPassed': rec['pbc_pass'],
            'Status': 'PENDING', 'Density_kg_m3': base['density_kg_m3'],
            'Mass_kg_m2': base['areal_mass_kg_m2'],
        }
        for i, label in enumerate(VOIGT):
            row[f'C{label}_Analytical_MPa'] = float(A[i, i])
        if base['target_metrics'] is not None:
            row['Analytical_min_C_target_eigenvalue'] = base['target_metrics']['min_generalized_eigenvalue']
        if result_path.exists():
            try:
                raw = json.loads(result_path.read_text(encoding='utf-8'))
                if raw.get('status') != 'COMPLETED':
                    row['Status'] = 'FE_FAILED'
                    row['Error'] = raw.get('error', 'FE summary status is not COMPLETED')
                else:
                    FE = _matrix(raw.get('CH', raw.get('CH_energy')), name='FE stiffness')
                    if raw.get('identification_method') == CH_METHOD_DIRECT:
                        check = raw.get('direct_reaction_audit')
                        if not isinstance(check, dict):
                            raise ValueError('Missing direct-reaction work/reciprocity audit')
                        if float(check.get('max_energy_reaction_work_error', math.inf)) > 0.01:
                            raise ValueError('PBC reaction-to-strain-energy work mismatch exceeds 1%')
                        if float(check.get('macro_reaction_symmetry_error', math.inf)) > 0.02:
                            raise ValueError('Raw CH reciprocity error exceeds 2%')

                    if not np.all(np.diag(FE) > 0):
                        raise ValueError('Non-positive FE diagonal stiffness')
                    if 'pbc_audit' in raw and not raw['pbc_audit'].get('passed', False):
                        raise ValueError('FE summary contains failed PBC audit')
                    metrics = _fe_comparison(FE, A, target)
                    row.update({'Status': 'FE_COMPLETED',
                                'FE_vs_analytical_Frobenius_pct': 100 * metrics['fe_vs_analytical_relative_frobenius'],
                                'FE_vs_analytical_MaxAbs_MPa': metrics['fe_vs_analytical_max_abs_coefficient_MPa'],
                                'FE_reaction_energy_difference': raw.get('energy_vs_reaction_relative_error'),
                                'FE_positive_definite': bool(np.linalg.eigvalsh(FE).min() > 0)})
                    if target is not None:
                        row['FE_min_C_target_eigenvalue'] = metrics['fe_target_min_generalized_eigenvalue']
                        row['FE_target_scaled_difference'] = metrics['fe_target_target_scaled_difference']
                    for i, label in enumerate(VOIGT):
                        row[f'C{label}_FE_MPa'] = float(FE[i, i])
                        row[f'C{label}_FE_Difference_pct'] = metrics['fe_vs_analytical_diagonal_error_percent'][i]
                    fe_out = root / 'fe_matrices' / rec['case_id']; fe_out.mkdir(parents=True, exist_ok=True)
                    np.savetxt(fe_out / f'C_FE_h{str(rec["mesh_mm"]).replace(".","p")}_MPa.csv', FE, delimiter=',', fmt='%.14g')
            except Exception as exc:
                row['Status'] = 'FE_INVALID'
                row['Error'] = str(exc)
        rows.append(row)
    fields = list(dict.fromkeys(key for r in rows for key in r))
    with (root / 'fe_convergence_summary.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)
    _make_report(root, manifest, rows)
    _make_figures(root, rows, manifest)
    return rows


def _make_report(root, manifest, rows):
    completed = sum(r['Status'] == 'FE_COMPLETED' for r in rows)
    total = len(rows)
    lines = [
        '# RVE Studio — independent Abaqus FE validation campaign', '',
        f'**Completed FE mesh studies: {completed}/{total}.**', '',
        'Only completed, extracted Abaqus homogenisation results count as FE evidence. '
        'Analytical values and generated input decks are not FE validation.', '',
        '**Models:** Periodic S4 layered shells, ' + ('six direct macro-reaction strain cases' if manifest.get('identification_method') == CH_METHOD_DIRECT else '21 strain-energy cases') + '; exact lattice domain where selected.',
        '**Analytical reference:** cylindrical bending only. Full 6×6 internal order: 11, 22, 33, 23, 13, 12.',
        '**Shell shear naming:** Q44=13 and Q55=23 (not internal C44/C55).',
        '**Scope:** periodic shell homogenisation, not full 3D continuum solid FE or eigenbuckling validation.', '',
        '| Design | h (mm) | S4 elements | PBC audit | FE status |',
        '|---|---:|---:|---|---|',
    ]
    for r in rows:
        lines.append(f'| {r["CaseID"]} | {r["Mesh_mm"]:g} | {r["ElementsS4"]:,} | {"PASS" if r["PBCAuditPassed"] else "FAIL"} | {r["Status"]} |')
    if completed == 0:
        lines += ['', '**FE validation pending.** Abaqus was not run or no completed homogenisation result JSON was found.']
    else:
        lines += ['', '### Actual completed FE comparisons', '', '| Design | h (mm) | Relative FE/native full-C difference |', '|---|---:|---:|']
        for r in rows:
            if r['Status'] == 'FE_COMPLETED':
                lines.append(f'| {r["CaseID"]} | {r["Mesh_mm"]:g} | {r["FE_vs_analytical_Frobenius_pct"]:.3f}% |')
    lines += ['', '### Execution', '',
              '```bash',
              'python -m rvestudio.cli fe-campaign run --root "<campaign directory>" --abaqus-command abaqus',
              'python -m rvestudio.cli fe-campaign collect --root "<campaign directory>"',
              '```', '',
              'Use the CSV and generated plots to assess mesh convergence before interpreting FE/native differences.',
              'The 21-state procedure identifies a shell-RVE homogenised stiffness by energy equivalence; '
              'the full 3D *solid* FE and independent eigenbuckling studies remain separate work.', '']
    (root / 'FE_VALIDATION_REPORT.md').write_text('\n'.join(lines), encoding='utf-8')


def _make_figures(root, rows, manifest):
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    figs = root / 'figures'; figs.mkdir(exist_ok=True)
    cases = sorted({r['CaseID'] for r in rows})
    target = np.asarray(manifest['target'], float) if manifest['target'] is not None else None
    for case in cases:
        a = json.loads((root / 'analytical' / case / 'analytical_reference.json').read_text())
        C = np.asarray(a['CH'], float)
        fig, ax = plt.subplots(figsize=(7, 4))
        x = np.arange(6)
        ax.bar(x-0.18, np.diag(C), width=0.36, label='Analytical cylindrical')
        if target is not None:
            ax.bar(x+0.18, np.diag(target), width=0.36, label='Target C')
        ax.set_xticks(x, [f'C{k}' for k in VOIGT]); ax.set_ylabel('Stiffness (MPa)')
        ax.set_title(f'{case}: analytical / target before FE solution')
        ax.legend(); fig.tight_layout(); fig.savefig(figs / f'{case}_analytical_vs_target.png', dpi=160)
        plt.close(fig)
        found = [r for r in rows if r['CaseID'] == case and r['Status'] == 'FE_COMPLETED']
        project = root / 'projects' / f'{case}.rvejson'
        if project.exists():
            geometry, _, _ = load_project(project)
            fig, ax = plt.subplots(figsize=(5, 5))
            for k, wall in enumerate(geometry.elements):
                v = geometry.nodes[[wall.node_i, wall.node_j]]
                ax.plot(v[:, 0], v[:, 1], linewidth=2.6)
                ax.text(float(v[:, 0].mean()), float(v[:, 1].mean()), f'W{k+1}', fontsize=9)
            ax.scatter(geometry.nodes[:, 0], geometry.nodes[:, 1], s=18)
            o = np.asarray(geometry.nodes[0]); a1 = np.asarray(geometry.a1); a2 = np.asarray(geometry.a2)
            boundary = np.array([o, o+a1, o+a1+a2, o+a2, o])
            ax.plot(boundary[:, 0], boundary[:, 1], '--', color='black', alpha=.6, label='Primitive lattice')
            ax.set_aspect('equal', adjustable='box'); ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')
            ax.set_title(f'{case}: RVE wall midsurfaces'); ax.legend()
            fig.tight_layout(); fig.savefig(figs / f'{case}_RVE_geometry.png', dpi=180); plt.close(fig)
        if found:
            found.sort(key=lambda r: r['Mesh_mm'], reverse=True)
            fig, ax = plt.subplots(figsize=(7, 4))
            for label in ('11', '33', '44', '55', '66'):
                ax.plot([r['Mesh_mm'] for r in found], [r[f'C{label}_FE_MPa'] for r in found], marker='o', label=f'C{label} FE')
                ax.axhline(float(C[VOIGT.index(label), VOIGT.index(label)]), linestyle=':', alpha=.35)
            ax.invert_xaxis(); ax.set_xlabel('Shell element size (mm)'); ax.set_ylabel('Effective stiffness (MPa)')
            ax.set_title(f'{case}: mesh convergence (Abaqus)'); ax.legend(ncol=2, fontsize=8)
            fig.tight_layout(); fig.savefig(figs / f'{case}_FE_mesh_convergence.png', dpi=160);plt.close(fig)
