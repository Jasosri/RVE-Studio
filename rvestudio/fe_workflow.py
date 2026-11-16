"""FE verification study: one folder, the built-in check and Abaqus, one report.

A study holds a set of RVE projects and mesh sizes.  Every FE model uses the
FE RVE of the project (automatic, or chosen by the user) and one of two boundary-value
problems; every FE value is compared with the analytical model:

* ``pbc`` -- strain-energy homogenisation on the periodic RVE (top and bottom periodic
  through the height, U3+ - U3- = E33 h).  21 S4 jobs per mesh: six unit strains and the
  15 pairs; C_ii = 2 W_i / (V e^2), C_ij = (W_ij - W_i - W_j) / (V e^2) from ALLSE give the
  full 6x6 C, from which follow the Q matrix (sigma33 = 0) and the engineering constants.
  The built-in check gets the same C from the mutual energies u_i^T K u_j / V of six states.
* ``dbc`` -- finite-height core, ten S4 jobs per mesh.  Six cases of thesis
  Table 3.4 (top/bottom U3 free) give the shell terms Q11, Q22, Q12, Q66, Q44, Q55
  (sigma33 = 0); four cases with the same side faces and U3 = 0 on the top and
  bottom faces give C1111, C2222, C1122, C1212 (eps33 = 0); C1313 and C2323 come
  from the Q44 and Q55 cases.

``check`` solves both models with the built-in shell element in seconds;
``prepare -> run -> collect`` does the same with Abaqus.  Both write long-format
result rows (``verification.result_rows``) into ``results/`` and refresh one HTML
report and one Excel workbook.  Every step reads its state from disk, so a study
can be prepared here, run on a licensed workstation or cluster (from Python or from
the generated ``run_all_jobs`` scripts, without Python), copied back and collected.

Nothing here invents FE numbers: an Abaqus result counts only when an ODB extractor
has written ``case_results.json`` with ``status=COMPLETED``.
"""
from __future__ import annotations

import csv
import json
import math
import os
import shutil
import subprocess
from dataclasses import dataclass, asdict, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable, Sequence

import numpy as np

from .io import load_project, save_project
from .examples_data import paper_hexagonal_case
from .abaqus import _run_shell_command, _safe_name, detect_abaqus_executable
from .abaqus_thesis_dbc import (DBC_CASES, ThesisDBCSettings, prepare_thesis_dbc, reconstruct_dbc_terms,
                                FACE_WALL_SECTIONS)
from .abaqus_homogenization import (_periodic_files_from_directory, PeriodicHomogenizationSettings,
                                    run_periodic_homogenization, PBC_DOMAIN_PUBLISHED, CH_METHOD_ENERGY)
from .verification import (analytical_reference, paper_configuration_of, builtin_check, result_rows, pbc_values,
                           SHELL_TERMS, SOLID_TERMS, GROUPS, SOLVERS, best_values, pct)

STUDY_FILE = 'fe_study.json'
SCHEMA = 'rve_fe_study_v1'
FORMAT = 3  # 1.13: 21-state strain-energy PBC (full C, Q, engineering constants); rows with units
CHECK_DBC = 'dbc'
CHECK_PBC = 'pbc'
CHECKS = (CHECK_DBC, CHECK_PBC)
MODEL_OF = {CHECK_DBC: 'DBC', CHECK_PBC: 'PBC'}

# Job states, read from disk.
PENDING, RUNNING, SOLVED, DONE, FAILED = 'PENDING', 'RUNNING', 'SOLVED', 'DONE', 'FAILED'

CHECK_DESCRIPTIONS = {
    CHECK_DBC: ('Finite-height core, side faces as thesis Table 3.4. Top/bottom U3 free: shell terms '
                'Q11, Q22, Q12, Q66, Q44, Q55. U3 = 0 on top and bottom: 3D terms C1111, C2222, C1122, C1212.'),
    CHECK_PBC: ('Strain energy on the periodic RVE (top and bottom periodic through the height): the full 6x6 C, '
                'the Q matrix (sigma33 = 0) and the engineering constants. Abaqus: 21 load cases per mesh.'),
}
JOBS_PER_MESH = {CHECK_DBC: len(DBC_CASES), CHECK_PBC: 21}
BUILTIN_RESULTS = 'builtin_results.csv'
ABAQUS_RESULTS = 'abaqus_results.csv'


@dataclass(slots=True)
class FEStudySettings:
    meshes: list[float] = field(default_factory=lambda: [10.0])
    checks: list[str] = field(default_factory=lambda: [CHECK_DBC, CHECK_PBC])
    cpus: int = 1
    # Fixed by the method; kept in the study file so that every deck is reproducible.
    pbc_domain: str = 'published-rve'
    pbc_repeats: tuple[int, int] = (1, 1)
    pbc_identification: str = CH_METHOD_ENERGY
    pbc_strain: float = 1e-3
    dbc_displacement: float = 1.0
    dbc_face_walls: str = 'half-thickness'

    def validate(self) -> None:
        m = [float(h) for h in self.meshes]
        if not m or any((not math.isfinite(h)) or h <= 0 for h in m) or len(set(m)) != len(m):
            raise ValueError('Mesh sizes must be distinct, positive and finite.')
        bad = [c for c in self.checks if c not in CHECKS]
        if bad or not self.checks:
            raise ValueError(f'Checks must be a non-empty subset of {CHECKS}; got {self.checks}')
        if self.pbc_domain != 'published-rve' or self.pbc_identification != CH_METHOD_ENERGY:
            raise ValueError('The FE verification study uses the FE RVE of each design and the 21-state '
                             'strain-energy PBC.')
        if int(self.cpus) < 1:
            raise ValueError('CPU count must be at least 1')
        if self.dbc_face_walls not in FACE_WALL_SECTIONS:
            raise ValueError(f'DBC face-wall treatment must be one of {FACE_WALL_SECTIONS}')


# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------

def _mesh_tag(h: float) -> str:
    return 'h' + format(float(h), '.8g').replace('.', 'p')


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def _read_json(path: Path) -> dict | None:
    try:
        return json.loads(Path(path).read_text(encoding='utf-8'))
    except Exception:
        return None


def _write_json(path: Path, data) -> None:
    Path(path).write_text(json.dumps(data, indent=2), encoding='utf-8')


def _write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for r in rows for k in r)) or ['empty']
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def read_rows(path: str | Path) -> list[dict]:
    """Result rows from a CSV, numbers as floats."""
    path = Path(path)
    if not path.exists():
        return []
    out = []
    with path.open(encoding='utf-8') as f:
        for r in csv.DictReader(f):
            row = {}
            for k, v in r.items():
                if k in ('mesh_mm', 'FE', 'analytical') or k.endswith('_MPa') or k.endswith('_pct') \
                        or k.startswith('FE_') or k in ('extrapolated', 'p', 'gci_fine_pct', 'asymptotic_ratio'):
                    try:
                        row[k] = float(v)
                    except (TypeError, ValueError):
                        row[k] = math.nan
                elif k == 'negligible':
                    row[k] = str(v).strip().lower() == 'true'
                else:
                    row[k] = v
            out.append(row)
    return out


def parse_config_list(text: str | Iterable[int] | None) -> list[int]:
    """Parse '1-14', '1,3,7' or '1-3,9' into a sorted list of paper configurations."""
    if text is None:
        return []
    if not isinstance(text, str):
        out = sorted({int(x) for x in text})
    else:
        out = set()
        for part in text.replace(' ', '').split(','):
            if not part:
                continue
            if '-' in part:
                a, b = part.split('-', 1)
                out.update(range(int(a), int(b) + 1))
            else:
                out.add(int(part))
        out = sorted(out)
    if any(not (1 <= c <= 14) for c in out):
        raise ValueError('Paper configurations must lie in 1..14')
    return out


def write_paper_config_projects(dest: str | Path, configs: Iterable[int]) -> list[Path]:
    """Write the benchmark hexagonal configurations as ordinary .rvejson projects."""
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    paths = []
    for c in parse_config_list(configs):
        model, mats, lams = paper_hexagonal_case(c)
        p = dest / f'paper_config_{c:02d}.rvejson'
        save_project(p, model, mats, lams)
        paths.append(p)
    return paths


def _case_state(case_dir: Path, job: str) -> tuple[str, str]:
    """Infer one Abaqus job's state from the files it leaves behind."""
    res = _read_json(case_dir / 'case_results.json')
    if res is not None:
        st = str(res.get('status', '')).upper()
        if st == 'COMPLETED':
            return DONE, ''
        return FAILED, str(res.get('error', 'extractor reported an error'))[:300]
    if (case_dir / f'{job}.lck').exists():
        return RUNNING, ''
    sta = case_dir / f'{job}.sta'
    if sta.exists():
        try:
            txt = sta.read_text(encoding='utf-8', errors='replace')
        except Exception:
            txt = ''
        if 'HAS NOT BEEN COMPLETED' in txt.upper():
            return FAILED, 'Abaqus .sta reports the analysis was not completed'
    if (case_dir / f'{job}.odb').exists():
        return SOLVED, 'ODB present; results not yet extracted'
    return PENDING, ''


# --------------------------------------------------------------------------------------
# Prepare
# --------------------------------------------------------------------------------------

def _reference(model, mats, lams) -> dict:
    ana = analytical_reference(model, mats, lams)
    return {**ana, 'paper_config': paper_configuration_of(model, mats, lams)}


def add_projects(root: str | Path, projects: Sequence[str | Path] = (), *,
                 paper_configs: str | Iterable[int] | None = None, settings: FEStudySettings | None = None,
                 progress: Callable[[str], None] | None = None) -> dict:
    """Create the study folder if needed and add projects (analytical reference only, no decks)."""
    s = settings or FEStudySettings()
    s.validate()
    root = Path(root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    say = progress or (lambda msg: None)
    existing = _read_json(root / STUDY_FILE)
    if existing is not None:
        if existing.get('schema') != SCHEMA:
            raise ValueError(f'{root} holds an unrecognised study file')
        problems = stale_runs(root, existing)
        if problems:
            raise ValueError(STALE_STUDY_MESSAGE)
        old = existing['settings']
        for key in ('pbc_identification', 'pbc_strain', 'dbc_displacement', 'dbc_face_walls'):
            if old.get(key) != asdict(s).get(key):
                raise ValueError(f'Cannot extend this study with a different {key!r} '
                                 f'({old.get(key)!r} -> {asdict(s).get(key)!r}); start a new study folder.')
        s.meshes = sorted(set(map(float, old['meshes'])) | set(map(float, s.meshes)), reverse=True)
        s.checks = [c for c in CHECKS if c in set(old['checks']) | set(s.checks)]
        study = existing
    else:
        study = {'schema': SCHEMA, 'format': FORMAT, 'created_utc': _now(), 'projects': {}, 'runs': []}
    study['settings'] = asdict(s)
    study['settings']['pbc_repeats'] = list(s.pbc_repeats)
    study['updated_utc'] = _now()
    study['descriptions'] = CHECK_DESCRIPTIONS
    proj_dir = root / 'projects'
    proj_dir.mkdir(exist_ok=True)
    sources = [Path(p).expanduser().resolve() for p in projects]
    if paper_configs:
        sources += write_paper_config_projects(proj_dir, paper_configs)
    for src in sources:
        base = _safe_name(src.stem, 'P')
        pid, k = base, 1
        # A different design under a name already in the study gets its own id (never a stale reference).
        while pid in study['projects'] and src != proj_dir / f'{pid}.rvejson' and \
                (proj_dir / f'{pid}.rvejson').read_bytes() != src.read_bytes():
            k += 1
            pid = f'{base}_{k}'
        dst = proj_dir / f'{pid}.rvejson'
        if src != dst:
            shutil.copy2(src, dst)
        if pid in study['projects']:
            continue
        model, mats, lams = load_project(dst)
        say(f'Analytical reference: {pid}')
        study['projects'][pid] = {'project': str(dst.relative_to(root)), 'name': model.name,
                                  'analytical': _reference(model, mats, lams)}
    if not study['projects']:
        raise ValueError('Add at least one project or paper configuration.')
    _write_json(root / STUDY_FILE, study)
    return study


def prepare_study(root: str | Path, projects: Sequence[str | Path] = (), *,
                  paper_configs: str | Iterable[int] | None = None,
                  settings: FEStudySettings | None = None,
                  progress: Callable[[str], None] | None = None) -> dict:
    """Add projects and write the Abaqus decks of every missing project/mesh/model.  Never runs Abaqus."""
    study = add_projects(root, projects, paper_configs=paper_configs, settings=settings, progress=progress)
    root = Path(root).expanduser().resolve()
    say = progress or (lambda msg: None)
    s = study['settings']
    have = {(r['check'], r['project'], float(r['mesh'])) for r in study['runs']}
    for pid, info in study['projects'].items():
        model, mats, lams = load_project(root / info['project'])
        for h in s['meshes']:
            for check in s['checks']:
                if (check, pid, float(h)) in have:
                    continue
                run = {'check': check, 'project': pid, 'mesh': float(h), 'run_dir': None,
                       'prepared': False, 'error': ''}
                name = f'{pid}_{_mesh_tag(h)}'
                try:
                    if check == CHECK_DBC:
                        say(f'DBC decks: {name}')
                        target = root / 'dbc' / name
                        if not (target / 'thesis_dbc_snapshot.json').exists():
                            if target.exists():
                                shutil.rmtree(target)
                            prepare_thesis_dbc(root / 'dbc', model, mats, lams,
                                               ThesisDBCSettings(mesh_size=h, displacement=s['dbc_displacement'],
                                                                 cpus=s['cpus'], face_wall_section=s['dbc_face_walls']),
                                               run_name=name)
                        run['run_dir'] = str(target.relative_to(root))
                    else:
                        say(f'PBC decks: {name}')
                        from .abaqus_homogenization import prepare_periodic_homogenization
                        ps = PeriodicHomogenizationSettings(
                            mesh_size=h, strain_amplitude=s['pbc_strain'], cpus=s['cpus'],
                            identification_method=CH_METHOD_ENERGY, domain_mode=PBC_DOMAIN_PUBLISHED,
                            strict_pbc_audit=True)
                        files = prepare_periodic_homogenization(root / 'pbc', model, mats, lams, ps,
                                                                run_name=name, resume=True)
                        snap = _read_json(files.snapshot_file) or {}
                        meta = (snap.get('cases') or [{}])[0].get('deck_meta', {})
                        if not meta.get('pbc_audit', {}).get('passed', False):
                            raise ValueError('periodic boundary audit failed')
                        run['run_dir'] = str(files.run_dir.relative_to(root))
                        run['pbc_settings'] = asdict(ps)
                    run['prepared'] = True
                except Exception as exc:  # keep going: one awkward geometry must not block the study
                    run['error'] = f'{type(exc).__name__}: {exc}'
                study['runs'].append(run)
    study['runs'].sort(key=lambda r: (r['project'], CHECKS.index(r['check']), -r['mesh']))
    _write_json(root / STUDY_FILE, study)
    write_run_scripts(root)
    return study


STALE_STUDY_MESSAGE = ('This study folder was prepared with an older RVE Studio. Since 1.13 the PBC model is a '
                       '21-state strain-energy homogenisation (full C, Q and engineering constants) and results carry '
                       'units; before 1.11.1 the Abaqus wall orientation put the ply 1-axis along the core height. '
                       'Prepare a new study folder (or rename this one) with RVE Studio 1.13.')


def _deck_orientation_stale(run_dir: Path) -> bool:
    """True for decks whose *ORIENTATION has no '3, 0.' line (Abaqus then projected the ply axes about local 1)."""
    for inp in sorted(Path(run_dir).rglob('*.inp'))[:1]:
        lines = [ln.strip() for ln in inp.read_text(errors='ignore').splitlines() if not ln.startswith('**')]
        for i, ln in enumerate(lines):
            if ln.upper().startswith('*ORIENTATION'):
                return not (i + 2 < len(lines) and lines[i + 2].replace(' ', '') == '3,0.')
    return False


def stale_runs(root: Path, study: dict) -> list[str]:
    """Run folders prepared before 1.11 (no U3-restrained DBC cases, six-state PBC) or 1.9.1 (one-cell window)."""
    out = []
    for r in study.get('runs', []):
        if not r.get('run_dir'):
            continue
        if _deck_orientation_stale(Path(root) / r['run_dir']):
            out.append(r['run_dir'])
            continue
        if r.get('check') == CHECK_DBC:
            snap = _read_json(Path(root) / r['run_dir'] / 'thesis_dbc_snapshot.json') or {}
            dom = snap.get('domain', {})
            names = {c.get('name') for c in snap.get('cases', [])}
            if snap and ('primitive_cells' not in dom or not set(DBC_CASES) <= names):
                out.append(r['run_dir'])
        elif r.get('check') == CHECK_PBC:
            if (r.get('pbc_settings') or {}).get('identification_method') != CH_METHOD_ENERGY:
                out.append(r['run_dir'])
    return out


stale_dbc_runs = stale_runs  # name used by 1.9.1-1.10


def load_study(root: str | Path) -> tuple[Path, dict]:
    root = Path(root).expanduser().resolve()
    study = _read_json(root / STUDY_FILE)
    if study is None or study.get('schema') != SCHEMA:
        raise FileNotFoundError(f'No FE study found in {root} (expected {STUDY_FILE}).')
    if int(study.get('format', 0)) < FORMAT:
        raise ValueError(STALE_STUDY_MESSAGE)
    return root, study


# --------------------------------------------------------------------------------------
# Jobs and status
# --------------------------------------------------------------------------------------

def _run_cases(root: Path, run: dict) -> list[dict]:
    """List the individual Abaqus jobs of one prepared run."""
    if not run.get('prepared'):
        return []
    rdir = root / run['run_dir']
    out = []
    if run['check'] == CHECK_DBC:
        for i, case in enumerate(DBC_CASES, 1):
            cdir = rdir / f'{i:02d}_{case}'
            job = _safe_name(f'DBC_{case}', 'DBC')
            out.append({'case': case, 'dir': cdir, 'job': job, 'extractor': 'extract_dbc_odb.py'})
    else:
        snap = _read_json(rdir / 'homogenization_snapshot.json') or {}
        for c in _periodic_files_from_directory(rdir, snap).cases:
            out.append({'case': c.name, 'dir': c.case_dir, 'job': c.job_name, 'extractor': c.extractor_file.name})
    return out


def study_jobs(root: str | Path) -> list[dict]:
    """Every Abaqus job in the study with its state, read fresh from disk."""
    root, study = load_study(root)
    rows = []
    for run in study['runs']:
        if not run.get('prepared'):
            rows.append({'check': run['check'], 'project': run['project'], 'mesh': run['mesh'], 'case': '-',
                         'job': '-', 'dir': '', 'state': 'NOT PREPARED', 'note': run.get('error', '')})
            continue
        for c in _run_cases(root, run):
            state, note = _case_state(c['dir'], c['job'])
            rows.append({'check': run['check'], 'project': run['project'], 'mesh': run['mesh'],
                         'case': c['case'], 'job': c['job'], 'dir': str(c['dir'].relative_to(root)),
                         'state': state, 'note': note})
    return rows


def study_status(root: str | Path) -> dict:
    jobs = study_jobs(root)
    counts: dict[str, int] = {}
    for j in jobs:
        counts[j['state']] = counts.get(j['state'], 0) + 1
    return {'jobs': jobs, 'counts': counts, 'total': len(jobs)}


# --------------------------------------------------------------------------------------
# Abaqus environment check
# --------------------------------------------------------------------------------------

def abaqus_doctor(abaqus_command: str | None = None, timeout: float = 120.0) -> dict:
    """Check that an Abaqus command exists and answers (``abaqus information=release``)."""
    cmd = (abaqus_command or '').strip() or detect_abaqus_executable() or ''
    report = {'command': cmd, 'found': False, 'responds': False, 'release': '', 'messages': []}
    if not cmd:
        report['messages'].append('No Abaqus command found on PATH. Enter the full path to abaqus.bat '
                                  '(Windows) or the abaqus launcher (Linux), e.g. C:\\SIMULIA\\Commands\\abaqus.bat.')
        return report
    resolved = shutil.which(cmd) or (cmd if Path(cmd).exists() else None)
    report['found'] = bool(resolved)
    report['resolved_path'] = resolved or ''
    if not resolved:
        report['messages'].append(f'"{cmd}" is not on PATH and is not an existing file.')
        return report
    try:
        if os.name == 'nt':
            p = subprocess.run(subprocess.list2cmdline([cmd, 'information=release']), shell=True,
                               capture_output=True, text=True, timeout=timeout)
        else:
            p = subprocess.run([cmd, 'information=release'], capture_output=True, text=True, timeout=timeout)
        text = (p.stdout or '') + (p.stderr or '')
        report['responds'] = p.returncode == 0 and bool(text.strip())
        report['release'] = '\n'.join(l for l in text.splitlines() if l.strip())[:2000]
        if not report['responds']:
            report['messages'].append(f'Abaqus returned exit code {p.returncode}. Check the licence server.')
    except subprocess.TimeoutExpired:
        report['messages'].append(f'Abaqus did not answer within {timeout:g} s.')
    except Exception as exc:
        report['messages'].append(f'Could not start Abaqus: {exc}')
    if report['responds']:
        report['messages'].append('Abaqus responds. A full solve still needs a licence token at run time.')
    return report


# --------------------------------------------------------------------------------------
# Run
# --------------------------------------------------------------------------------------

def _solve_and_extract(exe: str, cdir: Path, job: str, extractor: str, cpus: int) -> tuple[bool, str]:
    with (cdir / 'abaqus_case.log').open('a', encoding='utf-8', errors='replace') as log:
        log.write(f'\n=== {_now()} solve {job} ===\n'); log.flush()
        odb = cdir / f'{job}.odb'
        if not odb.exists():
            for stale in cdir.glob(f'{job}.lck'):
                stale.unlink()
            rc = _run_shell_command([exe, f'job={job}', f'input={job}.inp', f'cpus={int(cpus)}', 'interactive'], cdir, log)
            if rc != 0 or not odb.exists():
                return False, f'Abaqus exit code {rc}; ODB {"present" if odb.exists() else "missing"}. See {cdir / "abaqus_case.log"}'
        log.write('\n--- extraction ---\n'); log.flush()
        _run_shell_command([exe, 'python', extractor, odb.name, 'case_results.json'], cdir, log)
    res = _read_json(cdir / 'case_results.json')
    if res is None:
        return False, 'Extractor wrote no case_results.json'
    if str(res.get('status', '')).upper() != 'COMPLETED':
        return False, str(res.get('error', 'extractor error'))[:300]
    return True, ''


def run_study(root: str | Path, abaqus_command: str, *, checks: Sequence[str] | None = None,
              projects: Sequence[str] | None = None, retry_failed: bool = False,
              progress: Callable[[int, int, dict], None] | None = None,
              should_stop: Callable[[], bool] | None = None) -> dict:
    """Solve every pending job, one at a time; a failure is recorded and the study moves on.

    Safe to interrupt and call again: finished jobs are skipped.  Collects at the end.
    """
    root, study = load_study(root)
    if stale_runs(root, study):
        raise ValueError(STALE_STUDY_MESSAGE)
    exe = str(abaqus_command or '').strip()
    if not exe:
        raise ValueError('Give the Abaqus command (for example "abaqus" or the full path to abaqus.bat).')
    cpus = int(study['settings'].get('cpus', 1))
    todo_states = {PENDING, SOLVED} | ({FAILED} if retry_failed else set())
    jobs = [j for j in study_jobs(root)
            if j['state'] in todo_states
            and (not checks or j['check'] in checks)
            and (not projects or j['project'] in projects)]
    done = failed = 0
    for i, j in enumerate(jobs, 1):
        if should_stop and should_stop():
            break
        if progress:
            progress(i, len(jobs), j)
        cdir = root / j['dir']
        if j['state'] == FAILED:
            for f in ('case_results.json', f"{j['job']}.odb"):
                try:
                    (cdir / f).unlink()
                except FileNotFoundError:
                    pass
        extractor = 'extract_dbc_odb.py' if j['check'] == CHECK_DBC else 'extract_hom_odb.py'
        ok, msg = _solve_and_extract(exe, cdir, j['job'], extractor, cpus)
        done += ok
        failed += (not ok)
        if not ok and not (cdir / 'case_results.json').exists():
            _write_json(cdir / 'case_results.json', {'status': 'ERROR', 'error': msg})
    summary = collect_study(root)
    return {'attempted': done + failed, 'completed': done, 'failed': failed, 'collect': summary['counts']}


def write_run_scripts(root: str | Path, abaqus_command: str = 'abaqus') -> list[Path]:
    """Write run_all_jobs.bat / .sh that solve and extract every job without Python.

    Use these on a cluster or workstation where RVE Studio is not installed; copy
    the folder back and run ``fe report``.  Finished jobs are skipped.
    """
    root, study = load_study(root)
    cpus = int(study['settings'].get('cpus', 1))
    bat = ['@echo off', 'REM RVE Studio FE study: solve and extract every job. Re-run to resume.',
           f'set ABQ={abaqus_command}', 'set ROOT=%~dp0', '']
    sh = ['#!/usr/bin/env bash', '# RVE Studio FE study: solve and extract every job. Re-run to resume.',
          f'ABQ="${{ABQ:-{abaqus_command}}}"', 'ROOT="$(cd "$(dirname "$0")" && pwd)"', '']
    for run in study['runs']:
        for c in _run_cases(root, run):
            rel = c['dir'].relative_to(root)
            job, ex = c['job'], c['extractor']
            win = str(rel).replace('/', '\\')
            bat += [f'if not exist "%ROOT%{win}\\case_results.json" (',
                    f'  pushd "%ROOT%{win}"',
                    f'  if not exist {job}.odb call %ABQ% job={job} input={job}.inp cpus={cpus} interactive',
                    f'  if exist {job}.odb call %ABQ% python {ex} {job}.odb case_results.json',
                    '  popd', ')']
            sh += [f'if [ ! -f "$ROOT/{rel.as_posix()}/case_results.json" ]; then',
                   f'  ( cd "$ROOT/{rel.as_posix()}" && {{ [ -f {job}.odb ] || "$ABQ" job={job} input={job}.inp cpus={cpus} interactive; }}'
                   f' && "$ABQ" python {ex} {job}.odb case_results.json )',
                   'fi']
    bat += ['', 'echo Finished. Now run: python -m rvestudio.cli fe report --root "%ROOT%."']
    sh += ['', 'echo "Finished. Now run: python -m rvestudio.cli fe report --root \\"$ROOT\\""']
    pb, ps = root / 'run_all_jobs.bat', root / 'run_all_jobs.sh'
    pb.write_text('\r\n'.join(bat) + '\r\n', encoding='utf-8')
    ps.write_text('\n'.join(sh) + '\n', encoding='utf-8')
    try:
        ps.chmod(0o755)
    except Exception:
        pass
    return [pb, ps]


# --------------------------------------------------------------------------------------
# Results
# --------------------------------------------------------------------------------------

def richardson(hs: Sequence[float], fs: Sequence[float]) -> dict:
    """Celik et al. (2008) grid-convergence estimate from the three finest meshes.

    Returns the observed order p, the extrapolated value and the fine-grid GCI
    (safety factor 1.25).  Oscillatory or non-converging triples are flagged and
    not extrapolated.
    """
    pts = sorted(zip(map(float, hs), map(float, fs)))[:3]  # finest first
    if len(pts) < 3:
        return {'status': 'need 3 meshes'}
    (h1, f1), (h2, f2), (h3, f3) = pts
    r21, r32 = h2 / h1, h3 / h2
    e21, e32 = f2 - f1, f3 - f2
    # Mesh-independent values (e.g. Poisson ratios when every stiffness scales alike, or exact antiplane
    # shear terms) differ only by round-off -- about 1e-12 after the 21-state energy subtraction.
    if max(abs(e21), abs(e32)) <= 1e-9 * max(abs(f1), 1e-300):
        return {'status': 'converged', 'p': math.nan, 'extrapolated': f1, 'gci_fine_pct': 0.0}
    change = 100.0 * max(abs(e21), abs(e32)) / max(abs(f1), 1e-300)
    ratio = e32 / e21
    if ratio <= 0:
        return {'status': 'oscillatory', 'p': math.nan, 'extrapolated': math.nan,
                'gci_fine_pct': math.nan, 'max_change_pct': change}
    s = 1.0
    p = abs(math.log(abs(ratio))) / math.log(r21)
    for _ in range(100):
        q = math.log((r21 ** p - s) / (r32 ** p - s))
        p_new = abs(math.log(abs(ratio)) + q) / math.log(r21)
        if abs(p_new - p) < 1e-10:
            p = p_new
            break
        p = p_new
    denom = r21 ** p - 1.0
    if not math.isfinite(p) or denom <= 0:
        return {'status': 'diverging', 'p': p, 'extrapolated': math.nan, 'gci_fine_pct': math.nan, 'max_change_pct': change}
    ext = (r21 ** p * f1 - f2) / denom
    gci = 1.25 * abs((f1 - f2) / f1) / denom * 100.0 if f1 else math.nan
    gci32 = 1.25 * abs((f2 - f3) / f2) / (r32 ** p - 1.0) * 100.0 if f2 and r32 ** p > 1 else math.nan
    # Roache's asymptotic-range check: GCI_32 / (r21^p GCI_21) = 1 when all three meshes are in the asymptotic range
    asym = gci32 / (r21 ** p * gci) if gci and math.isfinite(gci32) else math.nan
    return {'status': 'monotone', 'p': p, 'extrapolated': ext, 'gci_fine_pct': gci, 'gci_coarse_pct': gci32,
            'asymptotic_ratio': asym, 'max_change_pct': change}


MESH_TOLERANCE_PCT = 1.0      # default discretisation tolerance on a stiffness term (GCI of the finest mesh)


def mesh_verdict(rich: dict, f_fine: float | None = None, f_medium: float | None = None,
                 tol_pct: float = MESH_TOLERANCE_PCT) -> str:
    """Plain-language mesh-convergence verdict for one term from ``richardson`` (three meshes, r = 2).

    * converged            -- the meshes agree to round-off or differ by less than tol/100, or the triple is
                              monotone, the observed order is
                              plausible (0.5 <= p <= 4), the meshes are in the asymptotic range
                              (|GCI_32 / (r^p GCI_21) - 1| <= 0.15) and GCI_fine <= tol;
    * refine (GCI)         -- monotone and asymptotic, but GCI_fine > tol: add a finer mesh;
    * refine (not asymptotic) -- p implausible or the asymptotic check fails: the coarse mesh is too coarse;
    * converged (oscillating) / refine (oscillating) -- non-monotone triple, judged by the change between the
                              two finest meshes against tol;
    * refine (diverging)   -- the differences grow with refinement.
    """
    st = rich.get('status')
    if st == 'converged':
        return 'converged'
    # the three meshes differ by far less than the tolerance: converged whatever the observed order says
    if math.isfinite(rich.get('max_change_pct', math.nan)) and rich['max_change_pct'] <= 0.01 * tol_pct:
        return 'converged'
    if st == 'monotone':
        p, gci, asym = rich.get('p', math.nan), rich.get('gci_fine_pct', math.nan), rich.get('asymptotic_ratio', math.nan)
        if not (0.5 <= p <= 4.0) or not (math.isfinite(asym) and abs(asym - 1.0) <= 0.15):
            return 'refine (not asymptotic)'
        return 'converged' if gci <= tol_pct else 'refine (GCI)'
    if st == 'oscillatory':
        if f_fine is not None and f_medium is not None and f_fine:
            return 'converged (oscillating)' if 100.0 * abs(f_fine - f_medium) / abs(f_fine) <= tol_pct else 'refine (oscillating)'
        return 'refine (oscillating)'
    if st == 'diverging':
        return 'refine (diverging)'
    return 'need 3 meshes'


def convergence(rows: list[dict]) -> list[dict]:
    """Richardson/GCI per project, solver, model and term (needs three meshes)."""
    out = []
    keys = sorted({(r['project'], r['solver'], r['model'], r['group'], r['term']) for r in rows})
    for pid, solver, model, group, term in keys:
        sel = [r for r in rows if (r['project'], r['solver'], r['model'], r['term']) == (pid, solver, model, term)]
        if all(str(r.get('negligible', '')).lower() in ('true', '1') or r.get('negligible') is True for r in sel):
            continue  # zero to round-off: nothing to extrapolate
        pts = sorted((float(r['mesh_mm']), float(r['FE'])) for r in sel)
        rich = richardson([p[0] for p in pts], [p[1] for p in pts])
        ana = float(sel[0]['analytical'])
        row = {'project': pid, 'solver': solver, 'model': model, 'group': group, 'term': term,
               'unit': sel[0].get('unit', 'MPa'), 'meshes': len(pts), 'finest_mesh_mm': pts[0][0],
               'FE_finest': pts[0][1], 'analytical': ana,
               **{k: rich.get(k) for k in ('status', 'p', 'extrapolated', 'gci_fine_pct', 'asymptotic_ratio')}}
        row['verdict'] = mesh_verdict(rich, pts[0][1], pts[1][1] if len(pts) > 1 else None)
        ext = row['extrapolated']
        ref = ext if isinstance(ext, float) and math.isfinite(ext) and rich.get('status') == 'monotone' else pts[0][1]
        row['FE_reference'] = ref
        row['analytical_vs_FE_pct'] = pct(ana, ref)
        out.append(row)
    return out


def _collect_dbc(root: Path, study: dict) -> list[dict]:
    rows = []
    for run in study['runs']:
        if run['check'] != CHECK_DBC or not run.get('prepared'):
            continue
        rdir = root / run['run_dir']
        snap = _read_json(rdir / 'thesis_dbc_snapshot.json') or {}
        dom = snap.get('domain', {})
        energies = {}
        for c in _run_cases(root, run):
            res = _read_json(c['dir'] / 'case_results.json')
            if res and str(res.get('status', '')).upper() == 'COMPLETED':
                energies[c['case']] = float(res['strain_energy_ALLSE'])
        disp = float(snap.get('settings', {}).get('displacement', study['settings']['dbc_displacement']))
        terms = reconstruct_dbc_terms(energies, float(dom['p']), float(dom['q']), float(dom['height']), disp)
        if not terms:
            continue
        info = study['projects'][run['project']]
        rows += result_rows(run['project'], 'Abaqus', run['mesh'], {'DBC': terms}, info['analytical'])
    return rows


def _collect_pbc(root: Path, study: dict) -> tuple[list[dict], list[dict]]:
    rows, audits = [], []
    for run in study['runs']:
        if run['check'] != CHECK_PBC or not run.get('prepared'):
            continue
        rdir = root / run['run_dir']
        cases = _run_cases(root, run)
        if not cases or not all(_case_state(c['dir'], c['job'])[0] == DONE for c in cases):
            continue
        summary = _read_json(rdir / 'homogenized_results.json')
        if not summary or str(summary.get('status', '')).upper() != 'COMPLETED':
            snap = _read_json(rdir / 'homogenization_snapshot.json')
            files = _periodic_files_from_directory(rdir, snap)
            ps = PeriodicHomogenizationSettings(**run['pbc_settings'])
            # Every case already has a COMPLETED extraction, so no Abaqus call is made here.
            summary = run_periodic_homogenization(files, 'abaqus', ps, resume=True)
        info = study['projects'][run['project']]
        rows += result_rows(run['project'], 'Abaqus', run['mesh'], {'PBC': pbc_values(summary['CH'])},
                            info['analytical'])
        audits.append({'project': run['project'], 'mesh_mm': run['mesh'], 'cases': summary.get('case_count'),
                       'energy_vs_reaction_error': summary.get('energy_vs_reaction_relative_error', math.nan),
                       'reaction_reciprocity_error': summary.get('macro_reaction_symmetry_error', math.nan),
                       'positive_definite': summary.get('positive_definite')})
    return rows, audits


def _paper_table(rows: list[dict], conv: list[dict], study: dict) -> tuple[list[dict], str]:
    """Paper-ready table: Abaqus PBC (strain energy) Q11 ... Q55 (extrapolated if available) vs analytical."""
    out = []
    for pid in study['projects']:
        b = best_values(rows, pid, 'Abaqus', 'PBC', 'Q', conv)
        if not b:
            continue
        ana = study['projects'][pid]['analytical']['Q15']
        r = {'project': pid}
        for t in SHELL_TERMS:
            if t in b:
                r[f'{t}_FE'] = b[t][0]
                r[f'{t}_analytical'] = ana[t]
                r[f'{t}_diff_pct'] = pct(ana[t], b[t][0])
        out.append(r)
    lines = [r'\begin{tabular}{l' + 'rrr' * len(SHELL_TERMS) + '}', r'\hline',
             'RVE & ' + ' & '.join(rf'\multicolumn{{3}}{{c}}{{$Q_{{{t[1:]}}}$}}' for t in SHELL_TERMS) + r' \\',
             ' & ' + ' & '.join(['FE & Model & $\\Delta$\\%'] * len(SHELL_TERMS)) + r' \\', r'\hline']
    for r in out:
        cells = []
        for t in SHELL_TERMS:
            if f'{t}_FE' in r:
                cells += [f"{r[f'{t}_FE']:.2f}", f"{r[f'{t}_analytical']:.2f}", f"{r[f'{t}_diff_pct']:+.2f}"]
            else:
                cells += ['--', '--', '--']
        lines.append(r['project'].replace('_', r'\_') + ' & ' + ' & '.join(cells) + r' \\')
    lines += [r'\hline', r'\end{tabular}']
    return out, '\n'.join(lines) + '\n'


def load_results(root: str | Path) -> dict:
    """Everything the GUI and the reports show, read from the study folder."""
    root, study = load_study(root)
    res = root / 'results'
    rows = read_rows(res / BUILTIN_RESULTS) + read_rows(res / ABAQUS_RESULTS)
    return {'root': root, 'study': study, 'rows': rows, 'convergence': convergence(rows),
            'pbc_audit': read_rows(res / 'pbc_audit.csv') if (res / 'pbc_audit.csv').exists() else [],
            'status': study_status(root)}


def _refresh_report(root: Path, out: dict) -> dict:
    try:
        from .fe_report import write_study_report, write_study_workbook
        out['html_report'] = str(write_study_report(root))
        out['workbook'] = str(write_study_workbook(root))
    except Exception as exc:  # a report problem must never lose the CSV results
        out['report_warning'] = f'HTML/Excel report not written: {exc}'
    return out


def collect_study(root: str | Path) -> dict:
    """Read every finished Abaqus job and write the results, the HTML report and the workbook.

    ``results/abaqus_results.csv`` (every term, mesh and model), ``convergence.csv``,
    ``paper_table.csv/.tex``, ``pbc_audit.csv`` and ``job_status.csv``.  Pending jobs
    are listed, never filled in.
    """
    root, study = load_study(root)
    if stale_runs(root, study):
        raise ValueError(STALE_STUDY_MESSAGE)
    status = study_status(root)
    rows = _collect_dbc(root, study)
    prows, audits = _collect_pbc(root, study)
    rows += prows
    res = root / 'results'
    res.mkdir(exist_ok=True)
    _write_csv(res / 'job_status.csv', status['jobs'])
    _write_csv(res / ABAQUS_RESULTS, rows)
    _write_csv(res / 'pbc_audit.csv', audits)
    conv = convergence(rows + read_rows(res / BUILTIN_RESULTS))
    _write_csv(res / 'convergence.csv', conv)
    paper_rows, tex = _paper_table(rows, conv, study)
    _write_csv(res / 'paper_table.csv', paper_rows)
    (res / 'paper_table.tex').write_text(tex, encoding='utf-8')
    out = {'counts': status['counts'], 'total': status['total'], 'rows': len(rows),
           'dbc_rows': sum(r['model'] == 'DBC' for r in rows), 'pbc_rows': sum(r['model'] == 'PBC' for r in rows)}
    return _refresh_report(root, out)


def check_study(root: str | Path, meshes: Sequence[float] | None = None, *,
                progress: Callable[[int, int, str], None] | None = None,
                checks: Sequence[str] | None = None) -> dict:
    """Built-in shell FE of every project and mesh (no Abaqus): PBC (strain energy: full C, Q, engineering
    constants) and DBC (shell and 3D terms).

    Same FE RVE, node sets, displacements and energy formulae as the Abaqus decks.
    Meshes finer than 2.5 mm are skipped (memory).  Writes ``results/builtin_results.csv``,
    the displacement fields of the coarsest mesh (``results/builtin_fields``, for the
    deformed shapes) and refreshes the report and the workbook.
    """
    root, study = load_study(root)
    s = study['settings']
    hs = sorted({float(h) for h in (meshes or s['meshes']) if float(h) >= 2.5}, reverse=True)
    if not hs:
        raise ValueError('The built-in shell check needs at least one mesh size of 2.5 mm or coarser.')
    models = [MODEL_OF[c] for c in (checks or s['checks'])]
    res = root / 'results'
    fdir = res / 'builtin_fields'
    fdir.mkdir(parents=True, exist_ok=True)
    done = {(r['project'], float(r['mesh_mm'])) for r in read_rows(res / BUILTIN_RESULTS)}
    keep = [r for r in read_rows(res / BUILTIN_RESULTS)]
    todo = [(pid, h) for pid in study['projects'] for h in hs]
    rows = []
    for i, (pid, h) in enumerate(todo, 1):
        if progress:
            progress(i, len(todo), f'{pid} h={h:g}')
        info = study['projects'][pid]
        model, mats, lams = load_project(root / info['project'])
        want_fields = h == hs[0]
        fe = builtin_check(model, mats, lams, h, models=models, displacement=float(s['dbc_displacement']),
                           face_wall_section=s.get('dbc_face_walls', 'half-thickness'), fields=want_fields)
        rows += result_rows(pid, 'built-in', h, fe, info['analytical'])
        if want_fields:
            save_fields(fdir / f'{pid}.npz', fe.get('fields', {}), h)
    newkeys = {(r['project'], float(r['mesh_mm']), r['model']) for r in rows}
    rows = [r for r in keep if (r['project'], float(r['mesh_mm']), r['model']) not in newkeys] + rows
    rows.sort(key=lambda r: (r['project'], r['model'], -float(r['mesh_mm']), r['group'], r['term']))
    _write_csv(res / BUILTIN_RESULTS, rows)
    conv = convergence(rows + read_rows(res / ABAQUS_RESULTS))
    _write_csv(res / 'convergence.csv', conv)
    out = {'rows': len(rows), 'meshes': hs, 'file': str(res / BUILTIN_RESULTS), 'skipped_meshes':
           sorted({float(h) for h in (meshes or s['meshes'])} - set(hs)), 'previous': len(done)}
    return _refresh_report(root, out)


shell_check_study = check_study  # name used by 1.9.1-1.10


def save_fields(path: Path, fields: dict, mesh: float) -> None:
    """Compressed displacement fields (float32) of the built-in models for the deformed shapes."""
    data = {'mesh_mm': np.array(float(mesh))}
    for model, f in fields.items():
        data[f'{model}_coords'] = np.asarray(f['coords'], np.float32)
        data[f'{model}_elements'] = np.asarray(f['elements'], np.int32)
        for case, u in f['u'].items():
            data[f'{model}_u_{case}'] = np.asarray(u, np.float32)
    np.savez_compressed(path, **data)


def load_fields(path: str | Path) -> dict:
    """{'DBC': {'coords', 'elements', 'u': {case: N x 3}}, 'PBC': {...}, 'mesh_mm'} from save_fields."""
    path = Path(path)
    if not path.exists():
        return {}
    z = np.load(path)
    out: dict = {'mesh_mm': float(z['mesh_mm'])}
    for key in z.files:
        if key.endswith('_coords'):
            model = key[:-7]
            out[model] = {'coords': z[key], 'elements': z[f'{model}_elements'],
                          'u': {k[len(model) + 3:]: z[k] for k in z.files if k.startswith(f'{model}_u_')}}
    return out
