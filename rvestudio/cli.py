"""Command-line interface; also available as python -m rvestudio.cli."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from .api import benchmark, solve_example, solve_project, summary
from .database import generate_parametric_database, run_batch


def main(argv=None):
    p = argparse.ArgumentParser(prog='rvestudio', description='RVE Studio: homogenisation, verification and design of laminated-composite cellular cores')
    from . import __version__
    p.add_argument('--version', action='version', version=f'RVE Studio {__version__}')
    sub = p.add_subparsers(dest='cmd', required=True)
    b = sub.add_parser('benchmark', help='Solve the 14 benchmark designs with any wall model and compare with the reference values')
    b.add_argument('--tolerance', type=float, default=1.0, help='Allowed absolute percentage difference (default 1%%)')
    b.add_argument('--wall-beam', choices=('euler-bernoulli', 'timoshenko'), default='euler-bernoulli',
                   help='Wall beam theory: euler-bernoulli (original model, default) or timoshenko (wall shear flexibility)')
    b.add_argument('--bending', choices=('cylindrical', 'twist-corrected', 'anisotropic-strip', 'free-moment'), default='cylindrical',
                   help='Wall bending option')
    b.add_argument('--twist', type=float, nargs=3, metavar=('A', 'B', 'N'), help='Twist-release coefficients of the twist-corrected wall (default 2.0 2.15 0.73)')
    b.add_argument('--csv', type=Path, help='Export 84 paper-comparison rows as CSV')
    b.add_argument('--json', action='store_true', help='Emit summary as JSON')
    s = sub.add_parser('single', help='Solve one of the 14 benchmark designs')
    s.add_argument('--config', type=int, choices=range(1, 15), default=1)
    s.add_argument('--json', action='store_true', help='Emit machine-readable summary')
    s.add_argument('--output', type=Path, help='Save JSON summary to file')
    proj = sub.add_parser('solve-project', help='Solve an existing .rvejson/.json file')
    proj.add_argument('input', type=Path)
    proj.add_argument('--output', type=Path, help='Write JSON summary')
    proj.add_argument('--full', action='store_true', help='Include CH, Qplane and Qshear matrices')
    g = sub.add_parser('generate-db', help='Generate a compact parametric database')
    g.add_argument('output'); g.add_argument('--cases', type=int, default=5000); g.add_argument('--seed', type=int, default=20260906)
    g.add_argument('--wall-beam', choices=('euler-bernoulli', 'timoshenko'), default='euler-bernoulli',
                   help='Wall beam theory written to the WallShear column (default euler-bernoulli, the original model)')
    g.add_argument('--bending', choices=('cylindrical', 'twist-corrected', 'anisotropic-strip', 'free-moment'), default=None,
                   help='Wall bending option' + ' (default: the solver of each geometry definition)')
    g.add_argument('--twist', type=float, nargs=3, metavar=('A', 'B', 'N'), help='Twist-release coefficients of the twist-corrected wall (default 2.0 2.15 0.73)')
    r = sub.add_parser('batch', help='Solve a compact case database')
    r.add_argument('input'); r.add_argument('output'); r.add_argument('--workers', type=int, default=1); r.add_argument('--limit', type=int)
    r.add_argument('--buckling', action='store_true', help='Also calculate local-wall biaxial+shear buckling')
    r.add_argument('--buckling-stress', nargs=6, type=float, metavar=('S11','S22','S33','T23','T13','T12'), default=[-1.0, 0.0, -1.0, 0.0, 0.25, 0.0])
    r.add_argument('--buckling-m', type=int, default=8); r.add_argument('--buckling-n', type=int, default=8); r.add_argument('--buckling-sf', type=float, default=1.5)
    r.add_argument('--no-abd-condensation', action='store_true')
    fe = sub.add_parser('fe', help='FE design & validation: built-in shell FE (seconds) or Abaqus, one study folder and one report')
    fsub = fe.add_subparsers(dest='fe_action', required=True)

    def study_args(q, root_required=True):
        q.add_argument('--root', required=root_required, help='Study folder (created if needed)')
        q.add_argument('--project', action='append', default=[], help='.rvejson project (repeatable)')
        q.add_argument('--paper-configs', help='Benchmark configurations, e.g. "1-14" or "1,7,9"')
        q.add_argument('--meshes', nargs='+', type=float, help='Shell mesh sizes in mm (default 10)')
        q.add_argument('--models', nargs='+', choices=['pbc', 'dbc'],
                       help='pbc: strain energy, full C, Q and engineering constants (21 Abaqus cases per mesh); '
                            'dbc: shell and 3D terms (10 cases). Default both')
        q.add_argument('--cpus', type=int, default=1, help='CPUs per Abaqus job')
    study_args(fsub.add_parser('check', help='Built-in shell FE of every project and mesh, no Abaqus (seconds)'))
    study_args(fsub.add_parser('prepare', help='Write the Abaqus decks (no Abaqus needed)'))
    fr = fsub.add_parser('run', help='Prepare if needed, solve every pending Abaqus job (resumable) and report')
    study_args(fr)
    fr.add_argument('--abaqus-command', default='abaqus')
    fr.add_argument('--retry-failed', action='store_true')
    fs = fsub.add_parser('status', help='Abaqus job states, read from disk')
    fs.add_argument('--root', required=True)
    fs.add_argument('--all', action='store_true', help='List every job, not only the counts and problems')
    fp = fsub.add_parser('report', help='Collect finished Abaqus jobs; rewrite the HTML report and the workbook')
    fp.add_argument('--root', required=True)
    fw_ = fsub.add_parser('scripts', help='Rewrite run_all_jobs.bat/.sh (run Abaqus without Python)')
    fw_.add_argument('--root', required=True)
    fw_.add_argument('--abaqus-command', default='abaqus')
    rv = fsub.add_parser('rve', help='Show (or set) the FE RVE that every DBC and PBC model of a design uses')
    rsrc = rv.add_mutually_exclusive_group(required=True)
    rsrc.add_argument('--project', type=Path, help='RVE .rvejson project')
    rsrc.add_argument('--paper-config', type=int, choices=range(1, 15), help='Benchmark configuration 1-14')
    rv.add_argument('--png', type=Path, help='Save the plan view of the FE RVE')
    rv.add_argument('--cells', type=int, nargs=2, metavar=('N1', 'N2'),
                    help='Block size in rectangular periodic cells (centred on a complete cell)')
    rv.add_argument('--centre', type=float, nargs=2, metavar=('X1', 'X2'),
                    help='With --cells: place the block with its centre here (FE frame, x1 along a1, mm), '
                         'as a rectangle drawn over the replicated design')
    rv.add_argument('--free-centre', action='store_true',
                    help='With --centre: keep the centre exactly there (default: move it to the nearest cell centre, '
                         'i.e. the geometry centre, wherever possible)')
    rv.add_argument('--auto', action='store_true', help='Clear a design\'s own FE RVE choice (use the automatic rule)')
    rv.add_argument('--save', type=Path, metavar='RVEJSON',
                    help='Save the design with this FE RVE choice (FE studies then use it for DBC and PBC)')
    doc = fsub.add_parser('doctor', help='Check that Abaqus can be found and answers')
    doc.add_argument('--abaqus-command', default='')
    args = p.parse_args(argv)
    if args.cmd == 'fe':
        return _fe(args)
    if args.cmd == 'benchmark':
        out = benchmark(tolerance_pct=args.tolerance, wall_shear=args.wall_beam, bending=args.bending, twist_release=args.twist)
        if args.csv:
            args.csv.parent.mkdir(parents=True, exist_ok=True)
            out['detail'].to_csv(args.csv, index=False)
        if args.json:
            print(json.dumps({'reference': out['reference'], 'tolerance_pct': out['tolerance_pct'], 'passed': out['passed'],
                'within_tolerance': out['within_tolerance'], 'total': out['total'], 'max_abs_paper_pct': out['max_abs_paper_pct'],
                'wall_beam': out['wall_shear'], 'bending': out['bending'],
                'by_property': out['summary'].to_dict(orient='records')}, indent=2))
        else:
            print(f"Wall beam: {out['wall_shear']}; bending option: {out['bending']}")
            print(out['summary'].to_string(index=False))
            print(f"Benchmark designs vs reference values: {out['within_tolerance']}/{out['total']} entries within ±{out['tolerance_pct']:g}%; max difference {out['max_abs_paper_pct']:.4f}%; all pass={out['passed']}")
            if args.csv: print(f'Exported: {args.csv}')
    elif args.cmd in ('single', 'solve-project'):
        result = solve_example(args.config) if args.cmd == 'single' else solve_project(args.input)
        data = summary(result, include_ch=(args.cmd == 'solve-project' and args.full))
        payload = json.dumps(data, indent=2, allow_nan=False)
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(payload + '\n', encoding='utf-8')
            print(f'Saved: {args.output}')
        elif getattr(args, 'json', False) or args.cmd == 'solve-project':
            print(payload)
        else:
            print('CH (internal [11,22,33,23,13,12]):\n', result.CH)
            print('\nQplane [11,22,12]:\n', result.props.Qplane)
            print('\nQshear [13,23], paper Q44/Q55:\n', result.props.Qshear)
            rel = data['reliability']
            print(f"\nReliability: in-plane {rel['in_plane']}, out-of-plane {rel['out_of_plane']}")
            for m in rel['messages']:
                print('  - ' + m)
    elif args.cmd == 'generate-db':
        df = generate_parametric_database(args.output, args.cases, args.seed, wall_shear=args.wall_beam, bending=args.bending, twist_release=args.twist)
        print(f'Generated {len(df)} cases -> {args.output}')
    elif args.cmd == 'batch':
        buckling = None
        if args.buckling:
            buckling = {'control_mode':'stress','macro_values':args.buckling_stress, 'm_max':args.buckling_m, 'n_max':args.buckling_n,
                'safety_factor':args.buckling_sf, 'use_abd_condensation':not args.no_abd_condensation, 'convergence_check':True}
        df = run_batch(args.input, args.output, args.workers, args.limit,
            lambda i,n,row: print(f'[{i}/{n}] {row.get("CaseID")} {row.get("Status")}'), buckling_settings=buckling)
        print(df['Status'].value_counts())

def _print_table(cols, table):
    from .verification import signed
    def cell(c, x):
        if not isinstance(x, float):
            return str(x)
        return signed(x) if '(%)' in c else f'{x:.4g}'
    t = [[cell(c, x) for c, x in zip(cols, r)] for r in table]
    w = [max(len(c), *(len(r[k]) for r in t)) for k, c in enumerate(cols)]
    print('  '.join(c.rjust(n) for c, n in zip(cols, w)))
    print('  '.join('-' * n for n in w))
    for r in t:
        print('  '.join(x.rjust(n) for x, n in zip(r, w)))


def _print_matrix(title, labels, M, diff=False):
    from .verification import signed
    print(title)
    print('       ' + ''.join(f'{x:>11}' for x in labels))
    for i, r in enumerate(labels):
        cells = []
        for j in range(len(labels)):
            v = float(M[i][j])
            cells.append(f"{(signed(v) if diff else f'{v:.5g}') if v == v else ('≈0' if diff else '–'):>11}")
        print(f'{r:>7}' + ''.join(cells))


def _print_results(root):
    """Stiffness matrices, term tables (PBC: C, Q, engineering constants; DBC) and findings of every project."""
    from .fe_workflow import load_results
    from .verification import summary_table, findings, GROUP_TITLES, GROUPS, matrix_comparison
    R = load_results(root)
    for pid, info in R['study']['projects'].items():
        ana = info['analytical']
        print(f"\n=== {pid}" + (f" (benchmark configuration {ana['paper_config']})" if ana.get('paper_config') else ''))
        for g in ('C', 'Q'):
            mc = matrix_comparison(R['rows'], pid, g, ana, R['convergence'])
            if mc is not None:
                name = 'C (6x6, Voigt 11 22 33 23 13 12)' if g == 'C' else 'Q (sigma33 = 0; 11 22 12 13 23)'
                print()
                _print_matrix(f'{name}, analytical, MPa', mc['labels'], mc['analytical'])
                _print_matrix(f"{name}, FE ({mc['solver']} PBC, {mc['how']}), MPa", mc['labels'], mc['FE'])
                _print_matrix(f'{name}, analytical vs FE (%)', mc['labels'], mc['difference_pct'], diff=True)
        for g in GROUPS:
            cols, table, notes = summary_table(R['rows'], pid, g, ana, R['convergence'])
            if len(cols) > 3 and table:
                print(f"\n{GROUP_TITLES[g]} ({'; '.join(f'{k}: {v}' for k, v in notes.items())})")
                _print_table(cols, table)
        for f in findings(R['rows'], pid, ana, None, R['convergence']):
            print('-', f)
        conv = [c for c in R['convergence'] if c['project'] == pid and int(c.get('meshes', 0)) >= 3]
        if conv:
            from collections import Counter
            n = Counter(str(c.get('verdict', '')) for c in conv)
            print('- mesh convergence: ' + ', '.join(f'{v} {k}' for k, v in sorted(n.items())))
            worst = [c for c in conv if str(c.get('verdict', '')).startswith('refine')]
            for c in sorted(worst, key=lambda c: -float(c.get('gci_fine_pct') or 0))[:6]:
                print(f"    {c['solver']} {c['model']} {c['term']}: {c['verdict']}, p = {float(c.get('p') or float('nan')):.2f}, "
                      f"GCI = {float(c.get('gci_fine_pct') or float('nan')):.2f} %, asymptotic ratio = {float(c.get('asymptotic_ratio') or float('nan')):.2f}")
    print('\n(%) = 100 (analytical - FE) / FE; terms zero to round-off are not listed.')
    print('Mesh verdict: converged = monotone, 0.5 <= p <= 4, asymptotic ratio 1 +/- 0.15 and GCI_fine <= 1 %; '
          'otherwise add a finer mesh (halve the finest size).')


def _fe(args) -> int:
    from . import fe_workflow as fw
    a = args.fe_action
    if a == 'rve':
        from .io import load_project
        from .examples_data import paper_hexagonal_case
        from .fe_rve import select_fe_rve, fe_rve_summary, normalise_fe_rve_spec, fe_rve_layout, centred_window
        model, _mats, _lams = load_project(args.project) if args.project else paper_hexagonal_case(args.paper_config)
        if args.centre and not args.cells:
            raise SystemExit('--centre needs --cells N1 N2 (the size of the drawn rectangle).')
        if args.auto:
            model.fe_rve = None
        elif args.cells:
            if args.centre:
                model.fe_rve, note = centred_window(model, fe_rve_layout(model), args.cells, args.centre,
                                                    snap=not args.free_centre)
                print(note)
            else:
                model.fe_rve = normalise_fe_rve_spec({'mode': 'cells', 'cells': list(args.cells)})
        _dm, d = select_fe_rve(model)
        print(fe_rve_summary(d))
        if args.save:
            from .io import save_project
            save_project(args.save, model, _mats, _lams)
            print('Saved', args.save, 'with FE RVE:', d['rve'])
        print('Cell types:', ', '.join(f"{c['name']} ({c['area']:.6g} mm2)" for c in d['cell_types']))
        if args.png:
            from matplotlib.figure import Figure
            from .fe_plots import plot_rve_plan
            fig = Figure(figsize=(6, 5.5)); plot_rve_plan(fig.add_subplot(111), d, title=str(model.name))
            fig.savefig(args.png, dpi=160, bbox_inches='tight'); print('Wrote', args.png)
        return 0
    if a == 'doctor':
        rep = fw.abaqus_doctor(args.abaqus_command)
        print(json.dumps(rep, indent=2))
        return 0 if rep['responds'] else 1
    if a == 'status':
        stt = fw.study_status(args.root)
        print('Jobs:', stt['total'], ' '.join(f'{k}={v}' for k, v in sorted(stt['counts'].items())))
        for j in stt['jobs']:
            if args.all or j['state'] not in (fw.DONE, fw.PENDING):
                print(f"  {j['state']:<12} {j['check']} {j['project']} h={j['mesh']:g} {j['case']:<8} {j['note']}")
        return 0
    if a == 'scripts':
        for f in fw.write_run_scripts(args.root, args.abaqus_command):
            print('Wrote', f)
        return 0
    if a == 'report':
        out = fw.collect_study(args.root)
        _print_results(args.root)
        print('\nReport:', out.get('html_report', out.get('report_warning')), '\nWorkbook:', out.get('workbook', ''))
        return 0
    # check / prepare / run: create or extend the study first
    exists = (Path(args.root) / fw.STUDY_FILE).exists()
    old = fw.load_study(args.root)[1]['settings'] if exists else {}
    settings = fw.FEStudySettings(meshes=args.meshes or old.get('meshes') or [10.0],
                                  checks=args.models or old.get('checks') or ['dbc', 'pbc'],
                                  cpus=args.cpus if args.cpus != 1 else int(old.get('cpus', 1)))
    if args.project or args.paper_configs or not exists or args.meshes or args.models:
        if a == 'check':
            fw.add_projects(args.root, args.project, paper_configs=args.paper_configs, settings=settings, progress=print)
        else:
            study = fw.prepare_study(args.root, args.project, paper_configs=args.paper_configs, settings=settings,
                                     progress=print)
            bad = [r for r in study['runs'] if not r['prepared']]
            for r in bad:
                print(f"  NOT PREPARED {r['check']} {r['project']} h={r['mesh']:g}: {r['error']}")
    if a == 'check':
        out = fw.check_study(args.root, args.meshes,
                             progress=lambda i, n, t: print(f'[{i}/{n}] built-in shell FE {t}', flush=True))
        _print_results(args.root)
        print('\nReport:', out.get('html_report', out.get('report_warning')), '\nWorkbook:', out.get('workbook', ''))
        return 0
    if a == 'prepare':
        stt = fw.study_status(args.root)
        print(f"Abaqus decks ready: {stt['total']} jobs. Run them with: rvestudio fe run --root \"{args.root}\" "
              f"--abaqus-command abaqus (or run_all_jobs.bat / .sh in the study folder).")
        return 0
    out = fw.run_study(args.root, args.abaqus_command, retry_failed=args.retry_failed,
                       progress=lambda i, n, j: print(f"[{i}/{n}] {j['check']} {j['project']} h={j['mesh']:g} {j['case']}",
                                                      flush=True))
    print(json.dumps(out, indent=2))
    _print_results(args.root)
    print('\nReport:', Path(args.root).resolve() / 'FE_STUDY_REPORT.html')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
