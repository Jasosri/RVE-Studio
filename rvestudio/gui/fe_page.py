"""FE Design & Validation page: design the FE RVE and mesh, then validate the analytical model against shell FE.

Left, top to bottom: what to verify (projects, benchmark configurations, meshes,
models, study folder), **Validate now** (built-in shell FE, seconds, no licence) and
**Run with Abaqus** (prepare, run and collect in one go, resumable).  Right: the
results of the selected project -- the full stiffness matrix C, the Q matrix
(sigma33 = 0) and the engineering constants from the strain-energy PBC model next to
the analytical ones, the DBC terms, charts, its FE RVE, deformed shapes -- the Abaqus
jobs and a log.  Every FE value is compared with the analytical model only.

The page keeps no results of its own: everything shown is re-read from the study
folder (:mod:`rvestudio.fe_workflow`), so it survives restarts and runs done elsewhere.
"""
from __future__ import annotations

import math
import os
import tempfile
from pathlib import Path

from PySide6.QtCore import Qt, QThread, Signal, QTimer, QUrl
from PySide6.QtGui import QColor, QDesktopServices
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel, QPushButton, QLineEdit, QListWidget, QCheckBox,
    QComboBox, QSpinBox, QGroupBox, QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView,
    QProgressBar, QFileDialog, QMessageBox, QTabWidget, QSplitter, QPlainTextEdit, QScrollArea,
)

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from .. import fe_workflow as fw
from ..abaqus import detect_abaqus_executable
from ..io import save_project, load_project
from ..fe_rve import select_fe_rve, fe_rve_summary, describe_fe_rve_spec, normalise_fe_rve_spec
from ..fe_plots import plot_rve_plan, project_figure, deformed_figure
from ..fe_report import project_rve, abaqus_dbc_fields
from ..verification import (signed, summary_table, findings, GROUP_TITLES, BOUNDARY_CONDITIONS,
                            matrix_comparison)

DIFF_COLOURS = ((1.0, '#E5F4E4'), (5.0, '#F9EEDA'), (float('inf'), '#FBEBEB'))
STATE_COLOURS = {fw.DONE: '#D8F0DF', fw.FAILED: '#F8D7D7', 'NOT PREPARED': '#F8D7D7',
                 fw.RUNNING: '#D7E6FA', fw.SOLVED: '#FFF2CC', fw.PENDING: '#FFFFFF'}


def diff_colour(value) -> QColor | None:
    try:
        v = abs(float(value))
    except (TypeError, ValueError):
        return None
    if not math.isfinite(v):
        return None
    for lim, col in DIFF_COLOURS:
        if v <= lim:
            return QColor(col)
    return None


class _Worker(QThread):
    """Runs one study function off the GUI thread."""
    progress = Signal(int, int, str)
    message = Signal(str)
    done = Signal(object)
    failed = Signal(str)

    def __init__(self, fn, *args, **kwargs):
        super().__init__()
        self.fn, self.args, self.kwargs = fn, args, kwargs
        self.stop_requested = False

    def run(self):
        try:
            self.done.emit(self.fn(*self.args, **self.kwargs))
        except Exception as exc:  # shown to the user, never swallowed
            self.failed.emit(f'{type(exc).__name__}: {exc}')


def _heading(text: str, sub: str = '') -> QWidget:
    from .main_window import title_block   # same page header as every other workspace
    return title_block(text, sub)


def _row(*widgets, stretch_first: bool = False) -> QWidget:
    w = QWidget()
    lay = QHBoxLayout(w)
    lay.setContentsMargins(0, 0, 0, 0)
    for k, x in enumerate(widgets):
        lay.addWidget(x, 1 if (stretch_first and k == 0) else 0)
    if not stretch_first:
        lay.addStretch()
    return w


def _button(text: str, slot, *, primary: bool = False, tip: str = '') -> QPushButton:
    b = QPushButton(text)
    if primary:
        b.setObjectName('Primary')
        b.setMinimumHeight(38)
    if tip:
        b.setToolTip(tip)
    b.clicked.connect(slot)
    return b


def _table() -> QTableWidget:
    t = QTableWidget(0, 0)
    t.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
    t.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
    t.verticalHeader().setVisible(False)
    return t


def fill_table(table: QTableWidget, cols: list[str], rows: list[list]) -> None:
    """Numbers to 4 significant figures, '(%)' columns signed and colour-coded."""
    table.clear()
    table.setColumnCount(len(cols))
    table.setRowCount(len(rows))
    table.setHorizontalHeaderLabels(cols)
    for r, line in enumerate(rows):
        for c, val in enumerate(line):
            pc = cols[c].endswith('(%)')
            if isinstance(val, (int, float)) and not isinstance(val, bool):
                txt = (signed(val) if pc else f'{val:.4g}') if math.isfinite(float(val)) else '–'
            else:
                txt = str(val)
            it = QTableWidgetItem(txt)
            if c > 0:
                it.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            if pc:
                col = diff_colour(val)
                if col is not None:
                    it.setBackground(col)
            table.setItem(r, c, it)


def fill_matrix(table: QTableWidget, labels, M, *, diff: bool = False) -> None:
    """A square matrix with labelled rows and columns; difference matrices signed and colour-coded."""
    n = len(labels)
    table.clear()
    table.setRowCount(n)
    table.setColumnCount(n)
    table.setHorizontalHeaderLabels(list(labels))
    table.setVerticalHeaderLabels(list(labels))
    for i in range(n):
        for j in range(n):
            v = float(M[i][j])
            if diff:
                txt = signed(v) if math.isfinite(v) else '≈0'
            else:
                txt = f'{v:.5g}' if math.isfinite(v) else '–'
            it = QTableWidgetItem(txt)
            it.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            if diff:
                col = diff_colour(v)
                if col is not None:
                    it.setBackground(col)
                else:
                    it.setForeground(QColor('#8A94A6'))
            table.setItem(i, j, it)


class FEPage(QWidget):
    def __init__(self, main):
        super().__init__()
        self.main = main
        self.worker: _Worker | None = None
        self._results: dict | None = None
        outer = QVBoxLayout(self)
        outer.setContentsMargins(18, 15, 18, 15)
        outer.addWidget(_heading(
            'FE Design & Validation',
            'Design the FE model (FE RVE, mesh sizes, boundary conditions) and validate the analytical model against shell FE on the FE RVE of each design (automatic, or chosen or '
            'drawn with Define FE RVE…): in seconds with the built-in solver, or with Abaqus. The PBC model gives '
            'the full stiffness matrix C from strain energy, the Q matrix (σ33 = 0) and the engineering constants; '
            'the DBC model gives the shell terms (U3 free) and 3D terms (U3 = 0). Everything is compared with the '
            'analytical model.'))
        split = QSplitter(Qt.Orientation.Horizontal)
        outer.addWidget(split, 1)

        # ------------------------------------------------------------------ left
        inner = QWidget()
        left = QVBoxLayout(inner)
        left.setContentsMargins(0, 0, 8, 0)

        g1 = QGroupBox('1  Design the FE model: designs, FE RVE, meshes')
        g1.setStyleSheet('QGroupBox{font-weight:700;}')
        f1 = QFormLayout(g1)
        self.project_list = QListWidget()
        self.project_list.setMaximumHeight(90)
        self.project_list.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        f1.addRow('Designs to add', self.project_list)
        f1.addRow('', _row(_button('Current design', self.add_current, tip='The design open in the Single RVE Designer'),
                           _button('Project files…', self.add_files), _button('Remove', self.remove_selected)))
        f1.addRow('FE RVE', _row(_button('Define FE RVE…', self.define_rve,
                                         tip='Choose the rectangular FE RVE of the selected design: automatic, a block of '
                                             'n1 x n2 periodic cells, or a rectangle drawn over the replicated design. '
                                             'Every DBC and PBC model of that design uses it.'),
                                 _button('Automatic', self.reset_rve, tip='Use the automatic FE RVE')))
        self.paper_configs = QLineEdit('')
        self.paper_configs.setPlaceholderText('benchmark configurations, e.g. 1-14 or 1,7,9')
        f1.addRow('Benchmarks', self.paper_configs)
        self.meshes = QLineEdit('10')
        self.meshes.setToolTip('Shell element size in mm. Three meshes halving in size (10, 5, 2.5) give a '
                               'Richardson extrapolation; the built-in check skips meshes finer than 2.5 mm.')
        f1.addRow('Mesh sizes (mm)', self.meshes)
        self.use_pbc = QCheckBox('PBC, strain energy: full C, Q (σ33 = 0), engineering constants')
        self.use_pbc.setChecked(True)
        self.use_pbc.setToolTip(fw.CHECK_DESCRIPTIONS[fw.CHECK_PBC])
        self.use_dbc = QCheckBox('DBC: shell terms (U3 free) and 3D terms (U3 = 0)')
        self.use_dbc.setChecked(True)
        self.use_dbc.setToolTip(fw.CHECK_DESCRIPTIONS[fw.CHECK_DBC])
        f1.addRow('Models', self.use_pbc)
        f1.addRow('', self.use_dbc)
        self.root_edit = QLineEdit(str((Path.home() / 'RVE_Studio_FE' / 'Study_1').resolve()))
        self.root_edit.editingFinished.connect(self.reload)
        f1.addRow('Study folder', _row(self.root_edit, _button('Browse…', self.choose_root), stretch_first=True))
        self.study_label = QLabel('')
        self.study_label.setWordWrap(True)
        self.study_label.setStyleSheet('color:#5B6B82;')
        f1.addRow('', self.study_label)
        left.addWidget(g1)

        g2 = QGroupBox('2  Validate now (built-in shell FE, no Abaqus)')
        g2.setStyleSheet('QGroupBox{font-weight:700;}')
        f2 = QVBoxLayout(g2)
        self.check_btn = _button('Validate now', self.check_now, primary=True,
                                 tip='Built-in four-node laminated shell on the same FE RVE, node sets and boundary '
                                     'conditions as the Abaqus decks. Seconds per design.')
        f2.addWidget(self.check_btn)
        left.addWidget(g2)

        g3 = QGroupBox('3  Validate with Abaqus')
        g3.setStyleSheet('QGroupBox{font-weight:700;}')
        f3 = QFormLayout(g3)
        self.command = QLineEdit(detect_abaqus_executable() or 'abaqus')
        f3.addRow('Abaqus command', _row(self.command, _button('Detect', self.detect), _button('Test', self.doctor),
                                         stretch_first=True))
        self.cpus = QSpinBox()
        self.cpus.setRange(1, 128)
        self.cpus.setValue(max(1, min(4, os.cpu_count() or 1)))
        self.retry = QCheckBox('Retry failed jobs')
        f3.addRow('CPUs per job', _row(self.cpus, self.retry))
        self.run_btn = _button('Run with Abaqus', self.run_abaqus, primary=True,
                               tip='Writes the decks (if needed), solves every pending job one at a time and collects. '
                                   'Safe to stop and run again: finished jobs are kept.')
        self.stop_btn = _button('Stop after current job', self.stop)
        self.stop_btn.setEnabled(False)
        f3.addRow(_row(self.run_btn, self.stop_btn, stretch_first=True))
        f3.addRow(_button('Scripts for another machine…', self.write_scripts,
                          tip='run_all_jobs.bat / .sh solve every job without Python; copy the folder back and press '
                              'Reload results.'))
        left.addWidget(g3)

        self.progress = QProgressBar()
        self.progress.setRange(0, 1)
        self.progress.setValue(0)
        left.addWidget(self.progress)
        self.status = QLabel('Ready.')
        self.status.setWordWrap(True)
        left.addWidget(self.status)
        left.addWidget(_row(_button('Open report', self.open_report), _button('Open Excel', self.open_excel),
                            _button('Open folder', lambda: self._open(self.root())),
                            _button('Reload results', self.collect_and_reload)))
        left.addStretch()
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(inner)
        split.addWidget(scroll)

        # ------------------------------------------------------------------ right
        right = QWidget()
        rl = QVBoxLayout(right)
        rl.setContentsMargins(8, 0, 0, 0)
        self.project_combo = QComboBox()
        self.project_combo.currentIndexChanged.connect(self.show_project)
        rl.addWidget(_row(QLabel('Project'), self.project_combo))
        self.tabs = QTabWidget()
        rl.addWidget(self.tabs, 1)

        res = QWidget()
        rv = QVBoxLayout(res)
        self.findings = QLabel('Press Validate now, or Run with Abaqus, to see results here.')
        self.findings.setWordWrap(True)
        self.findings.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        rv.addWidget(self.findings)
        self.result_tabs = QTabWidget()
        rv.addWidget(self.result_tabs, 1)

        # Stiffness matrices: analytical | FE | difference, for C (6x6) and Q (5x5).
        mw = QWidget()
        ml = QVBoxLayout(mw)
        self.matrix_tables = {}
        self.matrix_titles = {}
        for g, name in (("C", "C, 6×6 (Voigt 11, 22, 33, 23, 13, 12), MPa"),
                        ("Q", "Q, σ33 = 0 (11, 22, 12, 13, 23), MPa")):
            t = QLabel(name)
            t.setStyleSheet('font-weight:600;')
            ml.addWidget(t)
            self.matrix_titles[g] = t
            row = QHBoxLayout()
            tabs3 = []
            for what in ('Analytical', 'FE', 'Analytical vs FE (%)'):
                box = QVBoxLayout()
                box.addWidget(QLabel(what))
                tb = _table()
                tb.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
                tb.verticalHeader().setVisible(True)
                box.addWidget(tb, 1)
                row.addLayout(box, 1)
                tabs3.append(tb)
            ml.addLayout(row, 1)
            self.matrix_tables[g] = tabs3
        self.result_tabs.addTab(mw, 'Stiffness matrices')

        # Term tables: C, Q and engineering constants (PBC); shell and 3D terms (DBC).
        self.term_tables = {}
        self.term_titles = {}
        for tab_name, groups in (('C, Q and engineering constants', ('C', 'Q', 'eng')),
                                 ('DBC terms', ('shell', '3d'))):
            tw = QWidget()
            tl = QVBoxLayout(tw)
            for g in groups:
                t = QLabel(GROUP_TITLES[g])
                t.setStyleSheet('font-weight:600;')
                t.setWordWrap(True)
                tl.addWidget(t)
                tb = _table()
                tl.addWidget(tb, 1)
                self.term_titles[g], self.term_tables[g] = t, tb
            sc = QScrollArea()
            sc.setWidgetResizable(True)
            sc.setWidget(tw)
            self.result_tabs.addTab(sc, tab_name)
        # Names kept for scripts and tests written against 1.11-1.12.
        self.shell_table, self.solid_table = self.term_tables['shell'], self.term_tables['3d']
        self.shell_title, self.solid_title = self.term_titles['shell'], self.term_titles['3d']
        note = QLabel('Differences in %: 100 (analytical − FE) / FE; green ≤ 1 %, amber ≤ 5 %, red > 5 %; terms that '
                      'are zero to round-off show ≈0. FE values are Richardson-extrapolated when three meshes converge, '
                      'else the finest mesh. C from strain energy; Q = C condensed to σ33 = 0; engineering constants '
                      'from S = C⁻¹.')
        note.setWordWrap(True)
        note.setStyleSheet('color:#5B6B82;')
        rv.addWidget(note)
        self.tabs.addTab(res, 'Results')

        ch = QWidget()
        cv = QVBoxLayout(ch)
        self.chart_fig = Figure(figsize=(11, 4.6))
        self.chart_canvas = FigureCanvas(self.chart_fig)
        cv.addWidget(self.chart_canvas, 1)
        self.tabs.addTab(ch, 'Charts')

        rvw = QWidget()
        rvl = QVBoxLayout(rvw)
        self.rve_fig = Figure(figsize=(6, 5))
        self.rve_canvas = FigureCanvas(self.rve_fig)
        rvl.addWidget(self.rve_canvas, 1)
        self.rve_label = QLabel('Every DBC and PBC model of a design uses this FE RVE.')
        self.rve_label.setWordWrap(True)
        rvl.addWidget(self.rve_label)
        bc = _table()
        fill_table(bc, ['model', 'side faces', 'top and bottom', 'terms'], [list(r) for r in BOUNDARY_CONDITIONS])
        bc.setMaximumHeight(130)
        rvl.addWidget(bc)
        self.tabs.addTab(rvw, 'FE RVE')

        dfw = QWidget()
        dfl = QVBoxLayout(dfw)
        self.def_source = QComboBox()
        for label, key in (('Built-in PBC (6 unit strains)', ('built-in', 'PBC')), ('Built-in DBC (10 cases)', ('built-in', 'DBC')),
                           ('Abaqus DBC', ('Abaqus', 'DBC'))):
            self.def_source.addItem(label, key)
        self.def_source.currentIndexChanged.connect(self.show_deformed)
        dfl.addWidget(_row(QLabel('Show'), self.def_source))
        self.def_fig = Figure(figsize=(11, 6))
        self.def_canvas = FigureCanvas(self.def_fig)
        dfl.addWidget(self.def_canvas, 1)
        self.tabs.addTab(dfw, 'Deformed shapes')

        jobs = QWidget()
        jl = QVBoxLayout(jobs)
        self.counts_label = QLabel('No Abaqus jobs prepared.')
        self.counts_label.setStyleSheet('font-weight:600;')
        self.job_filter = QComboBox()
        self.job_filter.addItems(['All', 'Problems', 'Pending', 'Done'])
        self.job_filter.currentIndexChanged.connect(self.refresh_jobs)
        jl.addWidget(_row(self.counts_label, QLabel('Show'), self.job_filter, stretch_first=True))
        self.jobs_table = _table()
        self.jobs_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.jobs_table.cellDoubleClicked.connect(self.open_job_folder)
        jl.addWidget(self.jobs_table, 1)
        hint = QLabel('Double-click a job to open its folder (input deck, log, .sta, .msg).')
        hint.setStyleSheet('color:#5B6B82;')
        jl.addWidget(hint)
        self.tabs.addTab(jobs, 'Abaqus jobs')

        lg = QWidget()
        ll = QVBoxLayout(lg)
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        ll.addWidget(self.log)
        self.tabs.addTab(lg, 'Log')
        split.addWidget(right)
        split.setSizes([470, 870])

        self.timer = QTimer(self)
        self.timer.setInterval(5000)
        self.timer.timeout.connect(self._tick)
        self.timer.start()
        self.reload()

    # ------------------------------------------------------------------ helpers
    def root(self) -> Path:
        return Path(self.root_edit.text().strip()).expanduser()

    def _say(self, text: str) -> None:
        self.log.appendPlainText(text)

    def _open(self, path: Path) -> None:
        if not Path(path).exists():
            QMessageBox.information(self, 'Not found', f'{path} does not exist yet.')
            return
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))

    def _busy(self, busy: bool, name: str = '') -> None:
        self.check_btn.setEnabled(not busy)
        self.run_btn.setEnabled(not busy)
        self.stop_btn.setEnabled(busy and name == 'run')

    def _tick(self) -> None:
        if self.isVisible() and self.worker is None and (self.root() / fw.STUDY_FILE).exists():
            self.refresh_jobs()

    def _inputs(self) -> dict:
        meshes = [float(p) for p in self.meshes.text().replace(';', ',').replace(' ', ',').split(',') if p.strip()]
        models = [c for c, box in ((fw.CHECK_DBC, self.use_dbc), (fw.CHECK_PBC, self.use_pbc)) if box.isChecked()]
        settings = fw.FEStudySettings(meshes=meshes, checks=models, cpus=self.cpus.value())
        settings.validate()
        configs = self.paper_configs.text().strip() or None
        if configs:
            fw.parse_config_list(configs)
        projects = [self._project_path(self.project_list.item(i)) for i in range(self.project_list.count())]
        if not projects and not configs and not (self.root() / fw.STUDY_FILE).exists():
            raise ValueError('Add a design (Current design or Project files…) or benchmark configurations.')
        return {'settings': settings, 'configs': configs, 'projects': projects}

    # ------------------------------------------------------------------ design list (path + FE RVE choice)
    PATH_ROLE = Qt.ItemDataRole.UserRole
    SPEC_ROLE = Qt.ItemDataRole.UserRole + 1
    FILE_SPEC_ROLE = Qt.ItemDataRole.UserRole + 2

    def _add_project_item(self, path: str, spec=None, *, from_file: bool = True) -> None:
        from PySide6.QtWidgets import QListWidgetItem
        if from_file:
            try:
                spec = load_project(path)[0].fe_rve
            except Exception:
                spec = None
        item = QListWidgetItem()
        item.setData(self.PATH_ROLE, str(path))
        item.setData(self.FILE_SPEC_ROLE, spec)
        self._set_item_spec(item, spec)
        self.project_list.addItem(item)

    def _set_item_spec(self, item, spec) -> None:
        item.setData(self.SPEC_ROLE, spec)
        path = Path(item.data(self.PATH_ROLE))
        item.setText(f'{path.name}   —   FE RVE: {describe_fe_rve_spec(spec)}')
        item.setToolTip(str(path))

    def _project_path(self, item) -> str:
        """The file to add to the study: the design itself, or a copy carrying a changed FE RVE choice."""
        path = item.data(self.PATH_ROLE) or item.text()
        spec = normalise_fe_rve_spec(item.data(self.SPEC_ROLE))
        if spec == normalise_fe_rve_spec(item.data(self.FILE_SPEC_ROLE)):
            return str(path)
        model, mats, lams = load_project(path)
        model.fe_rve = spec
        tag = 'auto' if spec is None else ('rve' + ('drawn' if spec['mode'] == 'window' else '')
                                           + f"{spec['cells'][0]}x{spec['cells'][1]}")
        out = Path(tempfile.mkdtemp(prefix='rve_fe_')) / f'{Path(path).stem}_{tag}.rvejson'
        save_project(out, model, mats, lams)
        return str(out)

    def _target_item(self):
        items = self.project_list.selectedItems()
        if not items and self.project_list.count() == 1:
            items = [self.project_list.item(0)]
        if len(items) != 1:
            QMessageBox.information(self, 'Define FE RVE', 'Select one design in the list (add it with Current design '
                                    'or Project files… first).')
            return None
        return items[0]

    def define_rve(self) -> None:
        item = self._target_item()
        if item is None:
            return
        try:
            model = load_project(item.data(self.PATH_ROLE) or item.text())[0]
            model.fe_rve = item.data(self.SPEC_ROLE)
            from .fe_rve_dialog import FERVEDialog
            dlg = FERVEDialog(model, self, title=Path(item.data(self.PATH_ROLE)).stem)
        except Exception as exc:
            QMessageBox.warning(self, 'Define FE RVE', f'Cannot build the FE RVE of this design:\n{exc}')
            return
        if dlg.exec():
            self._set_item_spec(item, dlg.spec())
            self._say(f'FE RVE of {Path(item.data(self.PATH_ROLE)).name}: {describe_fe_rve_spec(dlg.spec())}')

    def reset_rve(self) -> None:
        item = self._target_item()
        if item is not None:
            self._set_item_spec(item, None)

    # ------------------------------------------------------------------ inputs
    def add_current(self) -> None:
        try:
            model, mats, lams = self.main.single_page.parse_model()
        except Exception as exc:
            QMessageBox.warning(self, 'No valid design', f'The Single RVE Designer model is not valid:\n{exc}')
            return
        stem = ''.join(ch if ch.isalnum() or ch in '-_' else '_' for ch in (model.name or 'current_design')).strip('_')
        tmp = Path(tempfile.mkdtemp(prefix='rve_fe_')) / f'{stem or "current_design"}.rvejson'
        save_project(tmp, model, mats, lams)
        self._add_project_item(str(tmp), None, from_file=False)

    def add_files(self) -> None:
        files, _ = QFileDialog.getOpenFileNames(self, 'Add RVE projects', '', 'RVE project (*.rvejson *.json)')
        have = {self.project_list.item(i).data(self.PATH_ROLE) for i in range(self.project_list.count())}
        for f in files:
            if f not in have:
                self._add_project_item(f)

    def remove_selected(self) -> None:
        for item in self.project_list.selectedItems():
            self.project_list.takeItem(self.project_list.row(item))

    def choose_root(self) -> None:
        d = QFileDialog.getExistingDirectory(self, 'Study folder (new or existing)', str(self.root().parent))
        if d:
            self.root_edit.setText(d)
            self.reload()

    # ------------------------------------------------------------------ actions
    def check_now(self) -> None:
        try:
            x = self._inputs()
        except Exception as exc:
            QMessageBox.warning(self, 'Check the inputs', str(exc))
            return

        root = self.root()  # read the widgets here, never from the worker thread

        def job(progress):
            fw.add_projects(root, x['projects'], paper_configs=x['configs'], settings=x['settings'])
            return fw.check_study(root, x['settings'].meshes, progress=progress, checks=x['settings'].checks)
        self._say(f'Built-in check in {root} …')
        self._start('check', job)

    def run_abaqus(self) -> None:
        try:
            x = self._inputs()
        except Exception as exc:
            QMessageBox.warning(self, 'Check the inputs', str(exc))
            return
        cmd = self.command.text().strip()
        if not cmd:
            QMessageBox.warning(self, 'Abaqus command', 'Enter the Abaqus command (for example abaqus, or the full path to abaqus.bat).')
            return

        root, retry = self.root(), self.retry.isChecked()  # read the widgets here, never from the worker thread

        def job(progress, should_stop):
            fw.prepare_study(root, x['projects'], paper_configs=x['configs'], settings=x['settings'])
            return fw.run_study(root, cmd, retry_failed=retry, progress=progress, should_stop=should_stop)
        self._say(f'Abaqus study in {root} with "{cmd}" …')
        self._start('run', job)

    def stop(self) -> None:
        if self.worker is not None:
            self.worker.stop_requested = True
            self.stop_btn.setEnabled(False)
            self.status.setText('Stopping after the current job …')

    def detect(self) -> None:
        found = detect_abaqus_executable()
        if found:
            self.command.setText(found)
        else:
            QMessageBox.information(self, 'Abaqus not found', 'No Abaqus launcher was found on PATH. Enter the full path, '
                                    'for example C:\\SIMULIA\\Commands\\abaqus.bat.')

    def doctor(self) -> None:
        rep = fw.abaqus_doctor(self.command.text())
        text = '\n'.join(rep['messages']) + (('\n\n' + rep['release']) if rep.get('release') else '')
        self._say('Abaqus test:\n' + text)
        (QMessageBox.information if rep['responds'] else QMessageBox.warning)(self, 'Abaqus test', text or 'No output.')

    def write_scripts(self) -> None:
        try:
            x = self._inputs()
            fw.prepare_study(self.root(), x['projects'], paper_configs=x['configs'], settings=x['settings'])
            paths = fw.write_run_scripts(self.root(), self.command.text() or 'abaqus')
        except Exception as exc:
            QMessageBox.warning(self, 'Scripts', str(exc))
            return
        self.project_list.clear()
        self.reload()
        QMessageBox.information(self, 'Scripts written', 'Copy the whole study folder to the Abaqus machine and run:\n\n'
                                + '\n'.join(str(p) for p in paths)
                                + '\n\nThen copy it back, open it here and press Reload results.')

    def collect_and_reload(self) -> None:
        if (self.root() / fw.STUDY_FILE).exists():
            try:
                fw.collect_study(self.root())
            except Exception as exc:
                self._say(f'Reload: {exc}')
        self.reload()

    def open_report(self) -> None:
        self._open(self.root() / 'FE_STUDY_REPORT.html')

    def open_excel(self) -> None:
        self._open(self.root() / 'FE_STUDY_RESULTS.xlsx')

    # ------------------------------------------------------------------ worker plumbing
    def _start(self, name: str, fn) -> None:
        w = _Worker(lambda: None)
        w.setObjectName(name)
        if name == 'run':
            w.fn = lambda: fn(lambda i, n, j: w.progress.emit(i, n, f"Abaqus {j['check'].upper()} {j['project']} "
                                                                     f"h={j['mesh']:g} {j['case']}"),
                              lambda: w.stop_requested)
        else:
            w.fn = lambda: fn(lambda i, n, t: w.progress.emit(i, n, 'built-in ' + t))
        w.progress.connect(self._on_progress)
        w.message.connect(self._say)
        w.done.connect(lambda out, n=name: self._on_done(n, out))
        w.failed.connect(lambda msg, n=name: self._on_failed(n, msg))
        self.worker = w
        self._busy(True, name)
        self.progress.setRange(0, 0)
        w.start()

    def _on_progress(self, i: int, n: int, text: str) -> None:
        self.progress.setRange(0, n)
        self.progress.setValue(i - 1)
        self.status.setText(f'{i}/{n}: {text}')
        self._say(f'[{i}/{n}] {text}')

    def _on_done(self, name: str, out) -> None:
        self.worker = None
        self._busy(False)
        self.progress.setRange(0, 1)
        self.progress.setValue(1)
        if name == 'check':
            msg = (f"Built-in check finished: {out['rows']} results, meshes {', '.join(f'{h:g}' for h in out['meshes'])} mm."
                   + (f" Skipped (finer than 2.5 mm): {', '.join(f'{h:g}' for h in out['skipped_meshes'])}."
                      if out.get('skipped_meshes') else ''))
        else:
            msg = f"Abaqus: {out['completed']} jobs completed, {out['failed']} failed in this run."
        if isinstance(out, dict) and out.get('report_warning'):
            msg += ' ' + out['report_warning']
        self.status.setText(msg)
        self._say(msg)
        self.project_list.clear()
        self.reload()
        self.tabs.setCurrentIndex(0)

    def _on_failed(self, name: str, msg: str) -> None:
        self.worker = None
        self._busy(False)
        self.progress.setRange(0, 1)
        self.progress.setValue(0)
        self.status.setText(f'{"Check" if name == "check" else "Abaqus run"} failed: {msg}')
        self._say(self.status.text())
        QMessageBox.critical(self, 'FE design & validation', msg)

    # ------------------------------------------------------------------ views
    def reload(self) -> None:
        """Re-read the study folder and refresh every view."""
        try:
            self._results = fw.load_results(self.root())
        except Exception:
            self._results = None
            self.study_label.setText('New study: it is created when you press Validate now or Run with Abaqus.')
            self.project_combo.blockSignals(True)
            self.project_combo.clear()
            self.project_combo.blockSignals(False)
            self.refresh_jobs()
            return
        study = self._results['study']
        s = study['settings']
        n_built = len({r['project'] for r in self._results['rows'] if r['solver'] == 'built-in'})
        n_abq = len({r['project'] for r in self._results['rows'] if r['solver'] == 'Abaqus'})
        self.study_label.setText(
            f"{len(study['projects'])} designs, meshes {', '.join(f'{h:g}' for h in s['meshes'])} mm, models "
            f"{', '.join(c.upper() for c in s['checks'])}; results: built-in {n_built}, Abaqus {n_abq}.")
        current = self.project_combo.currentText()
        self.project_combo.blockSignals(True)
        self.project_combo.clear()
        self.project_combo.addItems(list(study['projects']))
        k = self.project_combo.findText(current)
        self.project_combo.setCurrentIndex(k if k >= 0 else 0)
        self.project_combo.blockSignals(False)
        self.show_project()
        self.refresh_jobs()

    def _pid(self) -> str:
        return self.project_combo.currentText()

    def show_project(self, *_):
        R = self._results
        pid = self._pid()
        if not R or not pid or pid not in R['study']['projects']:
            return
        ana = R['study']['projects'][pid]['analytical']
        tables = {}
        for g, table in self.term_tables.items():
            cols, rows, notes = summary_table(R['rows'], pid, g, ana, R['convergence'])
            fill_table(table, cols, rows)
            have = len(cols) > 3
            self.term_titles[g].setText(GROUP_TITLES[g] + (f"  —  FE values: {'; '.join(f'{k} {v}' for k, v in notes.items())}"
                                                          if have and notes else ('' if have else '  —  no FE results yet')))
            if have and rows:
                tables[g] = (cols, rows)
        for g, (ta, tf, td) in self.matrix_tables.items():
            mc = matrix_comparison(R['rows'], pid, g, ana, R['convergence'])
            if mc is None:
                for t in (ta, tf, td):
                    t.clear(); t.setRowCount(0); t.setColumnCount(0)
                continue
            fill_matrix(ta, mc['labels'], mc['analytical'])
            fill_matrix(tf, mc['labels'], mc['FE'])
            fill_matrix(td, mc['labels'], mc['difference_pct'], diff=True)
            self.matrix_titles[g].setText(self.matrix_titles[g].text().split('  —')[0]
                                          + f"  —  FE: {mc['solver']} PBC, {mc['how']}")
        rve = project_rve(R['root'], R['study'], pid)
        lines = findings(R['rows'], pid, ana, None, R['convergence'])
        head = f"<b>{pid}</b>" + (f" = benchmark configuration {ana['paper_config']}" if ana.get('paper_config') else '')
        self.findings.setText(head + '<ul>' + ''.join(f'<li>{t}</li>' for t in lines) + '</ul>')
        project_figure(rve, tables, self.chart_fig)
        self.chart_canvas.draw_idle()
        self.rve_fig.clear()
        if rve:
            plot_rve_plan(self.rve_fig.add_subplot(111), rve, title=pid)
            self.rve_label.setText(fe_rve_summary(rve))
        self.rve_canvas.draw_idle()
        self.show_deformed()

    def show_deformed(self, *_):
        R = self._results
        pid = self._pid()
        if not R or not pid:
            return
        solver, model = self.def_source.currentData() or ('built-in', 'DBC')
        fields = None
        if solver == 'built-in':
            fields = fw.load_fields(R['root'] / 'results' / 'builtin_fields' / f'{pid}.npz').get(model)
        else:
            for run in R['study'].get('runs', []):
                if run.get('project') == pid and run.get('check') == fw.CHECK_DBC and run.get('run_dir'):
                    fields = abaqus_dbc_fields(R['root'] / run['run_dir'])
                    if fields:
                        break
        deformed_figure(fields, model, self.def_fig, solver=solver)
        self.def_canvas.draw_idle()

    def refresh_jobs(self, *_):
        try:
            st = fw.study_status(self.root())
        except Exception:
            self.jobs_table.setRowCount(0)
            self.counts_label.setText('No Abaqus jobs prepared.')
            return
        c = st['counts']
        self.counts_label.setText(f"{st['total']} Abaqus jobs: " + (', '.join(f'{k.lower()} {v}' for k, v in sorted(c.items())) or 'none'))
        keep = {'All': None, 'Problems': {fw.FAILED, 'NOT PREPARED'}, 'Pending': {fw.PENDING, fw.SOLVED, fw.RUNNING},
                'Done': {fw.DONE}}[self.job_filter.currentText() or 'All']
        rows = [j for j in st['jobs'] if keep is None or j['state'] in keep]
        cols = ['Model', 'Project', 'Mesh (mm)', 'Case', 'State', 'Note']
        self.jobs_table.setColumnCount(len(cols))
        self.jobs_table.setHorizontalHeaderLabels(cols)
        self.jobs_table.setRowCount(len(rows))
        for r, j in enumerate(rows):
            for col, v in enumerate([j['check'].upper(), j['project'], f"{j['mesh']:g}", j['case'], j['state'], j['note']]):
                it = QTableWidgetItem(str(v))
                it.setData(Qt.ItemDataRole.UserRole, j['dir'])
                if col == 4:
                    it.setBackground(QColor(STATE_COLOURS.get(j['state'], '#FFFFFF')))
                self.jobs_table.setItem(r, col, it)

    def open_job_folder(self, row: int, _col: int) -> None:
        it = self.jobs_table.item(row, 0)
        rel = it.data(Qt.ItemDataRole.UserRole) if it else ''
        if rel:
            self._open(self.root() / rel)
