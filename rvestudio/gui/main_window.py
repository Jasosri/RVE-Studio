from __future__ import annotations

import html, json, math, os, re, sys, traceback
from pathlib import Path
import numpy as np
import pandas as pd

from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtGui import QAction, QColor, QDesktopServices, QKeySequence
from PySide6.QtWidgets import (
    QApplication,QMainWindow,QWidget,QVBoxLayout,QHBoxLayout,QGridLayout,QFormLayout,QLabel,QPushButton,
    QListWidget,QListWidgetItem,QStackedWidget,QFrame,QTableWidget,QTableWidgetItem,QComboBox,QDoubleSpinBox,
    QSpinBox,QLineEdit,QFileDialog,QMessageBox,QTabWidget,QPlainTextEdit,QProgressBar,QHeaderView,QCheckBox,
    QSplitter,QAbstractItemView,QGroupBox,QInputDialog,QDialog,QDialogButtonBox,QScrollArea
)
from PySide6.QtCore import QUrl
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from ..library import PACKAGE_ROOT,ensure_library_dirs,load_material_library,load_geometry_library,load_layup_library,save_json,material_mapping
from ..models import Material,Ply,Laminate,Element,RVEModel
from ..examples_data import paper_hexagonal_case
from ..homogenization import homogenize_rve
from ..reliability import assess as assess_reliability
from ..load_recovery import recover_wall_loads, MACRO_LABELS
from ..buckling import evaluate_wall_buckling
from ..laminate import analyze_laminate
from ..database import generate_parametric_database,run_batch,_read_table
from ..validation import run_all14_benchmark,PROPERTY_NAMES,PAPER_URL
from ..io import save_project,load_project
from ..geometry_analysis import analyse_geometry, summary_text
from ..parametric_geometry import make_scaled_xy_plugin_source
from ..plugins import discover_geometry_plugins
from ..postprocess import (enrich_results,fit_power_law,pareto_mask_min_x_max_y,directional_youngs_modulus,merge_input_metadata,
    ObjectiveSpec,ConstraintSpec,ManufacturingRules,optimize_batch_designs,robust_optimize_batch_designs,
    sensitivity_analysis,design_variable_columns,manufacturing_feasible_mask)
from ..geometry import BUILTIN_GEOMETRIES
from .. import wall_axes
from .sketch import RVESketchView
from ..guides import CITATIONS, GUIDES, guide_path
from .ui_helpers import TableEditorDialog, WallModelOptions, fit_text, hero_banner, picture_tile
from .geometry_wizard import GeometryDefinitionWizard, ParameterEditorDialog
from .fe_page import FEPage

COLORS={
    'nav':'#132238','nav2':'#1B2D49','accent':'#2F80ED','accent2':'#56CCF2','bg':'#F4F7FB','card':'#FFFFFF',
    'text':'#172033','muted':'#6B7280','green':'#1F9D72','orange':'#F2994A','red':'#D64545','purple':'#7A5AF8'
}

APP_QSS=f"""
QMainWindow {{ background:{COLORS['bg']}; }}
QWidget {{ font-family:'Segoe UI'; font-size:10pt; color:{COLORS['text']}; }}
QScrollArea, QScrollArea > QWidget > QWidget {{ background:{COLORS['bg']}; }}
QListWidget#Nav {{ background:qlineargradient(x1:0,y1:0,x2:0,y2:1,stop:0 {COLORS['nav']},stop:1 #0D1828); border:0; color:#DCE6F3; padding:12px 10px; outline:0; }}
QListWidget#Nav::item {{ height:40px; border-radius:8px; padding-left:8px; margin:1px 0; }}
QListWidget#Nav::item:selected {{ background:qlineargradient(x1:0,y1:0,x2:1,y2:0,stop:0 {COLORS['accent']},stop:1 #56A3F2); color:white; font-weight:600; }}
QListWidget#Nav::item:hover:!selected {{ background:{COLORS['nav2']}; }}
QFrame#Card, QGroupBox {{ background:white; border:1px solid #E3E8EF; border-radius:10px; }}
QGroupBox {{ margin-top:14px; padding:14px 8px 8px 8px; font-weight:600; }}
QGroupBox::title {{ subcontrol-origin:margin; subcontrol-position:top left; left:12px; padding:2px 8px; color:white;
                    background:qlineargradient(x1:0,y1:0,x2:1,y2:0,stop:0 {COLORS['accent']},stop:1 #56A3F2); border-radius:6px; }}
QFrame#PageHeader {{ background:qlineargradient(x1:0,y1:0,x2:1,y2:0,stop:0 #FFFFFF,stop:1 #EEF4FC); border:1px solid #E3E8EF;
                     border-left:6px solid {COLORS['accent']}; border-radius:10px; }}
QPushButton {{ background:#EEF3F9; border:1px solid #D7E0EA; border-radius:7px; padding:7px 14px; min-height:20px; }}
QPushButton:checked {{ background:#D6E6FB; border:1px solid {COLORS['accent']}; color:{COLORS['accent']}; font-weight:600; }}
QPushButton:hover {{ background:#E2EBF5; border-color:#B9CBE0; }}
QPushButton:pressed {{ background:#D3E0EF; }}
QPushButton:disabled {{ color:#9AA5B5; background:#F2F4F7; border-color:#E3E8EF; }}
QPushButton#Primary {{ background:qlineargradient(x1:0,y1:0,x2:1,y2:0,stop:0 {COLORS['accent']},stop:1 #4C95F0); color:white; border:0; font-weight:600; }}
QPushButton#Primary:hover {{ background:#1F6FD8; }}
QPushButton#Primary:disabled {{ background:#A9C7F2; color:white; }}
QPushButton#Success {{ background:{COLORS['green']}; color:white; border:0; font-weight:600; }}
QLineEdit,QSpinBox,QDoubleSpinBox,QComboBox,QPlainTextEdit,QTextEdit {{ background:white; border:1px solid #D7E0EA; border-radius:6px; padding:4px 6px; }}
QLineEdit:focus,QSpinBox:focus,QDoubleSpinBox:focus,QComboBox:focus,QPlainTextEdit:focus {{ border:1px solid {COLORS['accent']}; }}
QComboBox:hover,QLineEdit:hover,QSpinBox:hover,QDoubleSpinBox:hover {{ border-color:#B9CBE0; }}
QTableWidget,QTableView {{ background:white; border:1px solid #D7E0EA; border-radius:8px; gridline-color:#EDF1F6;
                           alternate-background-color:#F7FAFD; selection-background-color:#D6E6FB; selection-color:{COLORS['text']}; }}
QHeaderView::section {{ background:#EAF0F7; border:0; border-right:1px solid #D7E0EA; border-bottom:2px solid {COLORS['accent']}; padding:6px; font-weight:600; color:#24324A; }}
QTableCornerButton::section {{ background:#EAF0F7; border:0; }}
QTabWidget::pane {{ border:1px solid #D7E0EA; background:white; border-radius:8px; top:-1px; }}
QTabBar::tab {{ padding:8px 18px; background:#EAF0F7; border:1px solid #D7E0EA; border-bottom:0; border-top-left-radius:8px; border-top-right-radius:8px; margin-right:3px; color:#43506A; }}
QTabBar::tab:hover:!selected {{ background:#DDE7F3; }}
QTabBar::tab:selected {{ background:white; color:{COLORS['accent']}; font-weight:600; border-top:3px solid {COLORS['accent']}; }}
QProgressBar {{ border:1px solid #D7E0EA; border-radius:6px; text-align:center; background:white; }}
QProgressBar::chunk {{ background:qlineargradient(x1:0,y1:0,x2:1,y2:0,stop:0 {COLORS['green']},stop:1 {COLORS['accent']}); border-radius:5px; }}
QLabel#Banner {{ font-size:11pt; font-weight:600; padding:8px 12px; border-radius:8px; background:#EEF3F9; }}
QLabel#Banner[tone="good"] {{ background:#E5F6EF; color:#16704F; border:1px solid #BFE6D6; }}
QLabel#Banner[tone="warn"] {{ background:#FDF0E4; color:#A25A14; border:1px solid #F6D3B0; }}
QLabel#Banner[tone="info"] {{ background:#EAF2FE; color:#1F4E9A; border:1px solid #C9DCFA; }}
QScrollBar:vertical {{ background:transparent; width:10px; margin:2px; }}
QScrollBar::handle:vertical {{ background:#C5D1DF; border-radius:4px; min-height:30px; }}
QScrollBar:horizontal {{ background:transparent; height:10px; margin:2px; }}
QScrollBar::handle:horizontal {{ background:#C5D1DF; border-radius:4px; min-width:30px; }}
QScrollBar::add-line, QScrollBar::sub-line {{ width:0; height:0; }}
QSplitter::handle {{ background:#E3E8EF; }}
QStatusBar {{ background:white; border-top:1px solid #E3E8EF; color:{COLORS['muted']}; }}
QMenuBar {{ background:white; border-bottom:1px solid #E3E8EF; }}
QMenuBar::item:selected, QMenu::item:selected {{ background:#D6E6FB; color:{COLORS['text']}; }}
QToolTip {{ background:{COLORS['nav']}; color:white; border:0; padding:6px; border-radius:6px; }}
"""

# Accent colour of every workspace (navigation icon and page header bar).
PAGE_COLOURS={'Dashboard':'#2F80ED','Single RVE Designer':'#2F80ED','Laminate Tools':'#1F9D72','Material Library':'#0891B2',
              'Geometry Library':'#7A5AF8','Generate DB':'#7A5AF8','Generate Database':'#7A5AF8','Large Batch Analysis':'#F2994A','Results Explorer':'#F2994A',
              'Benchmark Designs':'#1F9D72','FE Design & Validation':'#D64545','Help':'#6B7280','About / Workflow':'#6B7280'}


def nav_icon(name:str):
    """Small rounded badge with the initials of a workspace, in the workspace colour."""
    from PySide6.QtGui import QIcon, QPainter, QPixmap, QFont
    pm=QPixmap(44,44); pm.fill(Qt.GlobalColor.transparent); p=QPainter(pm); p.setRenderHint(QPainter.RenderHint.Antialiasing)
    p.setBrush(QColor(PAGE_COLOURS.get(name,COLORS['accent']))); p.setPen(Qt.PenStyle.NoPen); p.drawRoundedRect(2,2,40,40,10,10)
    f=QFont('Segoe UI'); f.setPixelSize(17); f.setBold(True); p.setFont(f); p.setPen(QColor('white'))
    letters=''.join(w[0] for w in re.findall(r'[A-Za-z]+',name) if w[0].isupper())[:2] or name[:1]
    p.drawText(pm.rect(),Qt.AlignmentFlag.AlignCenter,letters); p.end(); return QIcon(pm)


def title_block(title:str,subtitle:str='') -> QWidget:
    """Page header: title and subtitle on a light gradient card with the workspace's accent bar."""
    w=QWidget(); outer=QVBoxLayout(w); outer.setContentsMargins(0,0,0,8)
    f=QFrame(); f.setObjectName('PageHeader'); colour=PAGE_COLOURS.get(title.split(' — ')[0],COLORS['accent'])
    f.setStyleSheet(f"QFrame#PageHeader{{border-left:6px solid {colour};}}")
    l=QVBoxLayout(f); l.setContentsMargins(16,10,16,10); l.setSpacing(2)
    t=QLabel(title); t.setStyleSheet(f'font-size:19pt;font-weight:700;color:{colour};background:transparent;border:0;'); l.addWidget(t)
    if subtitle:
        s=QLabel(subtitle); s.setWordWrap(True); s.setStyleSheet(f"color:{COLORS['muted']};font-size:10pt;background:transparent;border:0;"); l.addWidget(s)
    outer.addWidget(f); return w


def card(title:str,value:str,color:str|None=None)->QFrame:
    f=QFrame(); f.setObjectName('Card'); l=QVBoxLayout(f); a=QLabel(title.upper()); a.setStyleSheet(f"color:{COLORS['muted']};font-size:9pt;font-weight:600;"); b=QLabel(value); b.setStyleSheet(f"font-size:20pt;font-weight:700;color:{color or COLORS['text']};"); l.addWidget(a); l.addWidget(b); l.addStretch(); return f


def table_item(v):
    if v is None: return QTableWidgetItem('')
    if isinstance(v,float) and math.isnan(v): return QTableWidgetItem('')
    return QTableWidgetItem(str(v))


def matrix_text(M:np.ndarray,fmt='.6g')->str:
    return '\n'.join('   '.join(format(float(v),fmt) for v in row) for row in np.asarray(M))


def cell_text(table:QTableWidget,row:int,col:int)->str:
    widget=table.cellWidget(row,col)
    if isinstance(widget,QComboBox):
        return widget.currentText().strip()
    item=table.item(row,col)
    return item.text().strip() if item else ''


def set_combo_cell(table:QTableWidget,row:int,col:int,values,current=None,*,editable=False,callback=None)->QComboBox:
    combo=QComboBox(); combo.setEditable(editable)
    texts=[str(v) for v in values]
    combo.addItems(texts)
    cur='' if current is None else str(current)
    if cur and combo.findText(cur)<0 and editable:
        combo.addItem(cur)
    if cur:
        idx=combo.findText(cur)
        if idx>=0: combo.setCurrentIndex(idx)
        elif editable: combo.setEditText(cur)
    combo.activated.connect(lambda *_args,t=table,r=row,c=col:t.setCurrentCell(r,c))
    if callback is not None:
        combo.currentTextChanged.connect(callback)
    table.setCellWidget(row,col,combo)
    return combo


def fill_matrix_table(table:QTableWidget,matrix:np.ndarray,labels:list[str]|None=None,fmt='.8g'):
    arr=np.asarray(matrix,float); table.setRowCount(arr.shape[0]); table.setColumnCount(arr.shape[1])
    if labels and len(labels)==arr.shape[0]: table.setVerticalHeaderLabels(labels)
    if labels and len(labels)==arr.shape[1]: table.setHorizontalHeaderLabels(labels)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            item=table_item(format(float(arr[i,j]),fmt)); item.setTextAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignVCenter); table.setItem(i,j,item)
    table.resizeColumnsToContents()


class BatchThread(QThread):
    progress=Signal(int,int,str,str); done=Signal(object); failed=Signal(str)
    def __init__(self,input_file,output_file,workers,limit,buckling_settings=None): super().__init__(); self.input_file=input_file; self.output_file=output_file; self.workers=workers; self.limit=limit; self.buckling_settings=buckling_settings
    def run(self):
        try:
            df=run_batch(self.input_file,self.output_file,self.workers,self.limit,lambda i,n,r:self.progress.emit(i,n,str(r.get('CaseID','')),str(r.get('Status',''))),buckling_settings=self.buckling_settings)
            self.done.emit(df)
        except Exception: self.failed.emit(traceback.format_exc())


class DashboardPage(QWidget):
    TILES=[('rve3d.png','Periodic RVE designer','Sketch any periodic wall graph and solve the full 6 × 6 Cᴴ, Q matrix and engineering constants in milliseconds.','Open RVE designer',COLORS['accent'],('Single RVE Designer',None)),
           ('laminate.png','Laminate design','Ply stacks, A, B, D and S_eq, lamination parameters and the bend–twist ratio ψ of every wall.','Open laminate tools',COLORS['green'],('Laminate Tools',None)),
           ('deformed.png','FE design & validation','Periodic shell FE of the same RVE: full C, Q and engineering constants from strain energy, built-in or Abaqus.','Open FE design & validation',COLORS['red'],('FE Design & Validation',None)),
           ('database.png','Design databases','Seeded parametric studies and batch analysis of thousands of designs, with a reliability rating for each.','Generate a database',COLORS['purple'],('Generate DB',None)),
           ('pareto.png','Multi-objective optimisation','Pareto fronts, weighted and robust screening, manufacturing rules and buckling constraints.','Explore results',COLORS['orange'],('Results Explorer',None)),
           ('buckling.png','Wall loads & buckling','Recover the membrane resultants of every wall and screen wall buckling under a design load.','Open wall buckling','#0891B2',('Single RVE Designer',5))]

    def __init__(self,main):
        super().__init__(); self.main=main
        # No scroll area: the hero, the six tiles and the footer share the window height; the pictures shrink.
        v=QVBoxLayout(self); v.setContentsMargins(22,16,22,14); v.setSpacing(12)
        mats=load_material_library(); geoms=load_geometry_library(); lays=load_layup_library()
        chips=[('materials',str(sum(bool(x.get('active',True)) for x in mats)),'rgba(47,128,237,0.55)'),
               ('geometries',str(sum(bool(x.get('active',True)) for x in geoms)),'rgba(122,90,248,0.55)'),
               ('layups',str(sum(bool(x.get('active',True)) for x in lays)),'rgba(31,157,114,0.60)')]
        v.addWidget(hero_banner('RVE Studio — Composite cellular core design','Strain-energy homogenisation of laminated-composite cellular cores: periodic RVE design, laminate analysis, FE design and validation, design databases and multi-objective screening.',
                                [('Start designing an RVE',lambda:self._open('Single RVE Designer')),('FE design & validation',lambda:self._open('FE Design & Validation')),('Help',lambda:self._open('Help'))],chips),0)
        grid=QGridLayout(); grid.setSpacing(14)
        for k,(img,title,text,btn,col,(page,tab)) in enumerate(self.TILES):
            grid.addWidget(picture_tile(img,title,text,btn,col,lambda _=False,p=page,t=tab:self._open(p,t)),k//3,k%3)
        for c in range(3): grid.setColumnStretch(c,1)
        for r in range(2): grid.setRowStretch(r,1)
        v.addLayout(grid,1)
        lab=QLabel('<b>Research basis.</b> Sriharan et al., Engineering Structures 299 (2024) 117106 (Q44 = 1–3 shear, Q55 = 2–3 shear). '
                   'Every solve reports whether the wall model is inside its verified range (bend–twist index Π, χ_B, wall shear parameter β); '
                   'the 14 benchmark designs are under Benchmark Designs.')
        lab.setWordWrap(True); lab.setStyleSheet(f"color:{COLORS['muted']};font-size:9pt;"); v.addWidget(lab,0)

    def _open(self,page,tab=None):
        self.main.go_to(page)
        if page=='Single RVE Designer':   # a tile always opens its own tab (the designer tab unless another is named)
            self.main.single_page.tabs.setCurrentIndex(0 if tab is None else tab)
            self.main.single_page.tabs.currentWidget().setFocus()


class SingleRVEPage(QWidget):
    def __init__(self,main):
        super().__init__(); self.main=main; self.materials=material_mapping(); self.laminates=[]; self.model=None; self.result=None; self._syncing=False; self._history=[]; self._redo=[]; self._history_suspended=False
        v=QVBoxLayout(self); v.setContentsMargins(18,15,18,15); v.addWidget(title_block('Single RVE Designer','Sketch arbitrary periodic RVEs or edit numerical definitions. Wall lengths, orientations and lattice dimensions update automatically.'))
        self.tabs=QTabWidget(); v.addWidget(self.tabs,1); self._geometry_tab(); self._laminate_tab(); self._solver_tab(); self._results_tab(); self._wall_load_tab(); self._buckling_tab(); self.load_example_configuration(1)

    def _geometry_tab(self):
        w=QWidget(); l=QVBoxLayout(w); l.setContentsMargins(6,6,6,6)
        bar=QHBoxLayout()
        bar.addWidget(QLabel('Example'))
        self.example_config=QComboBox(); self.example_config.addItems([f'Configuration {i}' for i in range(1,15)]); bar.addWidget(self.example_config)
        load_example=QPushButton('Load example configuration'); load_example.clicked.connect(lambda:self.load_example_configuration(self.example_config.currentIndex()+1)); bar.addWidget(load_example)
        b=QPushButton('New sketch'); b.clicked.connect(self.new_sketch); bar.addWidget(b)
        b=QPushButton('Validate geometry'); b.clicked.connect(self.validate_geometry_action); bar.addWidget(b)
        hb=QPushButton('Help'); hb.clicked.connect(lambda:self.main.go_to('Help')); bar.addWidget(hb); bar.addStretch(); l.addLayout(bar)

        # Direct access to every geometry definition in the Geometry Library.
        # Parameters are shown as editable spin boxes so loading a geometry is
        # explicit rather than silently using an unknown size.
        gbar=QHBoxLayout(); gbar.addWidget(QLabel('Existing geometry'))
        self.geometry_selector=QComboBox(); self.geometry_selector.setMinimumWidth(230); gbar.addWidget(self.geometry_selector)
        self.geometry_parameter_labels=[]; self.geometry_parameter_spins=[]
        for _ in range(4):
            lab=QLabel(); spin=QDoubleSpinBox(); spin.setDecimals(6); spin.setRange(-1e12,1e12); spin.setMinimumWidth(105); lab.hide(); spin.hide()
            self.geometry_parameter_labels.append(lab); self.geometry_parameter_spins.append(spin); gbar.addWidget(lab); gbar.addWidget(spin)
        load_geom=QPushButton('Load selected geometry'); load_geom.clicked.connect(self.load_selected_geometry); gbar.addWidget(load_geom)
        refresh_geom=QPushButton('Refresh list'); refresh_geom.clicked.connect(self.refresh_geometry_choices); gbar.addWidget(refresh_geom)
        gbar.addStretch(); l.addLayout(gbar)
        self.geometry_selector.currentIndexChanged.connect(self._update_geometry_parameter_controls)
        self.refresh_geometry_choices()

        split=QSplitter(Qt.Orientation.Horizontal)
        left=QWidget(); ll=QVBoxLayout(left); ll.setContentsMargins(0,0,4,0)
        tools=QGroupBox('RVE sketch tools'); tg=QGridLayout(tools)
        self.mode_buttons={}
        for k,(label,mode) in enumerate([('Select / move','select'),('Node','node'),('Wall','wall'),('Delete','delete'),('Set a1','a1'),('Set a2','a2')]):
            b=QPushButton(label); b.setCheckable(True); b.clicked.connect(lambda checked,m=mode:self.set_sketch_mode(m)); self.mode_buttons[mode]=b; tg.addWidget(b,k//3,k%3)
        self.mode_buttons['select'].setChecked(True)
        self.grid_size=QDoubleSpinBox(); self.grid_size.setRange(1e-6,1e6); self.grid_size.setDecimals(4); self.grid_size.setValue(5.0); self.grid_size.setSuffix(' mm'); tg.addWidget(QLabel('Grid'),2,0); tg.addWidget(self.grid_size,2,1)
        self.snap=QCheckBox('Snap to grid'); self.snap.setChecked(True); tg.addWidget(self.snap,2,2)
        fit=QPushButton('Fit drawing'); fit.clicked.connect(lambda:self.sketch.fit_geometry()); tg.addWidget(fit,3,0)
        sync=QPushButton('Refresh from tables'); sync.clicked.connect(self.refresh_sketch_from_tables); tg.addWidget(sync,3,1,1,2)
        undo=QPushButton('Undo'); undo.clicked.connect(self.undo_sketch); tg.addWidget(undo,4,0)
        redo=QPushButton('Redo'); redo.clicked.connect(self.redo_sketch); tg.addWidget(redo,4,1)
        self.show_dims=QCheckBox('Show dimensions'); self.show_dims.toggled.connect(lambda checked:self.sketch.set_show_dimensions(checked) if hasattr(self,'sketch') else None); tg.addWidget(self.show_dims,4,2)
        self.show_dir=QCheckBox('Wall direction I → J (local 1̄)'); self.show_dir.setChecked(True); self.show_dir.setToolTip('Arrowheads point from Node I to Node J: the positive local 1̄ axis of each wall.'); self.show_dir.toggled.connect(lambda checked:self.sketch.set_show_wall_direction(checked) if hasattr(self,'sketch') else None); tg.addWidget(self.show_dir,6,0,1,2)
        wed=QPushButton('Wall editor…'); wed.clicked.connect(lambda:self.open_wall_editor()); tg.addWidget(wed,6,2)
        hcon=QPushButton('Horizontal'); hcon.clicked.connect(lambda:self.apply_wall_constraint('horizontal')); tg.addWidget(hcon,5,0)
        vcon=QPushButton('Vertical'); vcon.clicked.connect(lambda:self.apply_wall_constraint('vertical')); tg.addWidget(vcon,5,1)
        lcon=QPushButton('Set length…'); lcon.clicked.connect(lambda:self.apply_wall_constraint('length')); tg.addWidget(lcon,5,2)
        ll.addWidget(tools)
        self.sketch=RVESketchView(); self.sketch.setMinimumSize(520,380); ll.addWidget(self.sketch,1)
        self.sketch_status=QLabel('Select/drag nodes. Geometry values update automatically.'); self.sketch_status.setWordWrap(True); self.sketch_status.setStyleSheet(f"color:{COLORS['muted']};"); ll.addWidget(self.sketch_status)
        split.addWidget(left)

        right=QWidget(); rl=QVBoxLayout(right); rl.setContentsMargins(4,0,0,0)
        ib=QGroupBox('Geometry inspector'); il=QVBoxLayout(ib); self.geometry_state=QLabel('Geometry not checked'); self.geometry_state.setStyleSheet('font-weight:700;'); il.addWidget(self.geometry_state)
        self.detected_family=QLabel('Family: Custom / arbitrary periodic RVE'); self.detected_family.setWordWrap(True); il.addWidget(self.detected_family)
        self.detected_params=QTableWidget(0,3); self.detected_params.setHorizontalHeaderLabels(['Parameter','Value','Source walls']); self.detected_params.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch); self.detected_params.setMaximumHeight(135); self.detected_params.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers); self.detected_params.itemSelectionChanged.connect(self.highlight_detected_parameter); il.addWidget(self.detected_params)
        self.geometry_info=QPlainTextEdit(); self.geometry_info.setReadOnly(True); self.geometry_info.setMinimumWidth(370); il.addWidget(self.geometry_info,1); rl.addWidget(ib,1)
        f=QGroupBox('Periodic RVE unit cell'); form=QGridLayout(f); self.a1x=QDoubleSpinBox(); self.a1y=QDoubleSpinBox(); self.a2x=QDoubleSpinBox(); self.a2y=QDoubleSpinBox(); self.height=QDoubleSpinBox()
        for x in [self.a1x,self.a1y,self.a2x,self.a2y,self.height]: x.setRange(-1e6,1e6); x.setDecimals(6)
        self.height.setMinimum(1e-9); self.height.setSuffix(' mm')
        form.addWidget(QLabel('a1x'),0,0); form.addWidget(self.a1x,0,1); form.addWidget(QLabel('a1y'),0,2); form.addWidget(self.a1y,0,3); form.addWidget(QLabel('a2x'),1,0); form.addWidget(self.a2x,1,1); form.addWidget(QLabel('a2y'),1,2); form.addWidget(self.a2y,1,3); form.addWidget(QLabel('Core height h'),2,0); form.addWidget(self.height,2,1,1,3)
        pp=QPushButton('Populate detected periodic pairs'); pp.clicked.connect(self.populate_detected_pairs); form.addWidget(pp,3,0,1,4); rl.addWidget(f)
        split.addWidget(right); split.setSizes([900,440])
        # Drawing gets the full height; the Node/Wall/periodic-pair tables open in a resizable dialog editor.
        l.addWidget(split,1)
        tbar=QFrame(); tbar.setObjectName('Card'); th=QHBoxLayout(tbar); th.setContentsMargins(12,8,12,8)
        tlab=QLabel('Geometry tables'); tlab.setStyleSheet('font-weight:700;'); th.addWidget(tlab)
        self.table_summary=QLabel(''); self.table_summary.setStyleSheet(f"color:{COLORS['muted']};"); th.addWidget(self.table_summary,1)
        for k,(txt,tip) in enumerate((('Edit nodes…','Node table: node IDs and coordinates x1, x2'),('Edit walls…','Wall table: Node I → Node J, weight, stack ID, wall group'),('Edit periodic pairs…','Explicit periodic node pairs (master, slave, n1, n2)'))):
            eb=QPushButton(txt); eb.setToolTip(tip); eb.clicked.connect(lambda _=False,i=k:self.open_table_editor(i))
            if k==1: eb.setObjectName('Primary')
            th.addWidget(eb)
        wd=QPushButton('Wall editor (drawing)…'); wd.setToolTip("Large wall editor with a drawing of every wall's Node I → Node J (local 1̄) axis."); wd.clicked.connect(self.open_wall_editor); th.addWidget(wd)
        l.addWidget(tbar)

        raw=QTabWidget(); raw.setMinimumHeight(300)
        nbox=QWidget(); nl=QVBoxLayout(nbox); nbar=QHBoxLayout(); an=QPushButton('Add node row'); an.clicked.connect(self.add_node_row); dn=QPushButton('Delete selected node'); dn.clicked.connect(self.delete_node); nbar.addWidget(an); nbar.addWidget(dn); nbar.addStretch(); nl.addLayout(nbar); self.nodes=QTableWidget(0,3); self.nodes.setHorizontalHeaderLabels(['ID','x1 [mm]','x2 [mm]']); self.nodes.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch); nl.addWidget(self.nodes); raw.addTab(nbox,'Node table')
        wbox=QWidget(); wl=QVBoxLayout(wbox); wbar=QHBoxLayout(); aw=QPushButton('Add wall row'); aw.clicked.connect(self.add_wall_row); dw=QPushButton('Delete selected wall'); dw.clicked.connect(self.delete_wall); wbar.addWidget(aw); wbar.addWidget(dw)
        rv=QPushButton('Reverse I ↔ J'); rv.setToolTip('Swap Node I and Node J of the selected walls: reverses their local 1̄ (and 2̄) axes.'); rv.clicked.connect(self.reverse_selected_walls); wbar.addWidget(rv)
        we=QPushButton('Wall editor…'); we.setObjectName('Primary'); we.setToolTip('Large wall editor with a drawing of every wall\'s Node I → Node J (local 1̄) axis.'); we.clicked.connect(self.open_wall_editor); wbar.addWidget(we); wbar.addStretch(); wl.addLayout(wbar)
        wnote=QLabel(wall_axes.html_note()); wnote.setWordWrap(True); wnote.setTextFormat(Qt.TextFormat.RichText); wnote.setStyleSheet('background:#EEF4FB;border:1px solid #CFE0F5;border-radius:5px;padding:5px;'); wl.addWidget(wnote)
        self.walls=QTableWidget(0,8); self.walls.setHorizontalHeaderLabels(['ID','Node I (start)','Node J (end)','Weight','Stack ID','Length [mm]','Angle of 1̄ [deg]','Wall group']); self.walls.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        for c,key in ((1,'Node I'),(2,'Node J'),(3,'Weight'),(4,'Stack ID'),(5,'Length'),(6,'Angle'),(7,'Wall group')): self.walls.horizontalHeaderItem(c).setToolTip(wall_axes.COLUMN_TIPS[key])
        self.walls.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows); self.walls.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        wl.addWidget(self.walls,1); raw.addTab(wbox,'Wall table')
        pbox=QWidget(); pl=QVBoxLayout(pbox); pbar=QHBoxLayout(); ap=QPushButton('Add pair'); ap.clicked.connect(lambda:self.periodic.insertRow(self.periodic.rowCount())); dp=QPushButton('Delete pair'); dp.clicked.connect(self.delete_periodic); pbar.addWidget(ap); pbar.addWidget(dp); pbar.addStretch(); pl.addLayout(pbar); self.periodic=QTableWidget(0,4); self.periodic.setHorizontalHeaderLabels(['Master','Slave','n1','n2']); self.periodic.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch); pl.addWidget(self.periodic); raw.addTab(pbox,'Explicit periodic pairs')
        for t in (self.nodes,self.walls,self.periodic): t.horizontalHeader().setMinimumSectionSize(90)
        self.raw_tabs=raw; self.table_dialog=TableEditorDialog(raw,self,on_refresh=self.refresh_sketch_from_tables)
        for t in (self.nodes,self.walls,self.periodic):
            t.model().rowsInserted.connect(self._update_table_summary); t.model().rowsRemoved.connect(self._update_table_summary)

        self.sketch.geometryEdited.connect(self._sketch_changed); self.sketch.latticeVectorDefined.connect(self._set_lattice_vector); self.sketch.statusMessage.connect(self.sketch_status.setText)
        self.grid_size.valueChanged.connect(self.sketch.set_grid_spacing); self.snap.toggled.connect(self.sketch.set_snap_enabled)
        self.nodes.itemChanged.connect(self._tables_changed); self.walls.itemChanged.connect(self._tables_changed); self.periodic.itemChanged.connect(lambda *_:self.update_geometry_inspector())
        for x in [self.a1x,self.a1y,self.a2x,self.a2y,self.height]: x.valueChanged.connect(lambda *_:self.update_geometry_inspector())
        self.tabs.addTab(w,'RVE Sketch & Geometry')

    def _laminate_tab(self):
        w=QWidget(); l=QVBoxLayout(w)
        help=QLabel('Each wall selects a Stack ID in the Wall table. Use a different Stack ID for a wall when it needs its own material, ply angles or thickness. Stack ID, ply number, material and angle are drop-down selections; angle remains editable for non-standard values.')
        help.setWordWrap(True); l.addWidget(help)
        bar=QHBoxLayout()
        b=QPushButton('Add ply'); b.clicked.connect(self.add_ply_row); bar.addWidget(b)
        b=QPushButton('Delete ply'); b.clicked.connect(self.delete_ply); bar.addWidget(b)
        b=QPushButton('Unique stack for each wall'); b.clicked.connect(self.make_unique_wall_stacks); bar.addWidget(b)
        b=QPushButton('Copy'); b.clicked.connect(self.copy_ply_cells); bar.addWidget(b)
        b=QPushButton('Paste'); b.clicked.connect(self.paste_ply_cells); bar.addWidget(b)
        b=QPushButton('Auto wall groups'); b.clicked.connect(self.auto_wall_groups); bar.addWidget(b)
        b=QPushButton('Apply stack to group'); b.clicked.connect(self.apply_stack_to_group); bar.addWidget(b)
        bar.addStretch(); l.addLayout(bar)
        self.plies=QTableWidget(0,5); self.plies.setHorizontalHeaderLabels(['Stack ID','Ply','Material','Angle [deg]','Thickness [mm]']); self.plies.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch); self.plies.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection); l.addWidget(self.plies,1)
        copy_action=QAction('Copy',self.plies); copy_action.setShortcut(QKeySequence(QKeySequence.StandardKey.Copy)); copy_action.setShortcutContext(Qt.ShortcutContext.WidgetWithChildrenShortcut); copy_action.triggered.connect(self.copy_ply_cells); self.plies.addAction(copy_action)
        paste_action=QAction('Paste',self.plies); paste_action.setShortcut(QKeySequence(QKeySequence.StandardKey.Paste)); paste_action.setShortcutContext(Qt.ShortcutContext.WidgetWithChildrenShortcut); paste_action.triggered.connect(self.paste_ply_cells); self.plies.addAction(paste_action)
        self.tabs.addTab(w,'Wall laminates')

    def _solver_tab(self):
        from ..laminate import TWIST_RELEASE_DEFAULT
        w=QWidget(); l=QVBoxLayout(w); f=QGroupBox('Homogenization settings'); form=QFormLayout(f)
        self.solver=QComboBox(); self.solver.addItems(['Cylindrical','Twist-corrected','Anisotropic strip','Free moment'])
        self.solver.setToolTip('Cylindrical: original model (twist suppressed). Twist-corrected: D_b = D11(1 − ηψ²) with the calibrated twist release '
                               'η = (a + bψ²)λⁿ, λ = √(D66/A44)/L. Anisotropic strip / free moment: bounds.')
        self.tol=QDoubleSpinBox(); self.tol.setDecimals(12); self.tol.setRange(1e-14,1e-2); self.tol.setValue(1e-8)
        self.anchor=QSpinBox(); self.anchor.setRange(1,100000); self.anchor.setValue(1)
        self.wall_beam=QComboBox(); self.wall_beam.addItem('Euler–Bernoulli (default)','euler-bernoulli'); self.wall_beam.addItem('Timoshenko (with wall shear flexibility)','timoshenko')
        wf=QGroupBox('Wall model'); wform=QFormLayout(wf); wform.addRow('Wall beam',self.wall_beam); wform.addRow('Bending option',self.solver); l.addWidget(wf)
        # Editable twist-corrected model: eta = (a + b psi^2) lambda^n
        tg=QGroupBox('Twist-corrected wall: η = (a + b ψ²) λⁿ,  λ = √(D₆₆/A₄₄)/L'); tgl=QGridLayout(tg)
        self.twist_coeff=[]
        for k,(lab,val,rng,tip) in enumerate((('a',TWIST_RELEASE_DEFAULT[0],(0,50),'constant term'),('b',TWIST_RELEASE_DEFAULT[1],(-50,50),'ψ² term'),('n',TWIST_RELEASE_DEFAULT[2],(0,3),'exponent of the twist length λ'))):
            sp=QDoubleSpinBox(); sp.setDecimals(4); sp.setRange(*rng); sp.setSingleStep(0.01); sp.setValue(val); sp.setToolTip(tip)
            tgl.addWidget(QLabel(lab),0,2*k); tgl.addWidget(sp,0,2*k+1); self.twist_coeff.append(sp)
        rb=QPushButton('Restore calibrated values'); rb.clicked.connect(lambda:[sp.setValue(v) for sp,v in zip(self.twist_coeff,TWIST_RELEASE_DEFAULT)]); tgl.addWidget(rb,0,6)
        tn=QLabel('Calibrated on 320 walls with mesh-converged periodic shell FE (a = 2.0, b = 2.15, n = 0.73; within 2 % for |ψ| ≤ 0.85, β ≤ 0.10, '
                  '0.010 ≤ λ ≤ 0.100). Other values are used as given and the reliability rating reports them as user-defined.')
        tn.setWordWrap(True); tn.setStyleSheet(f"color:{COLORS['muted']};"); tgl.addWidget(tn,1,0,1,7); l.addWidget(tg); self.twist_group=tg
        self.solver.currentTextChanged.connect(lambda t:self.twist_group.setEnabled(t.lower().startswith('twist'))); self.twist_group.setEnabled(False)
        form.addRow('Periodic tolerance',self.tol); form.addRow('Anchor node',self.anchor); l.addWidget(f)
        note=QLabel('Cylindrical uses the paper-based wall formulation (D′11); Twist-corrected releases the calibrated part ηψ² of the bend–twist stiffness; '
                     'Anisotropic strip condenses extension–bending coupling ((D − BA⁻¹B)11, differs only for unsymmetric stacks); Free moment also releases twist (soft bound). '
                     'Euler–Bernoulli walls reproduce the original model (Sriharan et al. 2024); Timoshenko walls add shear '
                     'flexibility through β = 12D₁₁/(L²S_eq).'); note.setWordWrap(True); note.setStyleSheet(f"color:{COLORS['muted']};"); l.addWidget(note)
        b=QPushButton('Run'); b.setObjectName('Primary'); b.setMinimumHeight(44); b.clicked.connect(self.run_model); l.addWidget(b); self.status=QLabel('Ready'); self.status.setWordWrap(True); l.addWidget(self.status); l.addStretch(); self.tabs.addTab(w,'Solver')

    def twist_coefficients(self):
        """(a, b, n) from the solver tab, or None when they equal the calibrated defaults."""
        from ..laminate import TWIST_RELEASE_DEFAULT
        c=tuple(float(sp.value()) for sp in self.twist_coeff)
        return None if all(abs(x-y)<=1e-9 for x,y in zip(c,TWIST_RELEASE_DEFAULT)) else c

    def _results_tab(self):
        w=QWidget(); l=QVBoxLayout(w)
        split=QSplitter(Qt.Orientation.Horizontal)
        tables=QTabWidget()
        self.ch_table=QTableWidget(); self.ch_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers); tables.addTab(self.ch_table,'Cᴴ [MPa]')
        self.s_table=QTableWidget(); self.s_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers); tables.addTab(self.s_table,'S [1/MPa]')
        self.qplane_table=QTableWidget(); self.qplane_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers); tables.addTab(self.qplane_table,'In-plane Q [MPa]')
        self.qshear_table=QTableWidget(); self.qshear_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers); tables.addTab(self.qshear_table,'Shear Q [MPa]')
        self.eng_table=QTableWidget(0,3); self.eng_table.setHorizontalHeaderLabels(['Property','Value','Unit']); self.eng_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch); self.eng_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers); tables.addTab(self.eng_table,'Engineering properties')
        self.phys_table=QTableWidget(0,3); self.phys_table.setHorizontalHeaderLabels(['Quantity','Value','Unit']); self.phys_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch); self.phys_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers); tables.addTab(self.phys_table,'RVE / physical')
        split.addWidget(tables)
        preview=QWidget(); pv=QVBoxLayout(preview); pv.addWidget(QLabel('RVE geometry')); self.fig=Figure(figsize=(5,4)); self.canvas=FigureCanvas(self.fig); pv.addWidget(self.canvas,1); split.addWidget(preview); split.setSizes([850,500])
        l.addWidget(split,1); self.tabs.addTab(w,'Results')

    def _wall_load_tab(self):
        w=QWidget(); l=QVBoxLayout(w)
        note=QLabel('Recover the paper-based wall-local membrane resultants under an arbitrary combined macroscopic load. Compression is negative. Wall resultants are reported as N̄11, N̄33 and N̄13 in N/mm in each wall local (1̄,3̄) plate plane. The model is re-solved from the current geometry and laminate definition before recovery.')
        note.setWordWrap(True); l.addWidget(note)
        loadbox=QGroupBox('Macroscopic design load'); grid=QGridLayout(loadbox)
        self.wall_load_control=QComboBox(); self.wall_load_control.addItems(['Stress controlled [MPa]','Strain controlled [-]']); self.wall_load_control.currentTextChanged.connect(self._wall_load_control_changed)
        self.wall_load_preset=QComboBox(); self.wall_load_preset.addItems(['Custom','Example: σ11 + σ33 + τ13','Uniaxial σ11 compression','Uniaxial σ22 compression','Pure τ12 shear']); self.wall_load_preset.currentTextChanged.connect(self._apply_wall_load_preset)
        grid.addWidget(QLabel('Control'),0,0); grid.addWidget(self.wall_load_control,0,1); grid.addWidget(QLabel('Preset'),0,2); grid.addWidget(self.wall_load_preset,0,3)
        self.wall_load_spins=[]
        for i,lab in enumerate(MACRO_LABELS):
            box=QDoubleSpinBox(); box.setDecimals(8); box.setRange(-1e9,1e9); box.setSingleStep(0.1); box.setValue(0.0); self.wall_load_spins.append(box)
            grid.addWidget(QLabel(('σ' if i<3 else 'τ')+lab if i<3 else 'τ'+lab),1+(i//3)*2,i%3)
            grid.addWidget(box,2+(i//3)*2,i%3)
        self.wall_load_spins[0].setValue(-1.0); self.wall_load_spins[2].setValue(-1.0); self.wall_load_spins[4].setValue(0.25)
        buttons=QHBoxLayout(); run=QPushButton('Recover wall loads'); run.setObjectName('Primary'); run.clicked.connect(self.recover_current_wall_loads); exp=QPushButton('Export wall loads CSV'); exp.clicked.connect(self.export_wall_loads); buttons.addWidget(run); buttons.addWidget(exp); buttons.addStretch(); grid.addLayout(buttons,5,0,1,4)
        l.addWidget(loadbox)
        self.wall_load_macro=QTableWidget(6,3); self.wall_load_macro.setHorizontalHeaderLabels(['Component','Macro strain','Macro stress [MPa]']); self.wall_load_macro.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch); self.wall_load_macro.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        for r,lab in enumerate(MACRO_LABELS): self.wall_load_macro.setItem(r,0,table_item(lab))
        self.wall_load_validation=QTableWidget(0,3); self.wall_load_validation.setHorizontalHeaderLabels(['Validation check','Value','Status']); self.wall_load_validation.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch); self.wall_load_validation.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        top=QSplitter(Qt.Orientation.Horizontal); top.addWidget(self.wall_load_macro); top.addWidget(self.wall_load_validation); top.setSizes([500,700]); l.addWidget(top)
        self.wall_load_table=QTableWidget(0,12); self.wall_load_table.setHorizontalHeaderLabels(['Wall','Group','Stack','L [mm]','Angle [deg]','ε̄11','ε̄33','γ̄13','N̄11 [N/mm]','N̄33 [N/mm]','N̄13 [N/mm]','Umem [N mm]']); self.wall_load_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents); self.wall_load_table.horizontalHeader().setStretchLastSection(True); self.wall_load_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers); l.addWidget(self.wall_load_table,1)
        self.wall_load_status=QLabel('Run a load recovery to populate this table.'); self.wall_load_status.setWordWrap(True); l.addWidget(self.wall_load_status)
        self.wall_load_df=None; self.wall_load_recovery=None; self._wall_load_control_changed(); self.tabs.addTab(w,'Wall load recovery')

    def _buckling_tab(self):
        w=QWidget(); l=QVBoxLayout(w)
        note=QLabel('Linear local-wall buckling under the recovered combined biaxial membrane resultants N̄11, N̄33 and N̄13. The model uses a multi-term simply-supported CLPT Ritz/Navier eigenproblem. Compression is negative. For unsymmetric laminates the default D*=D−B A⁻¹B condensation is a reduced-order screening assumption; the reported factor should be validated against a higher-fidelity plate/FE model when membrane-bending coupling or edge restraint is important.')
        note.setWordWrap(True); l.addWidget(note)
        box=QGroupBox('Buckling settings'); g=QGridLayout(box)
        self.buck_m=QSpinBox(); self.buck_m.setRange(2,20); self.buck_m.setValue(8)
        self.buck_n=QSpinBox(); self.buck_n.setRange(2,20); self.buck_n.setValue(8)
        self.buck_sf=QDoubleSpinBox(); self.buck_sf.setRange(0.01,100); self.buck_sf.setDecimals(3); self.buck_sf.setValue(1.5)
        self.buck_abd=QCheckBox('Use ABD-condensed bending D*=D−B A⁻¹B'); self.buck_abd.setChecked(True)
        self.buck_conv=QCheckBox('Check modal convergence against (m−2)×(n−2)'); self.buck_conv.setChecked(True)
        g.addWidget(QLabel('m max'),0,0); g.addWidget(self.buck_m,0,1); g.addWidget(QLabel('n max'),0,2); g.addWidget(self.buck_n,0,3); g.addWidget(QLabel('Required factor'),0,4); g.addWidget(self.buck_sf,0,5)
        g.addWidget(self.buck_abd,1,0,1,3); g.addWidget(self.buck_conv,1,3,1,3)
        run=QPushButton('Evaluate buckling from recovered wall loads'); run.setObjectName('Primary'); run.clicked.connect(self.evaluate_current_buckling)
        rerun=QPushButton('Recover current loads + evaluate'); rerun.clicked.connect(self.recover_and_evaluate_buckling)
        exp=QPushButton('Export buckling CSV'); exp.clicked.connect(self.export_buckling)
        g.addWidget(run,2,0,1,2); g.addWidget(rerun,2,2,1,2); g.addWidget(exp,2,4,1,2); l.addWidget(box)
        self.buck_table=QTableWidget(0,16); self.buck_table.setHorizontalHeaderLabels(['Wall','Group','Stack','a=L [mm]','b=h [mm]','N̄11','N̄33','N̄13','λcr','Required','Margin','Pass','m','n','Convergence','Symmetry indicator']); self.buck_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents); self.buck_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers); l.addWidget(self.buck_table,1)
        self.buck_status=QLabel('Recover wall loads first, then evaluate buckling.'); self.buck_status.setWordWrap(True); l.addWidget(self.buck_status)
        self.buckling_df=None; self.tabs.addTab(w,'Wall buckling')

    def recover_and_evaluate_buckling(self):
        self.recover_current_wall_loads()
        if self.wall_load_recovery is not None:
            self.evaluate_current_buckling()

    def evaluate_current_buckling(self):
        if self.wall_load_recovery is None or self.result is None or self.model is None:
            QMessageBox.information(self,'Wall buckling','Recover wall loads first.'); return
        try:
            buck=evaluate_wall_buckling(self.model,self.result,self.wall_load_recovery,m_max=self.buck_m.value(),n_max=self.buck_n.value(),safety_factor=self.buck_sf.value(),use_abd_condensation=self.buck_abd.isChecked(),convergence_check=self.buck_conv.isChecked())
            self.buckling_summary=buck; self.buckling_df=pd.DataFrame(buck.walls); self.buck_table.setRowCount(len(buck.walls))
            cols=['wall','group','laminate','a_wall_length_mm','b_core_height_mm','N11_N_per_mm','N33_N_per_mm','N13_N_per_mm','lambda_cr','required_safety_factor','buckling_margin','pass','dominant_m','dominant_n','mode_convergence_relative_error','laminate_symmetry_indicator']
            for rr,row in enumerate(buck.walls):
                for cc,c in enumerate(cols):
                    v=row.get(c,'');
                    if isinstance(v,(float,np.floating)):
                        txt='∞' if math.isinf(float(v)) else (f'{v:.6e}' if c in ('mode_convergence_relative_error','laminate_symmetry_indicator') else f'{v:.8g}')
                    else: txt=str(v)
                    self.buck_table.setItem(rr,cc,table_item(txt))
            finite_conv=[w['mode_convergence_relative_error'] for w in buck.walls if math.isfinite(float(w['mode_convergence_relative_error']))]
            conv=max(finite_conv) if finite_conv else float('nan')
            state='PASS' if buck.passed else 'FAIL'
            self.buck_status.setText(f'{state}: minimum λcr = {buck.minimum_lambda:.6g}, critical wall = {buck.critical_wall}, required factor = {buck.safety_factor:.4g}. Maximum modal truncation change = {conv:.3%}.' if math.isfinite(conv) else f'{state}: minimum λcr = {buck.minimum_lambda:.6g}, critical wall = {buck.critical_wall}, required factor = {buck.safety_factor:.4g}.')
        except Exception as exc:
            self.buck_status.setText('Buckling evaluation failed.'); QMessageBox.critical(self,'Wall buckling',str(exc))

    def export_buckling(self):
        if self.buckling_df is None or self.buckling_df.empty:
            QMessageBox.information(self,'Wall buckling','Evaluate buckling first.'); return
        path,_=QFileDialog.getSaveFileName(self,'Export buckling results','wall_buckling.csv','CSV (*.csv);;Excel (*.xlsx)')
        if not path:return
        try:
            if path.lower().endswith('.xlsx'): self.buckling_df.to_excel(path,index=False)
            else:
                if not path.lower().endswith('.csv'): path += '.csv'
                self.buckling_df.to_csv(path,index=False)
        except Exception as exc: QMessageBox.critical(self,'Export buckling',str(exc))

    def _wall_load_control_changed(self,*_):
        stress=self.wall_load_control.currentText().startswith('Stress')
        for box in self.wall_load_spins:
            box.blockSignals(True)
            if stress:
                box.setRange(-1e9,1e9); box.setDecimals(8); box.setSingleStep(0.1); box.setSuffix(' MPa')
            else:
                box.setRange(-10.0,10.0); box.setDecimals(10); box.setSingleStep(0.0001); box.setSuffix('')
            box.blockSignals(False)

    def _apply_wall_load_preset(self,name):
        if name=='Custom': return
        self.wall_load_control.setCurrentText('Stress controlled [MPa]')
        vals={'Example: σ11 + σ33 + τ13':[-1.0,0,-1.0,0,0.25,0], 'Uniaxial σ11 compression':[-1.0,0,0,0,0,0], 'Uniaxial σ22 compression':[0,-1.0,0,0,0,0], 'Pure τ12 shear':[0,0,0,0,0,1.0]}.get(name,[0]*6)
        for box,val in zip(self.wall_load_spins,vals): box.setValue(val)

    def recover_current_wall_loads(self):
        try:
            s=self.current_geometry_summary(); errors=[i.message for i in s.issues if i.severity=='error']
            if errors: raise ValueError('Geometry validation failed:\n• '+'\n• '.join(errors))
            self.wall_load_status.setText('Solving current RVE and recovering wall loads…'); QApplication.processEvents()
            model,mats,lams=self.parse_model(); out=homogenize_rve(model,mats,lams); vals=[x.value() for x in self.wall_load_spins]
            if self.wall_load_control.currentText().startswith('Stress'):
                rec=recover_wall_loads(model,out,macro_stress=vals)
            else:
                rec=recover_wall_loads(model,out,macro_strain=vals)
            self.model=model; self.laminates=lams; self.result=out; self.wall_load_recovery=rec
            for r,lab in enumerate(MACRO_LABELS):
                self.wall_load_macro.setItem(r,0,table_item(lab)); self.wall_load_macro.setItem(r,1,table_item(f'{rec.macro_strain[r]:.10g}')); self.wall_load_macro.setItem(r,2,table_item(f'{rec.macro_stress[r]:.10g}'))
            cols=['wall','group','laminate','length_mm','angle_deg','eps11_local','eps33_local','gamma13_local','N11_N_per_mm','N33_N_per_mm','N13_N_per_mm','membrane_energy_Nmm']
            self.wall_load_df=pd.DataFrame(rec.walls)
            self.wall_load_table.setRowCount(len(rec.walls))
            for rr,row in enumerate(rec.walls):
                values=[row.get(c,'') for c in cols]
                for cc,v in enumerate(values):
                    txt=f'{v:.10g}' if isinstance(v,(float,np.floating)) and math.isfinite(float(v)) else str(v)
                    self.wall_load_table.setItem(rr,cc,table_item(txt))
            checks=[
                ('Global equilibrium residual',rec.validation['equilibrium_relative_residual'],1e-8),
                ('Macro/micro energy error',rec.validation['energy_relative_error'],1e-9),
                ('Macro work consistency error',rec.validation['macro_work_relative_error'],1e-9),
                ('Element energy sum error',rec.validation['element_energy_sum_relative_error'],1e-9),
                ('Membrane nodal-force identity',rec.validation['max_membrane_force_identity_error'],1e-9),
                ('N̄33 conjugate-force identity',rec.validation['max_N33_conjugate_identity_error'],1e-9),
                ('Periodic constraint residual',rec.validation['constraint_relative_residual'],1e-9),
            ]
            if 'stress_control_relative_error' in rec.validation: checks.append(('Stress-control reconstruction error',rec.validation['stress_control_relative_error'],1e-9))
            self.wall_load_validation.setRowCount(len(checks))
            passed=True
            for rr,(name,val,tol) in enumerate(checks):
                ok=math.isfinite(float(val)) and float(val)<=tol; passed=passed and ok
                for cc,text in enumerate([name,f'{val:.3e}','PASS' if ok else 'CHECK']): self.wall_load_validation.setItem(rr,cc,table_item(text))
            minN=min((min(w['N11_N_per_mm'],w['N33_N_per_mm']) for w in rec.walls),default=float('nan'))
            maxShear=max((abs(w['N13_N_per_mm']) for w in rec.walls),default=float('nan'))
            self.wall_load_status.setText(f"Recovery {'validated' if passed else 'completed with checks requiring attention'}. Most compressive local normal resultant = {minN:.6g} N/mm; maximum |N̄13| = {maxShear:.6g} N/mm. These are the quantities intended for the subsequent biaxial-compression + shear buckling model.")
        except Exception as exc:
            self.wall_load_status.setText('Wall load recovery failed.'); QMessageBox.critical(self,'Wall load recovery',str(exc))

    def export_wall_loads(self):
        if self.wall_load_df is None or self.wall_load_df.empty:
            QMessageBox.information(self,'Wall load recovery','Recover wall loads first.'); return
        path,_=QFileDialog.getSaveFileName(self,'Export wall loads','wall_loads.csv','CSV (*.csv);;Excel (*.xlsx)')
        if not path:return
        try:
            if path.lower().endswith('.xlsx'): self.wall_load_df.to_excel(path,index=False)
            else:
                if not path.lower().endswith('.csv'): path += '.csv'
                self.wall_load_df.to_csv(path,index=False)
        except Exception as exc: QMessageBox.critical(self,'Export wall loads',str(exc))

    def refresh_geometry_choices(self,select_id=None):
        """Reload the Geometry Library into the Single-RVE geometry selector."""
        if not hasattr(self,'geometry_selector'): return
        current=select_id or self.geometry_selector.currentData()
        try: rows=load_geometry_library()
        except Exception as exc:
            self.sketch_status.setText(f'Cannot read Geometry Library: {exc}') if hasattr(self,'sketch_status') else None
            return
        self._geometry_choice_rows={str(r.get('id','')).upper():dict(r) for r in rows if str(r.get('id','')).strip()}
        old=self.geometry_selector.blockSignals(True); self.geometry_selector.clear()
        try:
            for r in rows:
                gid=str(r.get('id','')).upper().strip()
                if not gid: continue
                suffix='' if r.get('active',True) else ' (inactive)'
                self.geometry_selector.addItem(f"{gid} — {r.get('name',gid)}{suffix}",gid)
            idx=self.geometry_selector.findData(str(current).upper()) if current else -1
            if idx<0 and self.geometry_selector.count(): idx=0
            if idx>=0:self.geometry_selector.setCurrentIndex(idx)
        finally:
            self.geometry_selector.blockSignals(old)
        self._update_geometry_parameter_controls()

    def _selected_geometry_definition(self):
        gid=str(self.geometry_selector.currentData() or '').upper() if hasattr(self,'geometry_selector') else ''
        return getattr(self,'_geometry_choice_rows',{}).get(gid)

    def _update_geometry_parameter_controls(self,*_):
        if not hasattr(self,'geometry_parameter_spins'): return
        gd=self._selected_geometry_definition() or {}
        params=list(gd.get('parameters',[]))[:4]
        for i,(lab,spin) in enumerate(zip(self.geometry_parameter_labels,self.geometry_parameter_spins)):
            if i>=len(params):
                lab.hide(); spin.hide(); continue
            p=params[i]; name=str(p.get('name',f'P{i+1}'))
            mn=float(p.get('min',0.0)); mx=float(p.get('max',mn+1.0))
            ref=p.get('reference',None)
            try:ref=float(ref) if ref is not None else 0.5*(mn+mx)
            except Exception:ref=0.5*(mn+mx)
            if not math.isfinite(mn):mn=-1e12
            if not math.isfinite(mx):mx=1e12
            if mx<mn:mn,mx=mx,mn
            ref=min(max(ref,mn),mx)
            lab.setText(name); lab.setToolTip(f'{name} range: {mn:.6g} to {mx:.6g}')
            spin.setDecimals(0 if name.upper()=='NODECOUNT' else 6); spin.setRange(mn,mx); spin.setValue(ref); spin.setSuffix('' if name.upper() in ('R','NODECOUNT','STRETCHX','STRETCHY','SHEAR') else ' mm')
            lab.show(); spin.show()

    def load_selected_geometry(self):
        """Generate and load the selected Geometry Library family into the sketch."""
        gd=self._selected_geometry_definition()
        if not gd:
            QMessageBox.information(self,'Existing geometry','Select a geometry from the drop-down first.'); return
        gid=str(gd.get('id','')).upper()
        try:
            registry=discover_geometry_plugins(PACKAGE_ROOT/'plugins')
            if gid not in registry: raise ValueError(f'Geometry generator {gid} was not found. If it is a custom geometry, check the plugins folder.')
            pvals={}
            for i,pdef in enumerate(gd.get('parameters',[])[:4]):
                pvals[str(pdef['name'])]=float(self.geometry_parameter_spins[i].value())
            # More than four parameters are uncommon, but remain supported using
            # their stored reference/mid-range values.
            for pdef in gd.get('parameters',[])[4:]:
                mn=float(pdef.get('min',0)); mx=float(pdef.get('max',mn)); pvals[str(pdef['name'])]=float(pdef.get('reference',0.5*(mn+mx)))
            geom=registry[gid](pvals)
            nodes=np.asarray(geom.nodes,float); conn=np.asarray(geom.connectivity,int)
            if nodes.ndim!=2 or nodes.shape[1]!=2 or len(nodes)==0: raise ValueError('Geometry generator returned invalid nodes.')
            if conn.ndim!=2 or conn.shape[1]!=2 or len(conn)==0: raise ValueError('Geometry generator returned invalid wall connectivity.')
            if len(getattr(geom,'wall_weights',[]))!=len(conn): raise ValueError('Geometry generator returned a wall-weight count that does not match the walls.')
            self._syncing=True
            try:
                self.nodes.setRowCount(len(nodes))
                for i,(x,y) in enumerate(nodes):
                    for j,val in enumerate((i+1,float(x),float(y))):self.nodes.setItem(i,j,table_item(val))
                self.walls.setRowCount(len(conn))
                groups=list(getattr(geom,'wall_groups',[])); weights=np.asarray(geom.wall_weights,float)
                for i,(ii,jj) in enumerate(conn):
                    group=groups[i] if i<len(groups) else f'wall-{i+1}'
                    self._set_wall_row(i,i+1,int(ii)+1,int(jj)+1,float(weights[i]),1,group=group)
                self.periodic.setRowCount(0)
                self.a1x.setValue(float(geom.a1[0])); self.a1y.setValue(float(geom.a1[1])); self.a2x.setValue(float(geom.a2[0])); self.a2y.setValue(float(geom.a2[1]))
            finally:self._syncing=False
            self.refresh_sketch_from_tables(); self.sketch.fit_geometry(); self._reset_history()
            ptxt=', '.join(f'{k}={v:.6g}' for k,v in pvals.items()) or 'fixed geometry'
            self.sketch_status.setText(f'Loaded {gid} from Geometry Library ({ptxt}).')
        except Exception as exc:
            QMessageBox.critical(self,'Cannot load geometry',str(exc))

    def set_sketch_mode(self,mode):
        for m,b in self.mode_buttons.items(): b.setChecked(m==mode)
        self.sketch.set_mode(mode)

    def add_node_row(self):
        r=self.nodes.rowCount(); self.nodes.insertRow(r); self.nodes.setItem(r,0,table_item(r+1)); self.nodes.setItem(r,1,table_item(0.0)); self.nodes.setItem(r,2,table_item(0.0)); self.refresh_sketch_from_tables()

    def _stack_choices(self):
        return [str(i) for i in range(1,51)]

    def _material_choices(self):
        names=list(self.materials.keys())
        for row in load_material_library():
            name=str(row.get('name','')).strip()
            if name and name not in names: names.append(name)
        return names

    def _set_wall_row(self,r,wall_id,ni,nj,weight,stack_id,length='',angle='',group=''):
        for c,v in enumerate([wall_id,ni,nj,weight]): self.walls.setItem(r,c,table_item(v))
        set_combo_cell(self.walls,r,4,self._stack_choices(),stack_id,callback=lambda *_:self._tables_changed())
        self.walls.setItem(r,5,table_item(length)); self.walls.setItem(r,6,table_item(angle)); self.walls.setItem(r,7,table_item(group))

    def add_wall_row(self):
        r=self.walls.rowCount(); self.walls.insertRow(r); self._set_wall_row(r,r+1,1,min(2,max(1,self.nodes.rowCount())),1.0,1); self.refresh_sketch_from_tables()

    def _set_ply_row(self,r,stack_id=1,ply_no=1,material=None,angle=0,thickness=0.5):
        material=material or (self._material_choices()[0] if self._material_choices() else '')
        set_combo_cell(self.plies,r,0,self._stack_choices(),stack_id,callback=lambda *_:self._laminate_table_changed())
        set_combo_cell(self.plies,r,1,[str(i) for i in range(1,51)],ply_no,callback=lambda *_:self._laminate_table_changed())
        set_combo_cell(self.plies,r,2,self._material_choices(),material,editable=True,callback=lambda *_:self._laminate_table_changed())
        angles=['-90','-75','-60','-45','-30','-15','0','15','30','45','60','75','90']
        set_combo_cell(self.plies,r,3,angles,angle,editable=True,callback=lambda *_:self._laminate_table_changed())
        self.plies.setItem(r,4,table_item(thickness))

    def add_ply_row(self):
        sid=1
        if self.plies.currentRow()>=0:
            try:sid=int(float(self._cell(self.plies,self.plies.currentRow(),0)))
            except Exception:pass
        ply_no=1+sum(1 for r in range(self.plies.rowCount()) if self._cell(self.plies,r,0)==str(sid))
        r=self.plies.rowCount(); self.plies.insertRow(r); self._set_ply_row(r,sid,ply_no)

    def _laminate_table_changed(self):
        if self._syncing:return
        self._refresh_wall_stack_combos()

    def _refresh_wall_stack_combos(self):
        for r in range(self.walls.rowCount()):
            current=self._cell(self.walls,r,4) or '1'
            old=self._syncing; self._syncing=True
            try:set_combo_cell(self.walls,r,4,self._stack_choices(),current,callback=lambda *_:self._tables_changed())
            finally:self._syncing=old

    def make_unique_wall_stacks(self):
        if self.walls.rowCount()==0:
            QMessageBox.information(self,'Wall laminates','Draw at least one wall first.'); return
        source_sid=1
        if self.plies.currentRow()>=0:
            try:source_sid=int(float(self._cell(self.plies,self.plies.currentRow(),0)))
            except Exception:pass
        source=[]
        for r in range(self.plies.rowCount()):
            try:
                if int(float(self._cell(self.plies,r,0)))==source_sid:
                    source.append((int(float(self._cell(self.plies,r,1))),self._cell(self.plies,r,2),float(self._cell(self.plies,r,3)),float(self._cell(self.plies,r,4))))
            except Exception:pass
        if not source:
            QMessageBox.warning(self,'Wall laminates',f'Stack {source_sid} has no plies to copy.'); return
        self._syncing=True
        try:
            self.plies.setRowCount(0)
            rr=0
            for wall in range(1,self.walls.rowCount()+1):
                set_combo_cell(self.walls,wall-1,4,self._stack_choices(),wall,callback=lambda *_:self._tables_changed())
                for ply_no,mat,ang,th in source:
                    self.plies.insertRow(rr); self._set_ply_row(rr,wall,ply_no,mat,ang,th); rr+=1
        finally:self._syncing=False
        self.refresh_sketch_from_tables()
        QMessageBox.information(self,'Wall laminates','Each wall now has its own stack. Material, angle and ply thickness can be edited independently for every wall.')

    def copy_ply_cells(self):
        indexes=self.plies.selectedIndexes()
        if not indexes and self.plies.currentRow()>=0 and self.plies.currentColumn()>=0:
            rows=[self.plies.currentRow()]; cols=[self.plies.currentColumn()]
        elif indexes:
            rows=range(min(i.row() for i in indexes),max(i.row() for i in indexes)+1); cols=range(min(i.column() for i in indexes),max(i.column() for i in indexes)+1)
        else:return
        QApplication.clipboard().setText('\n'.join('\t'.join(self._cell(self.plies,r,c) for c in cols) for r in rows))

    def _write_ply_cell(self,r,c,value):
        while r>=self.plies.rowCount():
            rr=self.plies.rowCount(); self.plies.insertRow(rr); self._set_ply_row(rr,1,rr+1)
        if c==0:set_combo_cell(self.plies,r,c,self._stack_choices(),value,callback=lambda *_:self._laminate_table_changed())
        elif c==1:set_combo_cell(self.plies,r,c,[str(i) for i in range(1,51)],value,callback=lambda *_:self._laminate_table_changed())
        elif c==2:set_combo_cell(self.plies,r,c,self._material_choices(),value,editable=True,callback=lambda *_:self._laminate_table_changed())
        elif c==3:set_combo_cell(self.plies,r,c,['-90','-75','-60','-45','-30','-15','0','15','30','45','60','75','90'],value,editable=True,callback=lambda *_:self._laminate_table_changed())
        elif c==4:self.plies.setItem(r,c,table_item(value))

    def paste_ply_cells(self):
        text=QApplication.clipboard().text()
        if not text.strip():return
        r0=max(0,self.plies.currentRow()); c0=max(0,self.plies.currentColumn())
        self._syncing=True
        try:
            for dr,line in enumerate(text.splitlines()):
                for dc,value in enumerate(line.split('\t')):
                    if c0+dc<self.plies.columnCount():self._write_ply_cell(r0+dr,c0+dc,value.strip())
        finally:self._syncing=False
        self._refresh_wall_stack_combos()

    def delete_node(self):
        r=self.nodes.currentRow()
        if r<0:return
        target=r+1; self.nodes.removeRow(r)
        wall_rows=[]
        for rr in range(self.walls.rowCount()):
            try:
                ni=int(float(self._cell(self.walls,rr,1))); nj=int(float(self._cell(self.walls,rr,2))); wt=float(self._cell(self.walls,rr,3) or 1); sid=int(float(self._cell(self.walls,rr,4) or 1)); grp=self._cell(self.walls,rr,7)
            except Exception: continue
            if target in (ni,nj): continue
            wall_rows.append([ni-1 if ni>target else ni,nj-1 if nj>target else nj,wt,sid,'','',grp])
        self._syncing=True
        try:
            for rr in range(self.nodes.rowCount()): self.nodes.setItem(rr,0,table_item(rr+1))
            self.walls.setRowCount(len(wall_rows))
            for rr,row in enumerate(wall_rows):self._set_wall_row(rr,rr+1,*row)
        finally:self._syncing=False
        self.refresh_sketch_from_tables()

    def _wall_table_rows(self):
        """Wall table as 1-based rows [node_i, node_j, weight, stack_id, group]."""
        rows=[]
        for r in range(self.walls.rowCount()):
            def num(c,cast,default):
                try:return cast(float(self._cell(self.walls,r,c)))
                except Exception:return default
            rows.append([num(1,int,0),num(2,int,0),num(3,float,1.0),num(4,int,1),self._cell(self.walls,r,7)])
        return rows

    def _write_wall_rows(self,rows):
        self._syncing=True
        try:
            self.walls.setRowCount(len(rows))
            for k,(ni,nj,wt,sid,grp) in enumerate(rows): self._set_wall_row(k,k+1,ni,nj,wt,sid,group=grp)
        finally:self._syncing=False
        self.refresh_sketch_from_tables(); self._commit_history()

    def reverse_selected_walls(self):
        rows=sorted({i.row() for i in self.walls.selectedIndexes()})
        if not rows:
            QMessageBox.information(self,'Reverse I ↔ J','Select one or more walls in the Wall table first (or click a wall in the sketch).'); return
        self._write_wall_rows(wall_axes.reverse_walls(self._wall_table_rows(),rows))
        from PySide6.QtCore import QItemSelectionModel
        sm=self.walls.selectionModel(); sm.clearSelection()
        for r in rows: sm.select(self.walls.model().index(r,0),QItemSelectionModel.SelectionFlag.Select|QItemSelectionModel.SelectionFlag.Rows)
        self.sketch_status.setText('Reversed '+', '.join(f'W{r+1}' for r in rows)+': their local 1̄ axes now point the other way (Node I ↔ Node J).')

    def open_table_editor(self,index:int=0):
        """Open the resizable Node / Wall / periodic-pair table editor at the given tab."""
        self.table_dialog.open_tab(index)

    def _update_table_summary(self,*_):
        if not hasattr(self,'table_summary'): return
        n=self.nodes.rowCount(); w=self.walls.rowCount(); p=self.periodic.rowCount()
        self.table_summary.setText(f"{n} node{'s' if n!=1 else ''} · {w} wall{'s' if w!=1 else ''} · "
                                   + (f"{p} explicit periodic pair{'s' if p!=1 else ''}" if p else 'periodic pairs detected automatically'))

    def open_wall_editor(self):
        from .wall_editor import WallEditorDialog
        try:
            nodes=[[float(self._cell(self.nodes,r,1)),float(self._cell(self.nodes,r,2))] for r in range(self.nodes.rowCount())]
        except Exception:
            QMessageBox.warning(self,'Wall editor','Complete the Node table first (every node needs x1 and x2).'); return
        dlg=WallEditorDialog(nodes,self._wall_table_rows(),[self.a1x.value(),self.a1y.value()],[self.a2x.value(),self.a2y.value()],self._stack_choices(),self)
        dlg.select_rows(sorted({i.row() for i in self.walls.selectedIndexes()}))
        if dlg.exec(): self._write_wall_rows(dlg.rows())

    def delete_wall(self):
        r=self.walls.currentRow()
        if r>=0:self.walls.removeRow(r); self._renumber_walls(); self.refresh_sketch_from_tables()

    def delete_periodic(self):
        r=self.periodic.currentRow()
        if r>=0:self.periodic.removeRow(r); self.update_geometry_inspector()

    def delete_ply(self):
        r=self.plies.currentRow()
        if r>=0:self.plies.removeRow(r); self._refresh_wall_stack_combos()

    def _renumber_walls(self):
        self._syncing=True
        try:
            for r in range(self.walls.rowCount()): self.walls.setItem(r,0,table_item(r+1))
        finally:self._syncing=False

    @staticmethod
    def _cell(table,row,col):
        return cell_text(table,row,col)

    def new_sketch(self):
        self._syncing=True
        try:
            self.nodes.setRowCount(0); self.walls.setRowCount(0); self.periodic.setRowCount(0); self.sketch.clear_geometry()
        finally:self._syncing=False
        self.result=None; self.model=None; self.update_geometry_inspector(); self.set_sketch_mode('node'); self._reset_history()

    def load_example_configuration(self,config=1):
        config=int(config); model,mats,lams=paper_hexagonal_case(config); self.materials=mats; self._syncing=True
        try:
            if hasattr(self,'example_config'):self.example_config.setCurrentIndex(max(0,min(13,config-1)))
            self.nodes.setRowCount(len(model.nodes))
            for i,(x,y) in enumerate(model.nodes):
                for j,v in enumerate([i+1,x,y]): self.nodes.setItem(i,j,table_item(v))
            self.walls.setRowCount(len(model.elements))
            for i,e in enumerate(model.elements):self._set_wall_row(i,i+1,e.node_i+1,e.node_j+1,e.weight,e.laminate+1,group=(model.wall_groups[i] if i<len(model.wall_groups) else ''))
            self.a1x.setValue(model.a1[0]); self.a1y.setValue(model.a1[1]); self.a2x.setValue(model.a2[0]); self.a2y.setValue(model.a2[1]); self.height.setValue(model.height); self.solver.setCurrentText('Cylindrical'); self.wall_beam.setCurrentIndex(0); self.periodic.setRowCount(0)
            rows=[]
            for sid,lam in enumerate(lams,1):
                for k,ply in enumerate(lam.plies,1): rows.append([sid,k,ply.material,ply.angle_deg,ply.thickness])
            self.plies.setRowCount(len(rows))
            for i,row in enumerate(rows):self._set_ply_row(i,*row)
        finally:self._syncing=False
        self.refresh_sketch_from_tables(); self.draw_geometry(model); self.update_geometry_inspector(); self._reset_history()

    def load_paper_example(self):
        # Compatibility alias for older callers/projects.
        self.load_example_configuration(1)

    def _geometry_rows(self):
        nodes=[]
        for r in range(self.nodes.rowCount()):
            x=self._cell(self.nodes,r,1); y=self._cell(self.nodes,r,2)
            if x=='' or y=='': raise ValueError(f'Node row {r+1} is incomplete')
            nodes.append([float(x),float(y)])
        elements=[]
        for r in range(self.walls.rowCount()):
            ni=self._cell(self.walls,r,1); nj=self._cell(self.walls,r,2)
            if ni=='' or nj=='': raise ValueError(f'Wall row {r+1} is incomplete')
            wt=float(self._cell(self.walls,r,3) or 1.0); elements.append(Element(int(float(ni))-1,int(float(nj))-1,0,wt))
        return np.asarray(nodes,float),elements

    def refresh_sketch_from_tables(self):
        if self._syncing:return
        try:
            nodes=[]
            for r in range(self.nodes.rowCount()): nodes.append((float(self._cell(self.nodes,r,1)),float(self._cell(self.nodes,r,2))))
            walls=[]
            for r in range(self.walls.rowCount()): walls.append((int(float(self._cell(self.walls,r,1))),int(float(self._cell(self.walls,r,2))),float(self._cell(self.walls,r,3) or 1),int(float(self._cell(self.walls,r,4) or 1))))
            self.sketch.set_geometry(nodes,walls)
        except Exception:
            pass
        self._refresh_wall_metrics(); self.update_geometry_inspector()

    def _tables_changed(self,*_):
        if self._syncing:return
        self.refresh_sketch_from_tables()

    def _sketch_changed(self):
        if self._syncing:return
        group_map={}
        for r in range(self.walls.rowCount()):
            try: group_map[tuple(sorted((int(float(self._cell(self.walls,r,1))),int(float(self._cell(self.walls,r,2))))))]=self._cell(self.walls,r,7)
            except Exception: pass
        nodes,walls=self.sketch.export_geometry(); self._syncing=True
        try:
            self.nodes.setRowCount(len(nodes))
            for i,(x,y) in enumerate(nodes):
                for c,v in enumerate([i+1,x,y]): self.nodes.setItem(i,c,table_item(f'{v:.12g}' if isinstance(v,float) else v))
            self.walls.setRowCount(len(walls))
            for i,(ni,nj,wt,sid) in enumerate(walls):
                grp=group_map.get(tuple(sorted((ni,nj))), '')
                self._set_wall_row(i,i+1,ni,nj,wt,sid,group=grp)
        finally:self._syncing=False
        self._refresh_wall_metrics(); self.update_geometry_inspector(); self._commit_history()

    def _refresh_wall_metrics(self):
        try:
            nodes,elements=self._geometry_rows(); s=analyse_geometry(nodes,elements,[self.a1x.value(),self.a1y.value()],[self.a2x.value(),self.a2y.value()],periodic_tol=self.tol.value() if hasattr(self,'tol') else 1e-8)
        except Exception:return
        self._syncing=True
        try:
            by_id={w.wall_id:w for w in s.walls}
            for r in range(self.walls.rowCount()):
                w=by_id.get(r+1)
                self.walls.setItem(r,5,table_item(f'{w.length:.8g}' if w else '')); self.walls.setItem(r,6,table_item(f'{w.angle_deg:.8g}' if w else ''))
        finally:self._syncing=False

    def _set_lattice_vector(self,kind,dx,dy):
        if kind=='a1': self.a1x.setValue(dx); self.a1y.setValue(dy)
        elif kind=='a2': self.a2x.setValue(dx); self.a2y.setValue(dy)
        self.update_geometry_inspector()

    def current_geometry_summary(self):
        nodes,elements=self._geometry_rows(); return analyse_geometry(nodes,elements,[self.a1x.value(),self.a1y.value()],[self.a2x.value(),self.a2y.value()],periodic_tol=self.tol.value() if hasattr(self,'tol') else 1e-8)

    def update_geometry_inspector(self):
        # Spin boxes emit valueChanged while examples/projects are being populated.
        # Defer geometry/scene work until the synchronized load is complete so the
        # inspector only sees a fully initialized model.
        if self._syncing or not hasattr(self,'geometry_info'):
            return
        self.sketch.set_lattice_vectors([self.a1x.value(),self.a1y.value()],[self.a2x.value(),self.a2y.value()])
        try:
            s=self.current_geometry_summary(); self.geometry_info.setPlainText(summary_text(s,self.height.value()))
            if s.detected_geometry_id:
                self.detected_family.setText(f'Family: {s.detected_geometry_id} — {s.detected_geometry_name}   ·   match {100*s.detection_confidence:.1f}%')
            else:self.detected_family.setText('Family: Custom / arbitrary periodic RVE')
            self.detected_params.setRowCount(len(s.detected_parameters))
            for rr,(name,value) in enumerate(s.detected_parameters.items()):
                src=s.parameter_sources.get(name,[]); unit='' if name=='R' else ' mm'
                for cc,text in enumerate([name,f'{value:.8g}{unit}',', '.join(f'W{i}' for i in src)]): self.detected_params.setItem(rr,cc,table_item(text))
                self.detected_params.item(rr,0).setData(Qt.ItemDataRole.UserRole,src)
            errors=sum(i.severity=='error' for i in s.issues); warnings=sum(i.severity=='warning' for i in s.issues)
            if errors:
                self.geometry_state.setText(f'✕ {errors} error(s), {warnings} warning(s)'); self.geometry_state.setStyleSheet(f"font-weight:700;color:{COLORS['red']};")
            elif warnings:
                self.geometry_state.setText(f'✓ Geometry usable — {warnings} warning(s)'); self.geometry_state.setStyleSheet(f"font-weight:700;color:{COLORS['orange']};")
            else:
                self.geometry_state.setText('✓ Geometry checks passed'); self.geometry_state.setStyleSheet(f"font-weight:700;color:{COLORS['green']};")
        except Exception as e:
            self.geometry_info.setPlainText(f'Geometry is incomplete.\n\n{e}'); self.geometry_state.setText('Geometry incomplete'); self.geometry_state.setStyleSheet(f"font-weight:700;color:{COLORS['muted']};")

    def highlight_detected_parameter(self):
        items=self.detected_params.selectedItems()
        if not items:return
        row=items[0].row(); item=self.detected_params.item(row,0); src=item.data(Qt.ItemDataRole.UserRole) if item else []
        self.sketch.highlight_walls(src or [])

    def _selected_wall_id(self):
        ids=self.sketch.selected_wall_ids()
        if ids:return ids[0]
        r=self.walls.currentRow(); return r+1 if r>=0 else None

    def apply_wall_constraint(self,kind):
        wid=self._selected_wall_id()
        if not wid:
            QMessageBox.information(self,'Sketch constraint','Select a wall in the sketch or Wall table first.'); return
        value=None
        if kind=='length':
            try:cur=float(self._cell(self.walls,wid-1,5) or 1.0)
            except Exception:cur=1.0
            value,ok=QInputDialog.getDouble(self,'Set wall length',f'Length of W{wid} [mm]:',cur,1e-9,1e9,6)
            if not ok:return
        if not self.sketch.constrain_wall(wid,kind,value): QMessageBox.warning(self,'Sketch constraint','Constraint could not be applied.')

    def _capture_sketch_state(self):
        try:
            nodes,walls=self.sketch.export_geometry()
            return (tuple((round(float(x),12),round(float(y),12)) for x,y in nodes), tuple((int(a),int(b),float(w),int(s)) for a,b,w,s in walls))
        except Exception:return None

    def _reset_history(self):
        state=self._capture_sketch_state(); self._history=[state] if state is not None else []; self._redo=[]

    def _commit_history(self):
        if self._history_suspended:return
        state=self._capture_sketch_state()
        if state is not None and (not self._history or state!=self._history[-1]): self._history.append(state); self._redo=[]
        if len(self._history)>60:self._history=self._history[-60:]

    def _restore_sketch_state(self,state):
        if not state:return
        nodes,walls=state; self._syncing=True
        try:self.sketch.set_geometry(nodes,walls)
        finally:self._syncing=False
        self._sketch_changed()

    def undo_sketch(self):
        if len(self._history)<2:return
        cur=self._history.pop(); self._redo.append(cur); state=self._history[-1]
        nodes,walls=state; self._history_suspended=True; self._syncing=True
        try:self.sketch.set_geometry(nodes,walls)
        finally:self._syncing=False
        try:self._sketch_changed()
        finally:self._history_suspended=False

    def redo_sketch(self):
        if not self._redo:return
        state=self._redo.pop(); nodes,walls=state; self._history_suspended=True; self._syncing=True
        try:self.sketch.set_geometry(nodes,walls)
        finally:self._syncing=False
        try:self._sketch_changed()
        finally:self._history_suspended=False
        if not self._history or self._history[-1]!=state:self._history.append(state)

    def auto_wall_groups(self):
        try:s=self.current_geometry_summary()
        except Exception as e: QMessageBox.warning(self,'Wall groups',str(e)); return
        groups={}
        if s.detected_geometry_id:
            for param,ids in s.parameter_sources.items():
                for wid in ids: groups[wid]=param
        for w in s.walls:
            if w.wall_id not in groups:
                # Direction family independent of wall orientation sign.
                a=(w.angle_deg%180.0); bucket=round(a/5.0)*5.0; groups[w.wall_id]=f'{bucket:g}° family'
        self._syncing=True
        try:
            for r in range(self.walls.rowCount()): self.walls.setItem(r,7,table_item(groups.get(r+1,'')))
        finally:self._syncing=False
        QMessageBox.information(self,'Wall groups','Wall groups were assigned from detected parameters where available, otherwise from wall-direction families.')

    def apply_stack_to_group(self):
        r=self.walls.currentRow()
        if r<0: QMessageBox.information(self,'Wall groups','Select a wall in the Wall table first.'); return
        grp=self._cell(self.walls,r,7); sid=self._cell(self.walls,r,4)
        if not grp: QMessageBox.information(self,'Wall groups','The selected wall has no group. Use Auto wall groups or enter a group name.'); return
        self._syncing=True
        try:
            for rr in range(self.walls.rowCount()):
                if self._cell(self.walls,rr,7)==grp:set_combo_cell(self.walls,rr,4,self._stack_choices(),sid,callback=lambda *_:self._tables_changed())
        finally:self._syncing=False
        QMessageBox.information(self,'Wall groups',f'Stack {sid} applied to group “{grp}”.')

    def validate_geometry_action(self):
        try:
            s=self.current_geometry_summary(); QMessageBox.information(self,'RVE geometry validation',summary_text(s,self.height.value()))
        except Exception as e: QMessageBox.warning(self,'RVE geometry validation',str(e))

    def populate_detected_pairs(self):
        try:s=self.current_geometry_summary()
        except Exception as e: QMessageBox.warning(self,'Periodic pairs',str(e)); return
        pairs=s.periodic_pairs if s.periodic_pairs is not None else np.empty((0,4))
        self._syncing=True
        try:
            self.periodic.setRowCount(len(pairs))
            for i,row in enumerate(pairs):
                for c,v in enumerate([int(row[0])+1,int(row[1])+1,int(row[2]),int(row[3])]): self.periodic.setItem(i,c,table_item(v))
        finally:self._syncing=False
        self.update_geometry_inspector(); self.sketch_status.setText(f'Added {len(pairs)} detected periodic pair(s).')

    def parse_model(self):
        nodes=[]
        for r in range(self.nodes.rowCount()): nodes.append([float(self._cell(self.nodes,r,1)),float(self._cell(self.nodes,r,2))])
        stack_rows={}
        for r in range(self.plies.rowCount()):
            sid=int(float(self._cell(self.plies,r,0))); mat=self._cell(self.plies,r,2); ang=float(self._cell(self.plies,r,3)); th=float(self._cell(self.plies,r,4)); stack_rows.setdefault(sid,[]).append((int(float(self._cell(self.plies,r,1))),Ply(mat,ang,th)))
        laminates=[]; sid_to_idx={}
        for sid in sorted(stack_rows): sid_to_idx[sid]=len(laminates); laminates.append(Laminate([x[1] for x in sorted(stack_rows[sid],key=lambda x:x[0])],f'Stack {sid}'))
        elements=[]
        for r in range(self.walls.rowCount()):
            ni=int(float(self._cell(self.walls,r,1)))-1; nj=int(float(self._cell(self.walls,r,2)))-1; wt=float(self._cell(self.walls,r,3)); sid=int(float(self._cell(self.walls,r,4)))
            if sid not in sid_to_idx: raise ValueError(f'Wall {r+1} references Stack ID {sid}, but that laminate stack is not defined.')
            elements.append(Element(ni,nj,sid_to_idx[sid],wt))
        pairs=[]
        for r in range(self.periodic.rowCount()):
            if not self._cell(self.periodic,r,0): continue
            pairs.append([int(float(self._cell(self.periodic,r,c)))-1 if c<2 else float(self._cell(self.periodic,r,c)) for c in range(4)])
        periodic_pairs=np.array(pairs,float) if pairs else None
        lib_materials=material_mapping(); custom={k:v for k,v in self.materials.items() if k not in lib_materials}; lib_materials.update(custom); self.materials=lib_materials
        method=self.solver.currentText().strip().lower()
        if method.startswith('free'):
            wall_type='anisotropicStrip'; strip_mode='freeMoment'
        elif method.startswith('twist'):
            wall_type='anisotropicStrip'; strip_mode='twistRelease'
        elif method.startswith('aniso'):
            wall_type='anisotropicStrip'; strip_mode='cylindrical'
        else:
            wall_type='paper'; strip_mode='cylindrical'
        model=RVEModel(np.array(nodes),elements,np.array([self.a1x.value(),self.a1y.value()]),np.array([self.a2x.value(),self.a2y.value()]),self.height.value(),periodic_pairs=periodic_pairs,periodic_tol=self.tol.value(),anchor_node=self.anchor.value()-1,wall_element_type=wall_type,strip_bending_mode=strip_mode,wall_groups=[self._cell(self.walls,r,7) for r in range(self.walls.rowCount())],wall_shear=self.wall_beam.currentData() or 'euler-bernoulli',twist_release=self.twist_coefficients())
        return model,self.materials,laminates

    def set_project(self,model,materials,laminates):
        self.materials=materials; self._syncing=True
        try:
            self.nodes.setRowCount(len(model.nodes))
            for i,(x,y) in enumerate(model.nodes):
                for j,v in enumerate([i+1,x,y]): self.nodes.setItem(i,j,table_item(v))
            self.walls.setRowCount(len(model.elements))
            for i,e in enumerate(model.elements):self._set_wall_row(i,i+1,e.node_i+1,e.node_j+1,e.weight,e.laminate+1,group=(model.wall_groups[i] if i<len(model.wall_groups) else ''))
            self.a1x.setValue(model.a1[0]); self.a1y.setValue(model.a1[1]); self.a2x.setValue(model.a2[0]); self.a2y.setValue(model.a2[1]); self.height.setValue(model.height)
            _aniso=model.wall_element_type.lower().replace('_','') in ('anisotropicstrip','anisotropic','fullanisotropic'); _sm=str(model.strip_bending_mode).lower().replace('_',''); method='Free moment' if _sm=='freemoment' else ('Twist-corrected' if _aniso and _sm in ('twistrelease','twistcorrected','twist') else ('Anisotropic strip' if _aniso else 'Cylindrical'))
            self.solver.setCurrentText(method); [sp.setValue(v) for sp,v in zip(self.twist_coeff,(getattr(model,'twist_release',None) or __import__('rvestudio.laminate',fromlist=['x']).TWIST_RELEASE_DEFAULT))]; k=self.wall_beam.findData('timoshenko' if str(getattr(model,'wall_shear','euler-bernoulli')).lower().startswith('timo') else 'euler-bernoulli'); self.wall_beam.setCurrentIndex(max(k,0)); self.tol.setValue(model.periodic_tol); self.anchor.setValue(model.anchor_node+1)
            self.periodic.setRowCount(0)
            if model.periodic_pairs is not None:
                pp=np.asarray(model.periodic_pairs); self.periodic.setRowCount(len(pp))
                for i,row in enumerate(pp):
                    for j,v in enumerate([row[0]+1,row[1]+1,row[2],row[3]]): self.periodic.setItem(i,j,table_item(v))
            rows=[]
            for sid,lam in enumerate(laminates,1):
                for k,ply in enumerate(lam.plies,1): rows.append([sid,k,ply.material,ply.angle_deg,ply.thickness])
            self.plies.setRowCount(len(rows))
            for i,row in enumerate(rows):self._set_ply_row(i,*row)
        finally:self._syncing=False
        self.refresh_sketch_from_tables(); self.draw_geometry(model); self.update_geometry_inspector(); self._reset_history()

    def run_model(self):
        try:
            s=self.current_geometry_summary(); errors=[i.message for i in s.issues if i.severity=='error']
            if errors:
                self.tabs.setCurrentIndex(0); raise ValueError('Geometry validation failed:\n• '+'\n• '.join(errors))
            self.status.setText('Solving…'); QApplication.processEvents(); model,mats,lams=self.parse_model(); out=homogenize_rve(model,mats,lams); self.model=model; self.laminates=lams; self.result=out; p=out.props
            labels6=['11','22','33','23','13','12']; fill_matrix_table(self.ch_table,out.CH,labels6); fill_matrix_table(self.s_table,out.S,labels6)
            fill_matrix_table(self.qplane_table,p.Qplane,['11','22','12']); fill_matrix_table(self.qshear_table,p.Qshear,['13 (Q44)','23 (Q55)'])
            eng=[('E1',p.E1,'MPa'),('E2',p.E2,'MPa'),('E3',p.E3,'MPa'),('G12',p.G12,'MPa'),('G13',p.G13,'MPa'),('G23',p.G23,'MPa'),('ν12',p.nu12,'—'),('ν21',p.nu21,'—'),('ν13',p.nu13,'—'),('ν31',p.nu31,'—'),('ν23',p.nu23,'—'),('ν32',p.nu32,'—')]
            self.eng_table.setRowCount(len(eng))
            for r,row in enumerate(eng):
                for c,v in enumerate([row[0],f'{row[1]:.9g}',row[2]]):self.eng_table.setItem(r,c,table_item(v))
            phys=[('RVE area',out.area,'mm²'),('RVE volume',out.volume,'mm³'),('Relative density',out.relative_density,'—'),('Mass',out.mass_kg,'kg'),('Average density',out.physical.average_density_kg_m3,'kg/m³'),('Solid volume',out.physical.solid_volume,'mm³'),('Void volume',out.physical.void_volume,'mm³')]
            self.phys_table.setRowCount(len(phys))
            for r,row in enumerate(phys):
                for c,v in enumerate([row[0],f'{row[1]:.9g}',row[2]]):self.phys_table.setItem(r,c,table_item(v))
            self.draw_geometry(model); rel=assess_reliability(out)
            self.status.setText(f"Solved successfully ({self.wall_beam.currentText().split(' (')[0]} walls, {self.solver.currentText().lower()} bending). Model reliability — {rel.short()}")
            self.status.setToolTip('\n'.join(rel.messages) if rel.messages else 'All walls are inside the verified range of the wall model (Theory guide, Section 9).'); self.tabs.setCurrentIndex(3)
        except Exception as e:self.status.setText('Failed'); QMessageBox.critical(self,'RVE solve failed',str(e))

    def draw_geometry(self,model):
        self.fig.clear(); ax=self.fig.add_subplot(111); nodes=model.nodes
        for k,e in enumerate(model.elements,1):
            a,b=nodes[e.node_i],nodes[e.node_j]; ax.plot([a[0],b[0]],[a[1],b[1]],'-o',lw=2); ax.text(*(0.5*(a+b)),str(k),fontsize=8)
        for i,p in enumerate(nodes,1): ax.text(p[0],p[1],f' {i}',fontsize=8)
        o=np.zeros(2); ax.arrow(*o,*model.a1,length_includes_head=True,head_width=max(np.linalg.norm(model.a1),1)*0.02,alpha=.5); ax.arrow(*o,*model.a2,length_includes_head=True,head_width=max(np.linalg.norm(model.a2),1)*0.02,alpha=.5); ax.set_aspect('equal','box'); ax.grid(True,alpha=.2); ax.set_xlabel('x1'); ax.set_ylabel('x2'); self.canvas.draw_idle()


class LaminateToolsPage(QWidget):
    def __init__(self,main):
        super().__init__(); self.main=main; self.materials=material_mapping(); v=QVBoxLayout(self); v.setContentsMargins(18,15,18,15); v.addWidget(title_block('Laminate Tools','Standalone ABD, material-invariant, lamination-parameter and strip analysis using the same routines as the RVE solver.'))
        b=QHBoxLayout(); add=QPushButton('Add ply'); add.clicked.connect(self.add_ply); b.addWidget(add); delete=QPushButton('Delete ply'); delete.clicked.connect(self.delete_ply); b.addWidget(delete); calc=QPushButton('Calculate laminate'); calc.setObjectName('Primary'); calc.clicked.connect(self.calculate); b.addWidget(calc); b.addStretch(); v.addLayout(b)
        sp=QSplitter(); self.table=QTableWidget(5,4); self.table.setHorizontalHeaderLabels(['Ply','Material','Angle [deg]','Thickness [mm]']); self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch); defaults=[45,-45,0,-45,45]
        for i,a in enumerate(defaults):self._set_row(i,i+1,'E-glass/Epoxy',a,0.6252703415)
        sp.addWidget(self.table); self.out=QPlainTextEdit(); self.out.setReadOnly(True); sp.addWidget(self.out); sp.setSizes([500,700]); v.addWidget(sp,1)

    def _set_row(self,r,ply=1,material='E-glass/Epoxy',angle=0,thickness=0.5):
        set_combo_cell(self.table,r,0,[str(i) for i in range(1,101)],ply)
        mats=list(self.materials.keys()); set_combo_cell(self.table,r,1,mats,material,editable=True)
        set_combo_cell(self.table,r,2,['-90','-75','-60','-45','-30','-15','0','15','30','45','60','75','90'],angle,editable=True)
        self.table.setItem(r,3,table_item(thickness))

    def add_ply(self):
        r=self.table.rowCount(); self.table.insertRow(r); self._set_row(r,r+1)

    def delete_ply(self):
        r=self.table.currentRow()
        if r>=0:self.table.removeRow(r)

    def calculate(self):
        try:
            plies=[]
            for r in range(self.table.rowCount()):
                material=cell_text(self.table,r,1)
                if not material:continue
                plies.append(Ply(material,float(cell_text(self.table,r,2)),float(cell_text(self.table,r,3))))
            out=analyze_laminate(Laminate(plies),self.materials); sec=out['section']; lp=out['lp']; text=f"A [N/mm]\n{matrix_text(sec.A)}\n\nB [N]\n{matrix_text(sec.B)}\n\nD classical [N mm]\n{matrix_text(sec.Doriginal)}\n\nD neutral [N mm]\n{matrix_text(sec.D)}\n\nNeutral offset d = {sec.d_neutral:.8g} mm\nSeq = {sec.Seq:.8g} N/mm\nShear correction = {sec.kappa_shear:.8g}\n\nxiA = {np.array2string(lp['xiA'],precision=8)}\nxiD = {np.array2string(lp['xiD'],precision=8)}\n\nDstrip cylindrical = {out['strip_cyl'].Dstrip:.8g}\nDstrip free moment = {out['strip_free'].Dstrip:.8g}\ncurvature ratio [k33/k11,k13/k11] = {np.array2string(out['strip_free'].curvature_ratio,precision=8)}\n\nMaterial invariants\n"
            for row in out['invariants']: text+=f"{row['material']}: U={np.array2string(row['U'],precision=8)}\n"
            self.out.setPlainText(text)
        except Exception as e: QMessageBox.critical(self,'Laminate analysis failed',str(e))


class MaterialLibraryPage(QWidget):
    HEADERS=['Active','Name','E1','E2','G12','nu12','Gth','Density kg/m3','Cost / kg (user units)','Source','Source URL']
    def __init__(self,main):
        super().__init__(); v=QVBoxLayout(self); v.setContentsMargins(18,15,18,15); v.addWidget(title_block('Material Library','Edit lamina elastic properties, transverse shear modulus, density and a user-defined cost per kg. Cost may be a currency/kg value or a normalized dummy cost index/kg.'))
        bar=QHBoxLayout()
        for text,fn in [('Add material',self.add),('Delete selected',self.delete),('Reload',self.reload),('Save library',self.save)]:
            b=QPushButton(text); b.clicked.connect(fn); bar.addWidget(b)
            if text=='Save library': b.setObjectName('Primary')
        bar.addStretch(); v.addLayout(bar); self.table=QTableWidget(0,len(self.HEADERS)); self.table.setHorizontalHeaderLabels(self.HEADERS); self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents); self.table.horizontalHeader().setStretchLastSection(True); v.addWidget(self.table,1); self.reload()

    def _set_row(self,i,r):
        set_combo_cell(self.table,i,0,['True','False'],'True' if r.get('active',True) else 'False')
        vals=[r.get('name',''),r.get('E1',''),r.get('E2',''),r.get('G12',''),r.get('nu12',''),r.get('Gth',''),r.get('density_kg_m3',''),r.get('cost_per_kg',''),r.get('source',''),r.get('source_url','')]
        for j,val in enumerate(vals,1):self.table.setItem(i,j,table_item(val))

    def reload(self):
        rows=load_material_library(); self.table.setRowCount(len(rows))
        for i,r in enumerate(rows):self._set_row(i,r)

    def add(self):
        r=self.table.rowCount(); self.table.insertRow(r); self._set_row(r,{'active':True,'name':f'New Material {r+1}','E1':1,'E2':1,'G12':1,'nu12':.3,'Gth':1})

    def delete(self):
        r=self.table.currentRow()
        if r>=0:self.table.removeRow(r)

    def save(self):
        try:
            rows=[]
            for r in range(self.table.rowCount()):
                txt=lambda c:cell_text(self.table,r,c)
                name=txt(1)
                if not name:continue
                density=float(txt(7)) if txt(7) else None; cost=float(txt(8)) if txt(8) else None
                rows.append({'active':txt(0).lower()=='true','name':name,'E1':float(txt(2)),'E2':float(txt(3)),'G12':float(txt(4)),'nu12':float(txt(5)),'Gth':float(txt(6) or txt(4)),'density_kg_m3':density,'cost_per_kg':cost,'source':txt(9),'source_url':txt(10)})
            save_json(PACKAGE_ROOT/'data'/'materials.json',rows); QMessageBox.information(self,'Material Library','Saved successfully.')
        except Exception as e: QMessageBox.critical(self,'Cannot save materials',str(e))


class GeometryLibraryPage(QWidget):
    HEADERS=['Active','ID','Name','Independent parameters','Sampling ranges','Weight','Method','Description']
    def __init__(self,main):
        super().__init__(); self.main=main; v=QVBoxLayout(self); v.setContentsMargins(18,15,18,15)
        v.addWidget(title_block('Geometry Library','Reusable geometry families. Parameter names are shown explicitly; no generic P1/P2 labels are used. Create new families with the Geometry Definition Wizard.'))
        guide=QLabel('<b>Workflow:</b> draw or load an RVE → Geometry Inspector identifies known parameters and their source walls → open Geometry Library → <b>Create from current sketch</b> → confirm independent variables and ranges → save. Unknown topologies can be parameterised by their drawing width and height without editing Python code.')
        guide.setWordWrap(True); guide.setStyleSheet(f"color:{COLORS['muted']};"); v.addWidget(guide)
        bar=QHBoxLayout()
        actions=[('Add definition',self.add),('Create from current sketch',self.capture_current_sketch),('Edit parameters',self.edit_parameters),('Parameter guide',self.parameter_guide),('Delete selected',self.delete),('Reload',self.reload),('Save library',self.save)]
        for text,fn in actions:
            b=QPushButton(text); b.clicked.connect(fn); bar.addWidget(b)
            if text=='Save library':b.setObjectName('Primary')
        op=QPushButton('Open plugins folder'); op.clicked.connect(lambda:QDesktopServices.openUrl(QUrl.fromLocalFile(str(PACKAGE_ROOT/'plugins')))); bar.addWidget(op); bar.addStretch(); v.addLayout(bar)
        self.table=QTableWidget(0,len(self.HEADERS)); self.table.setHorizontalHeaderLabels(self.HEADERS); self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents); self.table.horizontalHeader().setStretchLastSection(True)
        # Parameter diagram of the selected family: what each independent parameter controls.
        split=QSplitter(Qt.Orientation.Horizontal); split.addWidget(self.table)
        pv=QWidget(); pl=QVBoxLayout(pv); pl.setContentsMargins(8,0,0,0)
        head=QLabel('<b>Parameter diagram</b> — the selected family at its reference values (middle of each range). Coloured walls are set by the parameter of that colour; W1, W2… are the wall numbers used by the Geometry Inspector.'); head.setWordWrap(True); pl.addWidget(head)
        self.diagram_fig=Figure(figsize=(6.4,4.6)); self.diagram_canvas=FigureCanvas(self.diagram_fig); self.diagram_canvas.setMinimumHeight(320); pl.addWidget(self.diagram_canvas,1)
        self.diagram_caption=QLabel(''); self.diagram_caption.setWordWrap(True); self.diagram_caption.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse); pl.addWidget(self.diagram_caption)
        split.addWidget(pv); split.setSizes([760,560]); v.addWidget(split,1)
        self.table.currentCellChanged.connect(lambda *_:self.show_diagram())
        self.reload()

    def _row_entry(self,r):
        return {'id':cell_text(self.table,r,1).upper(),'name':cell_text(self.table,r,2),'parameters':self._parameters(r),'description':cell_text(self.table,r,7)}

    def show_diagram(self):
        """Draw the parameter diagram of the selected library row (never raises: a bad plug-in only shows a message)."""
        from ..geometry_diagram import entry_diagram, caption_html
        r=self.table.currentRow()
        self.diagram_fig.clear()
        if r<0 or not cell_text(self.table,r,1):
            self.diagram_caption.setText('Select a geometry row to see what its parameters control.'); self.diagram_canvas.draw_idle(); return
        entry=self._row_entry(r)
        try:
            _fig,info=entry_diagram(entry,discover_geometry_plugins(PACKAGE_ROOT/'plugins'),fig=self.diagram_fig,fontsize=9)
            self.diagram_caption.setText(caption_html(info))
        except Exception as exc:
            self.diagram_fig.clear(); self.diagram_caption.setText(f"<span style='color:#B03A2E'>No diagram for {entry['id']}: {exc}</span>")
        self.diagram_canvas.draw_idle()

    def _method_from_row(self,r):
        solver=str(r.get('solver','paper')).lower().replace('_',''); strip=str(r.get('strip_mode','cylindrical')).lower().replace('_','')
        return 'Free moment' if solver in ('anisotropicstrip','anisotropic','fullanisotropic') or strip=='freemoment' else 'Cylindrical'

    @staticmethod
    def _range_text(ps):
        return '; '.join(f"{p.get('name','')}: {float(p.get('min',0)):.6g} … {float(p.get('max',0)):.6g}" for p in ps)

    def _set_row(self,i,r):
        set_combo_cell(self.table,i,0,['True','False'],'True' if r.get('active',True) else 'False')
        self.table.setItem(i,1,table_item(r.get('id',''))); self.table.setItem(i,2,table_item(r.get('name','')))
        ps=[dict(p) for p in r.get('parameters',[])]
        pit=table_item(', '.join(p.get('name','') for p in ps) if ps else 'Fixed geometry'); pit.setData(Qt.ItemDataRole.UserRole,ps); pit.setFlags(pit.flags() & ~Qt.ItemFlag.ItemIsEditable); self.table.setItem(i,3,pit)
        rit=table_item(self._range_text(ps) if ps else '—'); rit.setFlags(rit.flags() & ~Qt.ItemFlag.ItemIsEditable); self.table.setItem(i,4,rit)
        self.table.setItem(i,5,table_item(r.get('weight',1)))
        set_combo_cell(self.table,i,6,['Cylindrical','Free moment'],self._method_from_row(r))
        self.table.setItem(i,7,table_item(r.get('description','')))

    def _parameters(self,row):
        it=self.table.item(row,3); data=it.data(Qt.ItemDataRole.UserRole) if it else []
        return [dict(p) for p in (data or [])]

    def _set_parameters(self,row,ps):
        it=table_item(', '.join(p.get('name','') for p in ps) if ps else 'Fixed geometry'); it.setData(Qt.ItemDataRole.UserRole,[dict(p) for p in ps]); it.setFlags(it.flags() & ~Qt.ItemFlag.ItemIsEditable); self.table.setItem(row,3,it); rit=table_item(self._range_text(ps) if ps else '—'); rit.setFlags(rit.flags() & ~Qt.ItemFlag.ItemIsEditable); self.table.setItem(row,4,rit)

    def reload(self):
        rows=load_geometry_library(); self.table.setRowCount(len(rows))
        for i,r in enumerate(rows):self._set_row(i,r)
        if hasattr(self.main,'single_page'): self.main.single_page.refresh_geometry_choices()
        if rows and self.table.currentRow()<0: self.table.setCurrentCell(0,1)
        if hasattr(self,'diagram_fig'): self.show_diagram()

    def add(self):
        r=self.table.rowCount(); self.table.insertRow(r); self._set_row(r,{'active':False,'id':f'CUSTOM{r+1}','name':'Custom periodic RVE','parameters':[],'weight':1,'description':'Use Create from current sketch for a guided parametric definition.'})

    def delete(self):
        r=self.table.currentRow()
        if r>=0:self.table.removeRow(r)

    def edit_parameters(self):
        r=self.table.currentRow()
        if r<0: QMessageBox.information(self,'Geometry parameters','Select a geometry row first.'); return
        dlg=ParameterEditorDialog(self._parameters(r),self)
        if dlg.exec(): self._set_parameters(r,dlg.parameters()); self.show_diagram()

    def parameter_guide(self):
        """Diagrams of the built-in families and of the current sketch, with the parameter evidence."""
        from ..geometry_diagram import builtin_guide_figure, draw_geometry_diagram
        from ..geometry import Geometry
        dlg=QDialog(self); dlg.setWindowTitle('Geometry parameters'); dlg.resize(1180,860); lay=QVBoxLayout(dlg)
        tabs=QTabWidget(); lay.addWidget(tabs,1)
        w1=QWidget(); l1=QVBoxLayout(w1); f1=builtin_guide_figure(Figure(figsize=(12,9))); l1.addWidget(FigureCanvas(f1),1)
        t1=QLabel('Coloured walls are set by the parameter of that colour; W1, W2… are wall numbers. '
                  'HEX — L1: wall normal to a1; L2: inclined walls. SQUARE — Lx, Ly: wall segments parallel and normal to a1 '
                  '(junction to cell edge, so the junction-to-junction wall is 2Lx, 2Ly). MIXED and TRI — Lh: reference walls '
                  'parallel to a1; R = Ld/Lh with Ld the inclined-wall length. The dashed outline is the unit cell (a1, a2), '
                  'centred on the geometry.'); t1.setWordWrap(True); l1.addWidget(t1)
        tabs.addTab(w1,'Built-in families')
        w2=QWidget(); l2=QVBoxLayout(w2); text=''
        try:
            s=self.main.single_page.current_geometry_summary(); nodes,walls=self.main.single_page.sketch.export_geometry()
            sp=self.main.single_page; a1=[sp.a1x.value(),sp.a1y.value()]; a2=[sp.a2x.value(),sp.a2y.value()]
            conn=np.array([[int(ni)-1,int(nj)-1] for ni,nj,_,_ in walls],int); wts=np.array([float(wt) for _,_,wt,_ in walls],float)
            geom=Geometry(s.detected_geometry_id or 'SKETCH','Current sketch',np.asarray(nodes,float),conn,np.asarray(a1,float),np.asarray(a2,float),[f'W{i+1}' for i in range(len(conn))],wts)
            vals=dict(s.detected_parameters) if s.detected_geometry_id else {}
            f2=Figure(figsize=(8,6)); draw_geometry_diagram(f2.add_subplot(111),geom,vals,title=f"Current sketch: {s.detected_geometry_id or 'custom'} — {s.detected_geometry_name}")
            f2.subplots_adjust(left=0.02,right=0.72,top=0.93,bottom=0.03); l2.addWidget(FigureCanvas(f2),1)
            if s.detected_geometry_id:
                text=(f"Recognised as {s.detected_geometry_id} ({100*s.detection_confidence:.1f}% match). "
                      +'; '.join(f"{k} = {v:.6g} from "+', '.join(f'W{i}' for i in s.parameter_sources.get(k,[])) for k,v in s.detected_parameters.items()))
            else:
                text='Custom / arbitrary periodic RVE. Create from current sketch makes a family whose independent drivers are the drawing width W and height H: node coordinates and lattice vectors are normalised to this sketch and scaled by W/H.'
        except Exception as exc:
            text=f'No valid sketch in the Single RVE Designer to draw ({exc}).'
        t2=QLabel(text); t2.setWordWrap(True); l2.addWidget(t2)
        tabs.addTab(w2,'Current sketch')
        bb=QDialogButtonBox(QDialogButtonBox.StandardButton.Close); bb.rejected.connect(dlg.reject); bb.accepted.connect(dlg.accept); lay.addWidget(bb)
        dlg.exec()

    def parameter_guide_text(self):
        """Former text-only guide (kept for reference and scripting)."""
        detected=''
        try:
            s=self.main.single_page.current_geometry_summary()
            if s.detected_geometry_id:
                params=', '.join(f'{k}={v:.6g}' for k,v in s.detected_parameters.items()); sources='; '.join(f"{k}: "+','.join(f'W{i}' for i in s.parameter_sources.get(k,[])) for k in s.detected_parameters)
                detected=f"\n\nCurrent sketch: {s.detected_geometry_id} — {s.detected_geometry_name} ({100*s.detection_confidence:.1f}% match)\nExtracted: {params}\nEvidence: {sources}"
            else:detected='\n\nCurrent sketch: custom / arbitrary periodic RVE. The wizard can use drawing width and height as explicit independent drivers.'
        except Exception:pass
        QMessageBox.information(self,'Geometry parameters',
            'Built-in families:\n\n'
            'HEX — L1: wall family normal to a1; L2: inclined-wall family.\n'
            'SQUARE — Lx: wall family parallel to a1; Ly: wall family normal to a1.\n'
            'MIXED — Lh: reference wall family parallel to a1; R = Ld/Lh.\n'
            'TRI — Lh: reference wall family parallel to a1; R = Ld/Lh.\n\n'
            'For a new topology, Create from current sketch opens a wizard. Its default custom parameterisation uses W (drawing width) and H (drawing height). Node coordinates and lattice vectors are normalized to the reference sketch and scaled by W/H. This makes the parameter meaning reproducible and visible in the library.'+detected)

    def capture_current_sketch(self):
        try:
            summary=self.main.single_page.current_geometry_summary(); errors=[i.message for i in summary.issues if i.severity=='error']
            if errors:raise ValueError('The current sketch is not ready to capture:\n• '+'\n• '.join(errors))
            nodes,walls=self.main.single_page.sketch.export_geometry()
        except Exception as e: QMessageBox.warning(self,'Create geometry from sketch',str(e)); return
        existing_ids=[cell_text(self.table,r,1).upper() for r in range(self.table.rowCount()) if cell_text(self.table,r,1)]
        dlg=GeometryDefinitionWizard(summary,self,existing_ids=existing_ids)
        if not dlg.exec():return
        try:definition=dlg.definition()
        except Exception as e: QMessageBox.warning(self,'Create geometry from sketch',str(e)); return
        gid=definition['id']; name=definition['name']; mode=definition['mode']; params=[{'name':p['name'],'min':p['min'],'max':p['max'],'reference':p.get('reference')} for p in definition['parameters']]
        existing={cell_text(self.table,r,1).upper():r for r in range(self.table.rowCount())}
        row=existing.get(gid)
        if mode=='builtin' and summary.detected_geometry_id:
            # Built-in generators already exist; the wizard only defines readable sampling ranges.
            gid=summary.detected_geometry_id
            row=existing.get(gid)
            entry={'active':True,'id':gid,'name':name,'parameters':params,'weight':1,'solver':'paper','strip_mode':'cylindrical','description':f'Built-in {summary.detected_geometry_name}; parameters confirmed by Geometry Definition Wizard.'}
        else:
            if row is not None:
                QMessageBox.warning(self,'Create geometry from sketch',f'Geometry ID {gid} already exists. Choose a different ID for a custom family.'); return
            pts=np.asarray(nodes,float); xmin,ymin=np.min(pts,axis=0); width=float(summary.width); height=float(summary.height_2d)
            if width<=1e-12 or height<=1e-12: QMessageBox.warning(self,'Create geometry from sketch','Custom width/height parameterisation requires a non-zero drawing width and height.'); return
            conn=[[int(ni)-1,int(nj)-1] for ni,nj,_,_ in walls]; weights=[float(wt) for _,_,wt,_ in walls]
            groups=[]
            for i in range(len(conn)):
                try:g=self.main.single_page._cell(self.main.single_page.walls,i,7)
                except Exception:g=''
                groups.append(g or f'wall-{i+1}')
            a1=[self.main.single_page.a1x.value(),self.main.single_page.a1y.value()]; a2=[self.main.single_page.a2x.value(),self.main.single_page.a2y.value()]
            # The first two wizard drivers are drawing width and height; names may be edited.
            pdefs=definition['parameters']; wname=next((p['name'] for p in pdefs if 'width' in p.get('driver','').lower()),pdefs[0]['name']); hname=next((p['name'] for p in pdefs if 'height' in p.get('driver','').lower()),pdefs[-1]['name'])
            plugin=PACKAGE_ROOT/'plugins'/f'{gid.lower()}_geometry.py'
            code=make_scaled_xy_plugin_source(
                geometry_id=gid,name=name,nodes=pts,connectivity=conn,a1=a1,a2=a2,
                wall_groups=groups,wall_weights=weights,width_parameter=wname,height_parameter=hname,
            )
            plugin.write_text(code,encoding='utf-8')
            # Validate the generated plug-in immediately.  Creation is not reported
            # as successful unless the file imports, generates geometry at the
            # reference values and passes the same geometry checks as the designer.
            try:
                registry=discover_geometry_plugins(PACKAGE_ROOT/'plugins')
                if gid not in registry: raise ValueError(f'Generated plug-in {gid} was not discovered.')
                refvals={p['name']:float(p.get('reference',0.5*(float(p['min'])+float(p['max'])))) for p in definition['parameters']}
                test_geom=registry[gid](refvals)
                test_elements=[Element(int(ii),int(jj),0,float(test_geom.wall_weights[k])) for k,(ii,jj) in enumerate(test_geom.connectivity)]
                check=analyse_geometry(np.asarray(test_geom.nodes,float),test_elements,np.asarray(test_geom.a1,float),np.asarray(test_geom.a2,float))
                errs=[i.message for i in check.issues if i.severity=='error']
                if errs: raise ValueError('Generated geometry did not validate: '+'; '.join(errs))
            except Exception as exc:
                try:plugin.unlink()
                except Exception:pass
                QMessageBox.critical(self,'Cannot create geometry',f'The generated geometry plug-in failed its validation check:\n\n{exc}')
                return
            entry={'active':True,'id':gid,'name':name,'parameters':params,'weight':1,'solver':'paper','strip_mode':'cylindrical','description':'Wizard-created custom family; independent drivers are drawing width and height.'}
        if row is None:
            row=self.table.rowCount(); self.table.insertRow(row)
        self._set_row(row,entry); self.table.setCurrentCell(row,1); self.save(show_message=False)
        self.main.single_page.refresh_geometry_choices(gid)
        QMessageBox.information(self,'Geometry created',f'{gid} is now in the Geometry Library and is available in the Existing geometry drop-down.\n\nIndependent parameters: '+(', '.join(p['name'] for p in params) or 'none')+'\nTheir ranges are visible directly in the library table.')

    def save(self,show_message=True):
        try:
            rows=[]
            for r in range(self.table.rowCount()):
                txt=lambda c:cell_text(self.table,r,c)
                if not txt(1):continue
                ps=self._parameters(r); method=txt(6).lower(); solver='anisotropicStrip' if method.startswith('free') else 'paper'; strip='freeMoment' if method.startswith('free') else 'cylindrical'
                rows.append({'active':txt(0).lower()=='true','id':txt(1).upper(),'name':txt(2),'parameters':ps,'weight':float(txt(5) or 1),'solver':solver,'strip_mode':strip,'description':txt(7)})
            save_json(PACKAGE_ROOT/'data'/'geometries.json',rows)
            if hasattr(self.main,'single_page'): self.main.single_page.refresh_geometry_choices()
            if show_message:QMessageBox.information(self,'Geometry Library','Saved successfully.')
        except Exception as e: QMessageBox.critical(self,'Cannot save geometries',str(e))


class JsonLibraryPage(QWidget):
    def __init__(self,main,kind):
        super().__init__(); self.main=main; self.kind=kind; self.path=PACKAGE_ROOT/'data'/f'{kind}.json'; v=QVBoxLayout(self); v.setContentsMargins(18,15,18,15); title={'materials':'Material Library','geometries':'Geometry Library','layups':'Layup Library'}[kind]; v.addWidget(title_block(title,'Editable JSON-backed library used by Generate DB. Changes are saved directly into the Python package data directory.'))
        b=QHBoxLayout(); load=QPushButton('Reload'); load.clicked.connect(self.reload); b.addWidget(load); save=QPushButton('Save'); save.setObjectName('Primary'); save.clicked.connect(self.save); b.addWidget(save); raw=QPushButton('Open JSON externally'); raw.clicked.connect(lambda:QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.path)))); b.addWidget(raw)
        if kind=='geometries':
            plug=QPushButton('Open geometry plugins folder'); plug.clicked.connect(lambda:QDesktopServices.openUrl(QUrl.fromLocalFile(str(PACKAGE_ROOT/'plugins')))); b.addWidget(plug)
        b.addStretch(); v.addLayout(b); self.edit=QPlainTextEdit(); self.edit.setStyleSheet('font-family:Consolas;'); v.addWidget(self.edit,1); self.reload()
    def reload(self): self.edit.setPlainText(self.path.read_text(encoding='utf-8'))
    def save(self):
        try: data=json.loads(self.edit.toPlainText()); save_json(self.path,data); QMessageBox.information(self,'Saved',f'{self.kind.title()} library saved.')
        except Exception as e: QMessageBox.critical(self,'Invalid JSON',str(e))


class GenerateDBPage(QWidget):
    def __init__(self,main):
        super().__init__(); v=QVBoxLayout(self); v.setContentsMargins(18,15,18,15); v.addWidget(title_block('Generate Database','Create large Latin-hypercube parametric databases from the active material, geometry and layup libraries. Different wall thicknesses are scaled to match target relative density.'))
        box=QGroupBox('Generation settings'); f=QFormLayout(box); self.n=QSpinBox(); self.n.setRange(1,2_000_000); self.n.setValue(5000); self.seed=QSpinBox(); self.seed.setRange(0,2_000_000_000); self.seed.setValue(20260906); self.rmin=QDoubleSpinBox(); self.rmin.setRange(.0001,.99); self.rmin.setValue(.03); self.rmax=QDoubleSpinBox(); self.rmax.setRange(.0001,.99); self.rmax.setValue(.15); self.output=QLineEdit(str(PACKAGE_ROOT/'data'/'RVE_5000_python.csv')); browse=QPushButton('Browse…'); browse.clicked.connect(self.browse); path=QWidget(); h=QHBoxLayout(path); h.setContentsMargins(0,0,0,0); h.addWidget(self.output); h.addWidget(browse); f.addRow('Cases',self.n); f.addRow('Seed',self.seed); f.addRow('rho* min',self.rmin); f.addRow('rho* max',self.rmax); f.addRow('Output',path); v.addWidget(box); self.wall_model=WallModelOptions('Wall model written to every case (beam theory and bending option)'); v.addWidget(self.wall_model); b=QPushButton('Generate DB'); b.setObjectName('Primary'); b.setMinimumHeight(44); b.clicked.connect(self.generate); v.addWidget(b); self.status=QPlainTextEdit(); self.status.setReadOnly(True); v.addWidget(self.status,1)
    def browse(self):
        f,_=QFileDialog.getSaveFileName(self,'Save database',self.output.text(),'CSV (*.csv);;Excel (*.xlsx);;Parquet (*.parquet)');
        if f:self.output.setText(f)
    def generate(self):
        try:
            self.status.setPlainText('Generating…'); QApplication.processEvents(); df=generate_parametric_database(self.output.text(),self.n.value(),self.seed.value(),(self.rmin.value(),self.rmax.value()),wall_shear=self.wall_model.wall_shear(),bending=self.wall_model.bending_option(),twist_release=self.wall_model.twist_release()); self.status.setPlainText(f'Generated {len(df):,} cases ({self.wall_model.describe()}).\n\nGeometry counts:\n{df.GeometryID.value_counts().to_string()}\n\nMaterial counts:\n{df.Material.value_counts().to_string()}\n\nSaved to:\n{self.output.text()}')
        except Exception as e: QMessageBox.critical(self,'Database generation failed',str(e))


class BatchPage(QWidget):
    def __init__(self,main):
        super().__init__(); self.thread=None; v=QVBoxLayout(self); v.setContentsMargins(18,15,18,15); v.addWidget(title_block('Large Batch Analysis','Solve compact Python databases sequentially or with multiple worker processes. Each failed RVE is isolated and reported without terminating the batch.'))
        box=QGroupBox('Batch files'); g=QGridLayout(box); self.input=QLineEdit(); self.output=QLineEdit(); bi=QPushButton('Browse input'); bo=QPushButton('Browse output'); bi.clicked.connect(self.browse_input); bo.clicked.connect(self.browse_output); g.addWidget(QLabel('Input'),0,0); g.addWidget(self.input,0,1); g.addWidget(bi,0,2); g.addWidget(QLabel('Output'),1,0); g.addWidget(self.output,1,1); g.addWidget(bo,1,2); self.workers=QSpinBox(); self.workers.setRange(1,max(1,os.cpu_count() or 1)); self.workers.setValue(1); self.limit=QSpinBox(); self.limit.setRange(0,2_000_000); self.limit.setSpecialValueText('All cases'); g.addWidget(QLabel('Workers'),2,0); g.addWidget(self.workers,2,1); g.addWidget(QLabel('Case limit'),3,0); g.addWidget(self.limit,3,1); v.addWidget(box)
        buck=QGroupBox('Optional local-wall buckling screening'); bg=QGridLayout(buck); self.batch_buck_enable=QCheckBox('Calculate local-wall buckling from recovered combined resultants for every solved RVE'); bg.addWidget(self.batch_buck_enable,0,0,1,6)
        self.batch_buck_spins=[]
        for i,lab in enumerate(MACRO_LABELS):
            sp=QDoubleSpinBox(); sp.setRange(-1e6,1e6); sp.setDecimals(6); sp.setSuffix(' MPa'); self.batch_buck_spins.append(sp); bg.addWidget(QLabel(('σ' if i<3 else 'τ')+lab),1+(i//3)*2,i%3); bg.addWidget(sp,2+(i//3)*2,i%3)
        self.batch_buck_spins[0].setValue(-1.0); self.batch_buck_spins[2].setValue(-1.0); self.batch_buck_spins[4].setValue(0.25)
        self.batch_buck_m=QSpinBox(); self.batch_buck_m.setRange(2,20); self.batch_buck_m.setValue(8); self.batch_buck_n=QSpinBox(); self.batch_buck_n.setRange(2,20); self.batch_buck_n.setValue(8); self.batch_buck_sf=QDoubleSpinBox(); self.batch_buck_sf.setRange(.01,100); self.batch_buck_sf.setValue(1.5); self.batch_buck_sf.setDecimals(3); self.batch_buck_abd=QCheckBox('ABD-condensed D*'); self.batch_buck_abd.setChecked(True)
        bg.addWidget(QLabel('m max'),5,0); bg.addWidget(self.batch_buck_m,5,1); bg.addWidget(QLabel('n max'),5,2); bg.addWidget(self.batch_buck_n,5,3); bg.addWidget(QLabel('Required λ'),5,4); bg.addWidget(self.batch_buck_sf,5,5); bg.addWidget(self.batch_buck_abd,6,0,1,3); v.addWidget(buck)
        self.run=QPushButton('Run batch'); self.run.setObjectName('Primary'); self.run.clicked.connect(self.start); v.addWidget(self.run); self.progress=QProgressBar(); v.addWidget(self.progress); self.log=QPlainTextEdit(); self.log.setReadOnly(True); v.addWidget(self.log,1)
    def browse_input(self):
        f,_=QFileDialog.getOpenFileName(self,'Open database','','Data (*.csv *.xlsx *.parquet)');
        if f: self.input.setText(f); p=Path(f); self.output.setText(str(p.with_name(p.stem+'_results.csv')))
    def browse_output(self):
        f,_=QFileDialog.getSaveFileName(self,'Save results',self.output.text(),'CSV (*.csv);;Excel (*.xlsx);;Parquet (*.parquet)');
        if f:self.output.setText(f)
    def start(self):
        if not self.input.text() or not self.output.text(): QMessageBox.warning(self,'Batch','Select input and output files first.'); return
        self.run.setEnabled(False); self.log.clear(); self.progress.setValue(0)
        buckling=None
        if self.batch_buck_enable.isChecked():
            buckling={'control_mode':'stress','macro_values':[x.value() for x in self.batch_buck_spins],'m_max':self.batch_buck_m.value(),'n_max':self.batch_buck_n.value(),'safety_factor':self.batch_buck_sf.value(),'use_abd_condensation':self.batch_buck_abd.isChecked(),'convergence_check':True}
        self.thread=BatchThread(self.input.text(),self.output.text(),self.workers.value(),self.limit.value() or None,buckling); self.thread.progress.connect(self.on_progress); self.thread.done.connect(self.on_done); self.thread.failed.connect(self.on_failed); self.thread.start()
    def on_progress(self,i,n,cid,status): self.progress.setMaximum(n); self.progress.setValue(i); self.log.appendPlainText(f'[{i}/{n}] {cid}  {status}')
    def on_done(self,df): self.run.setEnabled(True); self.log.appendPlainText('\nCompleted.\n'+df.Status.value_counts().to_string()); QMessageBox.information(self,'Batch completed',f'Results saved to\n{self.output.text()}')
    def on_failed(self,msg): self.run.setEnabled(True); self.log.appendPlainText(msg); QMessageBox.critical(self,'Batch failed',msg)


class ResultsPage(QWidget):
    """Batch-result post-processing with Ashby/Gibson-style and anisotropy views."""
    PRESETS={
        'Specific E1 vs density':('AverageDensityKgM3','SpecificE1_GPa_per_Mg_m3',True,True),
        'Specific E2 vs density':('AverageDensityKgM3','SpecificE2_GPa_per_Mg_m3',True,True),
        'Specific G12 vs density':('AverageDensityKgM3','SpecificG12_GPa_per_Mg_m3',True,True),
        'Gibson–Ashby: E1 vs relative density':('RelativeDensity','E1',True,True),
        'Gibson–Ashby: E2 vs relative density':('RelativeDensity','E2',True,True),
        'Gibson–Ashby: G12 vs relative density':('RelativeDensity','G12',True,True),
        'Geometry vs E anisotropy':('GeometryAspectRatio','E12_Anisotropy',False,False),
        'Geometry vs Q anisotropy':('GeometryAspectRatio','Q12_Anisotropy',False,False),
        'Relative density vs Q anisotropy':('RelativeDensity','Q12_Anisotropy',False,False),
        'Custom':(None,None,False,False),
    }
    def __init__(self,main):
        super().__init__(); self.main=main; self.df=None; self.raw_df=None
        v=QVBoxLayout(self); v.setContentsMargins(18,15,18,15)
        v.addWidget(title_block('Results Explorer','Post-process batch results with specific-stiffness maps, Gibson–Ashby-style scaling plots, geometry/stiffness anisotropy, directional stiffness and design ranking.'))
        top=QHBoxLayout(); self.path=QLineEdit(); self.path.setPlaceholderText('Batch results CSV / XLSX / Parquet')
        browse=QPushButton('Open results'); browse.clicked.connect(self.open); attach=QPushButton('Attach input database'); attach.clicked.connect(self.attach_input)
        export=QPushButton('Export derived data'); export.clicked.connect(self.export_derived)
        top.addWidget(self.path,1); top.addWidget(browse); top.addWidget(attach); top.addWidget(export); v.addLayout(top)
        self.info=QLabel('No result file loaded.'); self.info.setStyleSheet(f"color:{COLORS['muted']};"); v.addWidget(self.info)
        self.tabs=QTabWidget(); v.addWidget(self.tabs,1)
        self._build_scatter_tab(); self._build_anisotropy_tab(); self._build_direction_tab(); self._build_correlation_tab(); self._build_sensitivity_tab(); self._build_ranking_tab(); self._build_optimization_tab()

    def _build_scatter_tab(self):
        w=QWidget(); lay=QVBoxLayout(w); controls=QHBoxLayout()
        self.preset=QComboBox(); self.preset.addItems(list(self.PRESETS)); self.preset.currentTextChanged.connect(self.apply_preset)
        self.x=QComboBox(); self.y=QComboBox(); self.group=QComboBox(); self.group.addItem('None')
        self.logx=QCheckBox('log X'); self.logy=QCheckBox('log Y'); self.power=QCheckBox('Power-law fit'); self.pareto=QCheckBox('Pareto front')
        for c in (self.x,self.y,self.group): c.currentTextChanged.connect(self.plot_scatter)
        for c in (self.logx,self.logy,self.power,self.pareto): c.stateChanged.connect(self.plot_scatter)
        save=QPushButton('Save plot'); save.clicked.connect(lambda:self.save_figure(self.fig))
        controls.addWidget(QLabel('Preset')); controls.addWidget(self.preset); controls.addWidget(QLabel('X')); controls.addWidget(self.x); controls.addWidget(QLabel('Y')); controls.addWidget(self.y); controls.addWidget(QLabel('Group')); controls.addWidget(self.group)
        controls.addWidget(self.logx); controls.addWidget(self.logy); controls.addWidget(self.power); controls.addWidget(self.pareto); controls.addWidget(save)
        lay.addLayout(controls); self.scatter_note=QLabel(''); self.scatter_note.setStyleSheet(f"color:{COLORS['muted']};"); lay.addWidget(self.scatter_note)
        self.fig=Figure(); self.canvas=FigureCanvas(self.fig); lay.addWidget(self.canvas,1); self.tabs.addTab(w,'Ashby / scaling')

    def _build_anisotropy_tab(self):
        w=QWidget(); lay=QVBoxLayout(w); c=QHBoxLayout(); self.aniso_x=QComboBox(); self.aniso_y=QComboBox(); self.aniso_group=QComboBox(); self.aniso_group.addItem('None')
        for z in (self.aniso_x,self.aniso_y,self.aniso_group): z.currentTextChanged.connect(self.plot_anisotropy)
        save=QPushButton('Save plot'); save.clicked.connect(lambda:self.save_figure(self.aniso_fig))
        c.addWidget(QLabel('Geometry metric')); c.addWidget(self.aniso_x); c.addWidget(QLabel('Stiffness metric')); c.addWidget(self.aniso_y); c.addWidget(QLabel('Group')); c.addWidget(self.aniso_group); c.addStretch(); c.addWidget(save); lay.addLayout(c)
        self.aniso_fig=Figure(); self.aniso_canvas=FigureCanvas(self.aniso_fig); lay.addWidget(self.aniso_canvas,1); self.tabs.addTab(w,'Anisotropy')

    def _build_direction_tab(self):
        w=QWidget(); lay=QVBoxLayout(w); c=QHBoxLayout(); self.case=QComboBox(); self.case.currentTextChanged.connect(self.plot_directional); save=QPushButton('Save plot'); save.clicked.connect(lambda:self.save_figure(self.dir_fig)); c.addWidget(QLabel('Case')); c.addWidget(self.case); c.addStretch(); c.addWidget(save); lay.addLayout(c)
        self.dir_note=QLabel('In-plane E(θ) is reconstructed from E1, E2, G12 and ν12 using the orthotropic plane-stress compliance transformation.'); self.dir_note.setWordWrap(True); self.dir_note.setStyleSheet(f"color:{COLORS['muted']};"); lay.addWidget(self.dir_note)
        self.dir_fig=Figure(); self.dir_canvas=FigureCanvas(self.dir_fig); lay.addWidget(self.dir_canvas,1); self.tabs.addTab(w,'Directional E(θ)')

    def _build_correlation_tab(self):
        w=QWidget(); lay=QVBoxLayout(w); c=QHBoxLayout(); self.corr_set=QComboBox(); self.corr_set.addItems(['Core stiffness + density','Specific stiffness','Anisotropy metrics','All numeric (first 20)']); self.corr_set.currentTextChanged.connect(self.plot_correlation); save=QPushButton('Save plot'); save.clicked.connect(lambda:self.save_figure(self.corr_fig)); c.addWidget(QLabel('Variables')); c.addWidget(self.corr_set); c.addStretch(); c.addWidget(save); lay.addLayout(c)
        self.corr_fig=Figure(); self.corr_canvas=FigureCanvas(self.corr_fig); lay.addWidget(self.corr_canvas,1); self.tabs.addTab(w,'Correlation')

    def _build_sensitivity_tab(self):
        w=QWidget(); lay=QVBoxLayout(w)
        intro=QLabel('<b>Design sensitivity screening.</b> Numerical design variables are assessed using Pearson correlation, Spearman rank correlation and standardized multiple-regression coefficients. Geometry/material/layup categories use eta-squared. These are descriptive design-space measures, not causal derivatives.')
        intro.setWordWrap(True); intro.setStyleSheet(f"color:{COLORS['muted']};"); lay.addWidget(intro)
        c=QHBoxLayout(); self.sens_response=QComboBox(); self.sens_response.currentTextChanged.connect(self.run_sensitivity); self.sens_include_cat=QCheckBox('Include geometry/material/layup'); self.sens_include_cat.setChecked(True); self.sens_include_cat.stateChanged.connect(self.run_sensitivity); save=QPushButton('Save plot'); save.clicked.connect(lambda:self.save_figure(self.sens_fig)); export=QPushButton('Export table'); export.clicked.connect(self.export_sensitivity)
        c.addWidget(QLabel('Response')); c.addWidget(self.sens_response); c.addWidget(self.sens_include_cat); c.addStretch(); c.addWidget(export); c.addWidget(save); lay.addLayout(c)
        split=QSplitter(Qt.Orientation.Horizontal); self.sens_fig=Figure(); self.sens_canvas=FigureCanvas(self.sens_fig); split.addWidget(self.sens_canvas); self.sens_table=QTableWidget(0,0); self.sens_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents); split.addWidget(self.sens_table); split.setSizes([700,650]); lay.addWidget(split,1)
        self.sens_note=QLabel(''); self.sens_note.setWordWrap(True); self.sens_note.setStyleSheet(f"color:{COLORS['muted']};"); lay.addWidget(self.sens_note); self.sens_df=None; self.tabs.addTab(w,'Sensitivity')

    def run_sensitivity(self):
        if self.df is None or not hasattr(self,'sens_response') or not self.sens_response.currentText():return
        try:
            response=self.sens_response.currentText(); self.sens_df=sensitivity_analysis(self.df,response,include_categories=self.sens_include_cat.isChecked())
            self.sens_fig.clear(); ax=self.sens_fig.add_subplot(111)
            d=self.sens_df.head(18).iloc[::-1]
            if len(d):
                vals=pd.to_numeric(d.EffectStrength,errors='coerce'); ax.barh(range(len(d)),vals); ax.set_yticks(range(len(d))); ax.set_yticklabels(d.Variable.astype(str),fontsize=8); ax.set_xlabel('Sensitivity / effect strength (0–1 scale for rank or η²; |β| may exceed 1)'); ax.grid(True,axis='x',alpha=.22)
            else: ax.text(.5,.5,'No suitable design variables found. Attach the input database for older result files.',transform=ax.transAxes,ha='center')
            ax.set_title(f'Sensitivity of {response}'); self.sens_fig.tight_layout(); self.sens_canvas.draw_idle()
            cols=['Variable','VariableType','EffectStrength','SpearmanR','PearsonR','StandardizedBeta','EtaSquared','Samples']; self.sens_table.setRowCount(len(self.sens_df)); self.sens_table.setColumnCount(len(cols)); self.sens_table.setHorizontalHeaderLabels(cols)
            for i,(_,r) in enumerate(self.sens_df.iterrows()):
                for j,col in enumerate(cols):
                    v=r.get(col,''); txt=f'{v:.6g}' if isinstance(v,(float,np.floating)) and np.isfinite(v) else str(v); self.sens_table.setItem(i,j,QTableWidgetItem(txt))
            self.sens_note.setText('Large |Spearman ρ| indicates monotonic influence; standardized β is a joint linear-screening coefficient; categorical η² is the fraction of response variance associated with group means. Use these as screening indicators, especially when design variables are correlated.')
        except Exception as e:
            self.sens_note.setText(str(e))

    def export_sensitivity(self):
        if self.sens_df is None:return
        f,_=QFileDialog.getSaveFileName(self,'Export sensitivity','sensitivity.csv','CSV (*.csv);;Excel (*.xlsx)')
        if not f:return
        try:
            if Path(f).suffix.lower()=='.xlsx': self.sens_df.to_excel(f,index=False,sheet_name='Sensitivity')
            else:self.sens_df.to_csv(f,index=False)
        except Exception as e: QMessageBox.critical(self,'Cannot export sensitivity',str(e))

    def _build_ranking_tab(self):
        w=QWidget(); lay=QVBoxLayout(w); c=QHBoxLayout(); self.rank_metric=QComboBox(); self.rank_metric.currentTextChanged.connect(self.update_ranking); self.rank_n=QSpinBox(); self.rank_n.setRange(5,100); self.rank_n.setValue(20); self.rank_n.valueChanged.connect(self.update_ranking); c.addWidget(QLabel('Rank by')); c.addWidget(self.rank_metric); c.addWidget(QLabel('Top')); c.addWidget(self.rank_n); c.addStretch(); lay.addLayout(c)
        self.rank_table=QTableWidget(0,0); self.rank_table.setSelectionBehavior(QAbstractItemView.SelectRows); self.rank_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents); lay.addWidget(self.rank_table,1); self.tabs.addTab(w,'Top designs')

    def _build_optimization_tab(self):
        w=QWidget(); lay=QVBoxLayout(w)
        intro=QLabel('<b>Database-based multi-objective optimisation.</b> Select one or more Cᴴ/Qᴴ/engineering, mass, density or cost objectives. Optional robust k-sigma screening penalises assumed response uncertainty; manufacturing rules filter infeasible solved designs. The search remains over cases already analysed by RVE Studio.')
        intro.setWordWrap(True); intro.setStyleSheet(f"color:{COLORS['muted']};"); lay.addWidget(intro)
        bar=QHBoxLayout(); self.opt_preset=QComboBox(); self.opt_preset.addItems(['Custom','Maximise in-plane stiffness / minimise mass','Maximise stiffness / minimise cost','Balanced E1 + E2 / mass + cost'])
        self.opt_preset.currentTextChanged.connect(self.apply_optimization_preset); addo=QPushButton('Add objective'); addo.clicked.connect(self._add_objective_row); delo=QPushButton('Remove objective'); delo.clicked.connect(lambda:self._remove_current_row(self.obj_table)); addc=QPushButton('Add constraint'); addc.clicked.connect(self._add_constraint_row); delc=QPushButton('Remove constraint'); delc.clicked.connect(lambda:self._remove_current_row(self.con_table)); run=QPushButton('Run optimisation'); run.setObjectName('Primary'); run.clicked.connect(self.run_optimization)
        for x in [QLabel('Preset'),self.opt_preset,addo,delo,addc,delc,run]: bar.addWidget(x)
        bar.addStretch(); lay.addLayout(bar)

        settings=QGroupBox('Robustness and manufacturing screening'); sg=QGridLayout(settings)
        self.robust_enable=QCheckBox('Robust k-sigma optimisation'); self.robust_unc=QDoubleSpinBox(); self.robust_unc.setRange(0,100); self.robust_unc.setDecimals(2); self.robust_unc.setValue(5); self.robust_unc.setSuffix(' %'); self.robust_k=QDoubleSpinBox(); self.robust_k.setRange(0,5); self.robust_k.setDecimals(2); self.robust_k.setValue(1.0); self.robust_constraints=QCheckBox('Apply uncertainty to constraints'); self.robust_constraints.setChecked(True)
        sg.addWidget(self.robust_enable,0,0); sg.addWidget(QLabel('Assumed response uncertainty'),0,1); sg.addWidget(self.robust_unc,0,2); sg.addWidget(QLabel('k'),0,3); sg.addWidget(self.robust_k,0,4); sg.addWidget(self.robust_constraints,0,5)
        self.mfg_enable=QCheckBox('Manufacturing constraints'); self.mfg_min_t=QDoubleSpinBox(); self.mfg_max_t=QDoubleSpinBox(); self.mfg_ratio=QDoubleSpinBox(); self.mfg_aspect=QDoubleSpinBox(); self.mfg_rho=QDoubleSpinBox(); self.mfg_min_plies=QSpinBox(); self.mfg_max_plies=QSpinBox();
        for box in (self.mfg_min_t,self.mfg_max_t): box.setRange(0,1000); box.setDecimals(4); box.setSpecialValueText('Off'); box.setSuffix(' mm')
        self.mfg_ratio.setRange(0,1000); self.mfg_ratio.setDecimals(3); self.mfg_ratio.setSpecialValueText('Off'); self.mfg_aspect.setRange(0,1000); self.mfg_aspect.setDecimals(3); self.mfg_aspect.setSpecialValueText('Off'); self.mfg_rho.setRange(0,1); self.mfg_rho.setDecimals(4); self.mfg_rho.setSpecialValueText('Off'); self.mfg_min_plies.setRange(0,1000); self.mfg_min_plies.setSpecialValueText('Off'); self.mfg_max_plies.setRange(0,1000); self.mfg_max_plies.setSpecialValueText('Off'); self.mfg_balanced=QCheckBox('Balanced layup'); self.mfg_symmetric=QCheckBox('Symmetric layup'); self.mfg_angles=QLineEdit(); self.mfg_angles.setPlaceholderText('Allowed angles, e.g. 0,45,-45,90'); self.mfg_buckling=QDoubleSpinBox(); self.mfg_buckling.setRange(0,1000); self.mfg_buckling.setDecimals(3); self.mfg_buckling.setSpecialValueText('Off')
        sg.addWidget(self.mfg_enable,1,0); sg.addWidget(QLabel('min t'),1,1); sg.addWidget(self.mfg_min_t,1,2); sg.addWidget(QLabel('max t'),1,3); sg.addWidget(self.mfg_max_t,1,4); sg.addWidget(QLabel('max t ratio'),1,5); sg.addWidget(self.mfg_ratio,1,6)
        sg.addWidget(QLabel('max geometry aspect'),2,1); sg.addWidget(self.mfg_aspect,2,2); sg.addWidget(QLabel('max relative density'),2,3); sg.addWidget(self.mfg_rho,2,4); sg.addWidget(QLabel('plies min/max'),2,5); pp=QHBoxLayout(); pp.addWidget(self.mfg_min_plies); pp.addWidget(self.mfg_max_plies); sg.addLayout(pp,2,6); sg.addWidget(self.mfg_balanced,2,7); sg.addWidget(self.mfg_symmetric,2,8); sg.addWidget(QLabel('allowed ply angles'),3,1); sg.addWidget(self.mfg_angles,3,2,1,4); sg.addWidget(QLabel('min buckling factor λ'),3,6); sg.addWidget(self.mfg_buckling,3,7); lay.addWidget(settings)

        split=QSplitter(Qt.Orientation.Horizontal)
        og=QGroupBox('Objectives'); ol=QVBoxLayout(og); self.obj_table=QTableWidget(0,5); self.obj_table.setHorizontalHeaderLabels(['Use','Metric','Goal','Target','Weight']); self.obj_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents); self.obj_table.horizontalHeader().setStretchLastSection(True); ol.addWidget(self.obj_table); split.addWidget(og)
        cg=QGroupBox('Property constraints (optional)'); cl=QVBoxLayout(cg); self.con_table=QTableWidget(0,5); self.con_table.setHorizontalHeaderLabels(['Use','Metric','Relation','Value','Upper / tolerance']); self.con_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents); self.con_table.horizontalHeader().setStretchLastSection(True); cl.addWidget(self.con_table); split.addWidget(cg); split.setSizes([650,550]); lay.addWidget(split,1)
        self.opt_note=QLabel(''); self.opt_note.setWordWrap(True); self.opt_note.setStyleSheet(f"color:{COLORS['muted']};"); lay.addWidget(self.opt_note)
        lower=QSplitter(Qt.Orientation.Horizontal); self.opt_fig=Figure(); self.opt_canvas=FigureCanvas(self.opt_fig); lower.addWidget(self.opt_canvas)
        right=QWidget(); rl=QVBoxLayout(right); exp=QPushButton('Export optimisation results'); exp.clicked.connect(self.export_optimization); save=QPushButton('Save Pareto plot'); save.clicked.connect(lambda:self.save_figure(self.opt_fig)); rb=QHBoxLayout(); rb.addWidget(exp); rb.addWidget(save); rb.addStretch(); rl.addLayout(rb); self.opt_table=QTableWidget(0,0); self.opt_table.setSelectionBehavior(QAbstractItemView.SelectRows); self.opt_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents); rl.addWidget(self.opt_table,1); lower.addWidget(right); lower.setSizes([600,750]); lay.addWidget(lower,2)
        self.opt_df=None; self._robust_active=False; self._add_objective_row(metric='E1',goal='Maximise',weight=1.0); self._add_objective_row(metric='MassKg',goal='Minimise',weight=1.0); self.tabs.addTab(w,'Multi-objective optimisation')

    def _metric_items(self):
        if self.df is None:return []
        preferred=[]
        for c in self.df.columns:
            if pd.api.types.is_numeric_dtype(self.df[c]): preferred.append(c)
        return preferred

    def _add_objective_row(self,checked=False,metric='',goal='Maximise',target='',weight=1.0):
        # checked argument is ignored when Qt sends the clicked(bool) signal.
        r=self.obj_table.rowCount(); self.obj_table.insertRow(r); set_combo_cell(self.obj_table,r,0,['True','False'],'True')
        set_combo_cell(self.obj_table,r,1,self._metric_items(),metric,editable=False); set_combo_cell(self.obj_table,r,2,['Maximise','Minimise','Target'],'Maximise' if isinstance(goal,bool) else goal)
        self.obj_table.setItem(r,3,table_item(target)); self.obj_table.setItem(r,4,table_item(weight))

    def _add_constraint_row(self,checked=False,metric='',relation='>=',value='',upper=''):
        r=self.con_table.rowCount(); self.con_table.insertRow(r); set_combo_cell(self.con_table,r,0,['True','False'],'True'); set_combo_cell(self.con_table,r,1,self._metric_items(),metric,editable=False); set_combo_cell(self.con_table,r,2,['>=','<=','Between'],'>=' if isinstance(relation,bool) else relation); self.con_table.setItem(r,3,table_item(value)); self.con_table.setItem(r,4,table_item(upper))

    @staticmethod
    def _remove_current_row(table):
        r=table.currentRow()
        if r>=0: table.removeRow(r)

    def _refresh_optimizer_metrics(self):
        items=self._metric_items()
        for table,col in ((self.obj_table,1),(self.con_table,1)):
            for r in range(table.rowCount()):
                old=cell_text(table,r,col); set_combo_cell(table,r,col,items,old,editable=False)

    def apply_optimization_preset(self):
        if not hasattr(self,'obj_table') or self.df is None:return
        name=self.opt_preset.currentText(); self.obj_table.setRowCount(0); self.con_table.setRowCount(0)
        available=set(self._metric_items())
        def add(metric,goal,weight=1.0,target=''):
            if metric in available:self._add_objective_row(metric=metric,goal=goal,target=target,weight=weight)
        if name=='Maximise in-plane stiffness / minimise mass': add('E1','Maximise'); add('E2','Maximise'); add('MassKg','Minimise')
        elif name=='Maximise stiffness / minimise cost': add('E1','Maximise'); add('E2','Maximise'); add('TotalMaterialCost','Minimise')
        elif name=='Balanced E1 + E2 / mass + cost': add('E1','Maximise'); add('E2','Maximise'); add('MassKg','Minimise'); add('TotalMaterialCost','Minimise')
        else: add('E1','Maximise'); add('MassKg','Minimise')

    def _objective_specs(self):
        specs=[]
        for r in range(self.obj_table.rowCount()):
            if cell_text(self.obj_table,r,0).lower()!='true':continue
            metric=cell_text(self.obj_table,r,1); goal=cell_text(self.obj_table,r,2); target_txt=cell_text(self.obj_table,r,3); weight_txt=cell_text(self.obj_table,r,4)
            if not metric:continue
            target=float(target_txt) if target_txt else None; weight=float(weight_txt or 1.0)
            specs.append(ObjectiveSpec(metric,{'Maximise':'max','Minimise':'min','Target':'target'}[goal],target,weight))
        return specs

    def _constraint_specs(self):
        specs=[]
        for r in range(self.con_table.rowCount()):
            if cell_text(self.con_table,r,0).lower()!='true':continue
            metric=cell_text(self.con_table,r,1); relation=cell_text(self.con_table,r,2); value=cell_text(self.con_table,r,3); upper=cell_text(self.con_table,r,4)
            if metric and value: specs.append(ConstraintSpec(metric,relation,float(value),float(upper) if upper else None))
        return specs

    def _manufacturing_rules(self):
        if not self.mfg_enable.isChecked(): return None
        val=lambda w: float(w.value()) if float(w.value())>0 else None
        ival=lambda w: int(w.value()) if int(w.value())>0 else None
        angles=None
        if self.mfg_angles.text().strip():
            angles=tuple(float(x) for x in re.split(r'[/,;\s]+',self.mfg_angles.text().strip()) if x)
        buck_req=val(self.mfg_buckling)
        if buck_req is not None and (self.df is None or 'BucklingLambdaMin' not in self.df.columns):
            raise ValueError('A minimum buckling factor was requested, but the result table has no BucklingLambdaMin column. Re-run Large Batch Analysis with local-wall buckling screening enabled.')
        return ManufacturingRules(
            min_wall_thickness_mm=val(self.mfg_min_t), max_wall_thickness_mm=val(self.mfg_max_t),
            max_wall_thickness_ratio=val(self.mfg_ratio), max_geometry_aspect_ratio=val(self.mfg_aspect),
            max_relative_density=val(self.mfg_rho), min_plies=ival(self.mfg_min_plies), max_plies=ival(self.mfg_max_plies),
            require_balanced_layup=self.mfg_balanced.isChecked(), require_symmetric_layup=self.mfg_symmetric.isChecked(), allowed_ply_angles=angles, min_buckling_factor=buck_req)

    def run_optimization(self):
        if self.df is None: QMessageBox.warning(self,'Optimisation','Open batch results first.'); return
        try:
            objectives=self._objective_specs(); constraints=self._constraint_specs(); rules=self._manufacturing_rules()
            if not objectives: raise ValueError('Add at least one active objective.')
            self._robust_active=self.robust_enable.isChecked()
            if self._robust_active:
                self.opt_df=robust_optimize_batch_designs(self.df,objectives,constraints,self.robust_unc.value(),self.robust_k.value(),rules,self.robust_constraints.isChecked())
                fcol,pcol,scol='RobustOptimizationFeasible','RobustOptimizationPareto','RobustOptimizationScore'
                mode=f'robust k={self.robust_k.value():.3g}, assumed uncertainty={self.robust_unc.value():.3g}%'
            else:
                self.opt_df=optimize_batch_designs(self.df,objectives,constraints,manufacturing_rules=rules)
                fcol,pcol,scol='OptimizationFeasible','OptimizationPareto','OptimizationScore'; mode='nominal'
            feasible=int(self.opt_df[fcol].sum()); pareto=int(self.opt_df[pcol].sum())
            mfg=' with manufacturing screening' if rules is not None else ''
            self.opt_note.setText(f'{feasible:,} feasible cases | {pareto:,} Pareto-optimal cases | {mode}{mfg}. Robust mode uses a conservative k-sigma penalty based on the user-assumed relative response uncertainty; it is not a stochastic re-analysis of the solver.')
            self._plot_optimization(objectives); self._show_optimization_table(objectives)
        except Exception as e: QMessageBox.critical(self,'Optimisation failed',str(e))

    def _opt_columns(self):
        return ('RobustOptimizationFeasible','RobustOptimizationPareto','RobustOptimizationScore') if self._robust_active else ('OptimizationFeasible','OptimizationPareto','OptimizationScore')

    def _opt_metric_col(self,metric):
        c=f'RobustValue_{metric}'
        return c if self._robust_active and self.opt_df is not None and c in self.opt_df.columns else metric

    def _plot_optimization(self,objectives):
        self.opt_fig.clear(); ax=self.opt_fig.add_subplot(111); d=self.opt_df; fcol,pcol,scol=self._opt_columns()
        feasible=d[d[fcol]]; pareto=d[d[pcol]]
        if len(objectives)>=2:
            x0,y0=objectives[0].metric,objectives[1].metric; x,y=self._opt_metric_col(x0),self._opt_metric_col(y0); ax.scatter(pd.to_numeric(feasible[x],errors='coerce'),pd.to_numeric(feasible[y],errors='coerce'),s=18,alpha=.35,label='Feasible'); ax.scatter(pd.to_numeric(pareto[x],errors='coerce'),pd.to_numeric(pareto[y],errors='coerce'),s=38,label='Pareto'); ax.set_xlabel(('Robust ' if self._robust_active else '')+self._label(x0)); ax.set_ylabel(('Robust ' if self._robust_active else '')+self._label(y0))
        else:
            x0=objectives[0].metric; x=self._opt_metric_col(x0); xx=pd.to_numeric(feasible[x],errors='coerce'); yy=pd.to_numeric(feasible[scol],errors='coerce'); ax.scatter(xx,yy,s=20,alpha=.5,label='Feasible'); ax.scatter(pd.to_numeric(pareto[x],errors='coerce'),pd.to_numeric(pareto[scol],errors='coerce'),s=42,label='Pareto'); ax.set_xlabel(('Robust ' if self._robust_active else '')+self._label(x0)); ax.set_ylabel('Normalized trade-off score')
        if len(feasible):
            scores=pd.to_numeric(feasible[scol],errors='coerce'); best=feasible[scores==scores.max()]
            if len(best) and len(objectives)>=2:
                b=best.iloc[0]; ax.scatter([float(b[self._opt_metric_col(objectives[0].metric)])],[float(b[self._opt_metric_col(objectives[1].metric)])],marker='*',s=150,label='Highest weighted score')
        ax.grid(True,alpha=.22);
        if len(feasible): ax.legend(fontsize=8,loc='best')
        self.opt_fig.tight_layout(); self.opt_canvas.draw_idle()

    def _show_optimization_table(self,objectives):
        fcol,pcol,scol=self._opt_columns(); d=self.opt_df[self.opt_df[fcol]].copy(); d['_pareto_sort']=d[pcol].astype(int); d=d.sort_values(['_pareto_sort',scol],ascending=[False,False]).head(200)
        cols=[c for c in ['CaseID',pcol,scol,'GeometryID','Material','LayupID','MassKg','TotalMaterialCost','AverageDensityKgM3','MinWallThickness_mm','MaxWallThickness_mm','WallThicknessRatio','GeometryAspectRatio','BalancedLayup','SymmetricLayup'] if c in d.columns]
        for o in objectives:
            for c in [o.metric,self._opt_metric_col(o.metric)]:
                if c in d.columns and c not in cols: cols.append(c)
        self.opt_table.setRowCount(len(d)); self.opt_table.setColumnCount(len(cols)); self.opt_table.setHorizontalHeaderLabels(cols)
        for i,(_,row) in enumerate(d.iterrows()):
            for j,c in enumerate(cols):
                v=row[c]; txt=f'{v:.7g}' if isinstance(v,(float,np.floating)) and np.isfinite(v) else str(v); self.opt_table.setItem(i,j,QTableWidgetItem(txt))

    def export_optimization(self):
        if self.opt_df is None:return
        f,_=QFileDialog.getSaveFileName(self,'Export optimisation','multiobjective_optimisation.csv','CSV (*.csv);;Excel (*.xlsx)')
        if not f:return
        try:
            if Path(f).suffix.lower()=='.xlsx': self.opt_df.to_excel(f,index=False,sheet_name='Optimisation')
            else:self.opt_df.to_csv(f,index=False)
        except Exception as e: QMessageBox.critical(self,'Cannot export',str(e))

    def open(self):
        f,_=QFileDialog.getOpenFileName(self,'Open results','','Data (*.csv *.xlsx *.parquet)')
        if not f:return
        try:
            self.path.setText(f); self.raw_df=_read_table(f); self._set_dataframe(self.raw_df)
        except Exception as e: QMessageBox.critical(self,'Cannot open results',str(e))

    def attach_input(self):
        if self.raw_df is None: QMessageBox.warning(self,'Results Explorer','Open a result file first.'); return
        f,_=QFileDialog.getOpenFileName(self,'Attach original input database','','Data (*.csv *.xlsx *.parquet)')
        if not f:return
        try:
            merged=merge_input_metadata(self.raw_df,_read_table(f)); self.raw_df=merged; self._set_dataframe(merged); self.info.setText(self.info.text()+f'  |  metadata attached from {Path(f).name}')
        except Exception as e: QMessageBox.critical(self,'Cannot attach input',str(e))

    def _set_dataframe(self,df):
        d=df.copy()
        if 'Status' in d.columns: d=d[d.Status=='OK'].copy()
        self.df=enrich_results(d)
        numeric=[c for c in self.df.columns if pd.api.types.is_numeric_dtype(self.df[c])]
        categories=[]
        for c in ['GeometryID','GeometryType','Material','LayupID','LayupAngles']:
            if c in self.df.columns: categories.append(c)
        for combo in (self.x,self.y,self.aniso_x,self.aniso_y,self.rank_metric,self.sens_response): combo.clear()
        self.x.addItems(numeric); self.y.addItems(numeric)
        geom=[c for c in ['GeometryAspectRatio','GeometryLengthRatio','LengthRatio','LengthA_mm','LengthB_mm','Param1Value','Param2Value','RelativeDensity'] if c in numeric]
        anis=[c for c in ['E12_Anisotropy','Q12_Anisotropy','OutOfPlaneShear_Anisotropy','E1_over_E2','Q11_over_Q22','G13_over_G23','Q12_CouplingIndex'] if c in numeric]
        self.aniso_x.addItems(geom or numeric); self.aniso_y.addItems(anis or numeric)
        rank=[c for c in numeric if c.startswith('Specific') or c in ['MinInPlaneSpecificE_GPa_per_Mg_m3','E1','E2','G12','Q11','Q22']]; self.rank_metric.addItems(rank or numeric)
        responses=[c for c in ['E1','E2','E3','G12','G13','G23','Q11','Q22','Q12','Q66','Q44','Q55','MassKg','RelativeDensity','TotalMaterialCost'] if c in numeric]; self.sens_response.addItems(responses or numeric)
        for combo in (self.group,self.aniso_group): combo.blockSignals(True); combo.clear(); combo.addItem('None'); combo.addItems(categories); combo.blockSignals(False)
        self.case.blockSignals(True); self.case.clear(); self.case.addItems(self.df['CaseID'].astype(str).tolist() if 'CaseID' in self.df else [str(i) for i in self.df.index]); self.case.blockSignals(False)
        self._refresh_optimizer_metrics(); self.apply_optimization_preset()
        self.info.setText(f'{len(self.df):,} successful cases | {len(numeric):,} numeric/derived metrics | specific modulus unit: GPa/(Mg/m³) | multi-objective Pareto design selection available')
        self.preset.setCurrentText('Specific E1 vs density'); self.apply_preset(); self.plot_anisotropy(); self.plot_directional(); self.plot_correlation(); self.run_sensitivity(); self.update_ranking()

    def apply_preset(self):
        if self.df is None:return
        x,y,lx,ly=self.PRESETS.get(self.preset.currentText(),(None,None,False,False))
        if x in [self.x.itemText(i) for i in range(self.x.count())]: self.x.setCurrentText(x)
        if y in [self.y.itemText(i) for i in range(self.y.count())]: self.y.setCurrentText(y)
        self.logx.setChecked(lx); self.logy.setChecked(ly)
        self.power.setChecked(self.preset.currentText().startswith('Gibson'))
        self.plot_scatter()

    def _scatter_groups(self,d,group):
        if group and group!='None' and group in d.columns:
            return [(str(k),g) for k,g in d.groupby(group,dropna=False)]
        return [('',d)]

    def plot_scatter(self):
        if self.df is None or not self.x.currentText() or not self.y.currentText():return
        xname,yname=self.x.currentText(),self.y.currentText(); d=self.df
        self.fig.clear(); ax=self.fig.add_subplot(111)
        for label,g in self._scatter_groups(d,self.group.currentText()):
            xx=pd.to_numeric(g[xname],errors='coerce'); yy=pd.to_numeric(g[yname],errors='coerce'); mask=np.isfinite(xx)&np.isfinite(yy)
            if self.logx.isChecked(): mask &= xx>0
            if self.logy.isChecked(): mask &= yy>0
            ax.scatter(xx[mask],yy[mask],s=20,alpha=.68,label=label or None)
        note=[]
        if self.power.isChecked():
            fit=fit_power_law(d[xname],d[yname])
            if fit:
                xv=pd.to_numeric(d[xname],errors='coerce').to_numpy(float); pos=xv[np.isfinite(xv)&(xv>0)]
                if len(pos):
                    xs=np.geomspace(pos.min(),pos.max(),120); ax.plot(xs,fit.coefficient*xs**fit.exponent,linewidth=1.8,label=f'fit: y={fit.coefficient:.3g} x^{fit.exponent:.3f}')
                    note.append(f'Power-law exponent n={fit.exponent:.3f}, R²={fit.r_squared:.4f}, N={fit.n_points}')
        if self.pareto.isChecked():
            xx=pd.to_numeric(d[xname],errors='coerce').to_numpy(float); yy=pd.to_numeric(d[yname],errors='coerce').to_numpy(float); pm=pareto_mask_min_x_max_y(xx,yy)
            if pm.any():
                order=np.argsort(xx[pm]); ax.plot(xx[pm][order],yy[pm][order],marker='o',linewidth=1.4,label='Pareto: min X / max Y')
        if self.logx.isChecked(): ax.set_xscale('log')
        if self.logy.isChecked(): ax.set_yscale('log')
        ax.set_xlabel(self._label(xname)); ax.set_ylabel(self._label(yname)); ax.grid(True,which='both',alpha=.22)
        handles,labels=ax.get_legend_handles_labels()
        if labels: ax.legend(fontsize=8,loc='best')
        self.scatter_note.setText(' | '.join(note) if note else 'Ashby-style plots use effective stiffness and density already calculated by RVE Studio; no strength values are inferred.')
        self.fig.tight_layout(); self.canvas.draw_idle()

    def plot_anisotropy(self):
        if self.df is None or not self.aniso_x.currentText() or not self.aniso_y.currentText():return
        x,y=self.aniso_x.currentText(),self.aniso_y.currentText(); self.aniso_fig.clear(); ax=self.aniso_fig.add_subplot(111)
        for label,g in self._scatter_groups(self.df,self.aniso_group.currentText()):
            xx=pd.to_numeric(g[x],errors='coerce'); yy=pd.to_numeric(g[y],errors='coerce'); m=np.isfinite(xx)&np.isfinite(yy); ax.scatter(xx[m],yy[m],s=22,alpha=.7,label=label or None)
        ax.set_xlabel(self._label(x)); ax.set_ylabel(self._label(y)); ax.grid(True,alpha=.22); handles,labels=ax.get_legend_handles_labels();
        if labels: ax.legend(fontsize=8,loc='best')
        self.aniso_fig.tight_layout(); self.aniso_canvas.draw_idle()

    def plot_directional(self):
        if self.df is None or self.case.count()==0:return
        cid=self.case.currentText(); row=self.df[self.df.CaseID.astype(str)==cid].iloc[0] if 'CaseID' in self.df else self.df.iloc[int(cid)]
        self.dir_fig.clear(); ax=self.dir_fig.add_subplot(111,projection='polar')
        try:
            ang,E=directional_youngs_modulus(row.E1,row.E2,row.G12,row.nu12); ax.plot(np.deg2rad(ang),E,linewidth=2); ax.set_title(f'{cid} — in-plane Young\'s modulus E(θ)'); ax.grid(True,alpha=.3)
        except Exception as e: ax.text(.5,.5,str(e),transform=ax.transAxes,ha='center')
        self.dir_fig.tight_layout(); self.dir_canvas.draw_idle()

    def plot_correlation(self):
        if self.df is None:return
        mode=self.corr_set.currentText();
        if mode=='Core stiffness + density': cols=['RelativeDensity','AverageDensityKgM3','E1','E2','E3','G12','G13','G23','Q11','Q22','Q66']
        elif mode=='Specific stiffness': cols=[c for c in self.df.columns if c.startswith('Specific')]
        elif mode=='Anisotropy metrics': cols=['GeometryAspectRatio','E12_Anisotropy','Q12_Anisotropy','OutOfPlaneShear_Anisotropy','Q12_CouplingIndex']
        else: cols=[c for c in self.df.columns if pd.api.types.is_numeric_dtype(self.df[c])][:20]
        cols=[c for c in cols if c in self.df and pd.api.types.is_numeric_dtype(self.df[c])]
        self.corr_fig.clear(); ax=self.corr_fig.add_subplot(111)
        if len(cols)<2: ax.text(.5,.5,'Not enough variables',transform=ax.transAxes,ha='center')
        else:
            corr=self.df[cols].corr().to_numpy(); im=ax.imshow(corr,vmin=-1,vmax=1,aspect='auto'); ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols,rotation=60,ha='right',fontsize=7); ax.set_yticks(range(len(cols))); ax.set_yticklabels(cols,fontsize=7); self.corr_fig.colorbar(im,ax=ax,fraction=.025,pad=.02,label='Pearson r')
        self.corr_fig.tight_layout(); self.corr_canvas.draw_idle()

    def update_ranking(self):
        if self.df is None or not self.rank_metric.currentText():return
        metric=self.rank_metric.currentText(); d=self.df.copy(); d['_rank']=pd.to_numeric(d[metric],errors='coerce'); d=d[np.isfinite(d['_rank'])].sort_values('_rank',ascending=False).head(self.rank_n.value())
        cols=[c for c in ['CaseID','GeometryID','Material','LayupID','RelativeDensity','AverageDensityKgM3',metric,'GeometryAspectRatio','E12_Anisotropy'] if c in d.columns]
        self.rank_table.setRowCount(len(d)); self.rank_table.setColumnCount(len(cols)); self.rank_table.setHorizontalHeaderLabels(cols)
        for i,(_,r) in enumerate(d.iterrows()):
            for j,c in enumerate(cols):
                val=r[c]; txt=f'{val:.6g}' if isinstance(val,(float,np.floating)) and np.isfinite(val) else str(val); self.rank_table.setItem(i,j,QTableWidgetItem(txt))

    def export_derived(self):
        if self.df is None:return
        f,_=QFileDialog.getSaveFileName(self,'Export derived results','derived_batch_results.csv','CSV (*.csv);;Excel (*.xlsx)')
        if not f:return
        try:
            if Path(f).suffix.lower()=='.xlsx': self.df.to_excel(f,index=False,sheet_name='DerivedResults')
            else: self.df.to_csv(f,index=False)
        except Exception as e: QMessageBox.critical(self,'Cannot export',str(e))

    def save_figure(self,fig):
        f,_=QFileDialog.getSaveFileName(self,'Save figure','batch_postprocessing.png','PNG (*.png);;SVG (*.svg);;PDF (*.pdf)')
        if f: fig.savefig(f,dpi=300,bbox_inches='tight')

    @staticmethod
    def _label(name):
        labels={'RelativeDensity':'Relative density ρ*/ρs','AverageDensityKgM3':'Average density (kg/m³)','GeometryAspectRatio':'Geometry aspect ratio (≥1)','GeometryLengthRatio':'Geometry length ratio A/B','E12_Anisotropy':'In-plane E anisotropy (max/min)','Q12_Anisotropy':'In-plane Q anisotropy (max/min)','OutOfPlaneShear_Anisotropy':'Out-of-plane shear anisotropy (max/min)','Q12_CouplingIndex':'|Q12| / √(|Q11 Q22|)','MinInPlaneSpecificE_GPa_per_Mg_m3':'Minimum in-plane specific E [GPa/(Mg/m³)]','MassKg':'RVE mass (kg)','WeightN':'RVE weight (N)','TotalMaterialCost':'RVE material cost (user cost units)','EffectiveCostPerKg':'Effective material cost / kg','CostPerRVEVolume':'Material cost / RVE volume'}
        if name.startswith('Specific') and name.endswith('_GPa_per_Mg_m3'):
            prop=name[len('Specific'):-len('_GPa_per_Mg_m3')]; return f'Specific {prop} [GPa/(Mg/m³)]'
        return labels.get(name,name)


class ValidationPage(QWidget):
    """Benchmark designs: the 14 wall configurations of the source study solved with any wall model."""
    def __init__(self, main):
        super().__init__()
        self.last_benchmark = None
        v = QVBoxLayout(self)
        v.setContentsMargins(18, 15, 18, 15)
        v.addWidget(title_block('Benchmark Designs',
            '14 benchmark wall configurations of the hexagonal core (Sriharan et al. 2024, Engineering Structures 299, '
            '117106) × 6 stiffnesses. Solve them with any wall beam and bending option and compare with the '
            'reference values of the cylindrical Euler–Bernoulli model.'))
        top = QHBoxLayout()
        self.wall_model = WallModelOptions('Wall model for the benchmark')
        top.addWidget(self.wall_model, 3)
        box = QGroupBox('Run'); g = QGridLayout(box)
        g.addWidget(QLabel('Tolerance (%)'), 0, 0)
        self.tolerance = QDoubleSpinBox()
        self.tolerance.setRange(0.01, 100.0)
        self.tolerance.setDecimals(2)
        self.tolerance.setValue(1.0)
        self.tolerance.setToolTip('Maximum allowed absolute relative difference from the reference values; default ±1%.')
        g.addWidget(self.tolerance, 0, 1)
        run = QPushButton('Solve the 14 benchmark designs')
        run.setObjectName('Primary')
        run.clicked.connect(self.run)
        g.addWidget(run, 1, 0, 1, 2)
        self.export_btn = QPushButton('Export CSV…')
        self.export_btn.setEnabled(False)
        self.export_btn.clicked.connect(self.export_csv)
        ref = QPushButton('Read source paper')
        ref.setToolTip(PAPER_URL)
        ref.clicked.connect(lambda: QDesktopServices.openUrl(QUrl(PAPER_URL)))
        g.addWidget(self.export_btn, 2, 0); g.addWidget(ref, 2, 1)
        top.addWidget(box, 2)
        v.addLayout(top)
        controls2 = QHBoxLayout()
        self.filter_input = QLineEdit()
        self.filter_input.setPlaceholderText('Filter by configuration or stiffness (e.g. Q66 or 7)')
        self.filter_input.textChanged.connect(self._filter_rows)
        controls2.addWidget(self.filter_input, 1)
        self.only_outliers = QCheckBox('Show outside tolerance only')
        self.only_outliers.toggled.connect(self._filter_rows)
        controls2.addWidget(self.only_outliers)
        v.addLayout(controls2)
        self.summary = QLabel('Choose a wall model and solve the benchmark designs.')
        self.summary.setWordWrap(True)
        self.summary.setObjectName('Banner')
        v.addWidget(self.summary)
        self.table = QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(['Design', 'Stiffness', 'Calculated (MPa)', 'Reference (MPa)',
                                              'Difference %', 'Status'])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.verticalHeader().setVisible(False)
        v.addWidget(self.table, 1)
        foot = QLabel('Reference values: cylindrical Euler–Bernoulli model, DOI 10.1016/j.engstruct.2023.117106, '
            'Eq. (5), Tables 3–9 (Q44 = 1–3 shear, Q55 = 2–3 shear). With the default wall model RVE Studio reproduces them '
            'within 0.3 %; other wall models differ by design. Check a design against shell FE under '
            '"FE Design & Validation".')
        foot.setWordWrap(True)
        foot.setStyleSheet(f'color:{COLORS["muted"]};')
        v.addWidget(foot)

    def _filter_rows(self, *_):
        text = self.filter_input.text().strip().lower()
        for i in range(self.table.rowCount()):
            cfg = self.table.item(i, 0).text().lower()
            prop = self.table.item(i, 1).text().lower()
            status = self.table.item(i, 5).text()
            matches = not text or text in cfg or text in prop
            self.table.setRowHidden(i, not (matches and (not self.only_outliers.isChecked() or status == 'DIFFERS')))

    def run(self):
        try:
            QApplication.setOverrideCursor(Qt.WaitCursor)
            out = run_all14_benchmark(tolerance_pct=self.tolerance.value(), wall_shear=self.wall_model.wall_shear(),
                                      bending=self.wall_model.bending_option(),
                                      twist_release=self.wall_model.twist_release())
            self.last_benchmark = out
            d = out['detail']
            self.table.setRowCount(len(d))
            red = QColor(COLORS['orange']); green = QColor(COLORS['green'])
            for i, row in enumerate(d.itertuples(index=False)):
                status = 'WITHIN' if row.Within_tolerance else 'DIFFERS'
                vals = [int(row.Config), row.Property, f'{row.Calculated_MPa:.6f}',
                        f'{row.Paper_MPa:.6f}', f'{row.Error_pct:+.3f}', status]
                for j, value in enumerate(vals):
                    item = table_item(value)
                    if j == 5:
                        item.setForeground(green if status == 'WITHIN' else red)
                    self.table.setItem(i, j, item)
            self.export_btn.setEnabled(True)
            self._filter_rows()
            default = out['wall_shear'] == 'euler-bernoulli' and out['bending'] == 'cylindrical'
            msg = (f"{self.wall_model.describe()}: {out['within_tolerance']}/{out['total']} values within "
                   f"±{out['tolerance_pct']:g}% of the reference; largest difference {out['max_abs_paper_pct']:.3f}%. ")
            msg += ('Reference model reproduced.' if out['passed'] else
                    ('Some entries differ — review them.' if default else
                     'Differences show the effect of the chosen wall model on each stiffness.'))
            self.summary.setText(msg)
            self.summary.setProperty('tone', 'good' if out['passed'] else ('info' if not default else 'warn'))
            self.summary.style().unpolish(self.summary); self.summary.style().polish(self.summary)
        except Exception as exc:
            QMessageBox.critical(self, 'Benchmark failed', str(exc))
        finally:
            QApplication.restoreOverrideCursor()

    def export_csv(self):
        if self.last_benchmark is None: return
        path, _ = QFileDialog.getSaveFileName(self, 'Export benchmark results',
            'RVE_Studio_benchmark_designs.csv', 'CSV files (*.csv)')
        if not path: return
        try:
            self.last_benchmark['detail'].to_csv(path, index=False)
        except Exception as exc:
            QMessageBox.critical(self, 'CSV export failed', str(exc))



class HelpPage(QWidget):
    """Help: opens the user guide and the theory guide."""
    def __init__(self,main):
        super().__init__(); self.main=main
        v=QVBoxLayout(self); v.setContentsMargins(18,15,18,15)
        v.addWidget(title_block('Help','The user guide explains every workspace step by step; the theory guide gives the equations, the wall models, their verified limits and how to judge mesh convergence.'))
        row=QHBoxLayout(); row.setSpacing(16); self.guide_buttons={}
        for kind,img,colour,text in (('user','rve3d.png','#2F80ED','Installation, a first design, every workspace (RVE designer, libraries, databases, results explorer, benchmark designs, FE design & validation), command line, Python API, troubleshooting.'),
                                     ('theory','laminate.png','#1F9D72','Laminate and wall equations, the four bending options and two beam theories, periodic homogenisation, reliability limits, choosing a wall model, shell FE validation and mesh convergence.')):
            title=GUIDES[kind][2]
            tile=picture_tile(img,title,text,f'Open {title}',colour,lambda _=False,k=kind:self.open_guide(k))
            self.guide_buttons[kind]=tile.findChild(QPushButton)
            row.addWidget(tile,1)
        v.addLayout(row,1)

    def open_guide(self,kind):
        p=guide_path(kind)
        if p is None:
            QMessageBox.warning(self,'Guide not found',f'The {GUIDES[kind][2]} was not found in this installation.'); return False
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(p))); return True


class AboutPage(QWidget):
    def __init__(self,main):
        super().__init__()
        v=QVBoxLayout(self)
        v.setContentsMargins(24,20,24,20)
        v.addWidget(title_block('About / Workflow','RVE Studio · Version 1.0.0'))
        cites=''.join(f'<li>{html.escape(t)}' + (f' <a href="{u}">{u}</a>' if u else '') + '</li>' for t,u in CITATIONS)
        txt=QLabel(f'''<b>How to cite</b><br>
        RVE Studio builds on the papers below. Please cite them together, and the software release through its Zenodo DOI (CITATION.cff):
        <ol>{cites}</ol>
        <b>Workflow</b><br>
        Materials and ply layups → CLT/ABD and wall shear stiffness → energy-based wall elements → periodic RVE assembly → condensed C<sup>H</sup> → effective properties, reliability rating and buckling screening.<br><br>
        <b>Conventions</b><br>
        Internal 3D tensor uses engineering Voigt [11,22,33,23,13,12].
        The shell transverse matrix follows Eq. (5) of the 2024 paper: Q44 maps to 13 and Q55 maps to 23.<br><br>
        <b>Benchmark designs and validation</b><br>
        Benchmark Designs solves the 14 wall configurations of the 2024 paper with any wall model and compares them
        with its reference values (Euler–Bernoulli wall, cylindrical bending). FE Design &amp; Validation checks any design
        against periodic shell FE (built-in solver or Abaqus) with a mesh-convergence verdict.<br><br>
        <b>Licence</b><br>
        Open source, MIT licence (see LICENSE).''')
        txt.setOpenExternalLinks(True)
        txt.setWordWrap(True)
        txt.setStyleSheet('font-size:11pt;')
        v.addWidget(txt)
        v.addStretch()


class MainWindow(QMainWindow):
    def __init__(self):
        ensure_library_dirs()
        super().__init__(); self.setWindowTitle(f'RVE Studio {__import__("rvestudio").__version__} — Composite cellular core design'); self.resize(1500,900); self.setMinimumSize(1100,700)
        root=QWidget(); self.setCentralWidget(root); h=QHBoxLayout(root); h.setContentsMargins(0,0,0,0); h.setSpacing(0); self.nav=QListWidget(); self.nav.setObjectName('Nav'); self.nav.setFixedWidth(230); h.addWidget(self.nav); self.stack=QStackedWidget(); h.addWidget(self.stack,1)
        names=['Dashboard','Single RVE Designer','Laminate Tools','Material Library','Geometry Library','Generate DB','Large Batch Analysis','Results Explorer','Benchmark Designs','FE Design & Validation','Help','About / Workflow']; self.page_names=names
        self.single_page=SingleRVEPage(self); self.fe_page=FEPage(self); self.help_page=HelpPage(self); pages=[DashboardPage(self),self.single_page,LaminateToolsPage(self),MaterialLibraryPage(self),GeometryLibraryPage(self),GenerateDBPage(self),BatchPage(self),ResultsPage(self),ValidationPage(self),self.fe_page,self.help_page,AboutPage(self)]
        from PySide6.QtCore import QSize
        self.nav.setIconSize(QSize(24,24))
        for n,p in zip(names,pages):
            self.nav.addItem(QListWidgetItem(nav_icon(n),n))
            if isinstance(p,DashboardPage):
                self.stack.addWidget(p)
            else:  # a narrow window scrolls instead of squeezing buttons until their text is cut
                sa=QScrollArea(); sa.setWidgetResizable(True); sa.setFrameShape(QFrame.Shape.NoFrame); sa.setWidget(p); self.stack.addWidget(sa)
        self._text_fitted=False
        self.nav.currentRowChanged.connect(self.stack.setCurrentIndex); self.nav.setCurrentRow(0)
        file_menu=self.menuBar().addMenu('&File'); aopen=QAction('Open RVE project…',self); asave=QAction('Save RVE project…',self); aopen.triggered.connect(self.open_project); asave.triggered.connect(self.save_project); file_menu.addAction(aopen); file_menu.addAction(asave); file_menu.addSeparator(); aquit=QAction('Quit',self); aquit.triggered.connect(self.close); file_menu.addAction(aquit)
        help_menu=self.menuBar().addMenu('&Help'); auser=QAction('User Guide',self); auser.setShortcut('F1'); auser.triggered.connect(lambda:self.help_page.open_guide('user')); help_menu.addAction(auser); atheory=QAction('Theory Guide',self); atheory.triggered.connect(lambda:self.help_page.open_guide('theory')); help_menu.addAction(atheory)
        self.statusBar().showMessage('Ready — analytical homogenisation, FE design & validation')
    def page_at(self,index:int):
        """The workspace page at a stack index (pages other than the dashboard sit in a scroll area)."""
        w=self.stack.widget(index)
        return w.widget() if isinstance(w,QScrollArea) else w
    def current_page(self):
        return self.page_at(self.stack.currentIndex())
    def showEvent(self,event):
        if not self._text_fitted:
            fit_text(self); self._text_fitted=True
        super().showEvent(event)
    def go_to(self,name:str):
        """Select a workspace by its navigation label (robust to page reordering)."""
        if name in self.page_names: self.nav.setCurrentRow(self.page_names.index(name))
    def save_project(self):
        try:
            model,mats,lams=self.single_page.parse_model(); f,_=QFileDialog.getSaveFileName(self,'Save RVE project','','RVE project (*.rvejson *.json)')
            if f: save_project(f,model,mats,lams); self.statusBar().showMessage(f'Saved project: {f}')
        except Exception as e: QMessageBox.critical(self,'Save project failed',str(e))
    def open_project(self):
        f,_=QFileDialog.getOpenFileName(self,'Open RVE project','','RVE project (*.rvejson *.json)')
        if not f:return
        self.open_project_file(f)

    def open_project_file(self,f):
        try:
            model,mats,lams=load_project(f); self.single_page.set_project(model,mats,lams); self.nav.setCurrentRow(1); self.statusBar().showMessage(f'Loaded project: {f}')
        except Exception as e: QMessageBox.critical(self,'Open project failed',str(e))


def launch(project=None):
    """Start the application; a .rvejson path (argument, or the first command-line argument, as passed by the Windows
    file association) is opened in the Single RVE Designer."""
    app=QApplication.instance() or QApplication(sys.argv); app.setStyleSheet(APP_QSS)
    icon=Path(__file__).resolve().parent/'assets'/'rve_studio.ico'
    if icon.is_file():
        from PySide6.QtGui import QIcon
        app.setWindowIcon(QIcon(str(icon)))
    w=MainWindow(); w.show()
    project=project or next((a for a in sys.argv[1:] if a.lower().endswith(('.rvejson','.json')) and os.path.isfile(a)),None)
    if project:
        w.open_project_file(project)
    return app.exec()
