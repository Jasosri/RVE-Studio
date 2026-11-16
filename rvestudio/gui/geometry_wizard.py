from __future__ import annotations

import math
import re
from dataclasses import dataclass

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFormLayout, QGroupBox,
    QHBoxLayout, QLabel, QLineEdit, QMessageBox, QPushButton, QTableWidget,
    QTableWidgetItem, QVBoxLayout, QWidget, QHeaderView
)


@dataclass
class WizardParameter:
    name: str
    reference: float
    minimum: float
    maximum: float
    driver: str = ""
    source: str = ""


class ParameterEditorDialog(QDialog):
    """Edit the explicit independent parameter list for a geometry-library row."""
    def __init__(self, parameters: list[dict] | None = None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Geometry parameters")
        self.resize(620, 360)
        v=QVBoxLayout(self)
        note=QLabel("Each row is one independent design variable used by the geometry generator. The name must exactly match the variable read by generate(p).")
        note.setWordWrap(True); v.addWidget(note)
        bar=QHBoxLayout(); add=QPushButton("Add parameter"); rem=QPushButton("Delete selected"); bar.addWidget(add); bar.addWidget(rem); bar.addStretch(); v.addLayout(bar)
        self.table=QTableWidget(0,3); self.table.setHorizontalHeaderLabels(["Parameter name","Minimum","Maximum"]); self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch); v.addWidget(self.table,1)
        add.clicked.connect(self.add_row); rem.clicked.connect(self.delete_row)
        for p in parameters or []: self.add_row(p)
        buttons=QDialogButtonBox(QDialogButtonBox.StandardButton.Ok|QDialogButtonBox.StandardButton.Cancel); buttons.accepted.connect(self._accept); buttons.rejected.connect(self.reject); v.addWidget(buttons)

    def add_row(self,p=None):
        r=self.table.rowCount(); self.table.insertRow(r); p=p or {}
        for c,val in enumerate([p.get("name",f"P{r+1}"),p.get("min",0.0),p.get("max",1.0)]): self.table.setItem(r,c,QTableWidgetItem(str(val)))

    def delete_row(self):
        r=self.table.currentRow()
        if r>=0:self.table.removeRow(r)

    def parameters(self):
        out=[]; names=set()
        for r in range(self.table.rowCount()):
            name=(self.table.item(r,0).text() if self.table.item(r,0) else "").strip()
            if not name: continue
            if not re.match(r"^[A-Za-z_][A-Za-z0-9_]*$",name): raise ValueError(f"Invalid parameter name: {name}")
            if name in names: raise ValueError(f"Duplicate parameter name: {name}")
            names.add(name)
            mn=float(self.table.item(r,1).text()); mx=float(self.table.item(r,2).text())
            if not math.isfinite(mn) or not math.isfinite(mx) or mx < mn: raise ValueError(f"Invalid range for {name}")
            out.append({"name":name,"min":mn,"max":mx})
        return out

    def _accept(self):
        try:self.parameters()
        except Exception as e: QMessageBox.warning(self,"Geometry parameters",str(e)); return
        self.accept()


class GeometryDefinitionWizard(QDialog):
    """Guided creation of a geometry-library definition from the current sketch.

    The default action always creates a *new* reusable geometry family from the
    current sketch.  If the sketch is recognised as a built-in family, the user
    may instead choose the explicit "use detected built-in" mode to edit that
    family's sampling ranges.
    """
    def __init__(self, summary, parent=None, existing_ids=None):
        super().__init__(parent)
        self.summary=summary
        self.existing_ids={str(x).upper() for x in (existing_ids or [])}
        self._custom_id=self._next_custom_id()
        self.setWindowTitle("Geometry Definition Wizard")
        self.resize(760,540)
        v=QVBoxLayout(self)
        intro=QLabel(
            "Create a reusable geometry definition from the current RVE. "
            "By default this creates a new parametric family from the sketch. "
            "If the sketch matches a built-in family, you can explicitly choose "
            "the built-in mode to edit that existing family's parameter ranges."
        )
        intro.setWordWrap(True); v.addWidget(intro)

        identity=QGroupBox("1. Geometry identity"); f=QFormLayout(identity)
        self.id_edit=QLineEdit(self._custom_id)
        default_name=(f"Custom {summary.detected_geometry_name}" if summary.detected_geometry_id else "Custom periodic RVE")
        self.name_edit=QLineEdit(default_name)
        self.mode=QComboBox()
        # New custom family is deliberately first/default.  This avoids the old
        # behaviour where "Create" silently edited the detected built-in row.
        self.mode.addItem("Create NEW family from this sketch (width / height)", "scale_xy")
        if summary.detected_geometry_id:
            self.mode.addItem(f"Use existing detected family ({summary.detected_geometry_id})", "builtin")
        f.addRow("Library ID",self.id_edit); f.addRow("Name",self.name_edit); f.addRow("Definition method",self.mode); v.addWidget(identity)

        detection=QGroupBox("2. Recognition evidence"); dl=QVBoxLayout(detection)
        if summary.detected_geometry_id:
            txt=f"Detected {summary.detected_geometry_id} — {summary.detected_geometry_name} · confidence {summary.detection_confidence*100:.1f}%"
        else:
            txt="No built-in topology matched. The generic periodic solver can still use this RVE."
        dl.addWidget(QLabel(txt))
        src=[]
        for name,val in summary.detected_parameters.items():
            walls=", ".join(f"W{i}" for i in summary.parameter_sources.get(name,[])) or "geometry"
            src.append((f"{name} = {val:.6g}" if name=='R' else f"{name} = {val:.6g} mm")+f"   ← {walls}")
        ev=QLabel("\n".join(src) if src else "Custom geometry: width and height will be the initial independent dimensions.")
        ev.setWordWrap(True); dl.addWidget(ev); v.addWidget(detection)

        params=QGroupBox("3. Independent parameters and sampling bounds"); pl=QVBoxLayout(params)
        self.table=QTableWidget(0,6)
        self.table.setHorizontalHeaderLabels(["Name","Meaning / driver","Reference","Minimum","Maximum","Evidence"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        pl.addWidget(self.table); v.addWidget(params,1)

        # IMPORTANT: create the preview object before _populate().  In v1.3.0,
        # _populate() called _update_preview() first and caused the exact
        # AttributeError reported by users.
        preview=QGroupBox("4. Range preview"); pv=QVBoxLayout(preview)
        self.preview=QLabel(); self.preview.setWordWrap(True); pv.addWidget(self.preview)
        btn=QPushButton("Check parameter ranges"); btn.clicked.connect(self._update_preview); pv.addWidget(btn)
        v.addWidget(preview)

        self.mode.currentIndexChanged.connect(self._mode_changed)
        self._populate()

        buttons=QDialogButtonBox(QDialogButtonBox.StandardButton.Save|QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self._accept); buttons.rejected.connect(self.reject); v.addWidget(buttons)

    def _next_custom_id(self):
        base=(f"CUSTOM_{self.summary.detected_geometry_id}" if self.summary.detected_geometry_id else "CUSTOM_RVE")
        candidate=base; n=2
        while candidate.upper() in self.existing_ids:
            candidate=f"{base}_{n}"; n+=1
        return candidate

    def _mode_changed(self):
        mode=self.mode.currentData()
        if mode=="builtin" and self.summary.detected_geometry_id:
            self.id_edit.setText(self.summary.detected_geometry_id)
            self.name_edit.setText(self.summary.detected_geometry_name or self.summary.detected_geometry_id)
            self.id_edit.setReadOnly(True)
        else:
            self.id_edit.setReadOnly(False)
            current=self.id_edit.text().strip().upper()
            if not current or current==str(self.summary.detected_geometry_id or '').upper():
                self.id_edit.setText(self._custom_id)
            if not self.name_edit.text().strip() or self.name_edit.text().strip()==self.summary.detected_geometry_name:
                self.name_edit.setText(f"Custom {self.summary.detected_geometry_name}" if self.summary.detected_geometry_id else "Custom periodic RVE")
        self._populate()

    def _set_row(self,r,name,driver,ref,mn,mx,evidence=""):
        for c,val in enumerate([name,driver,f"{ref:.9g}",f"{mn:.9g}",f"{mx:.9g}",evidence]):
            it=QTableWidgetItem(str(val))
            if c in (1,2,5): it.setFlags(it.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(r,c,it)

    def _populate(self):
        self.table.setRowCount(0)
        mode=self.mode.currentData()
        if mode=="builtin" and self.summary.detected_geometry_id:
            for name,ref in self.summary.detected_parameters.items():
                if not math.isfinite(ref): continue
                r=self.table.rowCount(); self.table.insertRow(r)
                delta=0.30*abs(ref) if abs(ref)>1e-12 else 1.0
                if name=='R': mn=max(1e-6,ref-delta); mx=ref+delta; unit="ratio"
                else: mn=max(1e-9,ref-delta); mx=ref+delta; unit="wall length"
                evidence=", ".join(f"W{i}" for i in self.summary.parameter_sources.get(name,[]))
                self._set_row(r,name,unit,ref,mn,mx,evidence)
        else:
            for name,driver,ref in [("W","drawing width",self.summary.width),("H","drawing height",self.summary.height_2d)]:
                if not math.isfinite(ref) or ref<=1e-12: continue
                r=self.table.rowCount(); self.table.insertRow(r)
                self._set_row(r,name,driver,ref,0.7*ref,1.3*ref,"bounding dimension")
        self._update_preview()

    def parameters(self):
        out=[]; names=set()
        for r in range(self.table.rowCount()):
            name=self.table.item(r,0).text().strip(); driver=self.table.item(r,1).text().strip()
            ref=float(self.table.item(r,2).text()); mn=float(self.table.item(r,3).text()); mx=float(self.table.item(r,4).text())
            if not re.match(r"^[A-Za-z_][A-Za-z0-9_]*$",name): raise ValueError(f"Invalid parameter name: {name}")
            if name in names: raise ValueError(f"Duplicate parameter name: {name}")
            names.add(name)
            if not all(math.isfinite(v) for v in (ref,mn,mx)) or mn<=0 or mx<mn or not (mn<=ref<=mx):
                raise ValueError(f"Invalid range for {name}: require 0 < minimum <= reference <= maximum")
            evidence=self.table.item(r,5).text().strip() if self.table.item(r,5) else ""
            out.append(WizardParameter(name,ref,mn,mx,driver,evidence))
        if not out: raise ValueError("At least one independent parameter is required.")
        return out

    def _update_preview(self):
        # Defensive guard as well as the corrected construction order above.
        if not hasattr(self,"preview"): return
        try:
            ps=self.parameters()
            self.preview.setText(" · ".join(f"{p.name}: {p.minimum:.5g} ← {p.reference:.5g} → {p.maximum:.5g}" for p in ps))
        except Exception as e:
            self.preview.setText(str(e))

    def definition(self):
        gid=re.sub(r"[^A-Za-z0-9_]+","_",self.id_edit.text().strip()).upper().strip("_")
        if not gid: raise ValueError("A valid geometry ID is required.")
        mode=self.mode.currentData()
        if mode!="builtin" and gid in self.existing_ids:
            raise ValueError(f"Geometry ID {gid} already exists. Choose a new ID.")
        ps=self.parameters()
        return {
            "id":gid,"name":self.name_edit.text().strip() or gid,"mode":mode,
            "parameters":[{"name":p.name,"min":p.minimum,"max":p.maximum,"reference":p.reference,"driver":p.driver} for p in ps],
        }

    def _accept(self):
        try:self.definition()
        except Exception as e:
            QMessageBox.warning(self,"Geometry Definition Wizard",str(e)); return
        self.accept()
