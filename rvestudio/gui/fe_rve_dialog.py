"""Define FE RVE dialog: choose, or draw, the rectangular FE RVE of a design.

The design is shown replicated in the FE frame (x1 along a1).  Three choices:

* **Automatic** -- the rule of the source study (``fe_rve.select_fe_rve``), the default;
* **Block** -- n1 x n2 rectangular periodic cells centred on a complete cell;
* **Drawn** -- drag a rectangle over the replicated design.  Its size is rounded to
  whole rectangular periodic cells, so it is always exactly periodic and holds a whole
  number of primitive cells; its centre is where it was drawn (or the nearest cell
  centre, the default, which keeps the DBC model comparable with the automatic RVE).

The choice is returned by :meth:`FERVEDialog.spec` (``None`` = automatic) and is stored
on the design as ``RVEModel.fe_rve``, so every DBC and PBC model -- built-in and Abaqus --
uses it.  All geometry is done by :mod:`rvestudio.fe_rve`; this module only wires widgets.
"""
from __future__ import annotations

import copy

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox, QRadioButton, QButtonGroup, QSpinBox,
    QDoubleSpinBox, QComboBox, QCheckBox, QLabel, QPushButton, QDialogButtonBox, QWidget,
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT
from matplotlib.figure import Figure

from ..fe_rve import (MAX_RVE_CELLS, centred_window, fe_rve_layout, fe_rve_summary, normalise_fe_rve_spec,
                      select_fe_rve, snap_to_cell, spec_rectangle, window_from_drag)
from ..fe_plots import plot_rve_layout

AUTO, CELLS, WINDOW = "auto", "cells", "window"


class FERVEDialog(QDialog):
    def __init__(self, model, parent=None, *, title: str = ""):
        super().__init__(parent)
        self.setWindowTitle("Define FE RVE" + (f" — {title}" if title else ""))
        self.resize(1180, 720)
        self.model = copy.deepcopy(model)
        self.layout_data = fe_rve_layout(self.model)
        self._drag_start = None
        self._drag_rect = None
        self._busy = False
        # The spin boxes show the centre to 3 decimals; the exact value computed by snapping or
        # drawing is kept here so the rounding cannot cut a sliver off a wall.
        self._exact_centre = None
        self._ok_spec = None
        self._note = ""
        start = normalise_fe_rve_spec(getattr(model, "fe_rve", None))

        outer = QHBoxLayout(self)
        # ---------------------------------------------------------------- picture
        pic = QWidget()
        pv = QVBoxLayout(pic)
        pv.setContentsMargins(0, 0, 0, 0)
        self.fig = Figure(figsize=(7.5, 6.5))
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        pv.addWidget(self.toolbar)
        pv.addWidget(self.canvas, 1)
        hint = QLabel("Drag on the picture to draw the FE RVE (choose Drawn rectangle first, or just drag). "
                      "Grey lines: rectangular periodic cells; + : cell centres; dashed: automatic RVE. "
                      "Turn off zoom/pan in the toolbar to draw.")
        hint.setWordWrap(True)
        hint.setStyleSheet("color:#5B6B82;")
        pv.addWidget(hint)
        outer.addWidget(pic, 1)

        # ---------------------------------------------------------------- controls
        side = QWidget()
        side.setMinimumWidth(360)
        side.setMaximumWidth(420)
        sv = QVBoxLayout(side)
        box = QGroupBox("FE RVE")
        bv = QVBoxLayout(box)
        self.mode_group = QButtonGroup(self)
        self.rb_auto = QRadioButton("Automatic (rule of the source study)")
        self.rb_cells = QRadioButton("Block of rectangular cells, centred on a cell")
        self.rb_window = QRadioButton("Drawn rectangle")
        for k, rb in enumerate((self.rb_auto, self.rb_cells, self.rb_window)):
            self.mode_group.addButton(rb, k)
            bv.addWidget(rb)
        form = QFormLayout()
        rw, rh = self.layout_data["rect_width"], self.layout_data["rect_height"]
        self.n1 = QSpinBox()
        self.n2 = QSpinBox()
        for s in (self.n1, self.n2):
            s.setRange(1, MAX_RVE_CELLS)
        self.n1.setSuffix(f"  × {rw:.4g} mm")
        self.n2.setSuffix(f"  × {rh:.4g} mm")
        form.addRow("Cells along x1", self.n1)
        form.addRow("Cells along x2", self.n2)
        self.centre_cell = QComboBox()
        for c in self.layout_data["cell_types"]:
            self.centre_cell.addItem(f"{c['index']}: {c['name']} ({c['area']:.4g} mm²)", c["index"])
        form.addRow("Centre on", self.centre_cell)
        lo_x, hi_x, lo_y, hi_y = self.layout_data["extent"]
        self.cx = QDoubleSpinBox()
        self.cy = QDoubleSpinBox()
        for s, lo, hi in ((self.cx, lo_x, hi_x), (self.cy, lo_y, hi_y)):
            s.setRange(lo - 10 * (hi - lo), hi + 10 * (hi - lo))
            s.setDecimals(3)
            s.setSuffix(" mm")
            s.setSingleStep(max(rw, rh) / 20)
        form.addRow("Centre x1", self.cx)
        form.addRow("Centre x2", self.cy)
        self.snap = QCheckBox("Centre on the geometry (nearest cell centre)")
        self.snap.setChecked(True)
        self.snap.setToolTip("Keeps the faces cutting the walls at matching points, so the DBC model stays "
                             "comparable with the automatic RVE. PBC results do not depend on it.")
        form.addRow("", self.snap)
        self.snap_btn = QPushButton("Centre on nearest cell")
        self.snap_btn.clicked.connect(self._snap_now)
        form.addRow("", self.snap_btn)
        bv.addLayout(form)
        sv.addWidget(box)

        self.info = QLabel("")
        self.info.setWordWrap(True)
        self.info.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        sv.addWidget(self.info)
        note = QLabel("PBC terms do not depend on the RVE (it is periodic). DBC terms do: the automatic RVE is "
                      "that of the source study.")
        note.setWordWrap(True)
        note.setStyleSheet("color:#5B6B82;")
        sv.addWidget(note)
        sv.addStretch()
        self.buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        reset = self.buttons.addButton("Reset to automatic", QDialogButtonBox.ButtonRole.ResetRole)
        reset.clicked.connect(lambda: self.rb_auto.setChecked(True))
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        sv.addWidget(self.buttons)
        outer.addWidget(side)

        # ---------------------------------------------------------------- initial state
        auto_n = self.layout_data["auto_cells"]
        n1, n2 = (start["cells"] if start else auto_n)
        self.n1.setValue(int(n1))
        self.n2.setValue(int(n2))
        centre = (start or {}).get("centre") or self.layout_data["auto_centre"]
        self._show_centre(centre)
        if start and start["mode"] == CELLS:
            k = self.centre_cell.findData(int(start.get("centre_cell", 0)))
            self.centre_cell.setCurrentIndex(max(k, 0))
        if start and start["mode"] == WINDOW:
            # Keep a stored centre as it is, even if it is not a cell centre.
            self.snap.setChecked(snap_to_cell(self.layout_data, centre) == [float(centre[0]), float(centre[1])])
        {None: self.rb_auto, CELLS: self.rb_cells, WINDOW: self.rb_window}[start["mode"] if start else None].setChecked(True)

        for w in (self.n1, self.n2, self.cx, self.cy):
            w.valueChanged.connect(self.refresh)
        self.centre_cell.currentIndexChanged.connect(self.refresh)
        self.mode_group.idToggled.connect(lambda _id, on: on and self.refresh())
        self.snap.toggled.connect(lambda on: on and self._snap_now())
        self.canvas.mpl_connect("button_press_event", self._press)
        self.canvas.mpl_connect("motion_notify_event", self._motion)
        self.canvas.mpl_connect("button_release_event", self._release)
        self.refresh()

    # ---------------------------------------------------------------- state
    def mode(self) -> str:
        return {0: AUTO, 1: CELLS, 2: WINDOW}.get(self.mode_group.checkedId(), AUTO)

    def current_spec(self) -> dict | None:
        m = self.mode()
        if m == AUTO:
            return None
        cells = [self.n1.value(), self.n2.value()]
        if m == CELLS:
            return normalise_fe_rve_spec({"mode": CELLS, "cells": cells,
                                          "centre_cell": int(self.centre_cell.currentData() or 0)})
        return normalise_fe_rve_spec({"mode": WINDOW, "cells": cells, "centre": self._centre()})

    def _show_centre(self, centre) -> None:
        """Put a centre in the spin boxes and remember its exact value."""
        self.cx.setValue(float(centre[0]))
        self.cy.setValue(float(centre[1]))
        self._exact_centre = ([float(centre[0]), float(centre[1])], (self.cx.value(), self.cy.value()))

    def _centre(self) -> list[float]:
        """The exact centre, unless the user has since typed a different one."""
        shown = (self.cx.value(), self.cy.value())
        if self._exact_centre is not None and self._exact_centre[1] == shown:
            return list(self._exact_centre[0])
        return [shown[0], shown[1]]

    def spec(self) -> dict | None:
        """The accepted choice (None = automatic)."""
        return self._ok_spec

    def accept(self) -> None:
        try:
            self._ok_spec = self.current_spec()
        except ValueError:
            return
        super().accept()

    # ---------------------------------------------------------------- drawing
    def _set_quiet(self, n1=None, n2=None, centre=None) -> None:
        self._busy = True
        try:
            if n1 is not None:
                self.n1.setValue(int(n1))
            if n2 is not None:
                self.n2.setValue(int(n2))
            if centre is not None:
                self._show_centre(centre)
        finally:
            self._busy = False

    def _snap_now(self) -> None:
        spec, self._note = centred_window(self.model, self.layout_data, [self.n1.value(), self.n2.value()],
                                          self._centre(), snap=True)
        self._set_quiet(centre=spec["centre"])
        if self.mode() != WINDOW:
            self.rb_window.setChecked(True)
        else:
            self.refresh()

    def refresh(self, *_):
        if self._busy:
            return
        m = self.mode()
        for w in (self.n1, self.n2):
            w.setEnabled(m != AUTO)
        self.centre_cell.setEnabled(m == CELLS)
        for w in (self.cx, self.cy, self.snap, self.snap_btn):
            w.setEnabled(m == WINDOW)
        if m == AUTO:
            self._set_quiet(*self.layout_data["auto_cells"])
        ok, rve, rect = True, None, None
        try:
            spec = self.current_spec()
            trial = copy.deepcopy(self.model)
            trial.fe_rve = spec
            _blk, rve = select_fe_rve(trial)
            rect = spec_rectangle(self.layout_data, spec)
            warn = rve.get("selection", {}).get("warnings", [])
            body = fe_rve_summary({**rve, "selection": {**rve.get("selection", {}), "warnings": []}})
            note = f"<br><i>{self._note}</i>" if (self._note and self.mode() == WINDOW) else ""
            self.info.setText(f"<b>{rve['rve'][:1].upper()}{rve['rve'][1:]}</b><br>{body}{note}"
                              + "".join(f"<br><span style='color:#B03A2E'>⚠ {w}</span>" for w in warn))
        except Exception as exc:  # shown, never swallowed: the user can fix the input
            ok = False
            self.info.setText(f"<span style='color:#B03A2E'><b>Not a valid FE RVE:</b> {exc}</span>")
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(ok)
        self._draw(rect, rve)

    def _draw(self, rect, rve) -> None:
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        plot_rve_layout(ax, self.layout_data, rect=rect, rve=rve, drag=self._drag_rect)
        self.fig.tight_layout()
        self.canvas.draw_idle()

    def _press(self, ev) -> None:
        if ev.inaxes is None or ev.button != 1 or self.toolbar.mode:
            return
        self._drag_start = (ev.xdata, ev.ydata)

    def _motion(self, ev) -> None:
        if self._drag_start is None or ev.inaxes is None:
            return
        self._drag_rect = (*self._drag_start, ev.xdata, ev.ydata)
        ax = self.fig.axes[0] if self.fig.axes else None
        if ax is None:
            return
        # Light-weight feedback while dragging: only move the outline.
        if not hasattr(self, "_drag_patch") or self._drag_patch.axes is not ax:
            from matplotlib.patches import Rectangle
            self._drag_patch = Rectangle((0, 0), 0, 0, fc="#0072B2", alpha=0.08, ec="#0072B2", ls=":", zorder=7)
            ax.add_patch(self._drag_patch)
        x0, y0, x1, y1 = self._drag_rect
        self._drag_patch.set_bounds(min(x0, x1), min(y0, y1), abs(x1 - x0), abs(y1 - y0))
        self.canvas.draw_idle()

    def _release(self, ev) -> None:
        if self._drag_start is None:
            return
        x0, y0 = self._drag_start
        self._drag_start = None
        self._drag_rect = None
        x1, y1 = (ev.xdata, ev.ydata) if ev.inaxes is not None else (x0, y0)
        if abs(x1 - x0) < 1e-9 and abs(y1 - y0) < 1e-9:
            # A click moves the current rectangle's centre there.
            cells, centre = [self.n1.value(), self.n2.value()], [x0, y0]
        else:
            drawn = window_from_drag(self.layout_data, x0, y0, x1, y1, snap_centre=False)
            cells, centre = drawn["cells"], drawn["centre"]
        # Centred on the geometry (nearest cell centre) wherever possible.
        spec, self._note = centred_window(self.model, self.layout_data, cells, centre, snap=self.snap.isChecked())
        self.set_spec(spec)

    def set_spec(self, spec) -> None:
        """Apply a spec to the controls (used by drawing and by tests)."""
        s = normalise_fe_rve_spec(spec)
        if s is None:
            self.rb_auto.setChecked(True)
            self.refresh()
            return
        self._set_quiet(*s["cells"], centre=s.get("centre"))
        if s["mode"] == CELLS:
            k = self.centre_cell.findData(int(s.get("centre_cell", 0)))
            self._busy = True
            self.centre_cell.setCurrentIndex(max(k, 0))
            self._busy = False
        target = self.rb_cells if s["mode"] == CELLS else self.rb_window
        if target.isChecked():
            self.refresh()
        else:
            target.setChecked(True)
