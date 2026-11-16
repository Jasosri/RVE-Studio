"""Wall editor: a large pop-up for defining walls, with the Node I -> Node J local axis drawn.

The dialog edits plain data -- 1-based wall rows ``[node_i, node_j, weight, stack_id, group]``
against fixed node coordinates -- and returns the edited rows from :meth:`WallEditorDialog.rows`.
The Single RVE Designer copies them into its Wall table.  Geometry and drawing come from
:mod:`rvestudio.wall_axes`; this module only wires widgets.
"""
from __future__ import annotations

from PySide6.QtCore import Qt, QItemSelectionModel
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QTableWidget, QTableWidgetItem, QHeaderView,
    QAbstractItemView, QComboBox, QDialogButtonBox, QSplitter, QWidget,
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from ..wall_axes import COLUMN_TIPS, HTML, draw_wall_axes, html_note, wall_frames

COLS = ["ID", "Node I (start)", "Node J (end)", "Weight", "Stack ID", "Length [mm]", "Angle of 1̄ [deg]", "Wall group"]
TIPS = [None, COLUMN_TIPS["Node I"], COLUMN_TIPS["Node J"], COLUMN_TIPS["Weight"], COLUMN_TIPS["Stack ID"],
        COLUMN_TIPS["Length"], COLUMN_TIPS["Angle"], COLUMN_TIPS["Wall group"]]
C_ID, C_I, C_J, C_W, C_S, C_L, C_A, C_G = range(8)
BAD = QColor("#FBEBEB")


def _item(v, editable: bool = True) -> QTableWidgetItem:
    it = QTableWidgetItem("" if v is None else str(v))
    if not editable:
        it.setFlags(it.flags() & ~Qt.ItemFlag.ItemIsEditable)
        it.setForeground(QColor("#5B6B82"))
    return it


class WallEditorDialog(QDialog):
    def __init__(self, nodes, walls, a1=None, a2=None, stack_choices=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Wall editor — Node I → Node J defines the local 1̄ axis")
        self.resize(1320, 800)
        self.nodes = [tuple(map(float, p)) for p in nodes]
        self.a1, self.a2 = a1, a2
        self.stacks = [str(s) for s in (stack_choices or range(1, 51))]
        self._busy = False

        lay = QVBoxLayout(self)
        note = QLabel(html_note())
        note.setWordWrap(True)
        note.setTextFormat(Qt.TextFormat.RichText)
        note.setStyleSheet("background:#EEF4FB;border:1px solid #CFE0F5;border-radius:6px;padding:8px;")
        lay.addWidget(note)

        split = QSplitter(Qt.Orientation.Horizontal)
        left = QWidget()
        ll = QVBoxLayout(left)
        ll.setContentsMargins(0, 0, 6, 0)
        bar = QHBoxLayout()
        for text, slot, tip in (("Add wall", self.add_wall, "New wall from Node 1 to Node 2; edit the nodes in the table."),
                                ("Delete selected", self.delete_selected, ""),
                                ("Reverse I ↔ J", self.reverse_selected,
                                 "Swap Node I and Node J of the selected walls: reverses their local 1̄ (and 2̄) axes.")):
            b = QPushButton(text)
            b.setToolTip(tip)
            b.clicked.connect(slot)
            bar.addWidget(b)
        bar.addStretch()
        ll.addLayout(bar)
        self.table = QTableWidget(0, len(COLS))
        self.table.setHorizontalHeaderLabels(COLS)
        for c, tip in enumerate(TIPS):
            if tip:
                self.table.horizontalHeaderItem(c).setToolTip(tip)
        hh = self.table.horizontalHeader()
        hh.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        hh.setSectionResizeMode(C_ID, QHeaderView.ResizeMode.ResizeToContents)
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.table.setAlternatingRowColors(True)
        ll.addWidget(self.table, 1)
        self.status = QLabel("")
        self.status.setWordWrap(True)
        ll.addWidget(self.status)
        split.addWidget(left)

        right = QWidget()
        rl = QVBoxLayout(right)
        rl.setContentsMargins(6, 0, 0, 0)
        self.fig = Figure(figsize=(6.2, 6.2))
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setMinimumSize(420, 420)
        rl.addWidget(self.canvas, 1)
        hint = QLabel(f"Arrows point from Node I to Node J (positive {HTML[1]}); green arrows are the wall normals "
                      f"{HTML[2]}. Select rows to highlight walls; neighbouring cells are grey.")
        hint.setWordWrap(True)
        hint.setTextFormat(Qt.TextFormat.RichText)
        hint.setStyleSheet("color:#5B6B82;")
        rl.addWidget(hint)
        split.addWidget(right)
        split.setSizes([760, 560])
        lay.addWidget(split, 1)

        self.buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Apply to Wall table")
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        lay.addWidget(self.buttons)

        self._busy = True
        for w in walls:
            self._append(w)
        self._busy = False
        self.table.itemChanged.connect(self._changed)
        self.table.itemSelectionChanged.connect(self.refresh)
        self.refresh()

    # ------------------------------------------------------------------ rows
    def _append(self, w) -> None:
        w = list(w) + [None] * (5 - len(w))
        ni, nj, wt, sid, grp = w[:5]
        r = self.table.rowCount()
        self.table.insertRow(r)
        self.table.setItem(r, C_ID, _item(r + 1, False))
        self.table.setItem(r, C_I, _item(int(ni) if ni not in (None, "") else ""))
        self.table.setItem(r, C_J, _item(int(nj) if nj not in (None, "") else ""))
        self.table.setItem(r, C_W, _item(1.0 if wt in (None, "") else float(wt)))
        combo = QComboBox()
        combo.addItems(self.stacks)
        s = str(int(float(sid))) if sid not in (None, "") else self.stacks[0]
        if combo.findText(s) < 0:
            combo.addItem(s)
        combo.setCurrentText(s)
        self.table.setCellWidget(r, C_S, combo)
        self.table.setItem(r, C_L, _item("", False))
        self.table.setItem(r, C_A, _item("", False))
        self.table.setItem(r, C_G, _item(grp or ""))

    def _text(self, r: int, c: int) -> str:
        if c == C_S:
            w = self.table.cellWidget(r, c)
            return w.currentText().strip() if w is not None else ""
        it = self.table.item(r, c)
        return it.text().strip() if it else ""

    def _node_pairs(self) -> list[tuple]:
        out = []
        for r in range(self.table.rowCount()):
            try:
                out.append((int(float(self._text(r, C_I))), int(float(self._text(r, C_J)))))
            except ValueError:
                out.append((0, 0))
        return out

    def problems(self) -> list[str]:
        """Plain-language problems that stop the walls being applied (empty = fine)."""
        n = len(self.nodes)
        bad = []
        seen = {}
        for r, (i, j) in enumerate(self._node_pairs()):
            if not (1 <= i <= n and 1 <= j <= n):
                bad.append(f"W{r + 1}: Node I and Node J must be node numbers 1–{n}.")
            elif i == j:
                bad.append(f"W{r + 1}: Node I and Node J are the same node.")
            else:
                key = frozenset((i, j))
                if key in seen:
                    bad.append(f"W{r + 1} repeats W{seen[key] + 1} (same two nodes).")
                seen.setdefault(key, r)
            try:
                if float(self._text(r, C_W) or 1.0) <= 0:
                    bad.append(f"W{r + 1}: weight must be positive.")
            except ValueError:
                bad.append(f"W{r + 1}: weight must be a number.")
        return bad

    def rows(self) -> list[list]:
        """Edited walls, 1-based: [node_i, node_j, weight, stack_id, group]."""
        out = []
        for r, (i, j) in enumerate(self._node_pairs()):
            out.append([i, j, float(self._text(r, C_W) or 1.0), int(float(self._text(r, C_S) or 1)), self._text(r, C_G)])
        return out

    # ------------------------------------------------------------------ actions
    def select_rows(self, rows) -> None:
        """Select several rows at once (selectRow() would keep only the last)."""
        sm = self.table.selectionModel()
        sm.clearSelection()
        flags = QItemSelectionModel.SelectionFlag.Select | QItemSelectionModel.SelectionFlag.Rows
        for r in rows:
            if 0 <= int(r) < self.table.rowCount():
                sm.select(self.table.model().index(int(r), 0), flags)

    def selected_rows(self) -> list[int]:
        return sorted({i.row() for i in self.table.selectedIndexes()})

    def add_wall(self) -> None:
        self._busy = True
        self._append([1, min(2, max(1, len(self.nodes))), 1.0, self.stacks[0], ""])
        self._busy = False
        self.table.selectRow(self.table.rowCount() - 1)
        self.refresh()

    def delete_selected(self) -> None:
        self._busy = True
        for r in reversed(self.selected_rows()):
            self.table.removeRow(r)
        for r in range(self.table.rowCount()):
            self.table.setItem(r, C_ID, _item(r + 1, False))
        self._busy = False
        self.refresh()

    def reverse_selected(self) -> None:
        rows = self.selected_rows()
        if not rows:
            self.status.setText("Select the walls to reverse first.")
            return
        self._busy = True
        for r in rows:
            a, b = self._text(r, C_I), self._text(r, C_J)
            self.table.setItem(r, C_I, _item(b))
            self.table.setItem(r, C_J, _item(a))
        self._busy = False
        self.select_rows(rows)
        self.refresh()
        self.status.setText(self.status.text() + f"  Reversed {', '.join(f'W{r + 1}' for r in rows)}: their "
                            "local 1̄ axes now point the other way.")

    def _changed(self, *_):
        if not self._busy:
            self.refresh()

    def refresh(self) -> None:
        frames = wall_frames(self.nodes, self._node_pairs()) if self.nodes else []
        self._busy = True
        try:
            for r, f in enumerate(frames):
                ok = f.get("ok", False)
                self.table.item(r, C_L).setText(f"{f['length']:.6g}" if ok else "")
                self.table.item(r, C_A).setText(f"{f['angle_deg']:.6g}" if ok else "")
                for c in (C_I, C_J):
                    it = self.table.item(r, c)
                    if it is not None:
                        it.setBackground(QColor("white") if ok else BAD)
        finally:
            self._busy = False
        bad = self.problems()
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(not bad)
        self.status.setText(("<span style='color:#B03A2E'>" + "<br>".join(bad[:6]) + "</span>") if bad else
                            f"{self.table.rowCount()} walls. Angle of 1̄ = direction from Node I to Node J, "
                            "counter-clockwise from x1.")
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        if self.nodes:
            # Invalid rows are skipped by the drawing; wall numbers stay those of the table.
            draw_wall_axes(ax, self.nodes, self._node_pairs(), self.a1, self.a2, highlight=self.selected_rows())
        self.fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.14)
        self.canvas.draw_idle()

