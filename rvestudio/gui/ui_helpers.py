"""Shared widgets: resizable table editor dialog, dashboard picture tiles and text-fitting helpers."""
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSize, Qt
from PySide6.QtGui import QColor, QPixmap
from PySide6.QtWidgets import (
    QAbstractButton, QComboBox, QDialog, QDoubleSpinBox, QFrame, QGridLayout, QGroupBox, QHBoxLayout, QLabel,
    QPushButton, QSizePolicy, QTabBar, QTabWidget, QVBoxLayout, QWidget,
)

ASSETS = Path(__file__).with_name("assets")


def asset(name: str) -> Path:
    return ASSETS / name


def fit_text(root: QWidget) -> None:
    """Never clip or elide text: buttons are at least as wide as their text, tabs are not elided."""
    for b in root.findChildren(QAbstractButton):
        if isinstance(b, QPushButton) and b.text():
            b.ensurePolished()
            b.setMinimumWidth(b.sizeHint().width())
    for bar in root.findChildren(QTabBar):
        bar.setElideMode(Qt.TextElideMode.ElideNone)
        bar.setUsesScrollButtons(True)
        bar.setExpanding(False)


class TableEditorDialog(QDialog):
    """Non-modal, resizable window holding the Node table, Wall table and Explicit periodic pairs.

    The tables keep their identity (the designer, the sketch and the wall editor keep working on the same
    widgets), so edits here update the sketch live and the dialog can stay open beside the drawing.
    """

    TITLES = ("Node table", "Wall table", "Explicit periodic pairs")

    def __init__(self, tabs: QTabWidget, parent=None, on_refresh=None):
        super().__init__(parent)
        self.setWindowTitle("RVE geometry tables")
        self.setWindowFlags(self.windowFlags() | Qt.WindowType.Window | Qt.WindowType.WindowMaximizeButtonHint)
        self.setModal(False)
        self.resize(1180, 640)
        self.setMinimumSize(820, 460)
        self.tabs = tabs
        v = QVBoxLayout(self); v.setContentsMargins(14, 12, 14, 12)
        head = QLabel("<b>Edit the RVE numerically.</b> Changes update the sketch immediately; keep this window open "
                      "beside the drawing, or close it and reopen it from the designer.")
        head.setWordWrap(True); head.setStyleSheet("color:#43506A;")
        v.addWidget(head)
        v.addWidget(tabs, 1)
        bar = QHBoxLayout(); bar.addStretch()
        if on_refresh is not None:
            r = QPushButton("Refresh sketch from tables"); r.clicked.connect(on_refresh); bar.addWidget(r)
        close = QPushButton("Close"); close.setObjectName("Primary"); close.clicked.connect(self.close); bar.addWidget(close)
        v.addLayout(bar)

    def open_tab(self, index: int) -> None:
        self.tabs.setCurrentIndex(index)
        self.setWindowTitle(f"RVE geometry tables — {self.TITLES[index]}")
        self.show(); self.raise_(); self.activateWindow()

    def showEvent(self, event):
        fit_text(self)
        super().showEvent(event)


class ScaledPicture(QLabel):
    """Picture that keeps its aspect ratio and follows the width of its tile."""

    def __init__(self, path: Path, ratio: float = 0.625, parent=None, fill: bool = False):
        super().__init__(parent)
        self._pix = QPixmap(str(path)) if Path(path).exists() else QPixmap()
        self._ratio = ratio
        self._fill = fill   # fill: take whatever height the layout gives (no height-for-width), never force scrolling
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)
        if fill:
            self.setMinimumSize(120, 48)
            self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Ignored)
        else:
            self.setMinimumSize(220, int(220 * ratio))
            self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        self.setStyleSheet("background:transparent;border:0;")
        if self._pix.isNull():
            self.setText("(picture missing)")

    def sizeHint(self):
        return QSize(360, int(360 * self._ratio))

    def minimumSizeHint(self):
        # A label's default minimum is its pixmap: that would stop the window from shrinking again.
        return QSize(self.minimumWidth(), self.minimumHeight())

    def hasHeightForWidth(self):
        return not self._fill

    def heightForWidth(self, w):
        return int(w * self._ratio)

    def resizeEvent(self, event):
        if not self._pix.isNull():
            w = max(1, self.width()); h = max(1, self.height())
            self.setPixmap(self._pix.scaled(w, h, Qt.AspectRatioMode.KeepAspectRatio,
                                            Qt.TransformationMode.SmoothTransformation))
        super().resizeEvent(event)


def picture_tile(image: str, title: str, text: str, button: str, colour: str, callback) -> QFrame:
    """Dashboard tile: coloured top strip, picture (shrinks with the window), title, description and an action button."""
    f = QFrame(); f.setObjectName("Tile")
    f.setStyleSheet(f"QFrame#Tile{{background:white;border:1px solid #E3E8EF;border-top:5px solid {colour};"
                    f"border-radius:10px;}} QFrame#Tile:hover{{border:1px solid {colour};border-top:5px solid {colour};}}")
    v = QVBoxLayout(f); v.setContentsMargins(12, 8, 12, 10); v.setSpacing(4)
    v.addWidget(ScaledPicture(asset(image), fill=True), 1)
    t = QLabel(title); t.setStyleSheet(f"font-size:12pt;font-weight:700;color:{colour};border:0;")
    d = QLabel(text); d.setWordWrap(True); d.setStyleSheet("color:#43506A;border:0;font-size:9.5pt;")
    b = QPushButton(button); b.setCursor(Qt.CursorShape.PointingHandCursor)
    b.setStyleSheet(f"QPushButton{{background:{colour};color:white;border:0;border-radius:6px;padding:6px 14px;"
                    f"font-weight:600;}}QPushButton:hover{{background:{QColor(colour).darker(115).name()};}}")
    b.clicked.connect(callback)
    v.addWidget(t); v.addWidget(d)
    row = QHBoxLayout(); row.addWidget(b); row.addStretch(); v.addLayout(row)
    return f


def stat_card(title: str, value: str, colour: str, note: str = "") -> QFrame:
    f = QFrame(); f.setObjectName("Stat")
    f.setStyleSheet(f"QFrame#Stat{{background:qlineargradient(x1:0,y1:0,x2:1,y2:1,stop:0 {colour},stop:1 {QColor(colour).lighter(130).name()});"
                    f"border:0;border-radius:12px;}} QLabel{{color:white;background:transparent;}}")
    v = QVBoxLayout(f); v.setContentsMargins(16, 12, 16, 12)
    a = QLabel(title.upper()); a.setStyleSheet("font-size:9pt;font-weight:700;letter-spacing:1px;")
    b = QLabel(value); b.setStyleSheet("font-size:24pt;font-weight:800;")
    v.addWidget(a); v.addWidget(b)
    if note:
        c = QLabel(note); c.setWordWrap(True); c.setStyleSheet("font-size:9pt;"); v.addWidget(c)
    return f


def hero_banner(title: str, subtitle: str, buttons: list[tuple[str, object]],
                chips: list[tuple[str, str, str]] | None = None) -> QFrame:
    """Gradient banner with a title, subtitle, action buttons, optional statistic chips (label, value, colour)
    and the 3D RVE picture."""
    f = QFrame(); f.setObjectName("Hero")
    f.setStyleSheet("QFrame#Hero{background:qlineargradient(x1:0,y1:0,x2:1,y2:0,stop:0 #132238,stop:0.55 #1F4E9A,"
                    "stop:1 #2F80ED);border-radius:14px;} QLabel{color:white;background:transparent;}")
    h = QHBoxLayout(f); h.setContentsMargins(22, 14, 16, 14)
    left = QVBoxLayout(); left.setSpacing(4)
    t = QLabel(title); t.setStyleSheet("font-size:22pt;font-weight:800;")
    s = QLabel(subtitle); s.setWordWrap(True); s.setStyleSheet("font-size:10.5pt;color:#DCE6F3;")
    left.addWidget(t); left.addWidget(s)
    row = QHBoxLayout()
    for k, (label, cb) in enumerate(buttons):
        b = QPushButton(label); b.setCursor(Qt.CursorShape.PointingHandCursor)
        b.setStyleSheet("QPushButton{background:%s;color:%s;border:0;border-radius:7px;padding:8px 16px;font-weight:700;}"
                        "QPushButton:hover{background:%s;}" % (("white", "#1F4E9A", "#EAF0F7") if k == 0 else
                                                               ("rgba(255,255,255,0.16)", "white", "rgba(255,255,255,0.28)")))
        b.clicked.connect(cb); row.addWidget(b)
    row.addStretch()
    for label, value, colour in chips or []:
        c = QLabel(f"<span style='font-size:15pt;font-weight:800'>{value}</span><br>"
                   f"<span style='font-size:8.5pt;letter-spacing:1px'>{label.upper()}</span>")
        c.setAlignment(Qt.AlignmentFlag.AlignCenter)
        c.setStyleSheet(f"background:{colour};border-radius:9px;padding:4px 12px;")
        row.addWidget(c)
    left.addSpacing(4); left.addLayout(row)
    h.addLayout(left, 3)
    pic = ScaledPicture(asset("rve3d_banner.png"), ratio=0.5, fill=True)
    pic.setMinimumSize(200, 100); pic.setMaximumHeight(170)
    h.addWidget(pic, 2)
    return f


class WallModelOptions(QGroupBox):
    """Wall beam (Euler-Bernoulli / Timoshenko), bending option and editable twist-release coefficients.

    Used by the benchmark and the database generator so both run any wall model the solver offers."""

    BENDING = (("Cylindrical (default)", "cylindrical"),
               ("Twist-corrected (calibrated, recommended for |ψ| > 0.3)", "twist-corrected"),
               ("Anisotropic strip (bound)", "anisotropic strip"),
               ("Free moment (lower bound)", "free moment"))

    def __init__(self, title: str = "Wall model", parent=None):
        super().__init__(title, parent)
        from ..laminate import TWIST_RELEASE_DEFAULT
        self._default = tuple(TWIST_RELEASE_DEFAULT)
        g = QGridLayout(self); g.setHorizontalSpacing(10)
        self.beam = QComboBox()
        self.beam.addItem("Euler–Bernoulli (default)", "euler-bernoulli")
        self.beam.addItem("Timoshenko (with wall shear flexibility)", "timoshenko")
        self.bending = QComboBox()
        for label, key in self.BENDING:
            self.bending.addItem(label, key)
        g.addWidget(QLabel("Wall beam"), 0, 0); g.addWidget(self.beam, 0, 1, 1, 3)
        g.addWidget(QLabel("Bending option"), 1, 0); g.addWidget(self.bending, 1, 1, 1, 3)
        self.coeff = []
        row = QHBoxLayout()
        row.addWidget(QLabel("η = (a + b ψ²) λⁿ:"))
        for lab, val, lo, hi in (("a", self._default[0], 0.0, 50.0), ("b", self._default[1], -50.0, 50.0),
                                 ("n", self._default[2], 0.05, 3.0)):
            sp = QDoubleSpinBox(); sp.setDecimals(3); sp.setRange(lo, hi); sp.setSingleStep(0.05); sp.setValue(val)
            row.addWidget(QLabel(lab)); row.addWidget(sp); self.coeff.append(sp)
        self.restore = QPushButton("Calibrated values")
        self.restore.clicked.connect(lambda: [sp.setValue(v) for sp, v in zip(self.coeff, self._default)])
        row.addWidget(self.restore); row.addStretch()
        self.twist_row = QWidget(); self.twist_row.setLayout(row); row.setContentsMargins(0, 0, 0, 0)
        g.addWidget(self.twist_row, 2, 0, 1, 4)
        self.bending.currentIndexChanged.connect(self._sync); self._sync()

    def _sync(self, *_):
        on = self.bending.currentData() == "twist-corrected"
        for w in (*self.coeff, self.restore):
            w.setEnabled(on)

    def wall_shear(self) -> str:
        return self.beam.currentData() or "euler-bernoulli"

    def bending_option(self) -> str:
        return self.bending.currentData() or "cylindrical"

    def twist_release(self):
        """User coefficients, or None when the calibrated defaults (or a non-twist option) are selected."""
        if self.bending_option() != "twist-corrected":
            return None
        c = tuple(round(float(sp.value()), 6) for sp in self.coeff)
        return None if all(abs(x - y) < 1e-9 for x, y in zip(c, self._default)) else c

    def describe(self) -> str:
        beam = "Euler–Bernoulli" if self.wall_shear() == "euler-bernoulli" else "Timoshenko"
        t = self.twist_release()
        extra = f", η coefficients {t}" if t else ""
        return f"{beam} walls, {self.bending_option()} bending{extra}"
