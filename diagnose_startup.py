from __future__ import annotations

import importlib
import platform
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "RVE_STUDIO_DIAGNOSTICS.txt"


def check(name: str) -> tuple[bool, str]:
    try:
        m = importlib.import_module(name)
        ver = getattr(m, "__version__", "unknown")
        return True, f"PASS  {name:<36} version={ver}"
    except Exception as exc:
        return False, f"FAIL  {name:<36} {type(exc).__name__}: {exc}"


lines = [
    "RVE Studio startup diagnostics",
    "=" * 40,
    f"Platform:   {platform.platform()}",
    f"Python:     {sys.version}",
    f"Executable: {sys.executable}",
    f"Project:    {ROOT}",
    "",
]
for name in ["numpy", "scipy", "pandas", "matplotlib", "openpyxl", "PySide6"]:
    _, line = check(name)
    lines.append(line)

lines += ["", "Application import checks"]
for name in ["rvestudio", "rvestudio.geometry_analysis", "rvestudio.parametric_geometry", "rvestudio.gui.sketch", "rvestudio.gui.geometry_wizard", "rvestudio.gui.main_window"]:
    _, line = check(name)
    lines.append(line)

lines += ["", "Qt / MainWindow smoke check"]
try:
    from PySide6.QtWidgets import QApplication
    from rvestudio.gui.main_window import MainWindow
    app = QApplication.instance() or QApplication([])
    window = MainWindow()
    lines.append("PASS  QApplication creation")
    lines.append("PASS  MainWindow construction")
    lines.append(f"PASS  Window title: {window.windowTitle()}")
    window.close()
except Exception:
    lines.append("FAIL  Qt/MainWindow smoke check")
    lines.append(traceback.format_exc())

text = "\n".join(lines) + "\n"
OUT.write_text(text, encoding="utf-8")
print(text)
print(f"Saved to: {OUT}")
