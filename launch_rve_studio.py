from __future__ import annotations

import os
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent
# The installed application may sit in a read-only Program Files folder: log to the user's library folder then.
LOG = (Path.home() / ".rve_studio_v1" if getattr(sys, "frozen", False) else ROOT) / "RVE_Studio_STARTUP_ERROR.txt"


def _write_log(text: str) -> None:
    try:
        LOG.parent.mkdir(parents=True, exist_ok=True)
        LOG.write_text(text, encoding="utf-8")
    except Exception:
        pass


def _windows_message(title: str, message: str) -> None:
    if os.name != "nt":
        return
    try:
        import ctypes
        ctypes.windll.user32.MessageBoxW(0, message, title, 0x10)
    except Exception:
        pass


def _fail(exc: BaseException) -> int:
    tb = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
    py = f"Python: {sys.version}\nExecutable: {sys.executable}\nProject: {ROOT}\n\n"
    if isinstance(exc, ModuleNotFoundError) and exc.name == "PySide6":
        advice = (
            "PySide6 is not installed in the Python environment used to start RVE Studio.\n\n"
            "Run install_windows.bat first (or repair the .venv), then start again.\n"
            "If installation fails, run diagnose_startup.bat and inspect RVE_STUDIO_DIAGNOSTICS.txt.\n\n"
        )
    else:
        advice = (
            "RVE Studio encountered an error while starting.\n\n"
            "Run diagnose_startup.bat for a dependency/import check.\n"
            "The full traceback has been saved to RVE_Studio_STARTUP_ERROR.txt.\n\n"
        )
    text = py + advice + tb
    _write_log(text)
    print(text, file=sys.stderr)
    _windows_message("RVE Studio startup error", advice + f"Error: {exc}\n\nLog: {LOG}")
    return 1


def main() -> int:
    try:
        # Import the real Qt window here so both missing dependencies and GUI import
        # problems are caught by the same startup handler.
        from rvestudio.gui.main_window import launch
        return int(launch() or 0)
    except BaseException as exc:
        if isinstance(exc, KeyboardInterrupt):
            raise
        return _fail(exc)


if __name__ == "__main__":
    raise SystemExit(main())
