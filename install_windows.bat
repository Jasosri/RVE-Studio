@echo off
setlocal EnableExtensions
cd /d "%~dp0"

echo RVE Studio - Windows installer
echo ============================================
echo.

set "PY_LAUNCHER="
set "PY_VERSION="

where py >nul 2>nul
if not errorlevel 1 goto :choose_py_launcher

goto :try_python_command

:choose_py_launcher
for %%V in (3.14 3.13 3.12 3.11 3.10) do (
  py -%%V -c "import sys,struct; raise SystemExit(0 if (3,10) <= sys.version_info[:2] < (3,15) and struct.calcsize('P')*8 == 64 else 1)" >nul 2>nul
  if not errorlevel 1 (
    set "PY_LAUNCHER=py"
    set "PY_VERSION=-%%V"
    goto :python_found
  )
)

goto :try_python_command

:try_python_command
where python >nul 2>nul
if errorlevel 1 goto :no_python
python -c "import sys,struct; raise SystemExit(0 if (3,10) <= sys.version_info[:2] < (3,15) and struct.calcsize('P')*8 == 64 else 1)" >nul 2>nul
if errorlevel 1 goto :unsupported_python
set "PY_LAUNCHER=python"
set "PY_VERSION="

goto :python_found

:no_python
echo ERROR: No suitable Python installation was found.
echo Install 64-bit Python 3.10, 3.11, 3.12, 3.13 or 3.14 from python.org.
echo Enable "Add Python to PATH" or install the Python Launcher.
goto :fail

:unsupported_python
echo ERROR: The available Python is not supported for this GUI environment.
echo Use a 64-bit Python version from 3.10 through 3.14.
goto :fail

:python_found
echo Using:
%PY_LAUNCHER% %PY_VERSION% -c "import sys,platform; print(sys.executable); print(sys.version); print(platform.architecture()[0])"
if errorlevel 1 goto :fail

if exist ".venv\Scripts\python.exe" (
  ".venv\Scripts\python.exe" -c "import sys,struct; raise SystemExit(0 if (3,10) <= sys.version_info[:2] < (3,15) and struct.calcsize('P')*8 == 64 else 1)" >nul 2>nul
  if errorlevel 1 (
    echo Existing .venv uses an incompatible Python. Recreating it...
    rmdir /s /q ".venv"
  ) else (
    echo Existing compatible .venv found. Repairing/updating it.
  )
)

if not exist ".venv\Scripts\python.exe" (
  echo Creating virtual environment...
  %PY_LAUNCHER% %PY_VERSION% -m venv .venv
  if errorlevel 1 goto :fail
)

set "VENV_PY=.venv\Scripts\python.exe"
if not exist "%VENV_PY%" goto :fail

echo Upgrading pip/setuptools/wheel...
"%VENV_PY%" -m pip install --upgrade pip setuptools wheel
if errorlevel 1 goto :fail

echo Installing runtime dependencies...
"%VENV_PY%" -m pip install -r requirements.txt
if errorlevel 1 goto :fail

echo Installing RVE Studio package...
"%VENV_PY%" -m pip install -e . --no-deps
if errorlevel 1 goto :fail

echo Verifying GUI imports...
"%VENV_PY%" -c "import PySide6, matplotlib, rvestudio; from PySide6.QtWidgets import QApplication; from rvestudio.gui.main_window import MainWindow; app=QApplication.instance() or QApplication([]); w=MainWindow(); print('GUI construction: PASS'); print('PySide6:', PySide6.__version__); print('RVE Studio:', rvestudio.__version__); w.close()"
if errorlevel 1 goto :fail

echo.
echo Installation completed successfully.
echo Start with START_RVE_STUDIO.bat
pause
exit /b 0

:fail
echo.
echo Installation failed.
echo Run diagnose_startup.bat to create RVE_STUDIO_DIAGNOSTICS.txt.
pause
exit /b 1
