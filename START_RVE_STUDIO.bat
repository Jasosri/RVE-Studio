@echo off
setlocal EnableExtensions
cd /d "%~dp0"

if not exist ".venv\Scripts\python.exe" (
  echo RVE Studio virtual environment was not found.
  echo Running the installer now...
  echo.
  call install_windows.bat
  if errorlevel 1 goto :fail_no_venv
)

call :check_gui
if errorlevel 1 (
  echo GUI dependency check failed. Attempting one repair installation...
  echo.
  call install_windows.bat
  if errorlevel 1 goto :fail
  call :check_gui
  if errorlevel 1 goto :fail
)

".venv\Scripts\python.exe" launch_rve_studio.py
set "RC=%ERRORLEVEL%"
if not "%RC%"=="0" goto :fail
exit /b 0

:check_gui
".venv\Scripts\python.exe" -c "import PySide6, matplotlib, rvestudio; from PySide6.QtWidgets import QApplication; from rvestudio.gui.main_window import MainWindow; app=QApplication.instance() or QApplication([]); w=MainWindow(); w.close()" >nul 2>"RVE_Studio_STARTUP_ERROR.txt"
exit /b %ERRORLEVEL%

:fail_no_venv
echo.
echo Installation did not complete, so RVE Studio cannot start.
echo Run diagnose_startup.bat after installing a supported 64-bit Python.
pause
exit /b 1

:fail
echo.
echo RVE Studio did not start successfully.
echo Running diagnostics; the output will be saved in this folder.
echo.
".venv\Scripts\python.exe" diagnose_startup.py 2>>"RVE_Studio_STARTUP_ERROR.txt"
echo.
echo Please keep RVE_Studio_STARTUP_ERROR.txt and RVE_STUDIO_DIAGNOSTICS.txt.
pause
exit /b 1
