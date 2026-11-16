@echo off
REM RVE Studio - FE design & validation of the 14 benchmark designs.
REM 1. Built-in shell FE at 10 and 5 mm (no Abaqus, a few minutes).
REM 2. Abaqus at 10, 5 and 2.5 mm: PBC strain energy (full C, Q, engineering constants; 21 cases per mesh)
REM    and DBC (shell terms with U3 free, 3D terms with U3 = 0). Everything is compared with the analytical model.
REM Safe to re-run: finished jobs are kept, so you can stop and resume at any time.
setlocal
cd /d "%~dp0"
set "ROOT=%~dp0FE_Study_PaperA"
set "ABQ=abaqus"
if not "%~1"=="" set "ABQ=%~1"
set "PY=py"
if exist ".venv\Scripts\python.exe" set "PY=.venv\Scripts\python.exe"

"%PY%" -m rvestudio.cli fe check --root "%ROOT%" --paper-configs 1-14 --meshes 10 5
if errorlevel 1 goto :failed
"%PY%" -m rvestudio.cli fe doctor --abaqus-command "%ABQ%"
if errorlevel 1 (
  echo.
  echo Abaqus did not answer. The built-in results are in "%ROOT%\FE_STUDY_REPORT.html".
  echo To run Abaqus, start this from the Abaqus Command prompt or pass the full path:
  echo   RUN_FE_VERIFICATION_WINDOWS.bat "C:\SIMULIA\Commands\abaqus.bat"
  pause
  exit /b 1
)
"%PY%" -m rvestudio.cli fe run --root "%ROOT%" --meshes 10 5 2.5 --abaqus-command "%ABQ%" --cpus 4
if errorlevel 1 goto :failed
echo.
echo Done. Open "%ROOT%\FE_STUDY_REPORT.html" (or FE_STUDY_RESULTS.xlsx).
pause
exit /b 0

:failed
echo.
echo The study stopped with an error (see above). A study folder prepared with RVE Studio 1.10 or earlier
echo must be prepared again: rename or delete "%ROOT%" and run this script again.
pause
exit /b 1
