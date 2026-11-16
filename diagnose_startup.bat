@echo off
setlocal
cd /d "%~dp0"
set "PYTHON_EXE="
if exist ".venv\Scripts\python.exe" set "PYTHON_EXE=.venv\Scripts\python.exe"
if not defined PYTHON_EXE (
  where py >nul 2>nul && set "PYTHON_EXE=py"
)
if not defined PYTHON_EXE set "PYTHON_EXE=python"
%PYTHON_EXE% diagnose_startup.py
set "RC=%ERRORLEVEL%"
echo.
echo Diagnostics saved to RVE_STUDIO_DIAGNOSTICS.txt
pause
exit /b %RC%
