@echo off
REM Build the Windows application and its installer:
REM   dist\RVE_Studio\RVE_Studio.exe          (PyInstaller, runs without Python)
REM   dist\RVE_Studio_Setup_<version>.exe     (Inno Setup installer, if Inno Setup 6 is installed)
REM Needs 64-bit Python 3.10+ (python.org) and, for the installer, Inno Setup 6 (https://jrsoftware.org/isinfo.php).
setlocal
cd /d "%~dp0"
if not exist .venv\Scripts\python.exe (
  echo Creating virtual environment .venv ...
  python -m venv .venv || goto :fail
)
call .venv\Scripts\activate.bat
python -m pip install --upgrade pip || goto :fail
python -m pip install . pyinstaller || goto :fail
for /f "delims=" %%v in ('python -c "import rvestudio; print(rvestudio.__version__)"') do set VERSION=%%v
echo Building RVE Studio %VERSION% ...
pyinstaller --noconfirm --clean RVE_Studio.spec || goto :fail
echo Application: dist\RVE_Studio\RVE_Studio.exe

set ISCC=
for %%p in ("%ProgramFiles(x86)%\Inno Setup 6\ISCC.exe" "%ProgramFiles%\Inno Setup 6\ISCC.exe" "%LocalAppData%\Programs\Inno Setup 6\ISCC.exe") do if exist %%p set ISCC=%%p
if not defined ISCC (
  echo.
  echo Inno Setup 6 was not found, so only the application folder was built.
  echo Install Inno Setup 6 and run this file again to build the installer.
  pause
  exit /b 0
)
%ISCC% /DAppVersion=%VERSION% installer\RVE_Studio.iss || goto :fail
echo.
echo Installer: dist\RVE_Studio_Setup_%VERSION%.exe
pause
exit /b 0
:fail
echo Build failed.
pause
exit /b 1
