# -*- mode: python ; coding: utf-8 -*-
# PyInstaller spec for the Windows application (dist\RVE_Studio\RVE_Studio.exe).
# Build with build_windows_installer.bat, which also wraps it in an installer (installer\RVE_Studio.iss).
from PyInstaller.utils.hooks import collect_all

datas = [('rvestudio/defaults', 'rvestudio/defaults'),      # default libraries, copied to ~/.rve_studio_v1
         ('rvestudio/gui/assets', 'rvestudio/gui/assets'),
         ('rvestudio/docs', 'rvestudio/docs')]              # user and theory guides opened from Help
binaries = []
hiddenimports = ['rvestudio.gui.main_window', 'matplotlib.backends.backend_qtagg']
for pkg in ('PySide6', 'matplotlib'):
    d, b, h = collect_all(pkg)
    datas += d; binaries += b; hiddenimports += h

# Qt modules the application does not use: left out to keep the installer small.
UNUSED = ('WebEngine', 'Qt3D', 'Quick3D', 'QtQuick', 'QtQml', 'Multimedia', 'QtPdf', 'QtCharts', 'DataVisualization',
          'QtGraphs', 'QtLocation', 'QtPositioning', 'QtSensors', 'QtSerial', 'QtBluetooth', 'QtNfc', 'QtWebSockets',
          'QtWebChannel', 'QtWebView', 'QtTextToSpeech', 'QtRemoteObjects', 'QtScxml', 'QtSpatialAudio', 'designer',
          'translations', 'qml', 'examples')
keep = lambda item: not any(u.lower() in str(item[0]).replace('\\', '/').lower() for u in UNUSED)
datas = [x for x in datas if keep(x)]
binaries = [x for x in binaries if keep(x)]
hiddenimports = [h for h in hiddenimports if not any(u.lower() in h.lower() for u in UNUSED)]

a = Analysis(
    ['launch_rve_studio.py'],
    pathex=[],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['tkinter', 'pytest', 'IPython', 'PyQt5', 'PyQt6'],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='RVE_Studio',
    icon='installer/rve_studio.ico',
    version='installer/version_info.txt',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name='RVE_Studio',
)
