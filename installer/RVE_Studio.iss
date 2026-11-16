; Inno Setup script: wraps dist\RVE_Studio (PyInstaller) into RVE_Studio_Setup_<version>.exe
; Build: build_windows_installer.bat, or  ISCC /DAppVersion=1.0.0 installer\RVE_Studio.iss
#define AppName "RVE Studio"
#ifndef AppVersion
  #define AppVersion "1.0.0"
#endif

[Setup]
AppId={{8C3E5A7E-6B1F-4C3A-9D3E-2F5B7A1C9E41}
AppName={#AppName}
AppVersion={#AppVersion}
AppVerName={#AppName} {#AppVersion}
AppPublisher=RVE Studio
AppPublisherURL=https://github.com/Jasosri/RVE-Studio
AppSupportURL=https://github.com/Jasosri/RVE-Studio/issues
DefaultDirName={autopf}\RVE Studio
DefaultGroupName=RVE Studio
DisableProgramGroupPage=yes
LicenseFile=..\LICENSE
OutputDir=..\dist
OutputBaseFilename=RVE_Studio_Setup_{#AppVersion}
SetupIconFile=rve_studio.ico
UninstallDisplayIcon={app}\RVE_Studio.exe
UninstallDisplayName={#AppName} {#AppVersion}
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
; Install for the current user without administrator rights by default; the dialog offers "all users".
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
ChangesAssociations=yes

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"
Name: "associate"; Description: "Open .rvejson project files with RVE Studio"; GroupDescription: "File association:"

[Files]
Source: "..\dist\RVE_Studio\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "..\docs\pdf\RVE_Studio_User_Guide.pdf"; DestDir: "{app}\docs"; Flags: ignoreversion
Source: "..\docs\pdf\RVE_Studio_Theory_Guide.pdf"; DestDir: "{app}\docs"; Flags: ignoreversion
Source: "..\LICENSE"; DestDir: "{app}"; DestName: "LICENSE.txt"; Flags: ignoreversion

[Icons]
Name: "{group}\RVE Studio"; Filename: "{app}\RVE_Studio.exe"
Name: "{group}\RVE Studio User Guide"; Filename: "{app}\docs\RVE_Studio_User_Guide.pdf"
Name: "{group}\RVE Studio Theory Guide"; Filename: "{app}\docs\RVE_Studio_Theory_Guide.pdf"
Name: "{group}\Uninstall RVE Studio"; Filename: "{uninstallexe}"
Name: "{autodesktop}\RVE Studio"; Filename: "{app}\RVE_Studio.exe"; Tasks: desktopicon

[Registry]
Root: HKA; Subkey: "Software\Classes\.rvejson\OpenWithProgids"; ValueType: string; ValueName: "RVEStudio.Project"; ValueData: ""; Flags: uninsdeletevalue; Tasks: associate
Root: HKA; Subkey: "Software\Classes\RVEStudio.Project"; ValueType: string; ValueName: ""; ValueData: "RVE Studio project"; Flags: uninsdeletekey; Tasks: associate
Root: HKA; Subkey: "Software\Classes\RVEStudio.Project\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\RVE_Studio.exe,0"; Tasks: associate
Root: HKA; Subkey: "Software\Classes\RVEStudio.Project\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\RVE_Studio.exe"" ""%1"""; Tasks: associate

[Run]
Filename: "{app}\RVE_Studio.exe"; Description: "{cm:LaunchProgram,RVE Studio}"; Flags: nowait postinstall skipifsilent

[UninstallDelete]
; user libraries in %USERPROFILE%\.rve_studio_v1 are kept on purpose (they hold the user's own materials and plug-ins)
Type: filesandordirs; Name: "{app}\__pycache__"
