@echo off
setlocal

set "APPDIR=%LOCALAPPDATA%\GeneWorkbench"
set "EXENAME=GeneWorkbench-1.2.0.exe"

if not exist "%APPDIR%" mkdir "%APPDIR%"
if not exist "%APPDIR%\skill\gene-workbench" mkdir "%APPDIR%\skill\gene-workbench"

if exist "%APPDIR%\%EXENAME%" (
  fc /b "%~dp0%EXENAME%" "%APPDIR%\%EXENAME%" >nul 2>nul
  if not errorlevel 1 goto exe_ready
)
copy /Y "%~dp0%EXENAME%" "%APPDIR%\%EXENAME%" >nul || exit /b 11
:exe_ready

copy /Y "%~dp0primer3_core.exe" "%APPDIR%\primer3_core.exe" >nul || exit /b 12
copy /Y "%~dp0seqkit.exe" "%APPDIR%\seqkit.exe" >nul || exit /b 13
copy /Y "%~dp0configure-workbuddy.ps1" "%APPDIR%\configure-workbuddy.ps1" >nul || exit /b 14
copy /Y "%~dp0unconfigure-workbuddy.ps1" "%APPDIR%\unconfigure-workbuddy.ps1" >nul || exit /b 15
copy /Y "%~dp0SKILL.md" "%APPDIR%\skill\gene-workbench\SKILL.md" >nul || exit /b 16

(
  echo @echo off
  echo powershell -NoProfile -ExecutionPolicy Bypass -File "%%LOCALAPPDATA%%\GeneWorkbench\unconfigure-workbuddy.ps1"
  echo powershell -NoProfile -ExecutionPolicy Bypass -Command "$d=Join-Path $env:LOCALAPPDATA 'GeneWorkbench'; Get-CimInstance Win32_Process ^| Where-Object { $_.ExecutablePath -like ($d+'\\GeneWorkbench*.exe') } ^| ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }; Start-Sleep -Milliseconds 300; Remove-Item $d -Recurse -Force -ErrorAction SilentlyContinue"
) > "%APPDIR%\Uninstall-GeneWorkbench.cmd"

powershell -NoProfile -ExecutionPolicy Bypass -File "%APPDIR%\configure-workbuddy.ps1" -InstallDir "%APPDIR%" -ExecutableName "%EXENAME%" || exit /b 20

echo Gene Workbench 1.2.0 installed. Restart or reload WorkBuddy to load the new MCP tools.
exit /b 0
