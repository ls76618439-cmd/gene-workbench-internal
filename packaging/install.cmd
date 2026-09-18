@echo off
setlocal

set "APPDIR=%LOCALAPPDATA%\GeneWorkbench"

taskkill /IM GeneWorkbench.exe /F >nul 2>nul

if not exist "%APPDIR%" mkdir "%APPDIR%"
if not exist "%APPDIR%\skill\gene-workbench" mkdir "%APPDIR%\skill\gene-workbench"

copy /Y "%~dp0GeneWorkbench.exe" "%APPDIR%\GeneWorkbench.exe" >nul || exit /b 11
copy /Y "%~dp0seqkit.exe" "%APPDIR%\seqkit.exe" >nul || exit /b 12
copy /Y "%~dp0configure-workbuddy.ps1" "%APPDIR%\configure-workbuddy.ps1" >nul || exit /b 13
copy /Y "%~dp0unconfigure-workbuddy.ps1" "%APPDIR%\unconfigure-workbuddy.ps1" >nul || exit /b 14
copy /Y "%~dp0SKILL.md" "%APPDIR%\skill\gene-workbench\SKILL.md" >nul || exit /b 15

(
  echo @echo off
  echo taskkill /IM GeneWorkbench.exe /F ^>nul 2^>nul
  echo powershell -NoProfile -ExecutionPolicy Bypass -File "%%LOCALAPPDATA%%\GeneWorkbench\unconfigure-workbuddy.ps1"
  echo rmdir /S /Q "%%LOCALAPPDATA%%\GeneWorkbench"
) > "%APPDIR%\Uninstall-GeneWorkbench.cmd"

powershell -NoProfile -ExecutionPolicy Bypass -File "%APPDIR%\configure-workbuddy.ps1" -InstallDir "%APPDIR%" || exit /b 20

echo Gene Workbench installed. Restart or reload WorkBuddy to load the MCP tools.
exit /b 0
