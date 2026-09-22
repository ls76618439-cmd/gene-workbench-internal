param(
  [string]$Repo = "ls76618439-cmd/gene-workbench-internal"
)

$ErrorActionPreference = "Stop"

function Fail([string]$Message) {
  Write-Error $Message
  exit 1
}

function Get-Sha256([string]$Path) {
  $sha = [System.Security.Cryptography.SHA256]::Create()
  $stream = [System.IO.File]::OpenRead($Path)
  try {
    $hash = $sha.ComputeHash($stream)
  } finally {
    $stream.Dispose()
    $sha.Dispose()
  }
  return (-join ($hash | ForEach-Object { $_.ToString("x2") }))
}

if ($env:OS -ne "Windows_NT") {
  Fail "Gene Workbench Company Edition currently supports Windows only."
}

$repoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$manifestPath = Join-Path $repoRoot "release-manifest.json"
if (-not (Test-Path $manifestPath)) {
  Fail "release-manifest.json is missing."
}

$manifest = Get-Content $manifestPath -Raw | ConvertFrom-Json
$assetName = [string]$manifest.asset
$executableName = [string]$manifest.executable
if ([string]::IsNullOrWhiteSpace($executableName)) { $executableName = "GeneWorkbench.exe" }
$expectedSha = ([string]$manifest.sha256).ToLowerInvariant()

$tempRoot = Join-Path $env:TEMP ("GeneWorkbenchBootstrap-" + [guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory -Force -Path $tempRoot | Out-Null
$assetPath = Join-Path $tempRoot $assetName

$downloadUrl = "https://github.com/$Repo/releases/download/$($manifest.tag)/$assetName"
try {
  Invoke-WebRequest -Uri $downloadUrl -OutFile $assetPath -UseBasicParsing
} catch {
  Fail "Failed to download $assetName from public GitHub Release $($manifest.tag): $($_.Exception.Message)"
}
if (-not (Test-Path $assetPath)) {
  Fail "Downloaded asset was not found: $assetPath"
}

$actualSha = (Get-Sha256 $assetPath).ToLowerInvariant()
if ($actualSha -ne $expectedSha) {
  Fail "SHA256 mismatch. Expected $expectedSha but got $actualSha."
}

if ([System.IO.Path]::GetExtension($assetPath).ToLowerInvariant() -eq ".zip") {
  $extractDir = Join-Path $tempRoot ("package-" + [string]$manifest.version)
  if (Test-Path $extractDir) { Remove-Item -Recurse -Force $extractDir }
  New-Item -ItemType Directory -Force -Path $extractDir | Out-Null
  Add-Type -AssemblyName System.IO.Compression.FileSystem
  [System.IO.Compression.ZipFile]::ExtractToDirectory($assetPath, $extractDir)
  $installCmd = Join-Path $extractDir "install.cmd"
  if (-not (Test-Path $installCmd)) { Fail "install.cmd is missing from the release package." }
  & $installCmd
  if ($LASTEXITCODE -ne 0) { Fail "install.cmd exited with code $LASTEXITCODE." }
} else {
  $proc = Start-Process -FilePath $assetPath -ArgumentList "/Q" -Wait -PassThru
  if ($proc.ExitCode -ne 0) { Fail "Installer exited with code $($proc.ExitCode)." }
}

$appDir = Join-Path $env:LOCALAPPDATA "GeneWorkbench"
$installExe = Join-Path $appDir $executableName
$primer3Core = Join-Path $appDir "primer3_core.exe"
$skillPath = Join-Path $env:USERPROFILE ".workbuddy\skills\gene-workbench\SKILL.md"
$mcpPath = Join-Path $env:USERPROFILE ".workbuddy\mcp.json"

if (-not (Test-Path $installExe)) { Fail "Installed Gene Workbench executable was not found: $installExe" }
if (-not (Test-Path $primer3Core)) { Fail "Installed Primer3 core executable was not found: $primer3Core" }
if (-not (Test-Path $skillPath)) { Fail "Gene Workbench Skill was not found." }
if (-not (Test-Path $mcpPath)) { Fail "WorkBuddy mcp.json was not found." }

$cfg = $null
$commandMatches = $false
for ($i = 0; $i -lt 10; $i++) {
  $cfg = Get-Content $mcpPath -Raw | ConvertFrom-Json
  if ($cfg.mcpServers.'gene-workbench') {
    $actualCommand = [string]$cfg.mcpServers.'gene-workbench'.command
    if ($actualCommand -eq $installExe) {
      $commandMatches = $true
      break
    }
  }
  Start-Sleep -Milliseconds 300
}
if (-not $cfg.mcpServers.'gene-workbench') {
  Fail "gene-workbench MCP entry is missing from WorkBuddy mcp.json."
}
if (-not $commandMatches) {
  Fail "gene-workbench MCP command does not point to the installed executable."
}

Write-Output "GENE_WORKBENCH_INSTALL_OK"
Write-Output "VERSION=$($manifest.version)"
Write-Output "MCP=$mcpPath"
Write-Output "SKILL=$skillPath"
Write-Output "EXE=$installExe"
Write-Output "PRIMER3_CORE=$primer3Core"
Write-Output "Restart or reload WorkBuddy once, then call gene-workbench tool_status for runtime discovery verification."
Remove-Item -Recurse -Force $tempRoot -ErrorAction SilentlyContinue
