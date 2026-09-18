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

function Resolve-Gh {
  $cmd = Get-Command gh -ErrorAction SilentlyContinue
  if ($cmd) { return $cmd.Source }

  $candidate = Join-Path $env:ProgramFiles "GitHub CLI\gh.exe"
  if (Test-Path $candidate) { return $candidate }

  $winget = Get-Command winget -ErrorAction SilentlyContinue
  if (-not $winget) {
    Fail "GitHub CLI is missing and winget is unavailable. Install GitHub CLI, then rerun."
  }

  Write-Output "Installing GitHub CLI..."
  & winget install --id GitHub.cli --exact --silent --accept-package-agreements --accept-source-agreements
  if ($LASTEXITCODE -ne 0) {
    Fail "Failed to install GitHub CLI with winget."
  }

  if (Test-Path $candidate) { return $candidate }
  $cmd = Get-Command gh -ErrorAction SilentlyContinue
  if ($cmd) { return $cmd.Source }
  Fail "GitHub CLI installation completed but gh.exe could not be located."
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
$ghExe = Resolve-Gh

& $ghExe auth status *> $null
if ($LASTEXITCODE -ne 0) {
  Fail "GitHub is not authenticated for this private repository. Run 'gh auth login' once, then rerun."
}

$tempRoot = Join-Path $env:TEMP "GeneWorkbenchBootstrap"
New-Item -ItemType Directory -Force -Path $tempRoot | Out-Null
$assetPath = Join-Path $tempRoot $assetName
if (Test-Path $assetPath) { Remove-Item -Force $assetPath }

& $ghExe release download $manifest.tag --repo $Repo --pattern $assetName --dir $tempRoot --clobber
if ($LASTEXITCODE -ne 0 -or -not (Test-Path $assetPath)) {
  Fail "Failed to download $assetName from GitHub Release $($manifest.tag)."
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

$installExe = Join-Path $env:LOCALAPPDATA "GeneWorkbench\GeneWorkbench.exe"
$skillPath = Join-Path $env:USERPROFILE ".workbuddy\skills\gene-workbench\SKILL.md"
$mcpPath = Join-Path $env:USERPROFILE ".workbuddy\mcp.json"

if (-not (Test-Path $installExe)) { Fail "Installed Gene Workbench executable was not found: $installExe" }
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
Write-Output "Restart or reload WorkBuddy once, then call gene-workbench tool_status for runtime discovery verification."