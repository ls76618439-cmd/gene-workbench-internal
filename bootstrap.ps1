param(
  [string]$Repo = "ls76618439-cmd/gene-workbench-internal"
)

$ErrorActionPreference = "Stop"

function Fail([string]$Message) {
  Write-Error $Message
  exit 1
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
$expectedSha = ([string]$manifest.sha256).ToLowerInvariant()

$gh = Get-Command gh -ErrorAction SilentlyContinue
if (-not $gh) {
  Fail "GitHub CLI (gh) is required for automatic installation from this private repository."
}

& gh auth status *> $null
if ($LASTEXITCODE -ne 0) {
  Fail "GitHub CLI is not authenticated. Run 'gh auth login' once, then rerun this installer."
}

$tempRoot = Join-Path $env:TEMP "GeneWorkbenchBootstrap"
New-Item -ItemType Directory -Force -Path $tempRoot | Out-Null
$assetPath = Join-Path $tempRoot $assetName
if (Test-Path $assetPath) {
  Remove-Item -Force $assetPath
}

& gh release download $manifest.tag --repo $Repo --pattern $assetName --dir $tempRoot --clobber
if ($LASTEXITCODE -ne 0 -or -not (Test-Path $assetPath)) {
  Fail "Failed to download $assetName from GitHub Release $($manifest.tag)."
}

$actualSha = (Get-FileHash -Algorithm SHA256 $assetPath).Hash.ToLowerInvariant()
if ($actualSha -ne $expectedSha) {
  Fail "SHA256 mismatch. Expected $expectedSha but got $actualSha."
}

$proc = Start-Process -FilePath $assetPath -ArgumentList "/Q" -Wait -PassThru
if ($proc.ExitCode -ne 0) {
  Fail "Installer exited with code $($proc.ExitCode)."
}

$installExe = Join-Path $env:LOCALAPPDATA "GeneWorkbench\GeneWorkbench.exe"
$skillPath = Join-Path $env:USERPROFILE ".workbuddy\skills\gene-workbench\SKILL.md"
$mcpPath = Join-Path $env:USERPROFILE ".workbuddy\mcp.json"

if (-not (Test-Path $installExe)) { Fail "Installed GeneWorkbench.exe was not found." }
if (-not (Test-Path $skillPath)) { Fail "Gene Workbench Skill was not found." }
if (-not (Test-Path $mcpPath)) { Fail "WorkBuddy mcp.json was not found." }

$cfg = Get-Content $mcpPath -Raw | ConvertFrom-Json
if (-not $cfg.mcpServers.'gene-workbench') {
  Fail "gene-workbench MCP entry is missing from WorkBuddy mcp.json."
}
if ([string]$cfg.mcpServers.'gene-workbench'.command -ne $installExe) {
  Fail "gene-workbench MCP command does not point to the installed executable."
}

Write-Output "GENE_WORKBENCH_INSTALL_OK"
Write-Output "MCP=$mcpPath"
Write-Output "SKILL=$skillPath"
Write-Output "EXE=$installExe"
Write-Output "Restart or reload WorkBuddy once, then call gene-workbench tool_status for runtime discovery verification."
