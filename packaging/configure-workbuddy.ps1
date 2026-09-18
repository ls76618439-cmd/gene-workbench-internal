param(
  [Parameter(Mandatory=$true)][string]$InstallDir,
  [string]$ExecutableName = 'GeneWorkbench.exe'
)
$ErrorActionPreference = 'Stop'
$wb = Join-Path $env:USERPROFILE '.workbuddy'
New-Item -ItemType Directory -Force -Path $wb | Out-Null
$mcpPath = Join-Path $wb 'mcp.json'
if (Test-Path $mcpPath) {
  $raw = Get-Content $mcpPath -Raw
  if ([string]::IsNullOrWhiteSpace($raw)) { $cfg = [pscustomobject]@{ mcpServers = [pscustomobject]@{} } }
  else { $cfg = $raw | ConvertFrom-Json }
} else { $cfg = [pscustomobject]@{ mcpServers = [pscustomobject]@{} } }
if (-not $cfg.PSObject.Properties['mcpServers']) { $cfg | Add-Member -NotePropertyName mcpServers -NotePropertyValue ([pscustomobject]@{}) }
$exePath = Join-Path $InstallDir $ExecutableName
if (-not (Test-Path $exePath)) { throw "Gene Workbench executable not found: $exePath" }
$server = [pscustomobject]@{
  type = 'stdio'
  command = $exePath
  args = @()
  timeout = 120000
}
if ($cfg.mcpServers.PSObject.Properties['gene-workbench']) { $cfg.mcpServers.'gene-workbench' = $server }
else { $cfg.mcpServers | Add-Member -NotePropertyName 'gene-workbench' -NotePropertyValue $server }
$cfg | ConvertTo-Json -Depth 20 | Set-Content -Encoding UTF8 $mcpPath

$skillSrc = Join-Path $InstallDir 'skill\gene-workbench'
$skillDstRoot = Join-Path $wb 'skills'
$skillDst = Join-Path $skillDstRoot 'gene-workbench'
New-Item -ItemType Directory -Force -Path $skillDstRoot | Out-Null
if (Test-Path $skillDst) { Remove-Item -Recurse -Force $skillDst }
Copy-Item -Recurse -Force $skillSrc $skillDst
Write-Output "Gene Workbench registered with WorkBuddy: $exePath"
