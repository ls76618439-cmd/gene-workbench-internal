$ErrorActionPreference = 'SilentlyContinue'
$wb = Join-Path $env:USERPROFILE '.workbuddy'
$mcpPath = Join-Path $wb 'mcp.json'
if (Test-Path $mcpPath) {
  $raw = Get-Content $mcpPath -Raw
  if (-not [string]::IsNullOrWhiteSpace($raw)) {
    $cfg = $raw | ConvertFrom-Json
    if ($cfg.PSObject.Properties['mcpServers'] -and $cfg.mcpServers.PSObject.Properties['gene-workbench']) {
      $cfg.mcpServers.PSObject.Properties.Remove('gene-workbench')
      $cfg | ConvertTo-Json -Depth 20 | Set-Content -Encoding UTF8 $mcpPath
    }
  }
}
$skillDst = Join-Path $wb 'skills\gene-workbench'
if (Test-Path $skillDst) { Remove-Item -Recurse -Force $skillDst }
