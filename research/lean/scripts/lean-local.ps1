param([Parameter(ValueFromRemainingArguments=$true)][string[]]$LeanArgs)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path $PSScriptRoot -Parent
$searchPaths = @((Join-Path $projectRoot '.lake/build/lib/lean'))
$searchPaths += Get-ChildItem (Join-Path $projectRoot '.lake/packages') -Directory | ForEach-Object {
    Join-Path $_.FullName '.lake/build/lib/lean'
}
$env:LEAN_PATH = $searchPaths -join ';'
Push-Location $projectRoot
try { & lean @LeanArgs; $leanExit = $LASTEXITCODE } finally { Pop-Location }
exit $leanExit
