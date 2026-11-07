param(
    [Parameter(Mandatory = $true)][string]$Manifest,
    [string]$DependencyLeanRoot,
    [int]$MaxParts = 0,
    [int]$MemoryMB = 8000,
    [switch]$BuildAssemblies
)

$ErrorActionPreference = 'Stop'
$leanRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../../..')).Path
$manifestPath = (Resolve-Path $Manifest).Path
$data = Get-Content -LiteralPath $manifestPath -Raw | ConvertFrom-Json
if ($data.status -ne 'source_emitted') {
    throw "The manifest does not record a full source emission: $($data.status)"
}
if (-not $DependencyLeanRoot) { $DependencyLeanRoot = $leanRoot }
$depRoot = (Resolve-Path $DependencyLeanRoot).Path
$depPath = (& powershell -NoProfile -Command "Set-Location -LiteralPath '$depRoot'; lake env powershell -NoProfile -Command 'Write-Output `$env:LEAN_PATH'")
if ($LASTEXITCODE -ne 0 -or -not $depPath) { throw 'Cannot get Lean dependency path' }
$buildLib = Join-Path $leanRoot '.lake/build/lib/lean'
$env:LEAN_PATH = "$buildLib;$depPath"
$outDir = Join-Path $buildLib 'SpinCodes/UniformMixer/Generated/FullCover'
$logDir = Join-Path (Split-Path -Parent $manifestPath) 'logs'
New-Item -ItemType Directory -Path $outDir,$logDir -Force | Out-Null
$baseOlean = Join-Path $buildLib 'SpinCodes/UniformMixer/ScalableDensePackedLeaf.olean'
if (-not (Test-Path -LiteralPath $baseOlean)) { throw "Missing base checker: $baseOlean" }

function Invoke-LeanSource([string]$Name, [string]$Source, [int]$Limit) {
    $output = Join-Path $outDir "$Name.olean"
    $sourceInfo = Get-Item -LiteralPath $Source
    $baseInfo = Get-Item -LiteralPath $baseOlean
    if (Test-Path -LiteralPath $output) {
        $outputInfo = Get-Item -LiteralPath $output
        if ($outputInfo.LastWriteTimeUtc -ge $sourceInfo.LastWriteTimeUtc -and
            $outputInfo.LastWriteTimeUtc -ge $baseInfo.LastWriteTimeUtc) {
            return 'skipped'
        }
    }
    $temporary = "$output.tmp"
    $log = Join-Path $logDir "$Name.log"
    $result = & lean -j1 -M $Limit -T 10000 -R $leanRoot -o $temporary $Source 2>&1
    $code = $LASTEXITCODE
    $text = ($result | Out-String)
    Set-Content -LiteralPath $log -Value $text -Encoding utf8
    if ($code -ne 0 -or $text.Contains('sorryAx')) {
        if (Test-Path -LiteralPath $temporary) {
            Remove-Item -LiteralPath $temporary
        }
        throw "Lean failed for $Name (exit $code); inspect $log"
    }
    Move-Item -LiteralPath $temporary -Destination $output -Force
    return 'compiled'
}

$built = 0
foreach ($part in $data.parts) {
    $name = "Part$($part.id)"
    $source = Join-Path $repoRoot $part.source
    if (-not (Test-Path -LiteralPath $source)) { throw "Missing source: $source" }
    $status = Invoke-LeanSource $name $source $MemoryMB
    if ($status -eq 'compiled') {
        $built++
        Write-Output "$name checked ($built new this run)"
        if ($MaxParts -gt 0 -and $built -ge $MaxParts) { break }
    }
}

if ($BuildAssemblies) {
    foreach ($part in $data.parts) {
        $olean = Join-Path $outDir "Part$($part.id).olean"
        if (-not (Test-Path -LiteralPath $olean)) {
            throw "Assembly requested before all parts are checked: $olean"
        }
    }
    foreach ($n in 0..($data.segments - 1)) {
        $name = "Segment$n"
        $source = Join-Path $repoRoot $data.assembly_sources[$n]
        $status = Invoke-LeanSource $name $source $MemoryMB
        Write-Output "$name $status"
    }
    $source = Join-Path $repoRoot $data.global_source
    $status = Invoke-LeanSource 'ScalableDenseFullCover' $source $MemoryMB
    Write-Output "ScalableDenseFullCover $status"
}
