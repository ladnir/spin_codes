param(
    [Parameter(Mandatory = $true)][string]$Manifest,
    [string]$DependencyLeanRoot,
    [int]$MaxParts = 0,
    [int]$MemoryMB = 8000,
    [switch]$BuildAssemblies,
    [switch]$RefreshCheckerDependencies
)

$ErrorActionPreference = 'Stop'
$leanRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '../../..')).Path
function Resolve-ManifestSource([string]$Source) {
    if ([System.IO.Path]::IsPathRooted($Source)) { return $Source }
    return Join-Path $repoRoot $Source
}

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

# Refresh explicitly: shared-dependency worktrees may not have local Lake
# packages, and merely resuming a completed cover must not provision them.
if ($RefreshCheckerDependencies) {
    Push-Location -LiteralPath $leanRoot
    try {
        & lake build SpinCodes.UniformMixer.ScalableDensePackedLeaf SpinCodes.UniformMixer.ScalableDenseCoverCore
        if ($LASTEXITCODE -ne 0) { throw 'Cannot build dense checker dependencies' }
    } finally { Pop-Location }
}
$baseDependencies = @(
    $baseOlean,
    (Join-Path $buildLib 'SpinCodes/UniformMixer/ScalableDenseCoverCore.olean')
)
foreach ($dependency in $baseDependencies) {
    if (-not (Test-Path -LiteralPath $dependency)) {
        throw "Missing checker object: $dependency. Build checker dependencies or pass -RefreshCheckerDependencies with local Lake dependencies ready."
    }
}

function Invoke-LeanSource([string]$Name, [string]$Source, [int]$Limit, [string[]]$Dependencies = @()) {
    $output = Join-Path $outDir "$Name.olean"
    $sourceInfo = Get-Item -LiteralPath $Source
    $newestInput = $sourceInfo.LastWriteTimeUtc
    foreach ($dependency in ($baseDependencies + $Dependencies)) {
        $dependencyInfo = Get-Item -LiteralPath $dependency
        if ($dependencyInfo.LastWriteTimeUtc -gt $newestInput) {
            $newestInput = $dependencyInfo.LastWriteTimeUtc
        }
    }
    if (Test-Path -LiteralPath $output) {
        $outputInfo = Get-Item -LiteralPath $output
        if ($outputInfo.LastWriteTimeUtc -ge $newestInput) {
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
    $source = Resolve-ManifestSource $part.source
    if (-not (Test-Path -LiteralPath $source)) { throw "Missing source: $source" }
    $status = Invoke-LeanSource $name $source $MemoryMB
    if ($status -eq 'compiled') {
        $built++
        Write-Output "$name checked ($built new this run)"
        if ($MaxParts -gt 0 -and $built -ge $MaxParts) {
            if ($BuildAssemblies) { throw 'Part build limit reached; resume without -MaxParts before assembling' }
            break
        }
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
        $source = Resolve-ManifestSource $data.assembly_sources[$n]
        # Conservatively depend on every part; segments import subsets.
        $partDependencies = @($data.parts | ForEach-Object {
            Join-Path $outDir "Part$($_.id).olean"
        })
        $status = Invoke-LeanSource $name $source $MemoryMB $partDependencies
        Write-Output "$name $status"
    }
    $source = Resolve-ManifestSource $data.global_source
    $segmentDependencies = @(0..($data.segments - 1) | ForEach-Object {
        Join-Path $outDir "Segment$_.olean"
    })
    $status = Invoke-LeanSource 'ScalableDenseFullCover' $source $MemoryMB $segmentDependencies
    Write-Output "ScalableDenseFullCover $status"
}
