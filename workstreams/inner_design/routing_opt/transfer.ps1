param([ValidateSet('upload','download')][string]$Mode='upload',[switch]$Deployment)
$ErrorActionPreference='Stop'
$repo='C:/Users/peter/repo/permute_conv-github-bch'
$remote='/tmp/inner-mixer-p3fxxt'
$archive=Join-Path ([System.IO.Path]::GetTempPath()) ('routing-opt-'+[guid]::NewGuid().ToString()+'.tar.gz')
$psi=[System.Diagnostics.ProcessStartInfo]::new()
$psi.FileName='C:/Users/peter/OneDrive/tools/plink.exe'
$psi.UseShellExecute=$false
$psi.CreateNoWindow=$true
foreach($arg in @('-ssh','-batch','-P','9022','-l','prindal','-i','C:/Users/peter/OneDrive/keys/key_cat-9-3-2025.ppk','-hostkey','ssh-ed25519 255 SHA256:6ZRPOoZl4NGXghNQToErE0anQ0PkFggvYu7eBjXbK9o','peach48.devcore4.com')) {$psi.ArgumentList.Add($arg)}
if($Mode -eq 'upload') {
    $dirs=@(Get-ChildItem -LiteralPath "$repo/workstreams/inner_design/generated" -Directory | Where-Object Name -Like 'routeopt_*' | ForEach-Object {"workstreams/inner_design/generated/$($_.Name)"})
    if($Deployment) {$dirs+=@('workstreams/bare_bch_rm2sub/Spin.cpp','workstreams/bare_bch_rm2sub/WorkspaceRouting.h','workstreams/bare_bch_rm2sub/CMakeLists.txt')}
    & tar -czf $archive -C $repo workstreams/inner_design/routing_opt @dirs
    if($LASTEXITCODE) {throw 'tar failed'}
    $psi.ArgumentList.Add("tar xzf - -C $remote")
    $psi.RedirectStandardInput=$true
    $process=[System.Diagnostics.Process]::Start($psi)
    $stream=[System.IO.File]::OpenRead($archive)
    try {$stream.CopyTo($process.StandardInput.BaseStream)} finally {$stream.Dispose();$process.StandardInput.Close()}
} else {
    $psi.ArgumentList.Add("tar czf - -C $remote routing-opt")
    $psi.RedirectStandardOutput=$true
    $process=[System.Diagnostics.Process]::Start($psi)
    $stream=[System.IO.File]::Create($archive)
    try {$process.StandardOutput.BaseStream.CopyTo($stream)} finally {$stream.Dispose()}
}
$process.WaitForExit()
if($process.ExitCode) {throw "transfer failed $($process.ExitCode)"}
if($Mode -eq 'download') {
    $entries=& tar -tzf $archive
    if($entries | Where-Object {$_ -notmatch '^routing-opt/' -or $_ -match '\.\.'}) {throw 'Unexpected archive path'}
    $destination="$repo/workstreams/inner_design/routing_opt/measurements"
    New-Item -ItemType Directory -Force -Path $destination | Out-Null
    & tar -xzf $archive -C $destination --strip-components=1
    if($LASTEXITCODE) {throw 'extract failed'}
}
Remove-Item -LiteralPath $archive
