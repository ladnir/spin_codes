param([ValidateSet('upload','download')][string]$Mode='upload')
$ErrorActionPreference='Stop'
$remote='/tmp/hd-lambda-20261009-JkBHbA'
$archive=Join-Path ([System.IO.Path]::GetTempPath()) ('hd-lambda-'+[guid]::NewGuid().ToString()+'.tar.gz')
$psi=[System.Diagnostics.ProcessStartInfo]::new()
$psi.FileName='C:/Users/peter/OneDrive/tools/plink.exe'
$psi.UseShellExecute=$false
$psi.CreateNoWindow=$true
foreach($arg in @('-ssh','-batch','-P','9022','-l','prindal','-i','C:/Users/peter/OneDrive/keys/key_cat-9-3-2025.ppk','-hostkey','ssh-ed25519 255 SHA256:6ZRPOoZl4NGXghNQToErE0anQ0PkFggvYu7eBjXbK9o','peach48.devcore4.com')) {$psi.ArgumentList.Add($arg)}
if($Mode -eq 'upload') {
    & tar -czf $archive -C $PSScriptRoot HdGroupedOuter.h hd_grouped_forward.cpp build.sh bench.sh
    if($LASTEXITCODE) {throw 'tar failed'}
    $psi.ArgumentList.Add("tar xzf - -C $remote")
    $psi.RedirectStandardInput=$true
    $process=[System.Diagnostics.Process]::Start($psi)
    $stream=[System.IO.File]::OpenRead($archive)
    try {$stream.CopyTo($process.StandardInput.BaseStream)} finally {$stream.Dispose();$process.StandardInput.Close()}
} else {
    $psi.ArgumentList.Add("tar czf - -C $remote build.log checks.log matched.log hashes.log natural-build.log natural-checks.log sanitize-build.log sanitize-checks.log")
    $psi.RedirectStandardOutput=$true
    $process=[System.Diagnostics.Process]::Start($psi)
    $stream=[System.IO.File]::Create($archive)
    try {$process.StandardOutput.BaseStream.CopyTo($stream)} finally {$stream.Dispose()}
}
$process.WaitForExit()
if($process.ExitCode) {throw "transfer failed $($process.ExitCode)"}
if($Mode -eq 'download') {
    $entries=& tar -tzf $archive
    if($entries | Where-Object {$_ -notmatch '^[a-z-]+\.log$'}) {throw 'Unexpected archive path'}
    & tar -xzf $archive -C $PSScriptRoot
    if($LASTEXITCODE) {throw 'extract failed'}
}
# Remove only this fresh, explicitly named temporary archive.
Remove-Item -LiteralPath $archive
