param([string]$Archive = "$PSScriptRoot/source.tar.gz")
$taskTransfer = New-Object System.Diagnostics.ProcessStartInfo
$taskTransfer.FileName = 'C:/Users/peter/OneDrive/tools/plink.exe'
$taskTransfer.UseShellExecute = $false
$taskTransfer.RedirectStandardInput = $true
$taskTransfer.CreateNoWindow = $true
$taskTransfer.Arguments = '-ssh -batch -P 9022 -l prindal -i "C:/Users/peter/OneDrive/keys/key_cat-9-3-2025.ppk" -hostkey "ssh-ed25519 255 SHA256:6ZRPOoZl4NGXghNQToErE0anQ0PkFggvYu7eBjXbK9o" peach48.devcore4.com "tar -xz -C /tmp/spin-scheduling-20260919"'
$taskProcess = [System.Diagnostics.Process]::Start($taskTransfer)
$taskArchiveStream = [System.IO.File]::OpenRead((Resolve-Path $Archive))
try { $taskArchiveStream.CopyTo($taskProcess.StandardInput.BaseStream) }
finally { $taskArchiveStream.Dispose(); $taskProcess.StandardInput.Close() }
$taskProcess.WaitForExit()
if ($taskProcess.ExitCode -ne 0) { throw "Transfer failed: $($taskProcess.ExitCode)" }
