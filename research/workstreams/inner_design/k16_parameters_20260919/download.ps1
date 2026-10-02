param([string]$Destination = "$PSScriptRoot/measurements/remote.tar.gz")
# Verified-host streaming transfer; the installed WinSCP cannot read the v3 key.
$taskTransfer = New-Object System.Diagnostics.ProcessStartInfo
$taskTransfer.FileName = 'C:/Users/peter/OneDrive/tools/plink.exe'
$taskTransfer.UseShellExecute = $false
$taskTransfer.RedirectStandardOutput = $true
$taskTransfer.CreateNoWindow = $true
$taskTransfer.Arguments = '-ssh -batch -P 9022 -l prindal -i "C:/Users/peter/OneDrive/keys/key_cat-9-3-2025.ppk" -hostkey "ssh-ed25519 255 SHA256:6ZRPOoZl4NGXghNQToErE0anQ0PkFggvYu7eBjXbK9o" peach48.devcore4.com "tar -cz -C /tmp/spin-scheduling-20260919/workstreams/inner_design/k16_parameters_20260919 measurements"'
$taskProcess = [System.Diagnostics.Process]::Start($taskTransfer)
$taskArchiveStream = [System.IO.File]::Create($Destination)
try { $taskProcess.StandardOutput.BaseStream.CopyTo($taskArchiveStream) }
finally { $taskArchiveStream.Dispose() }
$taskProcess.WaitForExit()
if ($taskProcess.ExitCode -ne 0) { throw "Transfer failed: $($taskProcess.ExitCode)" }
