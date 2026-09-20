$ErrorActionPreference = 'Stop'
$q2BuildDirectory = Join-Path $PSScriptRoot 'native_build'
New-Item -ItemType Directory -Force -Path $q2BuildDirectory | Out-Null
$q2CompilerSetup = 'C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat'
if (-not (Test-Path -LiteralPath $q2CompilerSetup)) { throw "Missing compiler setup: $q2CompilerSetup" }
Push-Location -LiteralPath $q2BuildDirectory
try {
    @(
        "@call `"$q2CompilerSetup`" >nul",
        '@if errorlevel 1 exit /b 1',
        '@cl /nologo /O2 /EHsc /std:c++17 /fp:precise /LD ..\activation_q2_kernel.cpp /Fe:activation_q2_kernel.dll /link /Brepro',
        '@exit /b %errorlevel%'
    ) | Set-Content -LiteralPath 'build_q2.cmd' -Encoding ASCII
    & cmd.exe /d /c build_q2.cmd
    if ($LASTEXITCODE -ne 0) { throw 'Q2 kernel compilation failed' }
} finally { Pop-Location }
