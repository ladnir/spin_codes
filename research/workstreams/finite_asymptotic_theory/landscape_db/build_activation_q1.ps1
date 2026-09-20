$ErrorActionPreference = 'Stop'
$q1BuildDirectory = Join-Path $PSScriptRoot 'native_build'
New-Item -ItemType Directory -Force -Path $q1BuildDirectory | Out-Null
$q1CompilerSetup = 'C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat'
if (-not (Test-Path -LiteralPath $q1CompilerSetup)) { throw "Missing compiler setup: $q1CompilerSetup" }
Push-Location -LiteralPath $q1BuildDirectory
try {
    @(
        "@call `"$q1CompilerSetup`" >nul",
        '@if errorlevel 1 exit /b 1',
        '@cl /nologo /O2 /EHsc /std:c++17 /fp:precise /LD ..\activation_q1_kernel.cpp /Fe:activation_q1_kernel.dll /link /Brepro',
        '@exit /b %errorlevel%'
    ) | Set-Content -LiteralPath 'build_q1.cmd' -Encoding ASCII
    & cmd.exe /d /c build_q1.cmd
    if ($LASTEXITCODE -ne 0) { throw 'Q1 kernel compilation failed' }
} finally { Pop-Location }
