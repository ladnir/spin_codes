$ErrorActionPreference = 'Stop'
$refreshBuildDirectory = Join-Path $PSScriptRoot 'native_build'
New-Item -ItemType Directory -Force -Path $refreshBuildDirectory | Out-Null
$refreshCompilerSetup = 'C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat'
if (-not (Test-Path -LiteralPath $refreshCompilerSetup)) { throw "Missing compiler setup: $refreshCompilerSetup" }
Push-Location -LiteralPath $refreshBuildDirectory
try {
    @(
        "@call `"$refreshCompilerSetup`" >nul",
        '@if errorlevel 1 exit /b 1',
        '@cl /nologo /O2 /EHsc /std:c++17 /fp:precise /LD ..\activation_refresh_kernel.cpp /Fe:activation_refresh_kernel.dll /link /Brepro',
        '@exit /b %errorlevel%'
    ) | Set-Content -LiteralPath 'build_refresh.cmd' -Encoding ASCII
    & cmd.exe /d /c build_refresh.cmd
    if ($LASTEXITCODE -ne 0) { throw 'Uniform-refresh kernel compilation failed' }
} finally { Pop-Location }
