@echo off
call "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 exit /b 1
cl /nologo /O2 /std:c++20 /EHsc /Fo:packing_audit.obj /Fe:packing_audit.exe packing_audit.cpp
