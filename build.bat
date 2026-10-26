@echo off
rem ===========================================================================
rem Build script for CligenWeather.dll (64-bit Windows)
rem Supports MinGW-w64 (g++) or Microsoft Visual C++ (cl)
rem ===========================================================================

echo Building CligenWeather.dll for GoldSim...

where cl >nul 2>nul
if %ERRORLEVEL% equ 0 (
    echo Detected MSVC cl.exe compiler.
    cl /O2 /std:c++17 /EHsc /LD goldsim_external.cpp /Fe:CligenWeather.dll
    goto done
)

where g++ >nul 2>nul
if %ERRORLEVEL% equ 0 (
    echo Detected MinGW g++ compiler.
    g++ -O2 -std=c++17 -Wall -shared goldsim_external.cpp -o CligenWeather.dll -static -Wl,--add-stdcall-alias
    goto done
)

echo ERROR: Neither MSVC (cl) nor MinGW (g++) compiler was found in PATH.
echo Please run from a Visual Studio Developer Command Prompt or install MinGW-w64.
exit /b 1

:done
if %ERRORLEVEL% equ 0 (
    echo Successfully built CligenWeather.dll!
) else (
    echo Build failed with error %ERRORLEVEL%.
)
