@echo off
setlocal
cd /d "%~dp0"
set PYTHONUTF8=1
if "%~1"=="" (
    echo Cach dung: CHAY_CONSOLE.bat examples\02_bessel.json
    exit /b 2
)
if exist ".venv\Scripts\python.exe" (
    ".venv\Scripts\python.exe" exam_cli.py "%~1"
) else (
    py -3 exam_cli.py "%~1"
)
