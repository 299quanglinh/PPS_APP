@echo off
setlocal
cd /d "%~dp0"
set PYTHONUTF8=1
echo Dang khoi dong PPS Python, vui long cho...
if exist ".venv\Scripts\python.exe" (
    ".venv\Scripts\python.exe" main.py
) else (
    py -3 main.py
)
if errorlevel 1 (
    echo.
    echo Khong khoi dong duoc. Chay CAI_DAT.bat de cai thu vien.
    pause
)
