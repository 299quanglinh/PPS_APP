@echo off
setlocal
cd /d "%~dp0"
set PYTHONUTF8=1
py -3 -m venv .venv
if errorlevel 1 goto fail
".venv\Scripts\python.exe" -m pip install -r requirements.txt
if errorlevel 1 goto fail
echo Cai dat thanh cong. Mo CHAY_PPS.bat de chay.
pause
exit /b 0
:fail
echo Cai dat that bai. Can Python 3.11 tro len va ket noi Internet.
pause
exit /b 1
