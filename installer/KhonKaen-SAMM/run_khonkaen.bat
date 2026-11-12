@echo off
rem Khon Kaen litter experiment (Thailand) with SAMM in MONICA
rem Runs the five treatments CT, DP, GN, RS, TM (or only those given as arguments, e.g. run_khonkaen.bat GN RS).
rem Expects this folder inside the MONICA release, e.g. monica_3.7.4.samm\projects\KhonKaen-SAMM.
rem If it is somewhere else, set MONICA_HOME to the release folder (the one containing bin and monica-parameters).
setlocal
cd /d "%~dp0"
if not defined MONICA_HOME set "MONICA_HOME=%~dp0..\.."
set "MONICA_PARAMETERS=%MONICA_HOME%\monica-parameters"
set "MONICA_EXE=%MONICA_HOME%\bin\monica-run.exe"
if not exist "%MONICA_EXE%" (
  echo monica-run.exe not found at "%MONICA_EXE%".
  echo Copy this folder into the "projects" folder of the MONICA release or set MONICA_HOME.
  pause
  exit /b 1
)
if not exist output mkdir output
set "TREATMENTS=%*"
if "%TREATMENTS%"=="" set "TREATMENTS=CT DP GN RS TM"
for %%T in (%TREATMENTS%) do (
  echo Running treatment %%T ...
  "%MONICA_EXE%" -o "output\out-%%T.csv" "sim-%%T.json" > "output\log-%%T.txt" 2>&1
  if errorlevel 1 (echo   failed, see output\log-%%T.txt) else (echo   done: output\out-%%T.csv)
)
echo.
echo Finished. Optional plot against the measurements: Rscript plot_results.R
pause
