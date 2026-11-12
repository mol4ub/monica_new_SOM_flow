@echo off
rem Genetic algorithm calibration with the objective function of the SAMM paper (see README.txt).
rem   calibrate.bat [modifiers] [start parameter file]
rem Writes output\ga_<modifiers>_best.csv (best parameter set, updated every 10 generations) and _log.csv.
rem Takes about 20-60 minutes; uses all processor cores.
setlocal
cd /d "%~dp0"
set "MOD=%~1"
if "%MOD%"=="" set "MOD=KirschbaumMONICA"
set "START=%~2"
if "%START%"=="" set "START=params\SAMM_v1.1.csv"
if not exist output mkdir output
samm.exe climate_file.csv AllMeasurements.csv ga "output\ga_%MOD%" "%START%" --modifiers %MOD% --pop 200 --gens 800
pause
