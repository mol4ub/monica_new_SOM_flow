@echo off
rem Runs SAMM for all five Khon Kaen treatments (CT, DP, GN, RS, TM).
rem   run_samm.bat [parameter file] [modifiers]
rem   parameter file: default params\SAMM_v1.1.csv
rem   modifiers: SAMM (paper: 2^((T-20)/10) x pF scalar, default) or KirschbaumMONICA
rem              (Kirschbaum temperature x MONICA moisture function, both at most 1)
rem Results: output\daily_<name>.csv (all pools, every day, all treatments) and
rem          output\obs_<name>.csv (every measurement with the simulated value)
setlocal
cd /d "%~dp0"
set "PARAMS=%~1"
if "%PARAMS%"=="" set "PARAMS=params\SAMM_v1.1.csv"
set "MOD=%~2"
if "%MOD%"=="" set "MOD=SAMM"
for %%F in ("%PARAMS%") do set "NAME=%%~nF_%MOD%"
if not exist output mkdir output
samm.exe climate_file.csv AllMeasurements.csv eval "%PARAMS%" "output\daily_%NAME%.csv" "output\obs_%NAME%.csv" --modifiers %MOD%
if errorlevel 1 (echo SAMM failed) else (echo Written: output\daily_%NAME%.csv and output\obs_%NAME%.csv)
echo Optional plot: Rscript plot_results.R output\daily_%NAME%.csv output\obs_%NAME%.csv output\%NAME%.png
pause
