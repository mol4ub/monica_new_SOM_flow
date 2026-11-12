SAMM standalone - Khon Kaen litter experiment (all five treatments)
====================================================================

A compiled version of the SAMM C++ core used in MONICA (src/core/samm.cpp), set up exactly like the SAMM v1.1
R code (Samm_main.R): one 0-15 cm layer, HYDRUS soil drivers (climate_file.csv), the yearly litter additions,
the daily root inputs, litterbag pools, classical RK4 with daily steps. With params\SAMM_v1.1.csv and the SAMM
modifiers it reproduces Samm_main.R to about 1e-11 (relative), and runs all five treatments in about 0.1 s.

Difference to Samm_main.R v1.1: the litterbag of the litterbag experiment (added on simulation day 4414) is followed
past the next litter addition until its last measurement (day 4799); the R code empties and refills the litterbag pools 13 days before every addition, so its
last litterbag value is compared with the next year's new bag. (Litterbag EF of v1.1: 0.69 instead of 0.12.)

Run (Windows)
-------------
  run_samm.bat                                          SAMM v1.1 parameters, SAMM modifiers
  run_samm.bat params\SAMM_v1.1.csv KirschbaumMONICA    any parameter file, Kirschbaum x MONICA moisture
Linux: ./samm_linux climate_file.csv AllMeasurements.csv eval params/SAMM_v1.1.csv out_daily.csv out_obs.csv --modifiers SAMM

The program prints the log-likelihood of the run (objective function below).
Optional plot (needs R): Rscript plot_results.R output\daily_<name>.csv output\obs_<name>.csv output\<name>.png

Files
-----
samm.exe / samm_linux    the model (Windows 64 bit, static, no DLLs needed / Linux x86-64)
climate_file.csv         daily drivers from the paper: soil temperature at 7.5 cm, water content (theta),
                         leaching [mm], pF (all from HYDRUS 1D)
AllMeasurements.csv      all measurements (by "Day of simulation", day 1 = 1995-03-01)
params\*.csv             parameter files: name,value (the 21 calibrated parameters of Table 4 of the paper)
run_samm.bat             run one parameter set;  calibrate.bat  genetic algorithm calibration
plot_results.R           plot of a run against the measurements
source\                  source code (MPL 2.0); build: g++ -O3 -std=c++17 -pthread samm_calib.cpp samm.cpp -o samm

Rate modifiers (--modifiers)
----------------------------
SAMM              (paper) 2^((T-20)/10) x pF scalar (1 between pF 1.5 and 2.5)
KirschbaumMONICA  Kirschbaum (1995): exp(-3.764 + 0.204 T (1 - T/73.8)), maximum 1 at 36.9 C,
                  x MONICA's moisture function exp(-18 (theta/0.468 - 0.45)^2), maximum 1
                  (0.468 = MONICA's saturated water content of this soil; options --thetaSat, --moistOpt)
                  Rates calibrated with the SAMM modifiers do NOT fit with this one (corr_F about 5x smaller);
                  it needs its own parameter set (calibration in progress).

Outputs
-------
output\daily_<name>.csv  time (day of simulation; row t = state at the start of day t, as the R output), the 14
                         SAMM pools (kg C or N per ha, 0-15 cm), LitterBagSTR_C, LitterBagLAB_C, LitterBagLAB_N, TRT
output\obs_<name>.csv    TRT, day, type, measured, simulated for every measurement used in the objective
Quantities: SOC = MIC_C + MAO_C + AggSTR_C + AggLAB_C + AggMAO_C; aggregate C = AggSTR_C + AggLAB_C + AggMAO_C;
free MAOC = MAO_C; microbial N = MIC_N; litterbag C = LitterBagSTR_C + LitterBagLAB_C.

Objective function (as Laub et al. 2024, GMD 17:931, section 2.4)
------------------------------------------------------------------
Gaussian log-likelihood of litterbag C (n = 28), microbial N (271), SOC (135), aggregate C (30) and free MAOC (30),
no covariances, variance of each data type = median of its "Variance" column:
  logLik = sum over observations of  -0.5 ln(2 pi var_type) - 0.5 (obs - sim)^2 / var_type
SAMM v1.1 with the SAMM modifiers: logLik -4277.7.

Calibration (calibrate.bat)
---------------------------
Real-valued genetic algorithm over the 21 parameters: population 200, up to 800 generations (stops after 60 without
improvement), linear rank selection, 5% elitism, BLX-0.5 crossover (p 0.8), Gaussian mutation (p 0.1, step shrinking
from 25% to 2% of the parameter range). Search bounds are in source\samm_calib.cpp (pdefs); mu_max <= 0.95.
