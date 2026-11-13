Khon Kaen litter experiment (Thailand) with SAMM in MONICA
=========================================================

Bare fallow on a sandy soil in Khon Kaen, NE Thailand (16.33 N, 102.82 E), 1995-03-01 to 2019-12-29, with
yearly additions of 10 t DM/ha of four litters (Laub et al. 2024, GMD 17:931):

  CT  control, no litter        GN  groundnut stover
  DP  dipterocarp leaves        RS  rice straw
                                TM  tamarind leaves

Two set-ups:
  run_khonkaen.bat             SAMM v1.1 with aggregates, MONICA's temperature and moisture rate modifiers (scaled
                               to SAMM's range), the "blue line" set-up of documentation/samm.md, section 7.3
  run_khonkaen_kirschbaum.bat  the recalibrated variant SAMMKirschbaum (Kirschbaum temperature x MONICA moisture
                               function, no scaling; density-dependent microbial death), with the weed roots of
                               its calibration as inputs and hoeing at the litter additions;
                               documentation/samm.md, section 7.5

How to run
----------
1. The release contains this folder as projects\KhonKaen-SAMM. If you copy it elsewhere, it should read
   ...\monica_<version>\projects\KhonKaen-SAMM
   (if you keep it elsewhere, set MONICA_HOME to the release folder before running).
2. Double-click run_khonkaen.bat (or run_khonkaen_kirschbaum.bat), or run it from a command prompt. All five
   treatments take well under a minute. Single treatments: run_khonkaen.bat GN RS
3. Results: output\out-<treatment>.csv (Kirschbaum: output\out-<treatment>-Kirschbaum.csv), one row per day;
   log files output\log-<treatment>.txt.
4. Optional, needs R: "Rscript plot_results.R" (Kirschbaum: "Rscript plot_results.R Kirschbaum") in this folder
   plots all treatments against the measurements (output\khonkaen_results[_Kirschbaum].png), prints RMSE/EF/bias
   and writes output\matched[_Kirschbaum].csv.

Reference results of the same set-up (Linux build): documentation/samm.md, section 7.3, and
documentation/images/samm_khonkaen_all_measurements.png.

Files
-----
sim-<treatment>.json   simulation settings and the list of output columns
site.json              soil (0-15 cm: 82.4% sand, 5% clay, i.e. measured silt + clay 0.176; bulk density 1.45),
                       15 cm layers, SAMM switched on
                       ("SOM_Model": "SAMM"); the "SAMM" block holds the start state, c_SORP (0.049), the
                       constant daily input of weed roots (2.28 kg C/ha/d, C:N 143) and the parameters
                       (mineral N coupling of scenario D switched off, as in the R code)
crop-<treatment>.json  the yearly litter additions as organic fertilisations (C, N, lignin and polyphenols of
                       Table 1 of the paper); no crop
climate.csv            daily weather from NASA POWER
measurements.csv       measurements of the experiment (SAMM paper data), by day of simulation (day 1 = 1995-03-01)

Variant SAMMKirschbaum (run_khonkaen_kirschbaum.bat):
sim-<treatment>-Kirschbaum.json, site-Kirschbaum.json
                       "variant": "SAMMKirschbaum" (parameters of the recalibration, DREAM run 7), start state
                       with the root chemistry of the calibration, no constant daily root input, aggregate
                       destruction by hoeing 0.25%, N deposition 9 kg N/ha/yr
crop-<treatment>-Kirschbaum.json
                       litter additions, hoeing (Tillage 0.15 m) at each litter addition, and the weed roots of
                       the calibration's weed model as organic matter on the 248 weeding days (0.6 kg C/ha/d in
                       the control to 5.5 in groundnut, C:N 21.5, 28% structural); MONICA has no weeds here
Expected (documentation/samm.md, section 7.5): SOC RMSE 849 kg C/ha, EF 0.63. The residue treatments end higher
than in the calibration because MONICA drains less water below 15 cm than the HYDRUS drivers of the calibration.

Output columns (layer 1 = 0-15 cm)
----------------------------------
SAMM_<pool>_1   SAMM pools in kg C m-3 or kg N m-3. Multiply by 1500 for kg/ha in 0-15 cm.
                STR, LAB = free litter (structural, labile); LMW = low molecular weight (dissolved) C/N;
                MIC = microbial biomass; MAO = free mineral-associated OM;
                AggSTR, AggLAB, AggMAO = the same inside aggregates.
SAMM_CO2_1      heterotrophic respiration (kg C m-3 d-1);  SAMM_LeachC  LMW-C moved below 15 cm (kg C/ha/d)
STemp_1, Mois_1, Sat_1, pF_1   soil temperature, water content, saturation and pF of the layer
SOC_1           SOC (kg C/kg soil) as passed on to the rest of MONICA
NH4_1, NO3_1    mineral N (kg N m-3);  WaterFlux_1, Precip  water flux and precipitation (mm/d)

Quantities compared with the measurements (as in the SAMM paper):
  SOC          = MIC_C + MAO_C + AggSTR_C + AggLAB_C + AggMAO_C  (Kirschbaum: + LMW_C, as in its calibration)
  Soil N       = LMW_N + MIC_N + MAO_N + AggLAB_N + AggMAO_N
  Microbial C  = MIC_C, Microbial N = MIC_N  (Kirschbaum: MIC_N * 3.1 / 2.004, i.e. as reported with the
                 fitted ninhydrin factor)
  Aggregate C  = AggSTR_C + AggLAB_C + AggMAO_C (measured only in 2019)
  Free MOAC    = MAO_C (measured only in 2019)

To use SAMM's own temperature and moisture functions instead, set "rateModifiers": "SAMM" in the "SAMM" block
of site.json (default here: "MONICA"). For this sandy soil they decompose much more slowly (see the documentation).
