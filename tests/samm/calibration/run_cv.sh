#!/bin/bash
# leave-one-residue-out cross-validation of the run 7 set-up; DREAM until all R-hat < 1.5
# run from the repository root; samm_calib binary in $S, scratch directory $C (needs $C/dream7_* of run 7 for cv_merge.py)
S=${S:-.}; C=${C:-cv_scratch}; mkdir -p $C/cv
O="--lik ar1 --weeds --ndep 7 --nfix 2 --tillage --flexN --leachN"
D="tests/samm/data/climate_file.csv tests/samm/data/AllMeasurements.csv"
for T in DP GN RS TM; do
  P=$C/cv/lo_$T
  date +%T > $P.start
  $S/samm_calib $D dream tests/samm/calibration/prior_kirschbaum_ar1_weeds_till.csv tests/samm/calibration/results/dream7/start.csv ${P}_chains.csv $O --holdout $T --chains 4 --maxGen 150000 --thin 10 --rhat 1.5 --seed 211 > $P.out 2> $P.err
  # continued to R-hat < 1.1 with the same command plus --resume and --rhat 1.1
  python3 -I tests/samm/calibration/dream_summary.py ${P}_chains.csv $P > ${P}_summary.txt
  $S/samm_calib $D eval ${P}_best.csv ${P}_daily.csv ${P}_obs.csv $O > ${P}_eval.txt 2>&1
  $S/samm_calib $D band ${P}_posterior.csv ${P}_band.csv --n 200 --seed 3 $O 2>/dev/null
  date +%T > $P.end
done
