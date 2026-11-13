# Leave-one-residue-out cross-validation of the run 7 set-up

Four DREAM calibrations with the set-up of run 7 (`../dream7/README.md`). In each one a residue treatment (DP
dipterocarp, GN groundnut, RS rice straw or TM tamarind) is left out with `samm_calib --holdout <TRT>`: its
measurements and its weed constraint do not enter the likelihood. The control is always calibrated. Each fold then
predicts the treatment it left out. Script: `../../run_cv.sh`; merge and statistics: `../../cv_merge.py`.

Sampler: DREAM(ZS), 4 chains, seed 211, same start as run 7, stop when all R̂ < 1.5. Generations to convergence:
DP 11,000, GN 11,000, RS 13,000, TM 21,000 (5–11 min each on 4 cores).

## Held-out vs in-sample (the four residue treatments pooled)

In-sample = run 7 (all treatments calibrated, R̂ < 1.1), evaluated on the same four treatments.

| Data type | n | RMSE in-sample | RMSE held-out | EF in-sample | EF held-out | Bias held-out |
|---|---|---|---|---|---|---|
| Litterbag C | 28 | 732 | 795 | 0.62 | 0.55 | +177 |
| Microbial N (as reported) | 219 | 24.8 | 26.2 | 0.21 | 0.12 | −0.3 |
| SOC | 108 | 657 | 1039 | 0.69 | 0.22 | +60 |
| Aggregate C | 24 | 980 | 861 | −6.83 | −5.04 | +687 |
| Free MAOC | 24 | 152 | 311 | 0.81 | 0.20 | +213 |
| Microbial C | 195 | 94 | 102 | 0.02 | −0.15 | +18.5 |
| Soil N | 96 | 121 | 136 | 0.09 | −0.16 | +16 |

Per treatment: `heldout_vs_insample.txt`. Plot: `evaluation_heldout.png`.

## Findings

- **Litterbags, microbial C and N, soil N** are predicted for an unseen residue almost as well as in-sample.
- **SOC and free MAOC do not transfer between residues.** The pooled SOC EF drops from 0.69 to 0.22:
  - Without DP, SOC under DP is overestimated by 1.3 t C/ha (about 25%) and free MAOC by 0.5 t C/ha.
  - Without TM, SOC under TM is underestimated by 0.9 t C/ha (in-sample already −0.66 t C/ha).
  - GN and RS are predicted reasonably well (SOC bias −0.39 and +0.21 t C/ha).
- **DP constrains the structural / aggregate pathway.** Without DP the posterior moves to:
  - kAgg 0.030 (all data: 0.008);
  - NonMicAgg 4.6 (1.5);
  - aggfactSTR_C 2.5 (1.4);
  - k_STR 0.0145 (0.0098).

  DP is the most lignin-rich residue. The effect of litter quality on SOC is therefore carried by fitted aggregate
  and STR parameters rather than following from the litter chemistry, so it does not extrapolate to a residue of
  different quality.
- **Aggregate C is no worse held-out than in-sample** (EF −5.0 vs −6.8). The in-sample misfit is structural, not overfitting.

## Caveats

- R̂ < 1.5 and short chains (posterior from the second half, 5,500–10,500 generations): posterior intervals are
  probably too narrow; the best sets are more reliable than the bands.
- The chains start from the run 7 starting point, which came from calibrations on all data. This favours the
  held-out predictions somewhat.
- leachShareN is 0.36–0.48 in all folds, against 0.11 in run 7. This parameter is poorly identified, and it may also
  indicate a second mode.

## Files

| File | Content |
|---|---|
| `best_without_<TRT>.csv` | highest-likelihood parameter set of the fold without TRT |
| `posterior_summary_without_<TRT>.csv` | posterior quantiles and R̂ of that fold |
| `heldout_obs_vs_sim.csv` | measurements of each residue treatment with the prediction of the fold without it |
| `heldout_vs_insample.txt` | RMSE, EF and bias per data type and treatment, in-sample vs held-out |
| `evaluation_heldout.png` | simulated (held-out) vs measured, as Fig. 5 of Laub et al. (2024) |
