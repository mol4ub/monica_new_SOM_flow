# SAMM calibration, DREAM run 7 (current best set-up)

Bayesian calibration of SAMM against all data of the Khon Kaen long-term residue experiment (1995–2019, five
treatments: CT control, DP dipterocarp, GN groundnut, RS rice straw, TM tamarind), with `samm_calib` and DREAM(ZS).

## Set-up

| Item | Setting |
|---|---|
| Rate modifiers | Kirschbaum temperature × MONICA moisture, no scaling factor (mu_max ≤ 0.95) |
| Numerics | RK4 with adaptive sub-steps (`stepAdaptive`) |
| Data | 855 observations: litterbag C, microbial N, SOC (including LMW), aggregate C, free MAOC, microbial C, soil N (without LAB_N) |
| Microbial observation model | MBN as reported = MIC_N × 3.1 / fN_ninhydrin (fitted, one factor for all years); MBC with k_EC 0.33 |
| Likelihood | `ar1`: Kalman filter per treatment × data type; measurement variance = median variance of the type; model error SD = sigmaRatio × measurement SD (0–3), correlation exp(−Δt/τ), τ 1–365 d |
| Microbial death | density dependent, k_MIC · (MIC/MIC_ref)^(β−1), β prior 1.5, ≤ 2 |
| Root input | simple weed model (`--weeds`) instead of constant daily litter: logistic growth limited by temperature, water and mineral N; weeding monthly Apr–Sep, 15 Nov/Jan/Mar and at litter days; above-ground biomass removed, roots to STR/LAB |
| Weed constraint | mean rainy-season above-ground weed biomass at weeding, GN 1.5 ± 0.75, CT 0.3 ± 0.15 t DM/ha (field photos) |
| N inputs | deposition 7, free-living fixation 2 kg N ha⁻¹ yr⁻¹ |
| Tillage | share aggDestroy of aggregates broken up at each hoeing / litter incorporation |
| Fit-driven options | `--flexN` (weed shoot N 0.8–2.5%), `--leachN` (leachable share of mineral N); consider dropping for definitive runs |
| Prior | `../../prior_kirschbaum_ar1_weeds_till.csv` (truncated Gaussian) |
| Sampler | DREAM(ZS), 4 chains, thin 10, seed 111, stop when all R̂ < 1.1 |

Command (from the repository root; `samm_calib` built with
`g++ -O3 -march=native -std=c++17 -pthread -I src/core tests/samm/calibration/samm_calib.cpp src/core/samm.cpp -o samm_calib`):

```bash
O="--lik ar1 --weeds --ndep 7 --nfix 2 --tillage --flexN --leachN"
./samm_calib tests/samm/data/climate_file.csv tests/samm/data/AllMeasurements.csv dream tests/samm/calibration/prior_kirschbaum_ar1_weeds_till.csv tests/samm/calibration/results/dream7/start.csv chains.csv $O --chains 4 --maxGen 150000 --thin 10 --seed 111
```

Converged after 68,000 generations (all R̂ < 1.10; about 40 min on 4 cores). Posterior from the second half of the
chains. Best log-likelihood −5098.3.

## Fit (best parameter set)

| Data type | n | RMSE (kg/ha) | EF | Bias (sim − obs) |
|---|---|---|---|---|
| Litterbag C | 28 | 732 | 0.62 | +438 |
| Microbial N (as reported) | 271 | 24.5 | 0.23 | −3.0 |
| SOC | 135 | 634 | 0.80 | −75 |
| Aggregate C | 30 | 880 | −2.36 | +688 |
| Free MAOC | 30 | 177 | 0.81 | −49 |
| Microbial C | 241 | 93 | 0.05 | −7.9 |
| Soil N | 120 | 116 | 0.25 | −2.8 |

Weeds (best set): mean rainy-season above-ground biomass at weeding CT 0.06, DP 0.23, GN 0.74, RS 0.26, TM 0.49 t DM/ha;
root C input CT 0.6, DP 1.9, GN 5.2, RS 1.9, TM 3.6 kg C ha⁻¹ d⁻¹.

## Notable parameters (posterior median, 95% interval)

| Parameter | Median | 95% interval | Comment |
|---|---|---|---|
| CUE_LAB | 0.083 | 0.059–0.114 | very low: litter C loss largely bypasses microbial growth |
| betaMIC | 1.94 | 1.77–2.00 | at the upper bound |
| fN_ninhydrin | 2.13 | 2.01–2.35 | below the 3.1 used for the data |
| MicMinCN / MicMaxCN | 5.0 / 6.2 | | narrow prior |
| kAgg | 0.0080 | 0.0043–0.0132 | |
| aggDestroy | 0.012 | 0.001–0.050 | tillage effect rejected by the data |
| leachShareN | 0.11 | 0.004–0.73 | poorly identified |

Full table: `posterior_summary.csv`.

## Open issues

- Aggregate C is overestimated by about 1000 kg C/ha in all residue treatments (bias makes up 61% of the MSE), while the control fits.
- Litterbag C is overestimated in late, strongly decomposed bags (groundnut in particular).
- Microbial C and N: treatment differences are reproduced, but variation between sampling dates within a treatment is not (regression slope about 0.3).
- Weeds in the control are lower than the photo-based constraint (0.06 vs 0.3 t/ha).

## Files

| File | Content |
|---|---|
| `best.csv` | highest-likelihood parameter set (input for `samm_calib eval`) |
| `map.csv` | highest-posterior parameter set |
| `posterior_summary.csv` | mean, SD, quantiles and R̂ per parameter |
| `posterior_thinned.csv` | every 7th posterior sample (1943 samples; input for `samm_calib band`) |
| `start.csv` | starting point of the chains |
| `obs_vs_sim.csv` | measurements with the simulated values of the best set |
| `evaluation.png` | simulated vs measured per data type, as Fig. 5 of Laub et al. (2024); `../../plot_eval.R` |
| `plot.png`, `microbes_C.png`, `microbes_N.png`, `weeds.png` | time series; `scripts/` (written for the scratch directory of the run, file names dream7_*) |
