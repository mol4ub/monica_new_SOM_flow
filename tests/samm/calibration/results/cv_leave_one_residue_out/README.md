# Leave-one-residue-out cross-validation of the run 7 set-up

Four DREAM calibrations with the set-up of run 7 (`../dream7/README.md`). In each one a residue treatment (DP
dipterocarp, GN groundnut, RS rice straw or TM tamarind) is left out with `samm_calib --holdout <TRT>`: its
measurements and its weed constraint do not enter the likelihood. The control is always calibrated. Each fold then
predicts the treatment it left out. Script: `../../run_cv.sh`; merge and statistics: `../../cv_merge.py`.

## Sampling

- DREAM(ZS), 4 chains, seed 211, same start as run 7.
- First run to R̂ < 1.5 (DP and GN 11,000 generations, RS 13,000, TM 21,000), then continued with `--resume` to R̂ < 1.1:

| Fold | Generations | Max R̂ | Note |
|---|---|---|---|
| without DP | 49,000 | 1.099 | |
| without GN | 54,000 | 1.099 | |
| without RS | 93,680 | 2.0 | **not converged**, stopped: CUE_STR bimodal (see below) |
| without TM | 77,000 | 1.094 | |

Posterior from the second half of each chain. Without RS, one chain stays in a mode with CUE_STR ≈ 0, i.e.
structural litter respired without microbial growth, which fits better (log-likelihood about −4084, against about
−4103 in the mode with CUE_STR around 0.7). The best set of this fold lies in that mode (CUE_STR 0.0015). The prior
allows CUE_STR down to 0; a physically motivated lower bound (e.g. 0.1–0.2) would rule it out.

## Held-out vs in-sample (the four residue treatments pooled, Nash–Sutcliffe EF)

- In-sample: run 7 (all treatments calibrated, R̂ < 1.1).
- Held-out: each treatment predicted by the fold without it.
- "best": highest-likelihood parameter set; "median": posterior median prediction (the proper cross-validation prediction; used in `evaluation_heldout.png`).

| Data type | In-sample best | In-sample median | Held-out R̂<1.5 best | Held-out R̂<1.5 median | **Held-out R̂<1.1 best** | **Held-out R̂<1.1 median** |
|---|---|---|---|---|---|---|
| Litterbag C | 0.62 | 0.71 | 0.55 | 0.57 | 0.36 | **0.56** |
| Microbial N | 0.21 | 0.21 | 0.12 | 0.13 | 0.08 | **0.11** |
| SOC | 0.69 | 0.71 | 0.22 | 0.34 | −0.05 | **0.09** |
| Aggregate C | −6.83 | −4.77 | −5.04 | −3.68 | −12.54 | **−6.25** |
| Free MAOC | 0.81 | 0.77 | 0.20 | 0.55 | −2.58 | **−0.55** |
| Microbial C | 0.02 | 0.09 | −0.15 | −0.01 | −0.37 | **−0.18** |
| Soil N | 0.09 | −0.02 | −0.16 | −0.07 | −0.24 | **−0.11** |

SOC bias of the held-out posterior median (R̂ < 1.1), kg C/ha:

| Treatment | In-sample | Held-out |
|---|---|---|
| DP | +37 | +1306 |
| GN | −159 | −486 |
| RS | +136 | +703 |
| TM | −654 | −1025 |

Per data type and treatment (best set): `heldout_vs_insample_best_rhat1.1.txt`, `heldout_vs_insample_best_rhat1.5.txt`.

## Findings

- **Litterbags and microbial N are predicted for an unseen residue almost as well as in-sample** (litterbag EF 0.56 vs 0.71).
- **SOC, free MAOC and aggregate C do not transfer between residues:**
  - The SOC EF drops from 0.71 to 0.09.
  - Without its data, SOC is overestimated by 1.3 t C/ha under DP (about 25%) and by 0.7 t C/ha under RS.
  - It is underestimated by 1.0 t C/ha under TM.
- **The R̂ < 1.5 results were too optimistic.** With longer sampling the chains move further away from the start point
  (which came from calibrations on all data), and the held-out fit becomes worse. The R̂ < 1.1 results are the
  honest estimate.
- **The posterior median predicts better than the best set**, in particular for the held-out treatments (free MAOC EF
  −0.55 vs −2.58). The highest-likelihood set is tuned to the calibration data.
- **The residue-specific part of the model is carried by aggregate and STR parameters, not by the litter chemistry.**
  Without DP, the most lignin-rich residue, the posterior moves to kAgg 0.032 (run 7: 0.008), aggfactSTR_C 3.8 (1.4)
  and NonMicAgg 4.0 (1.5). Without RS, CUE_STR becomes bimodal. The route from litter quality through STR to
  aggregates and MAOC is the weak part of SAMM.

## Files

| File | Content |
|---|---|
| `best_without_<TRT>.csv` | highest-likelihood parameter set of the fold without TRT (R̂ < 1.1; RS not converged) |
| `posterior_summary_without_<TRT>.csv` | posterior quantiles and R̂ of that fold |
| `heldout_obs_vs_sim.csv` | measurements of each residue treatment with the prediction of the fold without it: `sim` posterior median, `sim_best` best set |
| `heldout_vs_insample_best_rhat1.1.txt`, `..._rhat1.5.txt` | RMSE, EF and bias per data type and treatment (best sets) |
| `evaluation_heldout.png` | held-out posterior median vs measured, as Fig. 5 of Laub et al. (2024) |
