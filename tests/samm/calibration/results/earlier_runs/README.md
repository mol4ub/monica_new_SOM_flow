# Earlier calibration runs (before the final run 7)

Best parameter sets (`*_best.csv`, last line = log-likelihood), posterior summaries of the DREAM runs
(`*_posterior_summary.csv`: mean, sd, quantiles, R̂) and progress logs (`*_log.txt`) of the calibration runs that led
to run 7 (`../dream7`). Kept for reference; the full chains were not kept. Log-likelihoods are only comparable
between runs with the same likelihood and data types. Context: `documentation/samm-handover.md`, section 6.4, and
`documentation/samm-issues.md`, issues 8 and 9.

All runs use Kirschbaum × MONICA moisture (no scaling, mu_max ≤ 0.95) unless noted.

| Run | Method / set-up | Best logLik | Notes |
|---|---|---|---|
| `bs` | GA benchmark with SAMM's own modifiers, started from SAMM v1.1 | −4277.7 | stopped early; v1.1 reference value |
| `ga1`, `ga2` | first GA runs, paper likelihood (5 data types, 494 observations) | −4101.7 | exploratory |
| `k1`, `k2` | GA, paper likelihood, wide bounds | −3722.5, −3860.3 | without plausible bounds the GA ran to implausible values (root input 8 kg C/ha/d, CUE_STR 0.99, k_MAO 0.01, MicMaxCN 15) |
| `b1` | GA with plausible bounds | −4035.3 | stopped after 60 generations without improvement |
| `b2` | GA with plausible bounds | −3772.6 | best GA set; start of the DREAM runs |
| `sir` | SIR as in the paper, informed Gaussian prior (paper posterior, 50% wider), 200,000 draws | −4563.1 | failed: weights collapse (ESS ≈ 1) |
| `dream1` | DREAM(ZS), 4 chains from b2, paper likelihood, `prior_kirschbaum.csv` | −3649.0 | converged (42,000 gen.); implausible parameters |
| `dream2` | seven data types (855 observations), correlated-error likelihood `ar1` | −5072.8 | converged (112,000 gen.); microbial C fits, litterbag poor |
| `dream3` | microbial C excluded, model error ≤ 3 × measurement SD, τ ≤ 365 d | – | stopped at 28,000 generations to add betaMIC (log only) |
| `dream4` | as dream3 + density-dependent death betaMIC (prior 1.5, ≤ 2) | −3659.1 | converged (52,000 gen.); betaMIC 1.67 |
| `dream5` | + fitted ninhydrin factor fN, microbial C back in, narrow microbial C:N prior (`prior_kirschbaum_ar1_fN.csv`) | −5103.6 | converged (45,000 gen.); fN 2.09 at its lower bound |
| `dream6` | + simple weed model instead of constant root input (`prior_kirschbaum_ar1_weeds.csv`) | −5108.3 | converged (110,000 gen.); EF SOC 0.81, MAOC 0.81, litterbag 0.63, aggregate C −2.26 |
| run 7 | + N deposition / fixation, tillage, flexN, leachN | −5098.3 | final, see `../dream7` |
