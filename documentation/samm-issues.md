# SAMM in MONICA: critical issues found and their status

Status as of 2026-11-12, branch `samm-som-module`. Details of the implementation are in [samm.md](samm.md).

Summary:

| # | Issue | Status |
|---|---|---|
| 1 | Negative microbial C: one-day RK4 step unstable on hot, wet days | **fixed** |
| 2 | R v1.1 code: no lower limit of the LMW carbon use efficiency | **fixed** in MONICA (R code unchanged) |
| 3 | R v1.1 code: litterbag refilled before its last measurement | **fixed** in the calibration tool (R code unchanged) |
| 4 | DAISY results changed by the SAMM / Laub-flow code | **fixed** (DAISY bit-identical) |
| 5 | Measurement dates 15 days off the climate calendar | **fixed** (matched by simulation day) |
| 6 | Spin-up crash, climate file parsing, litter chemistry table parsing | **fixed** |
| 7 | Aggregate formation constant while aggregate breakdown follows the rate modifier: free MAOC drained into aggregates in the dry season | **open** (experimental option, overcorrects) |
| 8 | Rates scaled to the mean of another rate modifier: effective `mu_max` above 1 d⁻¹ | **in progress** (recalibration with Kirschbaum × MONICA moisture) |
| 9 | Recalibration: parameters at their bounds, very fast aggregate turnover, microbial C too high | **open** |
| 10 | SAMM's pF moisture scalar with MONICA's pF of sandy soils | **open** (MONICA's moisture function used instead) |
| 11 | Silt + clay of the Khon Kaen soil: 0.10 in MONICA, 0.176 (measured) in the R code | **fixed** (sand 82.4%) |
| 12 | Total MAOC lower with MONICA's drivers than with HYDRUS (control) | **open** (cause not confirmed) |
| 13 | Soil N overestimated in the litter treatments | **open** |
| 14 | N mineralisation at temperate sites 2.8-3.5 × DAISY | **open** (needs calibration) |
| 15 | Windows release `3.7.4.samm` built before the sub-steps | **open** (needs a new build) |

---

## Fixed

### 1. Negative microbial C (numerical instability)

- **Symptom:** in the Khon Kaen groundnut treatment, microbial C fell below zero from 2019-04-22 to 2019-05-14
  (minimum −25 kg C ha⁻¹). SOC and microbial C showed short dips in other years too.
- **Cause:** not the litter addition but the first heavy rain of the season (32 mm) on soil at 37.7 °C. MONICA's rate
  modifier reached corr_F = 7.27 (median 1.05). LMW-N was turned over 2.5-3 times per day, more than a one-day
  RK4 step can handle (stable up to about 2.8 d⁻¹).
  - In the RK4 stages LMW-N went 36 → 2.9 → 63.5 → −129 kg N ha⁻¹.
  - With negative LMW-N the CUE of LMW became −53 (no lower limit, see 2), so microbial growth and, through it,
    aggregate formation became strongly negative.
  - SAMM v1.1 in R driven by MONICA's drivers stays just stable, but has the same weakness.
- **Fix** (commit `a1f820f`):
  - adaptive sub-steps: a day is split into n steps when a pool would lose more than its size per step, and redone
    with 2n steps if an RK4 stage would make a pool negative (`maxTurnoverPerStep`, `maxSubsteps`; output
    `SAMM_Substeps`);
  - needed on 0-49 of 9070 days per treatment at Khon Kaen, never with the HYDRUS drivers or at temperate sites.
- **Checks:** C++ core vs R unchanged (2.4e-14), DAISY and SAMM-on test suites identical, no negative pools.

### 2. No lower limit of the LMW carbon use efficiency (R code)

- `CUE_CN_LMW` in `Samm_main.R` is capped at 1 but has no lower limit; a negative LMW-N or LMW-C (possible in an
  RK4 stage) gives a negative CUE and negative microbial growth.
- **Fix in MONICA's C++ core** (commit `a1f820f`): CUE = 0 if LMW-C or LMW-N ≤ 0, and aggregate formation is never
  negative. Only acts on negative stage states, so results are unchanged otherwise. `Samm_main.R` is unchanged.

### 3. Litterbag refilled before its last measurement (R code)

- `Samm_main.R` v1.1 empties the litterbag pools 13 days before each litter addition and refills them. The last
  litterbag measurement (simulation day 4799, 385 days after the bag was added on day 4414) is therefore compared
  with the next year's new bag (1689-4709 instead of 37-1721 kg C ha⁻¹). The R code still carries commented-out
  second-year litterbag pools.
- **Fix in the calibration tool** (commit `a8a7672`): the bag added on day 4414 is followed until its last
  measurement. Litterbag EF of SAMM v1.1: 0.69 instead of 0.12 (paper: 0.80). `Samm_main.R` is unchanged.
- Note: day 4414 is April 2007 on the model calendar (Table 2 of the paper says 2004).

### 4. DAISY results changed

- A changed order of floating-point operations in the Laub et al. (2020) flow made default DAISY results differ in
  the last digits. **Fixed** by restoring the original expression order. DAISY (default) is bit-identical to the
  original MONICA in all 66 test-suite simulations, checked again after every SAMM change.

### 5. Measurement dates 15 days off

- The date column of `AllMeasurements.csv` is 15 days later than the climate calendar. **Fixed** by matching by
  "Day of simulation", as the R calibration does.

### 11. Silt + clay of the Khon Kaen soil

- MONICA took silt + clay = 1 − sand = 0.10 from the 90% sand of the paper (an older classification); the measured
  value, used in the R code, is 0.176. `MAO_C_max` was therefore 10.7 instead of 18.8 t C ha⁻¹.
- **Fixed:** the 0-15 cm layer of the Khon Kaen set-up now has 82.4% sand (5% clay), so MONICA computes 0.176 itself
  (`make_setup.py`, `installer/KhonKaen-SAMM/site.json`). This also changes MONICA's water retention: saturation stays
  0.468, the median water content rises from 0.102 to 0.120. Fit of MONICA (default modifiers): SOC EF 0.40 → 0.54,
  aggregate C RMSE 1452 → 1112, free MAOC and microbial N unchanged, soil N EF −6.05 → −8.14.

### 6. Smaller fixes during the implementation

- Spin-up crash (vector index): spin-up used the wrong profile object. Fixed.
- Test harness: trailing empty rows of the climate file. Fixed.
- Litter chemistry: forage, vegetable and energy-crop tables of Thiébeau et al. were parsed from the wrong columns.
  Fixed (columns detected from the statistics row).
- Output ids: new outputs are appended, so existing output ids do not shift.

---

## Open

### 7. Aggregate formation constant, aggregate breakdown follows the rate modifier

- **Symptom:** with MONICA's drivers, free MAOC drains into aggregates in the dry season and is released when the
  rains come; in the control it drops to 20-80 kg C ha⁻¹ (measured about 1500-1900). Aggregate C swings between
  about 1500 and 7000 kg C ha⁻¹. With the HYDRUS drivers both are nearly constant, as measured.
- **Cause:** physicochemical aggregate formation (`NonMicAgg`, 51 kg C ha⁻¹ d⁻¹ of MAOC with `aggfactMAO_C`) does not
  depend on corr_F, while aggregate breakdown (`kAgg · corr_F`) does. With HYDRUS, corr_F stays at 1.6-2.9, so the
  equilibrium is steady; with MONICA's moisture function corr_F falls to about 0.3 in the dry season, and the
  equilibrium aggregate pool would exceed all MAOC. MONICA's moisture function also slows decomposition in wet soil
  (bell shape around 45% of saturation).
- **Not a MONICA coupling error:** SAMM v1.1 in R run with MONICA's temperature, moisture, percolation and modifier
  reproduces MONICA (mean SOC difference 16-254 kg C ha⁻¹; control free MAOC at the end 18 vs 18 kg C ha⁻¹).
  It is a property of SAMM's equations that HYDRUS's nearly constant corr_F hides.
- **Tested:** `NonMicAggCorrFRef` (formation proportional to corr_F, commit `e8fe404`, off by default). In MONICA it
  improves SOC (EF 0.36 → 0.54), aggregate C (RMSE 1471 → 891) and free MAOC (755 → 612), but in R with HYDRUS it
  removes the seasonal dip of aggregate C that the 2019 data show (aggregate C RMSE 419 → 848): it overcorrects.
- **Proposed:** partial coupling `NonMicAgg · (corr_F/ref)^β`, β calibrated together with the other parameters.

### 8. Rates scaled to the mean of another rate modifier

- SAMM's rates were calibrated with SAMM's modifiers (2^((T−20)/10) × pF scalar, mean corr_F 1.94 at Khon Kaen).
  MONICA's modifiers were matched to that mean with one factor (`rateScaleMONICA` = 0.545; Kirschbaum option
  `rateScaleKirschbaum` = 4.90). The peaks then differ: effective `mu_max` up to 0.337 × 7.27 ≈ 2.4 d⁻¹ (MONICA)
  or ≈ 1.7 d⁻¹ (Kirschbaum), although it should be a maximum rate. SAMM's own function also exceeds 1 on hot days
  (effective `mu_max` about 1.2 d⁻¹).
- **In progress:** recalibration with Kirschbaum's temperature function × MONICA's moisture function (both ≤ 1, no
  scale factor), `mu_max` ≤ 0.95, genetic algorithm with the objective of the paper (Gaussian log-likelihood of
  litterbag C, microbial N, SOC, aggregate C and free MAOC, median variance per data type). Tool:
  `tests/samm/calibration/`.
  - Reference: SAMM v1.1 with its own modifiers, log-likelihood −4277.7.
  - Bounded run b1: −4035.3 (stopped early after 60 generations without improvement); run b2 −3995 at
    generation 457 and still improving; benchmark run with SAMM's modifiers and the same bounds pending.
- The recalibrated set is not yet in MONICA (planned: `"variant": "SAMMKirschbaum"`, scale factor 1).

### 9. Problems of the recalibration so far

- **Parameters at their bounds.** Without plausible bounds the GA ran to 8 kg C ha⁻¹ d⁻¹ of weed-root input,
  `CUE_STR` 0.99, `k_MAO` 0.01 and `MicMaxCN` 15. With bounds (root input ≤ 5, `CUE_STR` ≤ 0.7, `k_MAO` ≤ 0.0041,
  `MicMaxCN` ≤ 12), b1 still ends at or near the bounds of `MicMaxCN`, `aggfactSTR_C`, root input, `k_MAO` and
  `kAgg`. The calibration compensates for something structural (most likely issue 7).
- **Very fast aggregate turnover** (`kAgg` 0.74 d⁻¹, v1.1: 0.032): aggregates react to every rain event, so aggregate
  C and SOC swing by up to ±2000 kg C ha⁻¹ within a year; free MAOC still drops in the dry season.
- **Microbial C too high** (EF −18.7, v1.1 −3.8): only microbial N is in the objective (as in the paper); with
  `MicMaxCN` 12 the microbial C:N becomes too wide. Options: bound `MicMaxCN` to the measured range (about 4-8) or
  add microbial C to the objective.
- **Fit (b1 vs v1.1):** better for litterbag C (EF 0.80 vs 0.69), microbial N (0.18 vs −0.11), SOC (0.65 vs 0.63,
  control much better) and soil N (−0.33 vs −1.24); worse for aggregate C (−0.13 vs 0.24) and free MAOC
  (0.06 vs 0.73).
- **No uncertainty:** the GA gives one best set; the SIR step of the paper was not done (with about 500 data
  points the importance weights would collapse onto very few samples).

### 10. SAMM's pF moisture scalar with MONICA's pF

- For the Khon Kaen sand MONICA's van Genuchten-Vereecken retention curve gives a median pF of 2.7 and up to 7.7 in
  the dry season (HYDRUS about 1.9 all year). With SAMM's pF scalar the mean moisture factor is 0.53 instead of
  1.00, and SOC builds up about twice as fast as in R. MONICA's moisture function (relative saturation) is used
  instead (default `"rateModifiers": "MONICA"`, and in the recalibration). MONICA's pF for sandy soils is not fixed.

### 12. Total MAOC lower with MONICA's drivers

- In the control, free + aggregated MAOC stays at about 1750 kg C ha⁻¹ with MONICA's drivers but builds up to about
  2400 with HYDRUS. Reproduced by R with MONICA's drivers, so it comes from the drivers. Likely cause (not confirmed):
  microbes decline almost to zero in the control with HYDRUS, which slows MAOC decomposition; with MONICA's
  drivers they stay at about 60 kg C ha⁻¹ (closer to the measurements).

### 13. Soil N overestimated in the litter treatments

- Measured soil N stays roughly constant; SAMM accumulates N in the litter treatments, most with groundnut and
  tamarind, in R (EF −1.24) and more in MONICA (−6.05). The recalibration improves it (−0.33). Soil N was not used
  in the calibration of the paper.

### 14. N mineralisation at temperate sites

- With SAMM switched on in the Müncheberg test suite, SOC fits better than DAISY (RMSE 3.9-4.6 t in fertilised
  treatments vs 7.9-18.7 for DAISY), but net N mineralisation is 2.8-3.5 times DAISY's. SAMM is calibrated for one
  tropical site only; temperate sites need a calibration (Kirschbaum's function is also much more temperature
  sensitive below 20 °C than SAMM's Q10 = 2, which Khon Kaen cannot constrain).

### 15. Windows release

- The released Windows build `3.7.4.samm` predates the adaptive sub-steps and the Kirschbaum option; a new build is
  needed to use them on Windows. The standalone SAMM package (`tests/samm/calibration/standalone/`) is up to date.

---

## Not an issue

- **Litterbags in MONICA:** litterbags are a calibration device of the R code and are not simulated in MONICA.
- **Adaptive sub-steps and the R reference:** with the HYDRUS drivers no sub-steps are needed, so the comparison with
  the R code stays exact.
