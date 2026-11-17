# SAMM in MONICA: project hand-over

Complete record of the work on this repository up to 2026-11-17, written so that a new person or a new
Claude Code session can continue from the repository alone. Read this first, then:

- [samm.md](samm.md): the technical documentation of SAMM in MONICA (model, every code change, parameters,
  validation). Section numbers below refer to it.
- [samm-issues.md](samm-issues.md): all critical issues found, fixed and open.
- [laub-som-flow.md](laub-som-flow.md): the earlier change (alternative DAISY SOM flow of Laub et al. 2020).
- [tests/samm/calibration/results/dream7/README.md](../tests/samm/calibration/results/dream7/README.md) and
  [results/cv_leave_one_residue_out/README.md](../tests/samm/calibration/results/cv_leave_one_residue_out/README.md):
  the final calibration and its cross-validation.

---

## 1. What the project is

**MONICA** (Model for Nitrogen and Carbon in Agro-ecosystems, ZALF, upstream `zalf-rpm/monica`) is a C++ crop and
soil model. This repository is the fork **`mol4ub/monica_new_SOM_flow`**. Its goal is to give MONICA better soil
organic matter (SOM) models:

1. **Laub et al. (2020) SOM flow** (done, release `3.7.4.laub-som-flow`): an optional alternative structure of the
   DAISY SOM pools, switch `"__enable_laub_SOM_flow__"`. See [laub-som-flow.md](laub-som-flow.md).
2. **SAMM** (main work, branch `samm-som-module`): the microbial / aggregate SOM model *Soil Aggregation through
   Microbial Mediation* (Laub et al. 2024, GMD 17:931, plus corrigendum = SAMM v1.1) as a second SOM module
   beside DAISY, switch `"SOM_Model": "SAMM"`.
3. **Recalibration of SAMM** for rate modifiers that never exceed 1 (Kirschbaum temperature × MONICA moisture),
   with a fast compiled calibration tool, Bayesian calibration (DREAM) against all data of the Khon Kaen litter
   experiment, and cross-validation. Result: the parameter set `"variant": "SAMMKirschbaum"` in MONICA.

The project owner works on SAMM and on the Khon Kaen experiment (Thailand), knows the model and data in depth, and
decides the scientific questions; the role of Claude was implementation, calibration runs, analysis and
documentation.

**Khon Kaen experiment** (calibration data of SAMM): bare fallow on a sandy soil, NE Thailand (16.33 N, 102.82 E),
1995-2019, yearly additions of about 10 t DM/ha of four litters: CT control, DP dipterocarp leaves, GN groundnut
stover, RS rice straw, TM tamarind leaves. Hand hoeing at each litter addition (also in CT), weeding monthly in the
rainy season (spontaneous weeds, above-ground biomass removed, roots stay). Measured: SOC, soil N, microbial C and N
(fumigation-extraction; microbial N reported with the factor 3.1, microbial C with k_EC 0.33), litterbag C,
aggregate C and free MAOC (2019 only). Silt + clay of 0-15 cm: 0.176 (measured; the 10% of the paper is an older
classification).

---

## 2. Repository map (what matters for this project)

| Path | Content |
|---|---|
| `src/core/samm.h`, `samm.cpp` | **SAMM core**, no MONICA dependencies: `Params` (incl. `samm_v11()`, `sammNoAgg_v11()`, `samm_kirschbaum()`), pools, fluxes, derivatives, RK4, adaptive sub-steps (`stepAdaptive`), LMW transport between layers, balance check, rate modifier functions |
| `src/core/samm-litter-chemistry.h/.cpp` | built-in lignin/polyphenol table (generated from `data/samm/`) |
| `src/core/soilorganic.h/.cpp` | MONICA's SOM module; SAMM branch: `initSAMM`, `spinupSAMM`, `fo_SAMM` (daily step, rate modifier), `addOrganicMatterSAMM`, `applyTillage`, getters |
| `src/core/monica-parameters.h/.cpp` | `SOM_Model`, `SAMMSettings` (JSON object `"SAMM"`, `merge()` handles `"variant"`), chemistry fields of organic matter |
| `src/io/build-output.cpp` | outputs `SAMM_*`, `SAMM_CO2`, `SAMM_LeachC`, `SAMM_Substeps`, `pF` |
| `data/samm/` | litter chemistry with sources (`litter-chemistry.csv`), generator |
| `tests/samm/data/` | SAMM v1.1 R code (`Samm_main.R`), HYDRUS drivers of the paper (`climate_file.csv`), all measurements (`AllMeasurements.csv`) |
| `tests/samm/khonkaen_forced.cpp`, `run_ref.R`, `compare_to_r.py` | C++ core against the R code with the paper's drivers |
| `tests/samm/khonkaen_monica/` | Khon Kaen set-up for MONICA (generator `make_setup.py`, three variants), comparison and evaluation scripts |
| `tests/samm/calibration/` | calibration tool `samm_calib.cpp`, priors, helper scripts, `results/` (run 7, cross-validation, earlier runs), `standalone/` (Windows/Linux package of the tool) |
| `installer/KhonKaen-SAMM/` | ready-to-run Windows project (copied into the release zip as `projects/KhonKaen-SAMM`) |
| `.github/workflows/on-tag.yml` | Windows build + draft release (manual run with input `tag`) |
| `documentation/samm*.md`, `documentation/images/samm_*` | documentation |
| `CLAUDE.md` | short orientation for Claude Code sessions |

Other documents in the repository from earlier sessions: `how-to-build.md` (all ways to build MONICA),
`PERFORMANCE.md`, `QUIRK-FIXES.md` (three C++ fixes), `issues-to-look-at.md`, `tests/README.md` (test suites).

### Branches and releases

| Branch | Content |
|---|---|
| `master` | upstream MONICA 3.7.4 (+ test suite additions) |
| `claude/fervent-albattani-83j1r5` | Laub et al. (2020) SOM flow; tag `3.7.4.laub-som-flow` |
| **`samm-som-module`** | all SAMM work (based on `3.7.4.laub-som-flow`); **work and push here** |

| Release (GitHub, Windows zip) | Commit | Content |
|---|---|---|
| `3.7.4.laub-som-flow` | `22c928c` | SOM flow switch |
| `3.7.4.samm` | `a41835a` | SAMM v1.1 in MONICA, before adaptive sub-steps and Kirschbaum |
| `3.7.4.samm-kirschbaum` | `c66d7e4` | **current**: SAMMKirschbaum, adaptive sub-steps, SAMM documentation and `projects/KhonKaen-SAMM` in the zip |

Leftover "untagged" draft releases (created when publishing re-triggers the workflow) can be deleted.

---

## 3. How SAMM works in MONICA (short; details in samm.md)

- Switch: `"SoilOrganicParameters": {"SOM_Model": "SAMM", "SAMM": {...}}`. DAISY stays the default and is
  bit-identical to MONICA 3.7.4 when SAMM is off (all 66 simulations of `tests/suites.toml`).
- Pools per organic layer (down to `MaxMineralisationDepth`): STR, LAB, LMW, MIC, MAO and aggregate pools
  AggSTR, AggLAB, AggMAO (C, and N except STR). SOC passed to the rest of MONICA = all pools except free litter.
- Daily step: classical RK4 for all organic layers together (layers coupled by LMW transport with MONICA's water
  flux); **adaptive sub-steps** when a pool would turn over more than `maxTurnoverPerStep` (default 1) per step or
  an RK4 stage would go negative; safety net (LMW CUE ≥ 0, aggregate formation ≥ 0).
- Rate modifier `corr_F` (`"rateModifiers"`):
  - `"MONICA"` (default): MONICA's kaiteew temperature × moisture functions × `rateScaleMONICA` 0.545;
  - `"SAMM"`: the paper's 2^((T−20)/10) × pF scalar;
  - `"Kirschbaum"`: exp(−3.764 + 0.204·T·(1 − T/73.8)) × exp(−18 (θ/θsat − 0.45)²), **no scaling** (only with the
    recalibrated set).
- Parameter sets (`"variant"`): `"SAMM"` = v1.1 (default), `"SAMMnoAgg"`, `"SAMMKirschbaum"` (run 7, also sets
  `rateModifiers` to `"Kirschbaum"`). Single parameters in `"params"` override the variant.
- Mineral N ("scenario D"): microbes take N from LMW-N; gross mineralisation `f_min`, `f_turn` to NH4; optional
  immobilisation (`min_imm_frac`, `min_imm_k`) and desorption (`k_des`).
- Inputs: crop residues, roots (crop module callback with root flag) and organic fertilisers are split into STR
  (0.63 × lignin + 0.53 × polyphenols, as C share) and LAB, chemistry from the parameters or the built-in table;
  optional constant daily litter input (`dailyLitterC`). Tillage breaks up `tillageAggregateDestruction` of the
  aggregates in the tilled layers and mixes them.
- Start state: fixed fractions of the measured SOC (`init_*`), or optional spin-up.

**Parameter set `SAMMKirschbaum`** (DREAM run 7 best set, rounded): k_STR 0.007702, k_LAB 0.1176, k_MIC 0.007573,
k_MAO 3.537e-05, mu_max 0.7057, kAgg 0.009948, K_M_MIC_C 2.568, m_MIC 9.721e-05, K_LMW_MAO 0.002656, CUE_STR 0.5584,
CUE_LAB 0.07513, MicMinCN 4.971, MicMaxCN 6.247, f_MICMAOM 0.07545, pc_STR_LAB 2.091, aggfactSTR_C 1.395,
aggfactMAO_C 0.1470, NonMicAgg 1.738, betaMIC 1.967, f_min 0.6575, f_turn 0.8632, min_imm_frac 0, min_imm_k 0,
k_des 0; defaults c_SORP 0.049, CUE_LMW 0.6, MIC_ref 100. Khon Kaen settings: tillageAggregateDestruction 0.002539,
root / start-litter C:N 21.49, structural share 0.283. Full precision: `tests/samm/calibration/results/dream7/best.csv`.

---

## 4. Building and running (cloud container or Linux)

Everything below runs from the repository root. The cloud containers have `nix`; `/tmp/.../scratchpad`
directories are lost when a container is replaced, so rebuild there and keep results in the repository.

### MONICA (Linux, nix)

```bash
nix build --extra-experimental-features 'nix-command flakes' '.?submodules=1' -o result-monica
export MONICA_PARAMETERS=$PWD/monica-parameters
```

Takes several minutes (dependencies come from the nix cache). Binary: `result-monica/bin/monica-run`. Uncommitted
changes to tracked files are included ("dirty" build); new files must be `git add`-ed to be seen by the flake.
Other ways (pixi, CMake, Windows): `how-to-build.md`.

### Calibration tool and C++ core tests (no MONICA needed, seconds)

```bash
g++ -O3 -march=native -std=c++17 -pthread -Wall -I src/core tests/samm/calibration/samm_calib.cpp src/core/samm.cpp -o samm_calib
g++ -O2 -std=c++17 -I src/core tests/samm/khonkaen_forced.cpp src/core/samm.cpp -o khonkaen_forced
```

Regression check after any change to the SAMM core (must print `logLik -5098.317`):

```bash
./samm_calib tests/samm/data/climate_file.csv tests/samm/data/AllMeasurements.csv eval tests/samm/calibration/results/dream7/best.csv daily.csv obs.csv --lik ar1 --weeds --ndep 7 --nfix 2 --tillage --flexN --leachN
```

### R (for plots and the R reference model)

```bash
nix build --impure --extra-experimental-features 'nix-command flakes' -o r-env --expr 'let pkgs = (builtins.getFlake "nixpkgs").legacyPackages.x86_64-linux; in pkgs.rWrapper.override { packages = with pkgs.rPackages; [ deSolve ]; }'
```

Then `r-env/bin/Rscript`. Plot scripts use base R only. R reference run: `cd tests/samm/data && ../../../r-env/bin/Rscript ../run_ref.R Samm_main.R R_out.csv`.

### Khon Kaen in MONICA

`tests/samm/khonkaen_monica/make_setup.py` writes three set-ups (already committed): `sim-<TRT>-SAMM.json` (SAMM
modifiers), `sim-<TRT>-MONICA.json` (v1.1 with MONICA's modifiers, "blue line" of samm.md 7.3) and
`sim-<TRT>-Kirschbaum.json` (SAMMKirschbaum with the weed roots of run 7, hoeing, N deposition 9; samm.md 7.5).
Run each in its own directory (output `out.csv` in the working directory):

```bash
mkdir -p runs/GN && cp tests/samm/khonkaen_monica/{sim-GN-Kirschbaum.json,crop-GN-Kirschbaum.json,site-Kirschbaum.json,climate.csv} runs/GN/
(cd runs/GN && ../../result-monica/bin/monica-run sim-GN-Kirschbaum.json)
python3 tests/samm/khonkaen_monica/eval_kirschbaum.py runs obs.csv      # needs runs/<TRT>/out.csv for all five
```

Evaluation plot (Fig. 5 style of the paper): `Rscript tests/samm/calibration/plot_eval.R <dir> tests/samm/data <run> out.png ["caption"]`
reads `<dir>/<run>_obs.csv` (columns TRT,day,type,obs,sim) and optionally `<dir>/<run>_band.csv`.

### Test suites

`python3 tests/run.py --monica <monica-run> --engine <label> --out <dir>`, compare with `tests/compare.py`, see
`tests/README.md`. DAISY must stay bit-identical; SAMM changes can be checked by running old and new builds on the
same set-up and comparing `out.csv` byte for byte (done for every change so far).

### Windows release

GitHub Actions workflow `on-tag.yml`, manual run (workflow_dispatch) on branch `samm-som-module` with input
`tag` (a new, unused tag name, e.g. `3.7.4.samm-<name>`). Takes about 10 min, creates a **draft** release with
`monica_win64_<tag>.zip` (executables, parameters, SAMM documentation, `projects/KhonKaen-SAMM`). The owner
publishes it in the browser (Releases page, logged in; draft links only work for people with write access).
Publishing creates the tag, which starts the workflow again; that second run fails or leaves an "untagged" draft
that can be deleted. In cloud sessions without `gh`, use the GitHub MCP tool `actions_run_trigger`
(`run_workflow`, `workflow_id: on-tag.yml`, `ref: samm-som-module`, `inputs: {"tag": "..."}`).

### Standalone SAMM package (calibration tool for Windows + Linux)

`tests/samm/calibration/standalone/build_package.sh <out_dir>` cross-compiles `samm_calib` with MinGW-w64 from
nixpkgs (static, with mcfgthread), adds data, parameter files, batch files, plot script and sources. The last
package (zip, built 2026-11-11/12) still has only SAMM v1.1 parameters and predates the weed model; **an update with
the run 7 parameters is open** (the owner postponed it until the calibration was final; it now is).

---

## 5. The calibration tool `samm_calib`

`tests/samm/calibration/samm_calib.cpp`: one 0-15 cm layer, exactly the set-up of `Samm_main.R` (start state,
litter additions, litterbag pools, classical RK4 with daily steps and the 4th stage on the next day's drivers),
HYDRUS drivers of `climate_file.csv`. All five treatments in about 0.4 s; 4 threads. Header comment lists all options.

Modes: `eval` (one parameter set → daily pools, measured vs simulated, weed root inputs), `ga`, `sir`, `dream`
(with `--resume`), `band` (posterior intervals from a posterior file).

Options used for the final set-up (run 7): `--lik ar1 --weeds --ndep 7 --nfix 2 --tillage --flexN --leachN`
(default modifiers = Kirschbaum × MONICA moisture with θsat 0.468 and optimum 0.45; `--modifiers SAMM` for the
paper's). Others: `--exclude "<type>"`, `--holdout <TRT>` (cross-validation), `--chains`, `--maxGen`, `--thin`,
`--rhat`, `--seed`, `--threads`, `--pop`, `--gens`, `--n`.

- **Data types**: Litterbag C, Microbial N, SOC (incl. LMW-C), Aggregate C, Free MOAC, Microbial C, Soil N
  (LMW_N + MIC_N + MAO_N + AggLAB_N + AggMAO_N, i.e. without free LAB_N, consistent with SOC).
- **Likelihood `paper`**: Gaussian, median measurement variance per type, five types (as the paper).
- **Likelihood `ar1`**: per treatment × data type a Kalman filter: measurement error (median variance r of the type)
  plus model error with variance (sigmaRatio)²·r (sigmaRatio 0-3) and exponential temporal correlation
  exp(−Δt/τ) (τ 1-365 d). sigmaRatio and τ per type are calibrated.
- **Observation model**: reported microbial N = MIC_N × 3.1 / fN_ninhydrin (fN fitted, one value for all years;
  the data were converted with 3.1, Amato & Ladd 1988); microbial C with k_EC 0.33 is taken as is.
- **Weed model (`--weeds`)**: replaces the constant daily root input. Logistic growth
  dW = rmax·fT·fW·fN·W·(1 − W/Wmax), fT Kirschbaum, fW linear between θ 0.06 and 0.15, fN = minN/(minN + KN);
  own mineral N pool fed by SAMM's gross mineralisation (f_min, f_turn calibrated, k_des 0, no immobilisation),
  N deposition + fixation, drained by weed uptake and leaching; weeding on the 15th Apr-Sep, 15 Nov/Jan/Mar and at
  litter days; above-ground biomass removed, roots (C 0.42 of DM, C:N and STR share = Daily_Litter_CN /
  Daily_LitterSTR_C) added to STR/LAB; restart from 20 kg DM/ha. Soft constraint from field photos: mean
  rainy-season above-ground weed biomass at weeding GN 1.5 ± 0.75, CT 0.3 ± 0.15 t DM/ha.
- **Tillage (`--tillage`)**: share aggDestroy of the aggregates broken up at every litter addition.
- **`--flexN`** (weed shoot N 0.8-2.5% by N availability) and **`--leachN`** (leachable share of mineral N) were
  added to improve the fit; they are fit-driven, not physically justified (consider dropping them).
- **DREAM(ZS)**: 4 chains, archive of 10·d start points (jitter 0.1 prior SD around the start file), archive
  addition every 10 generations, γ = 2.38/√(2δd*) with randomised crossover (1/3, 2/3, 1), 20% γ = 1 jumps, 10%
  snooker updates, folding at bounds; truncated Gaussian priors (`prior_*.csv`: name, mean, sd, lo, hi);
  Gelman-Rubin R̂ on the second half every 1000 generations, stop when all R̂ < `--rhat` (generation ≥ 10000).
  Summary: `dream_summary.py chains.csv prefix` → `_summary.csv`, `_best.csv`, `_map.csv`, `_posterior.csv`.
- **Speed**: about 28 s per 1000 DREAM generations (4 chains × 1 run) on 4 vCPUs; run 7 took 40 min.

---

## 6. Chronology of the work and what was learned

### 6.1 Laub et al. (2020) SOM flow (branch `claude/fervent-albattani-83j1r5`)
Optional switch for the alternative DAISY pool structure; off = bit-identical; release `3.7.4.laub-som-flow`.

### 6.2 SAMM in MONICA (from commit `e0a9551`)
- C++ core written from the paper and the R code; reproduces `Samm_main.R` to about 1e-11 relative with the
  paper's drivers (`khonkaen_forced`, samm.md 7.1).
- Integrated into `SoilOrganic` as an alternative module; litter chemistry researched for all MONICA residues, roots
  and manures (`data/samm/`); Khon Kaen set-up for MONICA; outputs; documentation.
- Silt + clay at Khon Kaen set to the measured 0.176 by reducing sand to 82.4% (no extra parameter; the owner's
  explicit wish).
- Comparison MONICA vs R (samm.md 7.3): with MONICA's own temperature/moisture/water flux, SAMM v1.1 in MONICA
  differs from R mainly through the drivers; R run with MONICA's drivers reproduces MONICA.

### 6.3 Problems found and fixed (samm-issues.md)
- **Negative microbial C** (GN 2019): not the litter but a hot first heavy rain gave corr_F 7.27 and an unstable
  daily RK4 step (LMW-N negative in stage 4 → negative CUE). Fixed with adaptive sub-steps + safety net.
- **R v1.1 litterbag**: the R code refills the litterbag pools before the last measurement of the litterbag
  experiment; the calibration tool follows the day-4414 bag (litterbag EF of v1.1 0.12 → 0.69).
- **Scaled rate modifiers**: MONICA's (×0.545) and Kirschbaum (×4.90) modifiers scaled to SAMM's mean made the
  effective mu_max exceed 1 d⁻¹. Hence the recalibration with unscaled Kirschbaum × MONICA moisture; the scale for
  Kirschbaum was removed from MONICA.
- Aggregate formation (`NonMicAgg`) independent of corr_F while breakdown depends on it: drains free MAOC into
  aggregates in the dry season with MONICA's drivers. Coupling formation to corr_F (`NonMicAggCorrFRef`) overcorrected
  and was removed again. Still open (issue 7).

### 6.4 Recalibration (Kirschbaum × MONICA moisture, mu_max ≤ 0.95)

Log-likelihoods are only comparable within the same likelihood and data set. Best sets, posterior summaries and
logs of all these runs: `tests/samm/calibration/results/earlier_runs/` (README there).

| Step | Set-up | Result |
|---|---|---|
| reference | SAMM v1.1, SAMM modifiers, paper likelihood | −4277.7 |
| GA b1, b2 | genetic algorithm (linear rank, BLX-0.5, Gaussian mutation, 5% elitism), paper likelihood, plausible bounds | b2 −3772.6; parameters at bounds |
| SIR | as in the paper, informed Gaussian prior (paper posterior, 50% wider) | failed: weights collapse, ESS ≈ 1 with ~850 observations, 45 parameters |
| DREAM run 1 | paper likelihood, 4 chains from b2 | −3649, implausible parameters |
| run 2 | seven data types, `ar1` likelihood | microbial C fits, litterbag poor |
| run 4 | + density-dependent microbial death betaMIC (prior 1.5, ≤ 2) | −3659 (other data set), betaMIC 1.67 |
| run 5 | + fitted ninhydrin factor, microbial C back in, narrow microbial C:N prior | −5104, fN 2.09 at its bound |
| run 6 | + weed model | −5108; EF SOC 0.81, MAOC 0.81, litterbag 0.63, aggregate C −2.26 |
| **run 7** | + N deposition 7, fixation 2, tillage, flexN, leachN; `prior_kirschbaum_ar1_weeds_till.csv` | **−5098.3**, converged at 68,000 generations (R̂ < 1.1) |

Run 7 fit (best set): litterbag C RMSE 732 / EF 0.62; microbial N 24.5 / 0.23; SOC 634 / 0.80; aggregate C 880 /
−2.36; free MAOC 177 / 0.81; microbial C 93 / 0.05; soil N 116 / 0.25. Notable: CUE_LAB 0.08 (litter loss bypasses
microbes), betaMIC at its upper bound, fN 2.0 instead of 3.1, aggDestroy ≈ 0 (tillage effect rejected), weeds in CT
below the photo constraint. Files and plots: `results/dream7/`.

### 6.5 Leave-one-residue-out cross-validation (`results/cv_leave_one_residue_out/`)
Four DREAM calibrations with the run 7 set-up, each without one residue treatment (`--holdout`), continued with
`--resume` to R̂ < 1.1 (RS not converged: CUE_STR bimodal, one mode near 0). Posterior-median prediction of the
held-out treatment: litterbags and microbial N transfer (litterbag EF 0.56 vs 0.71 in-sample), **SOC and free MAOC do
not** (SOC EF 0.09 vs 0.71; DP +1.3, RS +0.7, TM −1.0 t C/ha). Without DP the aggregate parameters change strongly.
Conclusion: the pathway litter quality → STR → aggregates / MAOC is the structural weak point of SAMM; CUE lower
bounds of 0.1-0.2 are advisable. The R̂ < 1.5 results were too optimistic (chains had not left the start region).

### 6.6 SAMMKirschbaum in MONICA (samm.md 7.5)
Variant added to MONICA; Khon Kaen variant Kirschbaum with run 7's weed roots (`results/dream7/weed_roots.csv`,
written by `samm_calib eval ... obs.csv roots.csv`) as organic matter on 248 weeding days per treatment. Control
reproduced; residue treatments end 0.5-1.6 t C/ha higher (SOC EF 0.63 vs 0.80). Cause confirmed: MONICA drains
781 mm/yr below 15 cm vs 1146 mm/yr in the HYDRUS drivers → less LMW leaching → more microbial growth and aggregates
(the tool driven by MONICA's drivers reproduces MONICA). The owner plans to calibrate the MONICA version
independently later.

---

## 7. Open points and possible next steps

1. Calibrate SAMM with MONICA's own drivers (the tool accepts any driver file in the format of
   `climate_file.csv`: date;day;soil temperature 7.5 cm;θ;drainage below 15 cm [mm];pF), or check MONICA's water
   balance at Khon Kaen first (issue 16).
2. Structural work on litter quality → STR → aggregates / MAOC (cross-validation); CUE lower bounds; aggregate
   formation vs corr_F (issue 7).
3. Clean-up decided with the owner: keep MONICA lean (done: `rateScaleKirschbaum`, `NonMicAggCorrFRef` removed);
   `samm_calib` keeps all its options. Possibly drop the fit-driven `--flexN` / `--leachN` in definitive runs.
4. Update the standalone package with the run 7 parameters (section 4).
5. Steady-state initialisation of SAMM (discussed, not implemented; current options: fixed fractions or spin-up).
6. Temperate sites: N mineralisation with SAMM is 2.8-3.5 × DAISY (issue 14); needs calibration.
7. Whether `SAMMKirschbaum` should become the default variant (currently v1.1 stays the default; existing set-ups
   unchanged).

---

## 8. Working conventions with the owner

- Push only to `samm-som-module`; no pull requests unless asked. Commit messages: imperative summary line,
  body with details; the repository's history shows the style.
- Before committing changes to MONICA: rebuild, check that existing set-ups give byte-identical output (DAISY and
  the SAMM v1.1 set-ups), run the `samm_calib` regression (`logLik -5098.317`).
- Do not introduce parameters the owner did not ask for; prefer justified mechanisms over fit-driven additions and
  say which is which. Ask before large new directions ("let's not overdo things").
- Report results with numbers (RMSE, EF, bias per data type), plots in the Fig. 5 style of the paper
  (`plot_eval.R`), and state caveats (convergence, start points, what is fitted).
- Documentation in English, plain and concrete; keep `samm.md`, `samm-issues.md` and this file up to date.
- Microbial data in plots: either as model pools with the observation model applied (microbial N as reported), or as
  extractable C (MBC × 0.33) and ninhydrin-N (MBN / 3.1), as the owner asked.
- Literature: the SAMM paper pages (gmd.copernicus.org) and zenodo were not reachable from the container; the owner
  uploads PDFs when needed.
