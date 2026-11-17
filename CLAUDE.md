# Notes for Claude Code sessions

This is `mol4ub/monica_new_SOM_flow`, a fork of MONICA (C++ crop/soil model, ZALF) in which the soil organic
matter model **SAMM** (Laub et al. 2024, GMD 17:931, v1.1) was added as an alternative to DAISY and recalibrated.

**Start by reading `documentation/samm-handover.md`**: what was done, how everything works, how to build and run,
results, open points and the conventions agreed with the owner. Technical details: `documentation/samm.md`;
issues: `documentation/samm-issues.md`.

## Essentials

- Work on and push to branch **`samm-som-module`** only (`git push -u origin samm-som-module`). No pull requests
  unless the owner asks.
- DAISY is the default SOM model and must stay bit-identical. SAMM is switched on with
  `"SOM_Model": "SAMM"`; parameter sets `"variant"`: `"SAMM"` (v1.1, default), `"SAMMnoAgg"`, `"SAMMKirschbaum"`
  (recalibrated, DREAM run 7).
- SAMM core: `src/core/samm.h/.cpp` (no MONICA dependencies). MONICA coupling: `src/core/soilorganic.cpp`
  (`fo_SAMM`, `initSAMM`, `addOrganicMatterSAMM`, `applyTillage`), parameters: `src/core/monica-parameters.*`
  (`SAMMSettings`).
- Calibration tool: `tests/samm/calibration/samm_calib.cpp` (GA, SIR, DREAM, eval, band; Khon Kaen, HYDRUS
  drivers). Final results: `tests/samm/calibration/results/`.

## Build and check

```bash
# MONICA (several minutes); new files must be git-added to be seen by the flake
nix build --extra-experimental-features 'nix-command flakes' '.?submodules=1' -o result-monica
export MONICA_PARAMETERS=$PWD/monica-parameters
# calibration tool (seconds)
g++ -O3 -march=native -std=c++17 -pthread -Wall -I src/core tests/samm/calibration/samm_calib.cpp src/core/samm.cpp -o samm_calib
# regression of the SAMM core: must print logLik -5098.317
./samm_calib tests/samm/data/climate_file.csv tests/samm/data/AllMeasurements.csv eval tests/samm/calibration/results/dream7/best.csv daily.csv obs.csv --lik ar1 --weeds --ndep 7 --nfix 2 --tillage --flexN --leachN
```

Before committing MONICA changes, compare `out.csv` of old and new builds on the Khon Kaen set-ups
(`tests/samm/khonkaen_monica/`) and, for anything touching DAISY, the test suites (`tests/README.md`).

Windows release: GitHub workflow `on-tag.yml`, manual run on `samm-som-module` with a new `tag` input; creates a
draft release that the owner publishes (see handover, section 4).

## Conventions

- Keep MONICA lean; do not add parameters the owner did not ask for; separate justified mechanisms from
  fit-driven additions and say which is which.
- Report fits with RMSE / EF / bias per data type and Fig. 5-style plots (`tests/samm/calibration/plot_eval.R`).
- Keep `samm.md`, `samm-issues.md` and `samm-handover.md` up to date with every change.
- Container scratch directories are temporary: commit results worth keeping.
