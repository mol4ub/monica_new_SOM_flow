# MONICA 3.7.4.laub-som-flow: running with and without the Laub et al. (2020) SOM flow

Hand-over notes for running the release `3.7.4.laub-som-flow` of MONICA with
the original DAISY soil organic matter (SOM) structure and with the alternative
structure from Laub et al. (2020), and comparing the two.

- Release: https://github.com/mol4ub/monica_new_SOM_flow/releases/tag/3.7.4.laub-som-flow
- Tag `3.7.4.laub-som-flow` = commit `22c928c` (branch `claude/fervent-albattani-83j1r5`)
- Paper: Laub, M. et al. (2020), DRIFTS band areas as measured pool size proxy to
  reduce parameter uncertainty in soil organic matter models, Biogeosciences 17,
  1393–1413, https://doi.org/10.5194/bg-17-1393-2020 (Fig. 1 = original DAISY,
  Fig. 7 = alternative structure)

## What changed

One new boolean in the soil organic parameters:

```json
"__enable_laub_SOM_flow__": true
```

It is **off by default** (`false`). With it off, every result is bit-identical to
MONICA 3.7.4 (checked with all 66 simulations of `tests/suites.toml`).

| | Original DAISY, Fig. 1 (switch off) | Laub et al. 2020, Fig. 7 (switch on) |
|---|---|---|
| Decomposed SOM_fast | fraction `f` → SOM_slow; rest → SMB_slow + CO₂ | all → SMB_slow + CO₂; **no SOM_fast → SOM_slow transfer** |
| Dead SMB carbon entering SOM | all → SOM_fast | fraction `f` → SOM_slow, `1 − f` → SOM_fast |
| Initial slow share of SOM (steady state) | `1 / (1 + k_slow / (f · k_fast))` | `1 / (1 + (1 − f) · k_slow / (f · k_fast))` |

- `f` is the existing parameter `PartSOM_Fast_to_SOM_Slow` (the paper's humification
  efficiency `f_SOM_slow`) in both structures. No new parameter. `k_slow`/`k_fast` are
  `SOM_SlowDecCoeffStandard`/`SOM_FastDecCoeffStandard`.
- With the same `f`, both structures give the same steady-state carbon flow into SOM_slow
  and SMB_slow. Only the SOM_fast pool is smaller, by `(1 − f)`.
- Dead SMB carbon that does not enter SOM is recycled into SMB_fast per pool (using
  `PartSMB_Slow_to_SOM_Fast` and `PartSMB_Fast_to_SOM_Fast` separately). With the default
  parameters, where both are 0.6, this is the same as the original.
- Code: `src/core/soilorganic.cpp` (constructor: initialisation; `fo_MIT()`: fluxes),
  `src/core/monica-parameters.{h,cpp}` (the switch).

## 1. Get MONICA

The release contains a **Windows** build only.

### Windows

Download `monica_win64_3.7.4.laub-som-flow.zip` from the release page and unzip it. It contains:

```
bin/monica-run.exe                        the model
monica-parameters/                        parameter library (MONICA_PARAMETERS points here)
monica-parameters/general/soil-organic.json
projects/Hohenfinow2-single-run/          example: sim.json, site.json, crop.json, climate.csv, run_monica.cmd
```

`projects/Hohenfinow2-single-run/run_monica.cmd` sets `MONICA_PARAMETERS=..\..\monica-parameters`
and runs `..\..\bin\monica-run.exe -o .\sim-out.csv sim.json`.

### Linux (and macOS) – build the tag with Nix

Tested; it takes about 2 minutes with the Nix binary cache:

```bash
git clone --branch 3.7.4.laub-som-flow --recurse-submodules https://github.com/mol4ub/monica_new_SOM_flow.git monica-laub
cd monica-laub
nix build --extra-experimental-features 'nix-command flakes' '.?submodules=1'
./result/bin/monica-run --version        # monica-run version 3.7.4.0
export MONICA_PARAMETERS="$PWD/monica-parameters"
```

The clone provides `monica-parameters/` (the Nix package contains only `bin/`) and the
example `installer/Hohenfinow2/` (`sim-min.json`, `site-min.json`, `crop-min.json`,
`climate-min.csv`). Without Nix, `pixi run build` builds into `_cmake_release/` (see
`how-to-build.md`).

`monica-run [options] path-to-sim-json`. Useful options: `-o FILE` (output file),
`-sd`/`-ed` (start/end date, ISO), `-d` (debug output). Relative paths in a `sim.json` are
relative to the sim file. `include-from-file` paths are relative to
`include-file-base-path`, which is normally `${MONICA_PARAMETERS}/`, unless absolute.

## 2. Run without the new flow (original DAISY)

Nothing to change: any existing setup runs the original structure.

```bash
cd installer/Hohenfinow2
../../result/bin/monica-run -o out-default.csv sim-min.json
```

## 3. Run with the new flow

The switch belongs in the soil organic parameters, which a `site.json` loads with
`"SoilOrganicParameters": ["include-from-file", "general/soil-organic.json"]`.

Don't replace that line with an object holding only the switch. The other parameters
would then fall back to the C++ defaults, which differ from `soil-organic.json`
(e.g. `SOM_SlowDecCoeffStandard` 4.3e-5 instead of 1.3e-5). Instead, either:

**A. A copy of `soil-organic.json` with the switch, and a site file pointing at it (recommended)**

```python
import json, os
P = os.environ["MONICA_PARAMETERS"]
so = json.load(open(f"{P}/general/soil-organic.json"))
so["__enable_laub_SOM_flow__"] = True
json.dump(so, open("soil-organic-laub.json", "w"), indent=2)

site = json.load(open("site-min.json"))
site["SoilOrganicParameters"] = ["include-from-file", os.path.abspath("soil-organic-laub.json")]
json.dump(site, open("site-laub.json", "w"), indent=2)

sim = json.load(open("sim-min.json"))
sim["site.json"] = "site-laub.json"
json.dump(sim, open("sim-laub.json", "w"), indent=2)
```

```bash
../../result/bin/monica-run -o out-laub.csv sim-laub.json
```

**B. A copy of the whole parameter folder**, with `"__enable_laub_SOM_flow__": true` added
to its `general/soil-organic.json`, and `MONICA_PARAMETERS` pointing at the copy. This is
handy for running whole test suites, e.g.
`python3 tests/run.py --monica ./result/bin/monica-run --params /path/to/params-laub --out results/laub`.

Use the same crop, site, climate and other parameters in both runs, so only the structure differs.

## 4. Check that the switch is on

MONICA ignores unknown keys silently, so a typo in `__enable_laub_SOM_flow__` silently
runs the original structure. MONICA doesn't log the switch either. Check the initial pools:
add SOM outputs to the sim's `output.events`, e.g.

```json
"events": ["daily", ["Date", ["SOMs", [1, 3]], ["SOMf", [1, 3]], ["SOC", [1, 3]], "Rh"]]
```

On the first day, the slow share `SOMs_1 / (SOMs_1 + SOMf_1)` must equal the formula of the
structure in the table above. With the default `soil-organic.json`
(`k_slow` = 1.3e-5, `k_fast` = 8.6e-5, `f` = 0.3) it is:

| | switch off | switch on |
|---|---|---|
| slow share of SOM on day 1 | **0.665** | **0.739** |
| Hohenfinow2 `sim-min.json`, 1991-01-01, layer 1 | SOMs 7.4271, SOMf 3.7426 kg C m⁻³ | SOMs 8.2571, SOMf 2.9125 kg C m⁻³ |

Total SOC at the start is the same in both. Only its split between the pools differs.

## 5. Outputs for comparing the two structures

| Output id | Unit | Meaning |
|---|---|---|
| `SOMs`, `SOMf` | kg C m⁻³ (per layer) | slow / fast soil organic matter pools |
| `SMBs`, `SMBf` | kg C m⁻³ (per layer) | slow / fast microbial biomass |
| `AOMs`, `AOMf` | kg C m⁻³ (per layer) | added organic matter (residues, manure) |
| `SOC` | kg C kg⁻¹ (per layer) | soil organic carbon content |
| `SOC-X-Y` | g C m⁻² (per layer, sum with `[1, 3, "SUM"]` for 0–30 cm) | SOC stock |
| `Rh` | kg C ha⁻¹ d⁻¹ | heterotrophic (decomposer) respiration |
| `NetNmin` | kg N ha⁻¹ d⁻¹ | net N mineralisation |

Layers are 10 cm. Organic matter turnover is simulated down to `MaxMineralisationDepth`
(0.4 m, i.e. layers 1–4). `tests/deviation.py results/default results/laub` summarises
how much every output changes between two result folders.

## 6. What to expect (MONICA default parameters, no recalibration)

- Initial SOM_slow about +11 %, SOM_fast about −22 % (the different steady-state split).
- Net N mineralisation higher, on average by about 4–9 % over a run (Hohenfinow2 and
  microlysimeter simulations).
- Müncheberg V140 long-term trial (6 treatments, 1963–2024): SOC in 2024 is 3–5 % lower
  with the new flow, and RMSE against measured SOC is lower in all 6 treatments. Both
  structures still overestimate SOC gains in fertilised and manured treatments, so most of
  that error comes from the parameters.
- Carbon is conserved in both structures (pool change + respiration balance on days without
  organic inputs).

A fair comparison may need recalibrating `SOM_SlowDecCoeffStandard`,
`SOM_FastDecCoeffStandard` and `PartSOM_Fast_to_SOM_Slow` for the new structure. Laub et
al. (2020) Table 5 gives Bayesian-calibrated Daisy values (normalised to 10 °C; MONICA's
temperature and moisture functions differ, so treat them as a starting range, not as
drop-in values): `k_fast` 2.29e-4 d⁻¹, `k_slow` 3.25e-5 d⁻¹, `f` 0.34.

## 7. Known limitations

- **Saved model states:** the switch isn't written into a serialised MONICA state (same as
  the existing `__enable_kaiteew_*__` switches). A run that loads a saved state
  (`serializedMonicaState.load`) runs the original structure unless the parameters are
  given again.
- **Not in `src/landcare-dss`:** that module has its own copy of the parameters.
- **Unchanged from MONICA 3.7.4, in both structures:** when a SOM pool would turn negative,
  `fo_MIT()` sets its change to `+pool` instead of `−pool`. It didn't occur in the tested
  simulations.
- **No Linux/macOS binary in the release:** build the tag as in section 1.
