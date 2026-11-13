"""Creates the MONICA set-up of the Khon Kaen litter experiment for SAMM (bare fallow, annual litter additions),
matching the SAMM v1.1 R code (tests/samm/data/Samm_main.R) as closely as MONICA allows.

usage: python3 make_setup.py [power.json]
  power.json: NASA POWER daily point data (only needed to (re)create climate.csv), downloaded with
  https://power.larc.nasa.gov/api/temporal/daily/point?parameters=T2M,T2M_MAX,T2M_MIN,PRECTOTCORR,
      ALLSKY_SFC_SW_DWN,RH2M,WS2M&community=AG&longitude=102.8167&latitude=16.3333&start=19950101&end=20201231&format=JSON

Writes climate.csv, site-<modifiers>.json, crop-<treatment>.json (litter parameters inline) and
sim-<treatment>-<modifiers>.json (modifiers: SAMM = SAMM's own temperature/moisture functions,
MONICA = MONICA's functions with the rate scaling) into this directory.

Variant "Kirschbaum" (site-Kirschbaum.json, crop-<treatment>-Kirschbaum.json, sim-<treatment>-Kirschbaum.json):
the recalibrated parameter set SAMMKirschbaum (DREAM run 7, tests/samm/calibration/results/dream7) as in its
calibration: start state with the root chemistry of run 7, no constant daily root input but the weed roots of run 7's
weed model (weed_roots.csv: added on each weeding day as organic matter), hoeing at each litter addition
(tillageAggregateDestruction = aggDestroy), N deposition 7 + free-living fixation 2 kg N/ha/yr.
"""
import datetime as dt
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# as in Samm_main.R: day 1 = 1995-03-01, litter added on these simulation days, run until day 9070
START = dt.date(1995, 3, 1)
ADD_DAYS = [1, 382, 759, 1125, 1490, 1874, 2231, 2596, 2961, 3323, 3688, 4049, 4414, 4779, 5144, 5509,
            5874, 6240, 6605, 6970, 7335, 7702, 8068, 8432, 8797, 9162]
LAST_DAY = 9070
# treatment: litter C [kg C ha-1], N [kg N ha-1], structural share of C (R code), C content (Table 1 of the paper)
TREATMENTS = {
    "CT": (0, 0, 0, 0.4),
    "DP": (4530, 57, 0.32, 0.453),
    "GN": (3880, 228, 0.13, 0.388),
    "RS": (3670, 47, 0.06, 0.367),
    "TM": (4300, 136, 0.17, 0.427),
}
DEPTH, BD, SICL = 0.15, 1450.0, 0.176
# DREAM run 7: weed root chemistry (Daily_Litter_CN, Daily_LitterSTR_C), aggregate destruction by hoeing (aggDestroy)
RUN7 = os.path.join(os.path.dirname(HERE), "calibration", "results", "dream7")
ROOT_C_CONTENT = 0.42
# start state of the R code (kg C ha-1 in 0-15 cm)
INIT = dict(LMW=100.0, MIC=100.0, MAO_total=1665.0, AggMAO=575.0, AggLitter=1400.0, AggLitterSTR=0.1,
            SystemC=4900.0, LitterSTR=0.16)


def date_of_day(d):
    return START + dt.timedelta(days=d - 1)


def write_climate(power_path):
    p = json.load(open(power_path))["properties"]["parameter"]
    days = sorted(p["T2M"])
    with open(os.path.join(HERE, "climate.csv"), "w") as f:
        f.write("iso-date,tavg,tmin,tmax,wind,globrad,precip,relhumid\n")
        f.write(",C_deg,C_deg,C_deg,m/s,MJ m-2 d-1,mm,%\n")
        for d in days:
            iso = f"{d[:4]}-{d[4:6]}-{d[6:]}"
            vals = [p[k][d] for k in ("T2M", "T2M_MIN", "T2M_MAX", "WS2M", "ALLSKY_SFC_SW_DWN", "PRECTOTCORR", "RH2M")]
            if any(v == -999 for v in vals):
                raise SystemExit(f"missing value on {iso}")
            f.write(iso + "," + ",".join(f"{v:g}" for v in vals) + "\n")
    print(f"climate.csv: {len(days)} days")


def run7_params():
    return {r.split(",")[0]: float(r.split(",")[1]) for r in open(os.path.join(RUN7, "best.csv")).read().split()[1:]}


def samm_settings(modifiers):
    soc = INIT["SystemC"] - (INIT["SystemC"] - INIT["LMW"] - INIT["MIC"] - INIT["MAO_total"] - INIT["AggLitter"])
    free_litter = INIT["SystemC"] - soc
    if modifiers == "Kirschbaum":
        x = run7_params()
        return {
            "variant": "SAMMKirschbaum",
            "dailyLitterC": 0.0,
            "init_fLMW": INIT["LMW"] / soc, "init_fMIC": INIT["MIC"] / soc,
            "init_fAggLitter": INIT["AggLitter"] / soc, "init_AggLitterSTR": INIT["AggLitterSTR"],
            "init_fMAOfree": (INIT["MAO_total"] - INIT["AggMAO"]) / INIT["MAO_total"],
            "init_freeLitter": free_litter / soc, "init_freeLitterSTR": x["Daily_LitterSTR_C"],
            "init_CN_MAO": 9, "init_CN_MIC": 10, "init_CN_LMW": 20, "init_CN_Litter": x["Daily_Litter_CN"],
            "cSorp": 0.049,
            "tillageAggregateDestruction": x["aggDestroy"],
        }, soc
    return {
        "variant": "SAMM",
        "rateModifiers": modifiers,
        # the R code's constant daily root input (weeds), its start state and sorption capacity
        "dailyLitterC": 2.28, "dailyLitterCN": 142.9, "dailyLitterSTR": 0.16,
        "init_fLMW": INIT["LMW"] / soc, "init_fMIC": INIT["MIC"] / soc,
        "init_fAggLitter": INIT["AggLitter"] / soc, "init_AggLitterSTR": INIT["AggLitterSTR"],
        "init_fMAOfree": (INIT["MAO_total"] - INIT["AggMAO"]) / INIT["MAO_total"],
        "init_freeLitter": free_litter / soc, "init_freeLitterSTR": INIT["LitterSTR"],
        "init_CN_MAO": 9, "init_CN_MIC": 10, "init_CN_LMW": 20, "init_CN_Litter": 142.9,
        "cSorp": 0.049,
        # pure SAMM v1.1 (no scenario D mineral N coupling), as the R code
        "params": {"f_min": 0, "f_turn": 0, "k_des": 0},
    }, soc


def write_site(modifiers):
    samm, soc = samm_settings(modifiers)
    soc_frac = soc / (DEPTH * 10000 * BD)  # kg C kg-1 soil
    site = {
        "SiteParameters": {
            "Latitude": 16.3333, "Slope": 0, "HeightNN": [175, "m"],
            "NDeposition": [9 if modifiers == "Kirschbaum" else 0, "kg N ha-1 y-1"],
            "SoilProfileParameters": [
                {"Thickness": [DEPTH, "m"], "SoilOrganicCarbon": [soc_frac * 100, "%"],
                 # measured silt + clay 0.176 (as in the SAMM R code); 90% sand of the paper is an older classification
                 "Sand": [round(1 - SICL, 3), "kg kg-1"], "Clay": [0.05, "kg kg-1"], "pH": 5.5, "CN": 10,
                 "SoilBulkDensity": [BD, "kg m-3"]},
                {"Thickness": [1.80, "m"], "SoilOrganicCarbon": [0.1, "%"],
                 "Sand": [0.88, "kg kg-1"], "Clay": [0.07, "kg kg-1"], "pH": 5.5, "CN": 10,
                 "SoilBulkDensity": [1500, "kg m-3"]},
            ],
        },
        "SoilTemperatureParameters": ["include-from-file", "general/soil-temperature.json"],
        "EnvironmentParameters": {"=": ["include-from-file", "general/environment.json"], "LeachingDepth": 1.5,
                                  "WindSpeedHeight": 2.0},
        "SoilOrganicParameters": {"=": ["include-from-file", "general/soil-organic.json"],
                                  "SOM_Model": "SAMM", "MaxMineralisationDepth": DEPTH, "SAMM": samm},
        "SoilTransportParameters": ["include-from-file", "general/soil-transport.json"],
        "SoilMoistureParameters": ["include-from-file", "general/soil-moisture.json"],
    }
    path = os.path.join(HERE, f"site-{modifiers}.json")
    json.dump(site, open(path, "w"), indent=2)


def organic_matter(oid, name, ccont, cn, fstr):
    """all C in one AOM pool with C:N cn, so that the organic N is C / cn; structural share fstr of the C
    (lignin = share * C content / 0.63, no polyphenols)"""
    return {"type": "OrganicFertiliserParameters", "id": oid, "name": name,
            "AOM_DryMatterContent": [1.0, "kg DM kg FM-1"], "CorgContent": [ccont, "kg C kg DM-1"],
            "AOM_NH4Content": [0.0, "kg N kg DM-1"], "AOM_NO3Content": [0.0, "kg N kg DM-1"],
            "AOM_CarbamidContent": 0.0,
            "PartAOM_to_AOM_Fast": [1.0, "kg kg-1"], "PartAOM_to_AOM_Slow": [0.0, "kg kg-1"],
            "CN_Ratio_AOM_Fast": cn, "CN_Ratio_AOM_Slow": 100,
            "AOM_FastDecCoeffStandard": [0.002, "d-1"], "AOM_SlowDecCoeffStandard": [0.0002, "d-1"],
            "PartAOM_Slow_to_SMB_Slow": [0.0, "kg kg-1"], "PartAOM_Slow_to_SMB_Fast": [1.0, "kg kg-1"],
            "NConcentration": 0.0,
            "SAMM_Lignin": [fstr * ccont / 0.63, "kg kg DM-1"], "SAMM_Polyphenols": [0.0, "kg kg DM-1"]}


def write_crops_kirschbaum():
    """litter additions, hoeing at each litter addition (all treatments) and the weed roots of run 7"""
    roots = [r.split(",") for r in open(os.path.join(RUN7, "weed_roots.csv")).read().split()[1:]]
    for trt, (c, n, fstr, ccont) in TREATMENTS.items():
        steps = []
        for d in ADD_DAYS:
            if d > LAST_DAY:
                continue
            if c > 0:
                steps.append((d, 0, {"type": "OrganicFertilization", "amount": [c / ccont, "kg FM ha-1"],
                                     "incorporation": True,
                                     "parameters": organic_matter(f"KK_{trt}", f"Khon Kaen litter {trt}", ccont, c / n, fstr)}))
            steps.append((d, 1, {"type": "Tillage", "depth": [DEPTH, "m"]}))
        for t, day, rc, rstr, rcn in roots:
            if t != trt or int(float(day)) > LAST_DAY:
                continue
            rc, rstr, rcn = float(rc), float(rstr), float(rcn)
            steps.append((int(float(day)), 2, {"type": "OrganicFertilization", "amount": [rc / ROOT_C_CONTENT, "kg FM ha-1"],
                                               "incorporation": True,
                                               "parameters": organic_matter("KK_weed_roots", "weed roots (run 7)",
                                                                            ROOT_C_CONTENT, rcn, rstr)}))
        steps.sort(key=lambda s: (s[0], s[1]))
        worksteps = [dict({"date": date_of_day(d).isoformat()}, **w) for d, _, w in steps]
        crop = {"crops": {}, "cropRotation": [{"worksteps": worksteps}], "CropParameters": ["include-from-file", "general/crop.json"]}
        json.dump(crop, open(os.path.join(HERE, f"crop-{trt}-Kirschbaum.json"), "w"), indent=1)


def write_organic_matter_and_crops():
    for trt, (c, n, fstr, ccont) in TREATMENTS.items():
        worksteps = []
        if c > 0:
            # all C in one AOM pool with the litter's C:N, so that the organic N is the litter N;
            # structural share as in the R code: lignin = share * C content / 0.63, no polyphenols
            om = {"type": "OrganicFertiliserParameters", "id": f"KK_{trt}", "name": f"Khon Kaen litter {trt}",
                  "AOM_DryMatterContent": [1.0, "kg DM kg FM-1"], "CorgContent": [ccont, "kg C kg DM-1"],
                  "AOM_NH4Content": [0.0, "kg N kg DM-1"], "AOM_NO3Content": [0.0, "kg N kg DM-1"],
                  "AOM_CarbamidContent": 0.0,
                  "PartAOM_to_AOM_Fast": [1.0, "kg kg-1"], "PartAOM_to_AOM_Slow": [0.0, "kg kg-1"],
                  "CN_Ratio_AOM_Fast": c / n, "CN_Ratio_AOM_Slow": 100,
                  "AOM_FastDecCoeffStandard": [0.002, "d-1"], "AOM_SlowDecCoeffStandard": [0.0002, "d-1"],
                  "PartAOM_Slow_to_SMB_Slow": [0.0, "kg kg-1"], "PartAOM_Slow_to_SMB_Fast": [1.0, "kg kg-1"],
                  "NConcentration": 0.0,
                  "SAMM_Lignin": [fstr * ccont / 0.63, "kg kg DM-1"], "SAMM_Polyphenols": [0.0, "kg kg DM-1"]}
            for d in ADD_DAYS:
                if d > LAST_DAY:
                    continue
                worksteps.append({"date": date_of_day(d).isoformat(), "type": "OrganicFertilization",
                                  "amount": [c / ccont, "kg FM ha-1"], "incorporation": True,
                                  "parameters": om})
        else:
            # MONICA needs a work step in the rotation; a zero-amount application does nothing
            worksteps.append({"date": date_of_day(1).isoformat(), "type": "OrganicFertilization",
                              "amount": [0.0, "kg FM ha-1"], "incorporation": True,
                              "parameters": ["include-from-file", "organic-fertilisers/CAM.json"]})
        crop = {"crops": {}, "cropRotation": [{"worksteps": worksteps}], "CropParameters": ["include-from-file", "general/crop.json"]}
        json.dump(crop, open(os.path.join(HERE, f"crop-{trt}.json"), "w"), indent=2)


def write_sims():
    pools = ["STR_C", "LAB_C", "LMW_C", "MIC_C", "MAO_C", "AggSTR_C", "AggLAB_C", "AggMAO_C",
             "LAB_N", "LMW_N", "MIC_N", "MAO_N", "AggLAB_N", "AggMAO_N"]
    for modifiers in ("SAMM", "MONICA", "Kirschbaum"):
        for trt in TREATMENTS:
            sim = {
                "crop.json": f"crop-{trt}-Kirschbaum.json" if modifiers == "Kirschbaum" else f"crop-{trt}.json", "site.json": f"site-{modifiers}.json", "climate.csv": "climate.csv",
                "climate.csv-options": {"no-of-climate-file-header-lines": 2, "csv-separator": ",",
                                        "start-date": date_of_day(1).isoformat(),
                                        "end-date": date_of_day(LAST_DAY).isoformat()},
                "debug?": False, "include-file-base-path": "${MONICA_PARAMETERS}/",
                "NumberOfLayers": 13, "LayerThickness": [DEPTH, "m"],
                "UseSecondaryYields": False, "NitrogenResponseOn": True, "WaterDeficitResponseOn": True,
                "EmergenceMoistureControlOn": False, "EmergenceFloodingControlOn": False,
                "UseAutomaticIrrigation": False, "UseNMinMineralFertilisingMethod": False,
                "output": {"write-file?": True, "path-to-output": "./", "file-name": "out.csv",
                           "csv-options": {"include-header-row": True, "include-units-row": True,
                                           "include-aggregation-rows": False, "csv-separator": ","},
                           "events": ["daily", ["Date"] + [[f"SAMM_{p}", 1] for p in pools]
                                      + [["STemp", 1], ["Mois", 1], ["Sat", 1], ["SAMM_CO2", 1], ["pF", 1], "SAMM_LeachC",
                                         ["SOC", 1], ["NH4", 1], ["NO3", 1], ["WaterFlux", 1], "Precip"]]},
            }
            json.dump(sim, open(os.path.join(HERE, f"sim-{trt}-{modifiers}.json"), "w"), indent=2)
        write_site(modifiers)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        write_climate(sys.argv[1])
    write_organic_matter_and_crops()
    write_crops_kirschbaum()
    write_sims()
    print("set-up written to", HERE)
