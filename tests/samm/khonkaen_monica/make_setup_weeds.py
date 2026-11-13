"""Khon Kaen with spontaneous weed vegetation simulated by MONICA instead of SAMM's constant daily root input.

The weeds are a generic tropical grass: species crops/Bacharia.json (C3, as defined there; switched to C4 it did
not grow with the available parameters), cultivar of perennial rye grass
without vernalisation/daylength requirement (crops/rye-grass/backup.json). They are sown again after every weeding and
removed by "harvesting" with all above-ground biomass exported (above-ground biomass was taken off the plots); living
roots stay in the soil (MONICA adds them at harvest; dead roots are added daily). Weeding: monthly in the rainy season
(April-September) and every second month in the dry season (Laub et al. 2022), and at every litter incorporation.
SAMM's mineral N coupling is switched on (f_min, f_turn of scenario D, no desorption), so mineralised N feeds the weeds.

usage: python3 make_setup_weeds.py params.csv out_dir
  params.csv: calibrated SAMM parameters (name,value as written by samm_calib); writes for every treatment
  sim-<TRT>-weeds.json (A) and sim-<TRT>-const.json (B: same parameters, constant daily root input, no weeds)."""
import copy, datetime as dt, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
TRTS = ["CT", "DP", "GN", "RS", "TM"]
SAMM_PARAMS = ["k_STR", "k_LAB", "k_MIC", "k_MAO", "mu_max", "kAgg", "K_M_MIC_C", "m_MIC", "K_LMW_MAO", "CUE_STR",
               "CUE_LAB", "MicMinCN", "MicMaxCN", "f_MICMAOM", "pc_STR_LAB", "aggfactSTR_C", "aggfactMAO_C", "NonMicAgg",
               "betaMIC"]
WEED = {
    "is-winter-crop": False,
    "cropParams": {
        "species": ["include-from-file", "crops/Bacharia.json"],  # as defined (C3); switching to C4 did not grow
        "cultivar": {"=": ["include-from-file", "crops/rye-grass/backup.json"],
                     "OrganIdsForPrimaryYield": [
                         {"organId": 2, "type": "YieldComponent", "yieldDryMatter": 1.0, "yieldPercentage": 1.0},
                         {"organId": 3, "type": "YieldComponent", "yieldDryMatter": 1.0, "yieldPercentage": 1.0}]},
    },
    "residueParams": ["include-from-file", "crop-residues/Bacharia.json"],
}


def find_samm(o):
    if isinstance(o, dict):
        if "rateModifiers" in o: return o
        for v in o.values():
            r = find_samm(v)
            if r is not None: return r
    elif isinstance(o, list):
        for v in o:
            r = find_samm(v)
            if r is not None: return r


def weeding_dates(start, end):
    d = []
    for y in range(start.year, end.year + 1):
        for m in range(1, 13):
            if 4 <= m <= 9 or m in (11, 1, 3):   # monthly in the rainy season, every second month in the dry season
                x = dt.date(y, m, 15)
                if start < x < end: d.append(x)
    return d


def main(params_csv, out):
    os.makedirs(out, exist_ok=True)
    p = dict(l.strip().split(",") for l in open(params_csv) if l.strip() and not l.startswith("name"))
    p = {k: float(v) for k, v in p.items()}
    site = json.load(open(os.path.join(HERE, "site-MONICA.json")))
    samm = find_samm(site)
    samm["rateModifiers"] = "Kirschbaum"
    samm["rateScaleKirschbaum"] = 1.0
    samm["params"] = {k: p[k] for k in SAMM_PARAMS}
    samm["params"].update({"f_min": 0.0, "f_turn": 0.0, "k_des": 0.0})
    samm["dailyLitterC"] = p["Daily_Litter_C"]; samm["dailyLitterCN"] = p["Daily_Litter_CN"]
    samm["dailyLitterSTR"] = p["Daily_LitterSTR_C"]
    site_const = site
    site_weeds = copy.deepcopy(site)
    sw = find_samm(site_weeds)
    sw["dailyLitterC"] = 0.0
    sw["params"].update({"f_min": 0.5, "f_turn": 0.10})
    json.dump(site_const, open(os.path.join(out, "site-const.json"), "w"), indent=1)
    json.dump(site_weeds, open(os.path.join(out, "site-weeds.json"), "w"), indent=1)
    for t in TRTS:
        crop = json.load(open(os.path.join(HERE, f"crop-{t}.json")))
        json.dump(crop, open(os.path.join(out, f"crop-{t}-const.json"), "w"), indent=1)
        cw = copy.deepcopy(crop)
        cw["crops"]["WEED"] = WEED
        ws = cw["cropRotation"][0]["worksteps"]
        litter_days = sorted(dt.date.fromisoformat(w["date"]) for w in ws)
        start, end = dt.date(1995, 3, 1), dt.date(2019, 12, 29)
        kills = sorted(set(weeding_dates(start, end)) | {d for d in litter_days if d > start})
        new = list(ws) + [{"date": (start + dt.timedelta(days=1)).isoformat(), "type": "Sowing", "crop": ["ref", "crops", "WEED"]}]
        for d in kills:
            new.append({"date": d.isoformat(), "type": "Harvest", "exported": True})
            if d + dt.timedelta(days=1) < end:
                new.append({"date": (d + dt.timedelta(days=1)).isoformat(), "type": "Sowing", "crop": ["ref", "crops", "WEED"]})
        order = {"OrganicFertilization": 0, "Harvest": 1, "Sowing": 2}
        new.sort(key=lambda w: (w["date"], order.get(w["type"], 3)))
        cw["cropRotation"][0]["worksteps"] = new
        json.dump(cw, open(os.path.join(out, f"crop-{t}-weeds.json"), "w"), indent=1)
        for var in ("const", "weeds"):
            sim = json.load(open(os.path.join(HERE, f"sim-{t}-MONICA.json")))
            sim["crop.json"] = f"crop-{t}-{var}.json"; sim["site.json"] = f"site-{var}.json"
            sim["climate.csv"] = os.path.join(HERE, "climate.csv")
            ev = sim["output"]["events"]
            ev[1] = ev[1] + ["Crop", "AbBiom", ["OrgBiom", "Root"], "SAMM_Substeps"]
            ev += ["Harvest", ["Date", "Crop", "AbBiom", "Yield", ["OrgBiom", "Root"]]]
            json.dump(sim, open(os.path.join(out, f"sim-{t}-{var}.json"), "w"), indent=1)
    print(f"{len(kills)} weedings, set-up written to {out}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
