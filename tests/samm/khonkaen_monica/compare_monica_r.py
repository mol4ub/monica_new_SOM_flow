"""Compares SAMM in MONICA (Khon Kaen set-up of make_setup.py) with SAMM v1.1 in R and with the measurements.

usage: python3 compare_monica_r.py R_out.csv runs_dir out_dir
  R_out.csv: output of Samm_main.R (all treatments, columns time, pools, TRT)
  runs_dir:  directory with <TRT>-<SAMM|MONICA>/out.csv (MONICA daily output of the sim files)
  out_dir:   writes comparison.csv (daily series) and prints the statistics
"""
import csv
import datetime as dt
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")
POOLS = ["STR_C", "LAB_C", "LMW_C", "MIC_C", "MAO_C", "AggSTR_C", "AggLAB_C", "AggMAO_C",
         "LAB_N", "LMW_N", "MIC_N", "MAO_N", "AggLAB_N", "AggMAO_N"]
TRTS = ["CT", "DP", "GN", "RS", "TM"]
START = dt.date(1995, 3, 1)
KGHA = 0.15 * 10000  # kg m-3 -> kg ha-1 in the 15 cm layer


def soc(p):  # SOC as in the R code's plot and the paper: MIC + MAO + aggregate pools
    return p["MIC_C"] + p["MAO_C"] + p["AggSTR_C"] + p["AggLAB_C"] + p["AggMAO_C"]


def read_r(path):
    out = {}
    for r in csv.DictReader(open(path)):
        out[(r["TRT"], int(float(r["time"])))] = {p: float(r[p]) for p in POOLS}
    return out


def read_monica(path):
    """state at the end of each date -> R time index (R row t = state at the start of day t)"""
    lines = open(path).read().splitlines()
    k = lines.index('"daily"')
    hdr = lines[k + 1].split(",")
    res, drivers = {}, {}
    for line in lines[k + 3:]:
        f = line.split(",")
        if len(f) != len(hdr) or not f[0][:4].isdigit():
            continue
        row = dict(zip(hdr, f))
        day = (dt.date.fromisoformat(row["Date"]) - START).days + 1
        res[day + 1] = {p: float(row[f"SAMM_{p}_1"]) * KGHA for p in POOLS}
        drivers[day] = {"temp": float(row["STemp_1"]), "theta": float(row["Mois_1"]),
                        "co2": float(row["SAMM_CO2_1"]) * KGHA, "leachC": float(row["SAMM_LeachC"])}
    return res, drivers


def read_measurements():
    m = {}
    for r in csv.DictReader(open(os.path.join(DATA, "AllMeasurements.csv"), encoding="utf-8-sig"), delimiter=";"):
        m.setdefault((r["Treatment"], r["Measurement type"]), []).append((int(r["Day of simulation"]), float(r["Value"])))
    return m


def read_climate():
    rows = []
    with open(os.path.join(DATA, "climate_file.csv"), encoding="utf-8-sig") as f:
        next(f)
        for line in f:
            p = line.strip().split(";")
            if len(p) < 6 or not p[2]:
                break
            rows.append((float(p[2]), float(p[3]), float(p[5])))
    return rows


def stats(obs, sim):
    n = len(obs)
    mo = sum(obs) / n
    sse = sum((o - s) ** 2 for o, s in zip(obs, sim))
    sst = sum((o - mo) ** 2 for o in obs)
    return math.sqrt(sse / n), (1 - sse / sst) if sst > 0 else float("nan"), sum(s - o for o, s in zip(obs, sim)) / n


def main(r_path, runs, out_dir):
    R = read_r(r_path)
    meas = read_measurements()
    runs_data = {}
    for trt in TRTS:
        for mod in ("SAMM", "MONICA"):
            runs_data[(trt, mod)] = read_monica(os.path.join(runs, f"{trt}-{mod}", "out.csv"))

    # 1. drivers: MONICA's simulated soil temperature and water content vs the HYDRUS drivers of the R code
    clim = read_climate()
    _, drv = runs_data[("CT", "SAMM")]
    days = [d for d in sorted(drv) if d <= len(clim)]
    t_r = [clim[d - 1][0] for d in days]; t_m = [drv[d]["temp"] for d in days]
    w_r = [clim[d - 1][1] for d in days]; w_m = [drv[d]["theta"] for d in days]
    def corr(a, b):
        ma, mb = sum(a) / len(a), sum(b) / len(b)
        return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b))
    print("Drivers (bare fallow CT), R/HYDRUS 7.5 cm vs MONICA layer 0-15 cm, mean and correlation:")
    print(f"  soil temperature  R {sum(t_r)/len(t_r):.2f} C   MONICA {sum(t_m)/len(t_m):.2f} C   r {corr(t_r, t_m):.2f}")
    print(f"  water content     R {sum(w_r)/len(w_r):.3f}     MONICA {sum(w_m)/len(w_m):.3f}     r {corr(w_r, w_m):.2f}")

    # 2. SOC of MONICA vs R
    print("\nSOC (MIC+MAO+aggregates, kg C/ha 0-15 cm): start, end (day 9070) and daily RMSE against R")
    print("  trt | R end | MONICA(SAMM mod.) end, RMSE vs R | MONICA(MONICA mod.) end, RMSE vs R")
    rows = []
    for trt in TRTS:
        line = f"  {trt} | {soc(R[(trt, 9070)]):7.0f} |"
        for mod in ("SAMM", "MONICA"):
            M, _ = runs_data[(trt, mod)]
            common = [t for t in range(2, 9071) if t in M and (trt, t) in R]
            rmse = math.sqrt(sum((soc(M[t]) - soc(R[(trt, t)])) ** 2 for t in common) / len(common))
            line += f" {soc(M[9070]):7.0f}, {rmse:6.0f} |"
        print(line)
        for t in range(1, 9071):
            r = R[(trt, t)]
            row = {"TRT": trt, "time": t, "date": (START + dt.timedelta(days=t - 1)).isoformat(),
                   "SOC_R": soc(r), "MICN_R": r["MIC_N"], "AggC_R": r["AggSTR_C"] + r["AggLAB_C"] + r["AggMAO_C"],
                   "MAOC_R": r["MAO_C"]}
            for mod in ("SAMM", "MONICA"):
                M, _ = runs_data[(trt, mod)]
                m = M.get(t)
                row[f"SOC_{mod}"] = soc(m) if m else ""
                row[f"MICN_{mod}"] = m["MIC_N"] if m else ""
                row[f"AggC_{mod}"] = (m["AggSTR_C"] + m["AggLAB_C"] + m["AggMAO_C"]) if m else ""
                row[f"MAOC_{mod}"] = m["MAO_C"] if m else ""
            rows.append(row)

    # 3. against the measurements, all treatments together
    print("\nAgainst the measurements (all treatments): RMSE, Nash-Sutcliffe EF, mean bias")
    targets = [("SOC", "SOC"), ("Microbial N", "MICN"), ("Aggregate C", "AggC"), ("Free MOAC", "MAOC")]
    lookup = {(r["TRT"], r["time"]): r for r in rows}
    for mtype, col in targets:
        obs, sims = [], {"R": [], "SAMM": [], "MONICA": []}
        for trt in TRTS:
            for day, val in meas.get((trt, mtype), []):
                t = max(day, 1)
                r = lookup.get((trt, t))
                if not r or r[f"{col}_SAMM"] == "" or r[f"{col}_MONICA"] == "":
                    continue
                obs.append(val)
                for k in sims:
                    sims[k].append(float(r[f"{col}_{k}"]))
        if not obs:
            continue
        out = f"  {mtype:12s} n={len(obs):3d}"
        for k, label in (("R", "R"), ("SAMM", "MONICA, SAMM mod."), ("MONICA", "MONICA, MONICA mod.")):
            rmse, ef, bias = stats(obs, sims[k])
            out += f" | {label}: {rmse:6.1f} / {ef:5.2f} / {bias:+6.1f}"
        print(out)

    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "comparison.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(out_dir, "measured_soc.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["TRT", "time", "SOC"])
        for trt in TRTS:
            for day, val in meas.get((trt, "SOC"), []):
                w.writerow([trt, day, val])


if __name__ == "__main__":
    main(*sys.argv[1:4])
