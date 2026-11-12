"""Simulated vs measured for the measurement types of the Khon Kaen experiment: SAMM in MONICA (default MONICA
rate modifiers, sim-<TRT>-MONICA.json) and SAMM v1.1 in R for reference. Litterbag C is left out: litterbags are
a calibration device of the R code that MONICA does not simulate.

usage: python3 all_measurements.py R_out.csv runs_dir out_dir
  runs_dir: <TRT>-MONICA/out.csv of the MONICA runs; writes out_dir/all_series.csv (daily simulated values),
  out_dir/all_obs.csv (measurements with simulated values on the same day) and prints statistics.

Measurements are matched by simulation day, as in the calibration of the R model (day 1 = 1995-03-01 of the
climate file; the date column of AllMeasurements.csv is 15 days later).

Simulated quantities (0-15 cm, kg ha-1):
  SOC          MIC_C + MAO_C + AggSTR_C + AggLAB_C + AggMAO_C (as in the SAMM paper)
  Soil N       LMW_N + MIC_N + MAO_N + AggLAB_N + AggMAO_N (same pools as SOC; free litter excluded)
  Microbial C  MIC_C;  Microbial N  MIC_N
  Aggregate C  AggSTR_C + AggLAB_C + AggMAO_C;  Free MOAC  MAO_C
"""
import csv
import datetime as dt
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")
START = dt.date(1995, 3, 1)
TRTS = ["CT", "DP", "GN", "RS", "TM"]
KGHA = 0.15 * 10000
POOLS = ["STR_C", "LAB_C", "LMW_C", "MIC_C", "MAO_C", "AggSTR_C", "AggLAB_C", "AggMAO_C",
         "LAB_N", "LMW_N", "MIC_N", "MAO_N", "AggLAB_N", "AggMAO_N"]
TYPES = ["SOC", "Soil N", "Microbial C", "Microbial N", "Aggregate C", "Free MOAC"]


def quantities(p):
    return {
        "SOC": p["MIC_C"] + p["MAO_C"] + p["AggSTR_C"] + p["AggLAB_C"] + p["AggMAO_C"],
        "Soil N": p["LMW_N"] + p["MIC_N"] + p["MAO_N"] + p["AggLAB_N"] + p["AggMAO_N"],
        "Microbial C": p["MIC_C"], "Microbial N": p["MIC_N"],
        "Aggregate C": p["AggSTR_C"] + p["AggLAB_C"] + p["AggMAO_C"], "Free MOAC": p["MAO_C"],
    }


def read_monica(path):
    """pools by R time index (R row t = state at the start of day t; MONICA writes the state at the end of a day)"""
    lines = open(path).read().splitlines()
    k = lines.index('"daily"')
    hdr = lines[k + 1].split(",")
    pools = {}
    for line in lines[k + 3:]:
        f = line.split(",")
        if len(f) != len(hdr) or not f[0][:4].isdigit():
            continue
        row = dict(zip(hdr, f))
        day = (dt.date.fromisoformat(row["Date"]) - START).days + 1
        pools[day + 1] = {p: float(row[f"SAMM_{p}_1"]) * KGHA for p in POOLS}
    return pools


def read_r(path):
    return {(r["TRT"], int(float(r["time"]))): {k: float(r[k]) for k in POOLS} for r in csv.DictReader(open(path))}


def stats(obs, sim):
    n = len(obs)
    mo = sum(obs) / n
    sse = sum((o - s) ** 2 for o, s in zip(obs, sim))
    sst = sum((o - mo) ** 2 for o in obs)
    return math.sqrt(sse / n), (1 - sse / sst) if sst > 0 else float("nan"), sum(s - o for o, s in zip(obs, sim)) / n


def main(r_path, runs, out_dir):
    R = read_r(r_path)
    meas = [m for m in csv.DictReader(open(os.path.join(DATA, "AllMeasurements.csv"), encoding="utf-8-sig"), delimiter=";")
            if m["Measurement type"] in TYPES]
    os.makedirs(out_dir, exist_ok=True)
    series, obs_rows = [], []
    for trt in TRTS:
        pools = read_monica(os.path.join(runs, f"{trt}-MONICA", "out.csv"))
        for t in range(2, 9071):
            if t not in pools:
                continue
            q, qr = quantities(pools[t]), quantities(R[(trt, t)])
            date = (START + dt.timedelta(days=t - 1)).isoformat()
            for ty in TYPES:
                series.append({"TRT": trt, "type": ty, "time": t, "date": date, "MONICA": q[ty], "R": qr[ty]})
        for m in meas:
            if m["Treatment"] != trt:
                continue
            t = max(int(m["Day of simulation"]), 2)
            ty = m["Measurement type"]
            obs_rows.append({"TRT": trt, "type": ty, "time": t, "date": (START + dt.timedelta(days=t - 1)).isoformat(),
                             "obs": float(m["Value"]), "sd": math.sqrt(float(m["Variance"])),
                             "MONICA": quantities(pools[t])[ty] if t in pools else None,
                             "R": quantities(R[(trt, t)])[ty]})
    for name, rows in (("all_series.csv", series), ("all_obs.csv", obs_rows)):
        with open(os.path.join(out_dir, name), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    print("Simulated vs measured: n, RMSE (kg/ha), Nash-Sutcliffe EF, mean bias (simulated - measured, kg/ha)")
    print(f"{'type':12s} {'trt':4s} | {'MONICA (blue)':>30s} | {'R':>30s}")
    for ty in TYPES:
        for trt in TRTS + ["all"]:
            sel = [r for r in obs_rows if r["type"] == ty and (trt == "all" or r["TRT"] == trt)]
            out = f"{ty:12s} {trt:4s}"
            for model in ("MONICA", "R"):
                rows = [r for r in sel if r[model] is not None]
                st = stats([r["obs"] for r in rows], [r[model] for r in rows])
                out += f" | n={len(rows):3d} {st[0]:8.1f} {st[1]:7.2f} {st[2]:+8.1f}"
            print(out)


if __name__ == "__main__":
    main(*sys.argv[1:4])
