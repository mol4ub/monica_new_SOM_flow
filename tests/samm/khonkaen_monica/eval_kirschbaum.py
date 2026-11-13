"""Measured vs simulated for the MONICA runs of variant Kirschbaum (sim-<TRT>-Kirschbaum.json), with the same
definitions as the calibration of SAMMKirschbaum (tests/samm/calibration/samm_calib.cpp, simValue):
  SOC          LMW_C + MIC_C + MAO_C + AggSTR_C + AggLAB_C + AggMAO_C (free litter STR_C, LAB_C excluded)
  Soil N       LMW_N + MIC_N + MAO_N + AggLAB_N + AggMAO_N
  Microbial C  MIC_C;  Microbial N as reported: MIC_N * 3.1 / fN_ninhydrin (run 7)
  Aggregate C  AggSTR_C + AggLAB_C + AggMAO_C;  Free MOAC  MAO_C
Litterbag C is not simulated in MONICA.

usage: python3 eval_kirschbaum.py runs_dir out_obs.csv
  runs_dir/<TRT>/out.csv: MONICA output of sim-<TRT>-Kirschbaum.json; writes TRT,day,type,obs,sim (format of
  samm_calib eval, for tests/samm/calibration/plot_eval.R) and prints RMSE / EF / bias per data type"""
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from all_measurements import DATA, TRTS, read_monica, stats

RUN7 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "calibration", "results", "dream7")
FN = {r.split(",")[0]: float(r.split(",")[1]) for r in open(os.path.join(RUN7, "best.csv")).read().split()[1:]}["fN_ninhydrin"]
TYPES = ["Microbial N", "SOC", "Aggregate C", "Free MOAC", "Microbial C", "Soil N"]


def sim_value(p, ty):
    return {"SOC": p["LMW_C"] + p["MIC_C"] + p["MAO_C"] + p["AggSTR_C"] + p["AggLAB_C"] + p["AggMAO_C"],
            "Soil N": p["LMW_N"] + p["MIC_N"] + p["MAO_N"] + p["AggLAB_N"] + p["AggMAO_N"],
            "Microbial C": p["MIC_C"], "Microbial N": p["MIC_N"] * 3.1 / FN,
            "Aggregate C": p["AggSTR_C"] + p["AggLAB_C"] + p["AggMAO_C"], "Free MOAC": p["MAO_C"]}[ty]


def main(runs, out):
    meas = [m for m in csv.DictReader(open(os.path.join(DATA, "AllMeasurements.csv"), encoding="utf-8-sig"), delimiter=";")
            if m["Measurement type"] in TYPES]
    rows = []
    for trt in TRTS:
        pools = read_monica(os.path.join(runs, trt, "out.csv"))
        for m in meas:
            if m["Treatment"] != trt:
                continue
            day = int(m["Day of simulation"])
            t = max(day, 2)  # pools[t]: state at the start of day t (end of day t - 1)
            if t in pools:
                rows.append({"TRT": trt, "day": day, "type": m["Measurement type"], "obs": float(m["Value"]),
                             "sim": sim_value(pools[t], m["Measurement type"])})
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["TRT", "day", "type", "obs", "sim"])
        w.writeheader()
        w.writerows(rows)
    print(f"{'type':12s} {'n':>4s} {'RMSE':>8s} {'EF':>6s} {'bias':>8s}")
    for ty in TYPES:
        sel = [r for r in rows if r["type"] == ty]
        rmse, ef, bias = stats([r["obs"] for r in sel], [r["sim"] for r in sel])
        print(f"{ty:12s} {len(sel):4d} {rmse:8.1f} {ef:6.2f} {bias:+8.1f}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
