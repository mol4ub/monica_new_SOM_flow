"""Evaluation statistics per data type (as Table 5 of the SAMM paper) from the observation file written by
`samm_calib ... eval params.csv daily.csv obs.csv`, plus microbial C and soil N (not used in calibration) from
daily.csv.  usage: python3 evaluate.py AllMeasurements.csv daily.csv obs.csv [label]"""
import csv, math, sys, statistics


def stats(o, s):
    n = len(o); mo = sum(o) / n
    sse = sum((a - b) ** 2 for a, b in zip(o, s)); sst = sum((a - mo) ** 2 for a in o)
    return n, math.sqrt(sse / n), 1 - sse / sst, sum(b - a for a, b in zip(o, s)) / n


def main(meas_path, daily_path, obs_path, label=""):
    rows = list(csv.DictReader(open(obs_path)))
    print(f"{label}\n{'type':12s} {'n':>4s} {'RMSE':>8s} {'EF':>6s} {'bias':>8s}")
    for ty in ("Litterbag C", "Microbial N", "SOC", "Aggregate C", "Free MOAC"):
        sel = [r for r in rows if r["type"] == ty]
        n, rmse, ef, bias = stats([float(r["obs"]) for r in sel], [float(r["sim"]) for r in sel])
        print(f"{ty:12s} {n:4d} {rmse:8.1f} {ef:6.2f} {bias:+8.1f}")
    daily = {(r["TRT"], int(r["time"])): r for r in csv.DictReader(open(daily_path))}
    meas = list(csv.DictReader(open(meas_path, encoding="utf-8-sig"), delimiter=";"))
    for ty, f in (("Microbial C", lambda r: float(r["MIC_C"])),
                  ("Soil N", lambda r: sum(float(r[k]) for k in ("LMW_N", "MIC_N", "MAO_N", "AggLAB_N", "AggMAO_N")))):
        sel = [m for m in meas if m["Measurement type"] == ty]
        o = [float(m["Value"]) for m in sel]
        s = [f(daily[(m["Treatment"], max(1, int(m["Day of simulation"])))]) for m in sel]
        n, rmse, ef, bias = stats(o, s)
        print(f"{ty:12s} {n:4d} {rmse:8.1f} {ef:6.2f} {bias:+8.1f}  (not used in calibration)")


if __name__ == "__main__":
    main(*sys.argv[1:])
