"""Evaluation statistics per data type (as Table 5 of the SAMM paper) from the observation file written by
`samm_calib ... eval params.csv daily.csv obs.csv` (all seven data types; SOC includes LMW-C, soil N is the N of the
SOC pools).  usage: python3 evaluate.py AllMeasurements.csv daily.csv obs.csv [label]"""
import csv, math, sys, statistics


def stats(o, s):
    n = len(o); mo = sum(o) / n
    sse = sum((a - b) ** 2 for a, b in zip(o, s)); sst = sum((a - mo) ** 2 for a in o)
    return n, math.sqrt(sse / n), 1 - sse / sst, sum(b - a for a, b in zip(o, s)) / n


def main(meas_path, daily_path, obs_path, label=""):
    rows = list(csv.DictReader(open(obs_path)))
    used = ("Litterbag C", "Microbial N", "SOC", "Aggregate C", "Free MOAC")
    print(f"{label}\n{'type':12s} {'n':>4s} {'RMSE':>8s} {'EF':>6s} {'bias':>8s}")
    present = {r["type"] for r in rows}
    for ty in used + ("Microbial C", "Soil N"):
        if ty in present:
            sel = [r for r in rows if r["type"] == ty]
            n, rmse, ef, bias = stats([float(r["obs"]) for r in sel], [float(r["sim"]) for r in sel])
        else:  # older observation files: computed from the daily output
            daily = {(r["TRT"], int(r["time"])): r for r in csv.DictReader(open(daily_path))}
            meas = [m for m in csv.DictReader(open(meas_path, encoding="utf-8-sig"), delimiter=";") if m["Measurement type"] == ty]
            f = (lambda r: float(r["MIC_C"])) if ty == "Microbial C" else \
                (lambda r: sum(float(r[k]) for k in ("LMW_N", "MIC_N", "MAO_N", "AggLAB_N", "AggMAO_N")))
            n, rmse, ef, bias = stats([float(m["Value"]) for m in meas],
                                      [f(daily[(m["Treatment"], max(1, int(m["Day of simulation"])))]) for m in meas])
        print(f"{ty:12s} {n:4d} {rmse:8.1f} {ef:6.2f} {bias:+8.1f}")


if __name__ == "__main__":
    main(*sys.argv[1:])
