"""Compare SAMM pool time series of two runs (e.g. the R model and the C++ core), per treatment and pool.

usage: python3 compare_to_r.py R_out.csv cpp_out.csv
Both files need the columns time, the 14 SAMM pools and TRT (the R output of Samm_main.R has them).
"""
import csv
import sys

POOLS = ["STR_C", "LAB_C", "LMW_C", "MIC_C", "MAO_C", "AggSTR_C", "AggLAB_C", "AggMAO_C",
         "LAB_N", "LMW_N", "MIC_N", "MAO_N", "AggLAB_N", "AggMAO_N"]


def read(path):
    data = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            data[(row["TRT"], int(float(row["time"])))] = {p: float(row[p]) for p in POOLS}
    return data


def main(a_path, b_path):
    a, b = read(a_path), read(b_path)
    keys = sorted(set(a) & set(b))
    print(f"{len(keys)} common (treatment, day) rows")
    worst = 0.0
    for trt in sorted({k[0] for k in keys}):
        print(f"\n{trt}: pool | max abs diff (kg/ha) | max rel diff | at day | value R / C++")
        for p in POOLS:
            mx, mr, at = 0.0, 0.0, None
            for k in keys:
                if k[0] != trt:
                    continue
                va, vb = a[k][p], b[k][p]
                d = abs(va - vb)
                r = d / max(abs(va), 1e-9)
                if d > mx:
                    mx, at = d, k[1]
                mr = max(mr, r if abs(va) > 1e-6 else 0)
            worst = max(worst, mr)
            print(f"  {p:9s} {mx:12.3e} {mr:12.3e} {at!s:>6} {a[(trt, at)][p] if at else 0:.6g} / {b[(trt, at)][p] if at else 0:.6g}")
        soc = lambda d, k: sum(d[k][p] for p in ["MIC_C", "MAO_C", "AggSTR_C", "AggLAB_C", "AggMAO_C"])
        end = max(k[1] for k in keys if k[0] == trt)
        print(f"  SOC (MIC+MAO+Agg) at day {end}: R {soc(a, (trt, end)):.4f}  C++ {soc(b, (trt, end)):.4f}")
    print(f"\nlargest relative difference of any pool: {worst:.3e}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
