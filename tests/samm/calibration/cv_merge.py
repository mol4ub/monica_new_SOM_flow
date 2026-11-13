"""Merge the leave-one-residue-out folds: held-out predictions of each fold (cv_obs.csv, cv_band.csv) and the
run 7 in-sample values for the same treatments (ins_obs.csv, ins_band.csv); print RMSE / EF per type and treatment.
usage: python3 cv_merge.py calib_dir [fold_dir (default calib_dir/cv)]
  calib_dir holds dream7_obs.csv and dream7_band.csv of run 7; fold_dir the lo_<TRT>_* files of the folds"""
import csv, math, os, sys
C = sys.argv[1]; CV = sys.argv[2] if len(sys.argv) > 2 else os.path.join(C, "cv"); TRTS = ["DP", "GN", "RS", "TM"]

def rows(path): return list(csv.DictReader(open(path)))
def write(path, rs):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rs[0].keys())); w.writeheader(); w.writerows(rs)

cv_obs, cv_band = [], []
for t in TRTS:
    cv_obs += [r for r in rows(f"{CV}/lo_{t}_obs.csv") if r["TRT"] == t]
    cv_band += [r for r in rows(f"{CV}/lo_{t}_band.csv") if r["TRT"] == t]
write(f"{CV}/cv_obs.csv", cv_obs); write(f"{CV}/cv_band.csv", cv_band)
ins_obs = [r for r in rows(f"{C}/dream7_obs.csv") if r["TRT"] in TRTS]
write(f"{CV}/ins_obs.csv", ins_obs)
write(f"{CV}/ins_band.csv", [r for r in rows(f"{C}/dream7_band.csv") if r["TRT"] in TRTS])

def st(rs):
    o = [float(r["obs"]) for r in rs]; s = [float(r["sim"]) for r in rs]; n = len(o); mo = sum(o) / n
    sse = sum((a - b) ** 2 for a, b in zip(o, s)); sst = sum((a - mo) ** 2 for a in o)
    return math.sqrt(sse / n), 1 - sse / sst if sst > 0 else float("nan"), sum(b - a for a, b in zip(o, s)) / n
types = ["Litterbag C", "Microbial N", "SOC", "Aggregate C", "Free MOAC", "Microbial C", "Soil N"]
print(f"{'type':12s} {'trt':4s} {'n':>4s} | {'in-sample RMSE':>14s} {'EF':>6s} {'bias':>8s} | {'held-out RMSE':>13s} {'EF':>6s} {'bias':>8s}")
for ty in types:
    for t in TRTS + ["all"]:
        a = [r for r in ins_obs if r["type"] == ty and (t == "all" or r["TRT"] == t)]
        b = [r for r in cv_obs if r["type"] == ty and (t == "all" or r["TRT"] == t)]
        if not a: continue
        x, y = st(a), st(b)
        print(f"{ty:12s} {t:4s} {len(a):4d} | {x[0]:14.1f} {x[1]:6.2f} {x[2]:+8.1f} | {y[0]:13.1f} {y[1]:6.2f} {y[2]:+8.1f}")
