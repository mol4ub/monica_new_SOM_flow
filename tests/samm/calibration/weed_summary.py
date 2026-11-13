import csv, sys, statistics, datetime as dt
rows = list(csv.DictReader(open(sys.argv[1])))
START = dt.date(1995, 3, 1)
add = {1, 382, 759, 1125, 1490, 1874, 2231, 2596, 2961, 3323, 3688, 4049, 4414, 4779, 5144, 5509, 5874, 6240, 6605, 6970, 7335, 7702, 8068, 8432, 8797, 9162}
def weeding(t):
    d = START + dt.timedelta(days=t - 1)
    return (d.day == 15 and (4 <= d.month <= 9 or d.month in (11, 1, 3))) or (t in add and t > 1)
print(f"{'trt':4s} {'weed AGB at weeding (t DM/ha): rainy':38s} {'dry':>6s} {'root C input (kg C/ha/d)':>26s} {'mineral N mean':>15s} {'SOC end':>8s}")
res = {}
for t in ("CT", "DP", "GN", "RS", "TM"):
    r = [x for x in rows if x["TRT"] == t]
    rainy, dry, rootC = [], [], 0
    for x in r:
        tt = int(x["time"])
        if weeding(tt):
            m = (START + dt.timedelta(days=tt - 1)).month
            (rainy if 4 <= m <= 9 else dry).append(float(x["WeedAGB"]) / 1000)
            rootC += float(x["WeedRoot"]) * 0.42
    soc = sum(float(r[-1][k]) for k in ("LMW_C", "MIC_C", "MAO_C", "AggSTR_C", "AggLAB_C", "AggMAO_C"))
    res[t] = statistics.mean(rainy)
    print(f"{t:4s} {statistics.mean(rainy):8.2f} (max {max(rainy):4.2f}){'':19s} {statistics.mean(dry):6.2f} {rootC / 9070:20.2f} {statistics.mean(float(x['MinN']) for x in r):15.1f} {soc:8.0f}")
print(f"ratio of rainy-season weed biomass GN / CT: {res['GN'] / res['CT']:.1f}")
