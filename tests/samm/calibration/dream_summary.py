"""Summary of a DREAM run of samm_calib: posterior statistics from the second half of the chains (burn-in = first
half), Gelman-Rubin R-hat per parameter, acceptance proxy and the highest-posterior / highest-likelihood sets.

usage: python3 dream_summary.py chains.csv out_prefix
  writes <out_prefix>_summary.csv, <out_prefix>_best.csv (highest likelihood) and <out_prefix>_map.csv
  (highest posterior), <out_prefix>_posterior.csv (the second-half samples)"""
import csv, math, statistics, sys


def main(path, prefix):
    rows = list(csv.DictReader(open(path)))
    names = [k for k in rows[0] if k not in ("gen", "chain", "logLik", "logPost")]
    gmax = max(int(r["gen"]) for r in rows)
    post = [r for r in rows if int(r["gen"]) > gmax / 2]
    chains = sorted({r["chain"] for r in rows})
    print(f"{gmax} generations, {len(chains)} chains; posterior from generation {gmax // 2 + 1}: {len(post)} samples")
    with open(f"{prefix}_summary.csv", "w") as f:
        f.write("name,mean,sd,q025,q50,q975,rhat\n")
        print(f"{'parameter':18s} {'mean':>10s} {'sd':>9s} {'2.5%':>10s} {'50%':>10s} {'97.5%':>10s} {'R-hat':>6s}")
        for n in names:
            xs = sorted(float(r[n]) for r in post)
            by = [[float(r[n]) for r in post if r["chain"] == c] for c in chains]
            m = min(len(b) for b in by); by = [b[:m] for b in by]
            means = [statistics.fmean(b) for b in by]; W = statistics.fmean(statistics.variance(b) for b in by)
            B = m * statistics.variance(means) if len(by) > 1 else 0
            rhat = math.sqrt(((m - 1) / m * W + B / m) / W) if W > 0 else float("nan")
            q = lambda p: xs[min(len(xs) - 1, int(p * len(xs)))]
            mu, sd = statistics.fmean(xs), statistics.pstdev(xs)
            f.write(f"{n},{mu:.6g},{sd:.3g},{q(0.025):.6g},{q(0.5):.6g},{q(0.975):.6g},{rhat:.3f}\n")
            print(f"{n:18s} {mu:10.4g} {sd:9.3g} {q(0.025):10.4g} {q(0.5):10.4g} {q(0.975):10.4g} {rhat:6.3f}")
    ll = [float(r["logLik"]) for r in post]
    print(f"logLik of the posterior samples: median {statistics.median(ll):.1f}, range {min(ll):.1f} to {max(ll):.1f}")
    for key, suffix in (("logLik", "best"), ("logPost", "map")):
        b = max(rows, key=lambda r: float(r[key]))
        with open(f"{prefix}_{suffix}.csv", "w") as f:
            f.write("name,value\n")
            for n in names: f.write(f"{n},{b[n]}\n")
            f.write(f"logLik,{b['logLik']}\n")
        print(f"highest {key}: {float(b[key]):.2f} (logLik {float(b['logLik']):.2f}, generation {b['gen']}, chain {b['chain']})")
    with open(f"{prefix}_posterior.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(post)


if __name__ == "__main__":
    main(*sys.argv[1:3])
