"""Sampling importance resampling (SIR) as in the SAMM paper (Laub et al. 2024, section 2.4): importance weights
proportional to the likelihood (standardised by the mean likelihood), resampling of the prior samples without
replacement with the weights as sampling probabilities.

usage: python3 sir_resample.py samples.csv out_prefix [n_resample=200] [seed=1]
  samples.csv: written by `samm_calib ... sir prior.csv samples.csv` (parameters..., logLik)
  writes <out_prefix>_posterior.csv (resampled sets), <out_prefix>_best.csv (highest likelihood),
  <out_prefix>_summary.csv and prints the effective sample size."""
import csv, math, random, sys


def main(path, prefix, n_res=200, seed=1):
    n_res, rnd = int(n_res), random.Random(int(seed))
    rows = list(csv.reader(open(path)))
    names, data = rows[0], [[float(x) for x in r] for r in rows[1:]]
    ll = [r[-1] for r in data]
    valid = [i for i, l in enumerate(ll) if math.isfinite(l)]
    m = max(ll[i] for i in valid)
    w = {i: math.exp(ll[i] - m) for i in valid}           # proportional to the likelihood
    tot = sum(w.values())
    ess = tot ** 2 / sum(v * v for v in w.values())        # effective sample size of the importance weights
    print(f"{len(data)} samples, {len(valid)} valid runs; max logLik {m:.2f}; effective sample size {ess:.2f}")
    top = sorted(w.items(), key=lambda kv: -kv[1])[:10]
    print("largest normalised weights:", ", ".join(f"{v / tot:.3g}" for _, v in top))
    # resampling without replacement, weights as sampling probabilities (sequential draws). The weights are
    # recomputed in log space relative to the best remaining sample at every draw: with hundreds of data points
    # most weights underflow to 0 relative to the overall best, and the remaining probabilities would be lost.
    pool = sorted(valid, key=lambda i: -ll[i]); chosen = []
    for _ in range(min(n_res, len(pool))):
        mx = ll[pool[0]]
        ws = [math.exp(ll[i] - mx) for i in pool[:2000]]  # weights beyond the 2000 best are < exp(-745) anyway
        u = rnd.random() * sum(ws); acc = 0.0
        for j, v in enumerate(ws):
            acc += v
            if acc >= u: break
        chosen.append(pool.pop(j))
    eff_chosen = sum(1 for i in chosen if w[i] / tot > 1e-6)
    print(f"resampled {len(chosen)} sets; {eff_chosen} of them carry a normalised weight > 1e-6")
    with open(f"{prefix}_posterior.csv", "w", newline="") as f:
        cw = csv.writer(f); cw.writerow(names + ["weight"])
        for i in chosen: cw.writerow(data[i] + [w[i] / tot])
    best = max(valid, key=lambda i: ll[i])
    with open(f"{prefix}_best.csv", "w") as f:
        f.write("name,value\n")
        for n, v in zip(names[:-1], data[best][:-1]): f.write(f"{n},{v:.10g}\n")
        f.write(f"logLik,{ll[best]:.6f}\n")
    with open(f"{prefix}_summary.csv", "w") as f:
        f.write("name,best,weighted_mean,weighted_sd,q025,q50,q975\n")
        for k, n in enumerate(names[:-1]):
            xs = sorted(data[i][k] for i in chosen)
            wm = sum(w[i] * data[i][k] for i in valid) / tot
            wsd = math.sqrt(max(0.0, sum(w[i] * (data[i][k] - wm) ** 2 for i in valid) / tot))
            q = lambda p: xs[min(len(xs) - 1, int(p * len(xs)))]
            f.write(f"{n},{data[best][k]:.6g},{wm:.6g},{wsd:.3g},{q(0.025):.6g},{q(0.5):.6g},{q(0.975):.6g}\n")


if __name__ == "__main__":
    main(*sys.argv[1:])
