# src/abtest/quick_eval.py
import argparse, pandas as pd, numpy as np
from math import sqrt

def srm_check(a_n, b_n, tol=0.02):
    """Simple SRM: A and B should be ~50/50. Flag if z>2.58 (~1% level)."""
    total = a_n + b_n
    exp = total / 2
    z = (a_n - exp) / sqrt(total * 0.5 * 0.5)
    return abs(z) > 2.58, z

def conv_rate(x):
    return float(x.mean()) if len(x) else 0.0

def two_prop_z(a_s, a_n, b_s, b_n):
    """Z-test for difference in two proportions, + 95% CI (normal approx)."""
    p1, p2 = a_s/a_n, b_s/b_n
    p = (a_s + b_s) / (a_n + b_n)
    # variance for Bernoulli p is p(1-p)
    var = p * (1 - p)
    se = sqrt(var * (1/a_n + 1/b_n)) if 0 < p < 1 else 0.0
    if se == 0.0:
        diff = p2 - p1
        return diff, 0.0, (diff, diff)
    z = (p2 - p1) / se
    diff = p2 - p1
    ci = (diff - 1.96*se, diff + 1.96*se)
    return diff, z, ci

def bootstrap_uplift(a, b, iters=5000, seed=42):
    """Non-parametric CI for uplift (B-A)."""
    rng = np.random.default_rng(seed)
    a = np.asarray(a); b = np.asarray(b)
    nA, nB = len(a), len(b)
    diffs = np.empty(iters, dtype=float)
    for i in range(iters):
        diffs[i] = b[rng.integers(0, nB, nB)].mean() - a[rng.integers(0, nA, nA)].mean()
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return float(diffs.mean()), (float(lo), float(hi))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", required=True, help="CSV with columns: user_id, variant (A/B), converted (0/1)")
    ap.add_argument("--variant_col", default="variant")
    ap.add_argument("--metric_col", default="converted")
    args = ap.parse_args()

    df = pd.read_csv(args.file)
    df[args.metric_col] = df[args.metric_col].astype(float)

    A = df[df[args.variant_col].str.upper().eq("A")][args.metric_col]
    B = df[df[args.variant_col].str.upper().eq("B")][args.metric_col]

    a_n, b_n = len(A), len(B)
    a_s, b_s = A.sum(), B.sum()

    print(f"n_A={a_n}, cr_A={conv_rate(A):.4f} | n_B={b_n}, cr_B={conv_rate(B):.4f}")

    srm_fail, z_srm = srm_check(a_n, b_n)
    print(f"SRM check: {'FAIL' if srm_fail else 'ok'} (z={z_srm:.2f})")

    diff, z, ci = two_prop_z(a_s, a_n, b_s, b_n)
    print(f"uplift (B-A): {diff:.4f}  z={z:.2f}  95% CI={ci[0]:.4f}..{ci[1]:.4f}")

    boot_mean, boot_ci = bootstrap_uplift(A.values, B.values)
    print(f"bootstrap uplift: mean={boot_mean:.4f}  95% CI={boot_ci[0]:.4f}..{boot_ci[1]:.4f}")

if __name__ == "__main__":
    main()
