# src/abtest/power_planner.py
from __future__ import annotations
import argparse, math
import numpy as np
from scipy.stats import norm

def sample_size_2prop(p1: float, p2: float, alpha: float = 0.05, power: float = 0.8) -> int:
    """Two-sided normal-approx sample size per group for two proportions."""
    p1 = float(np.clip(p1, 1e-9, 1-1e-9))
    p2 = float(np.clip(p2, 1e-9, 1-1e-9))
    delta = abs(p2 - p1)
    if delta <= 0:
        return 0
    z_alpha = norm.ppf(1 - alpha/2)       # two-sided
    z_beta  = norm.ppf(power)
    pbar = 0.5*(p1 + p2)
    term1 = z_alpha * np.sqrt(2 * pbar * (1 - pbar))
    term2 = z_beta  * np.sqrt(p1 * (1 - p1) + p2 * (1 - p2))
    n = ((term1 + term2)**2) / (delta**2)
    return int(math.ceil(n))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline_cr", type=float, required=True, help="e.g. 0.05")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--mde_abs",  type=float, help="absolute lift, e.g. 0.01 for +1pp")
    g.add_argument("--mde_rel",  type=float, help="relative lift, e.g. 0.2 for +20%")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--power", type=float, default=0.8)
    ap.add_argument("--daily_users", type=int, default=40000, help="total users/day across both variants")
    ap.add_argument("--split", type=float, default=0.5, help="traffic share per variant (A and B)")
    args = ap.parse_args()

    p1 = args.baseline_cr
    p2 = p1 + args.mde_abs if args.mde_abs is not None else p1 * (1 + args.mde_rel)
    p2 = float(np.clip(p2, 1e-9, 1-1e-9))

    n_per_variant = sample_size_2prop(p1, p2, alpha=args.alpha, power=args.power)
    per_variant_per_day = max(1, int(args.daily_users * args.split))
    days = math.ceil(n_per_variant / per_variant_per_day)

    mde_txt = f"+{args.mde_abs:.4f} abs" if args.mde_abs is not None else f"+{args.mde_rel:.1%} rel"
    print(f"Baseline CR = {p1:.4f}, target = {p2:.4f} ({mde_txt})")
    print(f"alpha={args.alpha:.3f}, power={args.power:.2f}")
    print(f"Required per variant: {n_per_variant:,}")
    print(f"With {per_variant_per_day:,}/day per variant → ~{days} day(s)")

if __name__ == "__main__":
    main()
