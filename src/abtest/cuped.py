# src/abtest/cuped.py
import argparse
import numpy as np, pandas as pd

def cuped_theta(y: np.ndarray, x: np.ndarray) -> float:
    """theta = Cov(y, x) / Var(x)"""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    vx = x.var(ddof=0)
    if vx == 0:
        return 0.0
    cov = ((x - x.mean()) * (y - y.mean())).mean()
    return float(cov / vx)

def cuped_adjust(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    """y* = y - theta * (x - mean(x))"""
    theta = cuped_theta(y, x)
    return y - theta * (x - x.mean())

def unadjusted_diff_and_ci(a: np.ndarray, b: np.ndarray):
    """Standard two-proportion diff with 95% CI (normal approx)."""
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    p1, p2 = a.mean(), b.mean()
    n1, n2 = len(a), len(b)
    se = np.sqrt(p1*(1-p1)/max(1,n1) + p2*(1-p2)/max(1,n2))
    diff = p2 - p1
    ci = (diff - 1.96*se, diff + 1.96*se)
    return diff, se, ci

def adjusted_diff_and_ci(a_star: np.ndarray, b_star: np.ndarray):
    """Diff of means on adjusted metric with 95% CI."""
    a_star = np.asarray(a_star, dtype=float); b_star = np.asarray(b_star, dtype=float)
    n1, n2 = len(a_star), len(b_star)
    m1, m2 = a_star.mean(), b_star.mean()
    v1, v2 = a_star.var(ddof=1), b_star.var(ddof=1)
    se = np.sqrt(v1/max(1,n1) + v2/max(1,n2))
    diff = m2 - m1
    ci = (diff - 1.96*se, diff + 1.96*se)
    return diff, se, ci

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", required=True, help="CSV with columns: variant, converted, pre_exposure")
    ap.add_argument("--variant_col", default="variant")
    ap.add_argument("--metric_col", default="converted")
    ap.add_argument("--covariate_col", default="pre_exposure")
    args = ap.parse_args()

    df = pd.read_csv(args.file)
    A = df[df[args.variant_col].str.upper().eq("A")]
    B = df[df[args.variant_col].str.upper().eq("B")]

    yA = A[args.metric_col].to_numpy(dtype=float)
    yB = B[args.metric_col].to_numpy(dtype=float)
    xA = A[args.covariate_col].to_numpy(dtype=float)
    xB = B[args.covariate_col].to_numpy(dtype=float)

    # unadjusted (binary proportion)
    diff_u, se_u, ci_u = unadjusted_diff_and_ci(yA, yB)

    # CUPED-adjusted
    yA_star = cuped_adjust(yA, xA)
    yB_star = cuped_adjust(yB, xB)
    diff_c, se_c, ci_c = adjusted_diff_and_ci(yA_star, yB_star)

    reduction = 0.0 if se_u == 0 else (1 - se_c/se_u) * 100

    print(f"Unadjusted: diff={diff_u:.4f}  SE={se_u:.5f}  95% CI={ci_u[0]:.4f}..{ci_u[1]:.4f}")
    print(f"CUPED     : diff={diff_c:.4f}  SE={se_c:.5f}  95% CI={ci_c[0]:.4f}..{ci_c[1]:.4f}")
    print(f"SE reduction from CUPED: {reduction:.1f}%")

if __name__ == "__main__":
    main()
