# scripts/self_check_quick_eval.py
import numpy as np, pandas as pd
from src.abtest.quick_eval import srm_check, conv_rate, two_prop_z, bootstrap_uplift

df = pd.read_csv("data/ab_sample.csv")
A = df[df["variant"].str.upper().eq("A")]["converted"].to_numpy(dtype=float)
B = df[df["variant"].str.upper().eq("B")]["converted"].to_numpy(dtype=float)

print(f"n_A={len(A)}, n_B={len(B)}")
print(f"cr_A={A.mean():.4f}, cr_B={B.mean():.4f}")

# SRM check
srm_fail, z_srm = srm_check(len(A), len(B))
print(f"SRM: {'FAIL' if srm_fail else 'ok'} (z={z_srm:.2f})")

# Two-proportion Z and 95% CI (pooled SE)
diff, z, ci = two_prop_z(A.sum(), len(A), B.sum(), len(B))
print(f"uplift(B-A)={diff:.4f}  z={z:.2f}  95% CI={ci[0]:.4f}..{ci[1]:.4f}")

# Bootstrap CI
boot_mean, boot_ci = bootstrap_uplift(A, B, iters=2000, seed=0)
print(f"bootstrap mean={boot_mean:.4f}  95% CI={boot_ci[0]:.4f}..{boot_ci[1]:.4f}")

# Loose sanity checks for your synthetic file
assert 2500 <= len(A) <= 2550 and 2450 <= len(B) <= 2495
assert abs(A.mean() - 0.0534) < 0.005
assert abs(B.mean() - 0.0644) < 0.005
assert abs(z_srm) < 2.6           # SRM should be ok
assert ci[0] < 0 < ci[1]          # CI should cross 0 for this sample
print("self-check passed")
