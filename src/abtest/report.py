# src/abtest/report.py
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import norm

from src.abtest.quick_eval import srm_check, two_prop_z
from src.abtest.cuped import cuped_theta, cuped_adjust, adjusted_diff_and_ci, unadjusted_diff_and_ci
from src.abtest.power_planner import sample_size_2prop

def to_md_table(df: pd.DataFrame) -> str:
    try:
        return df.to_markdown(index=False)
    except Exception:
        return "```\n" + df.to_csv(index=False) + "\n```"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", default="data/ab_sample_pre.csv",
                    help="CSV with columns: user_id, variant(A/B), converted(0/1) [, pre_exposure] [, date/timestamp]")
    ap.add_argument("--out", default="reports/ab_baseline.md")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--power", type=float, default=0.80)
    ap.add_argument("--daily_users", type=int, default=40000, help="total traffic per day (A+B), for duration calc")
    args = ap.parse_args()

    df = pd.read_csv(args.file)
    df["variant"] = df["variant"].astype(str).str.upper()
    if "date" in df.columns:
        df["date"] = pd.to_datetime(df["date"])
    elif "timestamp" in df.columns:
        df["date"] = pd.to_datetime(df["timestamp"])
    else:
        df["date"] = pd.NaT

    A = df[df["variant"].eq("A")]
    B = df[df["variant"].eq("B")]

    nA, nB = len(A), len(B)
    pA, pB = A["converted"].mean(), B["converted"].mean()
    diff, z, ci = two_prop_z(A["converted"].sum(), nA, B["converted"].sum(), nB)
    pval = 2 * (1 - norm.cdf(abs(z)))
    rel = (pB/pA - 1.0) if pA > 0 else np.nan

    srm_fail, z_srm = srm_check(nA, nB)

    rows = pd.DataFrame([
        {"variant":"A", "n":nA, "conversions":int(A["converted"].sum()), "cr":pA},
        {"variant":"B", "n":nB, "conversions":int(B["converted"].sum()), "cr":pB},
    ])

    md = []
    md += ["# A/B Test Report", ""]
    md += [f"**File:** `{args.file}`  |  **alpha:** {args.alpha:.2f}  |  **power target:** {args.power:.2f}", ""]
    md += ["## Summary", ""]
    md += [f"- Users: A **{nA:,}**, B **{nB:,}**",
           f"- Conversion: A **{pA:.4f}**, B **{pB:.4f}**",
           f"- Uplift (B−A): **{diff:.4f}** abs  |  **{rel*100:.2f}%** rel" if pA>0 else f"- Uplift (B−A): **{diff:.4f}**",
           f"- z = **{z:.2f}**, p = **{pval:.4f}**, 95% CI (diff) = **[{ci[0]:.4f}, {ci[1]:.4f}]**",
           f"- SRM check: **{'FAIL' if srm_fail else 'OK'}** (z={z_srm:.2f})",
           ""]
    md += ["## Variant table", "", to_md_table(rows.assign(cr=rows.cr.map(lambda x: f"{x:.4f}"))), ""]

    # CUPED (if covariate exists)
    if "pre_exposure" in df.columns:
        thetaA = cuped_theta(A["converted"].to_numpy(), A["pre_exposure"].to_numpy())
        thetaB = cuped_theta(B["converted"].to_numpy(), B["pre_exposure"].to_numpy())
        yA_star = cuped_adjust(A["converted"].to_numpy(), A["pre_exposure"].to_numpy())
        yB_star = cuped_adjust(B["converted"].to_numpy(), B["pre_exposure"].to_numpy())
        diff_u, se_u, ci_u = unadjusted_diff_and_ci(A["converted"].to_numpy(), B["converted"].to_numpy())
        diff_c, se_c, ci_c = adjusted_diff_and_ci(yA_star, yB_star)
        red = 0.0 if se_u == 0 else (1 - se_c/se_u)*100

        md += ["## CUPED (variance reduction)", ""]
        md += [f"- theta (A) = **{thetaA:.3f}**, theta (B) = **{thetaB:.3f}**",
               f"- Unadjusted: diff **{diff_u:.4f}**, SE **{se_u:.5f}**, CI **[{ci_u[0]:.4f}, {ci_u[1]:.4f}]**",
               f"- CUPED     : diff **{diff_c:.4f}**, SE **{se_c:.5f}**, CI **[{ci_c[0]:.4f}, {ci_c[1]:.4f}]**",
               f"- SE reduction from CUPED: **{red:.1f}%**", ""]

    # Power suggestions (two common asks)
    n_rel20 = sample_size_2prop(pA, pA*(1+0.20), alpha=args.alpha, power=args.power)
    n_abs1p = sample_size_2prop(pA, pA+0.01,       alpha=args.alpha, power=args.power)
    per_variant_per_day = max(1, int(args.daily_users * 0.5))
    days_rel20 = int(np.ceil(n_rel20 / per_variant_per_day))
    days_abs1p = int(np.ceil(n_abs1p / per_variant_per_day))

    md += ["## Power / duration (estimates)", "",
           f"- Detect +20% relative lift: **{n_rel20:,}/variant**  (~{days_rel20} day(s) at {per_variant_per_day:,}/day/variant)",
           f"- Detect +1pp absolute lift: **{n_abs1p:,}/variant**  (~{days_abs1p} day(s) at {per_variant_per_day:,}/day/variant)",
           ""]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(md), encoding="utf-8")
    print(f"[report] wrote {out}")

if __name__ == "__main__":
    main()
