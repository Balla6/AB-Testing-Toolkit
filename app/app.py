# app/app.py
import pandas as pd, numpy as np, streamlit as st
from scipy.stats import norm
from src.abtest.quick_eval import srm_check, two_prop_z
from src.abtest.cuped import cuped_adjust, adjusted_diff_and_ci, unadjusted_diff_and_ci
from src.abtest.power_planner import sample_size_2prop

st.set_page_config(page_title="A/B Evaluator", layout="wide")
st.title("A/B Testing — Quick Evaluator")

uploaded = st.file_uploader("Upload CSV (columns: variant, converted[, pre_exposure])", type=["csv"])
df = pd.read_csv(uploaded) if uploaded else pd.read_csv("data/ab_sample_pre.csv")
df["variant"] = df["variant"].astype(str).str.upper()

A, B = df[df["variant"].eq("A")], df[df["variant"].eq("B")]
nA, nB = len(A), len(B)
pA, pB = A["converted"].mean(), B["converted"].mean()
diff, z, ci = two_prop_z(A["converted"].sum(), nA, B["converted"].sum(), nB)
pval = 2*(1-norm.cdf(abs(z)))
rel = (pB/pA - 1.0) if pA>0 else np.nan
srm_fail, z_srm = srm_check(nA, nB)

c1,c2,c3,c4 = st.columns(4)
c1.metric("Users A", f"{nA:,}")
c2.metric("Users B", f"{nB:,}")
c3.metric("CR A / B", f"{pA:.4f} / {pB:.4f}")
c4.metric("Uplift (B-A)", f"{diff:.4f}" + (f" ({rel*100:.2f}% rel)" if pA>0 else ""))

st.write(f"z = **{z:.2f}**, p = **{pval:.4f}**, 95% CI **[{ci[0]:.4f}, {ci[1]:.4f}]**")
st.write(f"SRM: **{'FAIL' if srm_fail else 'OK'}** (z={z_srm:.2f})")

if "pre_exposure" in df.columns:
    yA, yB = A["converted"].to_numpy(), B["converted"].to_numpy()
    diff_u, se_u, ci_u = unadjusted_diff_and_ci(yA, yB)
    yA_star = cuped_adjust(yA, A["pre_exposure"].to_numpy())
    yB_star = cuped_adjust(yB, B["pre_exposure"].to_numpy())
    diff_c, se_c, ci_c = adjusted_diff_and_ci(yA_star, yB_star)
    red = 0.0 if se_u == 0 else (1 - se_c/se_u)*100
    st.markdown("### CUPED")
    st.write(f"Unadj SE **{se_u:.5f}** → CUPED SE **{se_c:.5f}** (↓ {red:.1f}%)")
    st.write(f"CUPED CI **[{ci_c[0]:.4f}, {ci_c[1]:.4f}]**")

st.markdown("---")
st.markdown("### Power / Duration")
alpha = st.slider("alpha", 0.01, 0.10, 0.05, 0.01)
power = st.slider("power", 0.5, 0.95, 0.80, 0.05)
daily = st.number_input("daily users (total A+B)", value=40000, step=1000)
col1, col2 = st.columns(2)
with col1:
    rel_lift = st.slider("target lift (relative %)", 0.01, 0.50, 0.20, 0.01)
    n_rel = sample_size_2prop(pA, pA*(1+rel_lift), alpha=alpha, power=power)
    st.write(f"Detect +{rel_lift*100:.0f}% → **{n_rel:,}/variant** (~{int(np.ceil(n_rel/(daily/2)))} days)")
with col2:
    abs_pp = st.slider("target lift (absolute pp)", 0.001, 0.05, 0.010, 0.001)
    n_abs = sample_size_2prop(pA, pA+abs_pp, alpha=alpha, power=power)
    st.write(f"Detect +{abs_pp:.3f}pp → **{n_abs:,}/variant** (~{int(np.ceil(n_abs/(daily/2)))} days)")
