# app/quick_eval_app.py
import math
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st
from scipy.stats import norm

TARGET = "converted"          # 0/1
VARIANT = "variant"           # 'A' / 'B'
PRE = "pre_exposure"          # optional float column for CUPED

# ---------- helpers ----------
def read_csv_any(src) -> pd.DataFrame:
    df = pd.read_csv(src)
    # coerce
    df[VARIANT] = df[VARIANT].astype(str)
    df[TARGET] = pd.to_numeric(df[TARGET], errors="coerce").fillna(0).astype(int)
    if PRE in df.columns:
        df[PRE] = pd.to_numeric(df[PRE], errors="coerce")
    return df[[c for c in [ "user_id", VARIANT, TARGET, PRE] if c in df.columns]]

def counts(df):
    g = df.groupby(VARIANT, observed=True)[TARGET].agg(["count", "sum"])
    nA = int(g.loc["A","count"]);  nB = int(g.loc["B","count"])
    cA = int(g.loc["A","sum"]);    cB = int(g.loc["B","sum"])
    crA = cA / nA if nA else 0.0
    crB = cB / nB if nB else 0.0
    return nA, nB, cA, cB, crA, crB

def srm_z(nA, nB, expected=0.5):
    N = nA + nB
    if N == 0:
        return 0.0
    p_hat = nA / N
    se = math.sqrt(expected * (1 - expected) / N)
    return (p_hat - expected) / se if se > 0 else 0.0

def diff_ci_ztest(crA, crB, nA, nB, alpha=0.05):
    diff = crB - crA
    p_pool = (crA*nA + crB*nB) / (nA + nB)
    se = math.sqrt(p_pool * (1 - p_pool) * (1/nA + 1/nB))
    z = diff / se if se > 0 else 0.0
    zcrit = norm.ppf(1 - alpha/2)
    lo, hi = diff - zcrit*se, diff + zcrit*se
    p = 2 * (1 - norm.cdf(abs(z)))
    return diff, se, (lo, hi), z, p

def bootstrap_diff(yA, yB, iters=2000, seed=0):
    rng = np.random.default_rng(seed)
    nA, nB = len(yA), len(yB)
    diffs = np.empty(iters, dtype=float)
    for i in range(iters):
        mA = yA[rng.integers(0, nA, nA)].mean()
        mB = yB[rng.integers(0, nB, nB)].mean()
        diffs[i] = mB - mA
    return diffs.mean(), np.quantile(diffs, [0.025, 0.975])

def cuped_adjust(y, x):
    # y: 0/1 outcomes, x: pre-exposure (float)
    # return adjusted y', theta, variance
    if x is None or x.isna().all():
        return y.to_numpy().astype(float), 0.0, float(y.var(ddof=1))
    x_c = x - x.mean()
    cov = np.cov(x_c, y, ddof=1)[0,1]
    varx = x_c.var(ddof=1)
    theta = 0.0 if varx == 0 else cov / varx
    y_adj = y - theta * x_c
    return y_adj.to_numpy().astype(float), float(theta), float(y_adj.var(ddof=1))

def cuped_se_diff(dfA, dfB):
    yA_adj, thetaA, varA = cuped_adjust(dfA[TARGET], dfA[PRE] if PRE in dfA else None)
    yB_adj, thetaB, varB = cuped_adjust(dfB[TARGET], dfB[PRE] if PRE in dfB else None)
    mA, mB = yA_adj.mean(), yB_adj.mean()
    nA, nB = len(yA_adj), len(yB_adj)
    se = math.sqrt(varA/nA + varB/nB)
    return (mB - mA), se, thetaA, thetaB

def samplesize_2prop(pA, pB, alpha=0.05, power=0.8):
    # normal approx (two-sided), equal allocation
    z_alpha = norm.ppf(1 - alpha/2)
    z_beta = norm.ppf(power)
    pbar = 0.5 * (pA + pB)
    num = (z_alpha * math.sqrt(2*pbar*(1-pbar)) + z_beta * math.sqrt(pA*(1-pA) + pB*(1-pB)))**2
    den = (pB - pA)**2
    return math.ceil(num / den)

# ---------- UI ----------
st.set_page_config(page_title="A/B Testing — Quick Evaluator", layout="wide")
st.title("A/B Testing — Quick Evaluator")

# default file (if present)
default_path = Path("data/ab_sample_pre.csv")
uploaded = st.file_uploader("Upload CSV (columns: variant, converted[, pre_exposure])", type=["csv"])
if uploaded is not None:
    df = read_csv_any(uploaded)
elif default_path.exists():
    df = read_csv_any(default_path)
else:
    st.info("Upload a CSV to begin.")
    st.stop()

# counts
assert set(df[VARIANT].unique()) >= {"A","B"}, "CSV must contain variants A and B"
nA, nB, cA, cB, crA, crB = counts(df)
yA = df.loc[df[VARIANT]=="A", TARGET].to_numpy()
yB = df.loc[df[VARIANT]=="B", TARGET].to_numpy()
z_srm = srm_z(nA, nB)
srm_ok = abs(z_srm) < 3.0

# headline
st.metric("Users A", f"{nA:,}")
st.metric("Users B", f"{nB:,}")
st.metric("CR A / B", f"{crA:0.4f} / {crB:0.4f}")
diff, se, (lo, hi), z, p = diff_ci_ztest(crA, crB, nA, nB)
st.metric("Uplift (B-A)", f"{diff:0.4f} ({(diff/crA):.2%} rel)")
st.write(f"z = **{z:0.2f}**, p = **{p:0.4f}**, 95% CI **[{lo:0.4f}, {hi:0.4f}]**")

st.write(f"SRM: **{'OK' if srm_ok else 'POSSIBLE MISMATCH'}** (z={z_srm:0.2f})")

# CUPED
st.subheader("CUPED")
if PRE in df.columns:
    dfA = df[df[VARIANT]=="A"]
    dfB = df[df[VARIANT]=="B"]
    diff_u, se_u, (lo_u, hi_u), _, _ = diff_ci_ztest(crA, crB, nA, nB)
    diff_c, se_c, thetaA, thetaB = cuped_se_diff(dfA, dfB)
    zc = diff_c / se_c if se_c > 0 else 0.0
    zcrit = norm.ppf(0.975)
    lo_c, hi_c = diff_c - zcrit*se_c, diff_c + zcrit*se_c
    se_red = (1 - se_c/se_u) * 100 if se_u > 0 else 0.0
    st.write(f"Unadj SE {se_u:0.5f} → **CUPED SE {se_c:0.5f}** (↓ {se_red:0.1f}%)")
    st.write(f"CUPED θA={thetaA:0.3f}, θB={thetaB:0.3f}")
    st.write(f"CUPED CI **[{lo_c:0.4f}, {hi_c:0.4f}]**")
else:
    st.info("No `pre_exposure` column found — CUPED skipped.")

# Power / Duration
st.subheader("Power / Duration")
alpha = st.slider("alpha", 0.01, 0.10, 0.05, 0.01)
power = st.slider("power", 0.50, 0.99, 0.80, 0.01)
daily_users = st.number_input("daily users (total A+B)", min_value=1000, value=40000, step=1000)
rel = st.slider("target lift (relative %)", 0.01, 0.50, 0.20, 0.01)
abs_pp = st.slider("target lift (absolute pp)", 0.001, 0.05, 0.01, 0.001)

pA = crA if crA > 0 else 0.05
pB_rel = pA * (1 + rel)
pB_abs = pA + abs_pp

def show_req(name, pB):
    n = samplesize_2prop(pA, pB, alpha=alpha, power=power)
    per_day = daily_users / 2.0  # 50/50 split
    days = math.ceil(n / per_day) if per_day > 0 else float("inf")
    st.write(f"Detect {name} → **{n:,}/variant** (~{days} day(s) at {int(per_day):,}/day/variant)")

show_req(f"+{rel:.0%} relative", pB_rel)
show_req(f"+{abs_pp:.3f}pp absolute", pB_abs)
