# A/B Testing & Experimentation Toolkit

Tiny toolkit to evaluate A/B tests from a CSV, reduce variance with CUPED, and plan sample size/duration. Includes a Streamlit UI and CLI.

## Run in 60 seconds

```bash
pip install -r requirements.txt

# Evaluate a CSV (SRM + uplift + CI + bootstrap)
python -m src.abtest.quick_eval --file data/ab_sample_pre.csv

# Power / duration (e.g., +20% relative lift at α=0.05, power=0.80)
python -m src.abtest.power_planner --baseline_cr 0.05 --mde_rel 0.2 --daily_users 40000 --split 0.5

# One-page Markdown report -> reports/ab_baseline.md
python -m src.abtest.report --file data/ab_sample_pre.csv

# Streamlit app
streamlit run app/quick_eval_app.py
```
## What it does
- SRM check (z-test on allocation A vs B)
- Uplift (two-proportion z-test) with 95% CI and bootstrap CI
- CUPED variance reduction (per-variant θ, adjusted SE & CI)
- Power planner (required n/variant for relative or absolute lift)

## CSV schema
- variant — "A" or "B"
- converted — 0/1 primary metric
- pre_exposure — (optional) numeric pre-period covariate for CUPED

## Repo layout

```bash
app/                  # Streamlit UI
scripts/              # helpers (sample data, etc.)
src/abtest/           # quick_eval, CUPED, power planner, report
data/                 # sample CSVs (kept small)
reports/              # generated markdown reports
requirements.txt
```

## Notes for implementers
- SRM uses a binomial z-test on traffic split (target 50/50).
- CUPED fits θ per variant to avoid bias; reduces SE when pre-metric correlates with conversion.
- Methods are educational; validate against your company’s experimentation platform before production use.

If you like it, save and push:

```bash
git add README.md
git commit -m "Docs: clean Run-in-60s README"
git push
