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
