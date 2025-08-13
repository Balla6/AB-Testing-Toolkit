# A/B Testing & Experimentation Toolkit

Tiny toolkit to **evaluate A/B tests** quickly from a CSV, **reduce variance** with CUPED, and **plan sample size / duration**. Comes with a Streamlit UI and CLI.

## What’s inside
- **Quick evaluator**: SRM check, conversion rates, uplift, z-test + CIs
- **CUPED variance reduction**: per-variant θ, adjusted SE & CI
- **Power planner**: required n per variant for relative or absolute lift
- **CLI + Streamlit app**
- **Sample data**: `data/ab_sample_pre.csv`

## CSV schema
- `variant` — e.g., `"A"` / `"B"`
- `converted` — `0/1`
- `pre_exposure` *(optional)* — continuous pre-period metric for CUPED

## Quickstart
```bash
# 1) create & activate env (optional)
python -m venv .venv
# Windows
.venv\Scripts\activate
# macOS/Linux
# source .venv/bin/activate

# 2) install deps
pip install -r requirements.txt
