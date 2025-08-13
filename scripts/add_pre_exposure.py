# scripts/add_pre_exposure.py
import numpy as np, pandas as pd

df = pd.read_csv("data/ab_sample.csv")
rng = np.random.default_rng(0)

# make a pre-exposure covariate correlated with conversion
# (higher before-experiment activity -> slightly more likely to convert)
noise = rng.normal(0, 0.3, size=len(df))
df["pre_exposure"] = (0.2 + 0.6*df["converted"] + noise).clip(lower=0)

out = "data/ab_sample_pre.csv"
df.to_csv(out, index=False)
print(f"wrote {out} with column 'pre_exposure'")
