# A/B Test Report

**File:** `data/ab_sample_pre.csv`  |  **alpha:** 0.05  |  **power target:** 0.80

## Summary

- Users: A **2,530**, B **2,470**
- Conversion: A **0.0534**, B **0.0644**
- Uplift (B−A): **0.0110** abs  |  **20.64%** rel
- z = **1.65**, p = **0.0979**, 95% CI (diff) = **[-0.0020, 0.0241]**
- SRM check: **OK** (z=0.85)

## Variant table

| variant   |    n |   conversions |     cr |
|:----------|-----:|--------------:|-------:|
| A         | 2530 |           135 | 0.0534 |
| B         | 2470 |           159 | 0.0644 |

## CUPED (variance reduction)

- theta (A) = **0.388**, theta (B) = **0.421**
- Unadjusted: diff **0.0110**, SE **0.00666**, CI **[-0.0020, 0.0241]**
- CUPED     : diff **0.0110**, SE **0.00587**, CI **[-0.0005, 0.0225]**
- SE reduction from CUPED: **11.9%**

## Power / duration (estimates)

- Detect +20% relative lift: **7,615/variant**  (~1 day(s) at 20,000/day/variant)
- Detect +1pp absolute lift: **8,626/variant**  (~1 day(s) at 20,000/day/variant)
