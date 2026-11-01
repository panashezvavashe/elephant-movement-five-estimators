"""Table 1 (data quality) and wet/dry centroid shift on the v6 screened data.

Uses the fixes that enter combined2.R after Section B (speed_km_h <= 5, 4557 excluded,
3905 cut at 2023-08-08). Season-year centroid shift = distance between the centroid of all
wet-season fixes and all dry-season fixes (pooled over years, as in the original Table 5),
plus the mean of per-season-year shifts (wet W_Y vs the following dry D_Y+1 is not used;
W_Y vs D_Y, the dry season preceding it).

    python table1_centroid_v6.py
"""
import numpy as np
import pandas as pd

d = pd.read_csv("data/processed/all_elephants_screened_v6.csv")
d["ts"] = pd.to_datetime(d.timestamp, format="%m/%d/%Y %H:%M")
d = d[(d.speed_km_h <= 5) & (d.elephant_id != 4557)]
d = d[~((d.elephant_id == 3905) & (d.ts >= "2023-08-08"))].sort_values(["elephant_id", "ts"])
wet = d.ts.dt.month.isin([11, 12, 1, 2, 3, 4])
d["season"] = np.where(wet, "wet", "dry")
# season-year label: wet season Nov Y - Apr Y+1 -> Y; dry May-Oct Y -> Y
d["sy"] = np.where(wet & (d.ts.dt.month <= 4), d.ts.dt.year - 1, d.ts.dt.year)

rows = []
for i, g in d.groupby("elephant_id"):
    gap = g.ts.diff().dt.total_seconds() / 3600
    span = (g.ts.max() - g.ts.min()).days
    cw = g[g.season == "wet"][["x_utm", "y_utm"]].mean()
    cd = g[g.season == "dry"][["x_utm", "y_utm"]].mean()
    shifts = []
    for y, gy in g.groupby("sy"):
        a, b = gy[gy.season == "wet"], gy[gy.season == "dry"]
        if len(a) >= 100 and len(b) >= 100:
            shifts.append(np.hypot(*(a[["x_utm", "y_utm"]].mean() - b[["x_utm", "y_utm"]].mean())) / 1000)
    rows.append(dict(
        ID=i, first_fix=g.ts.min().strftime("%d %b %Y"), last_fix=g.ts.max().strftime("%d %b %Y"),
        months=round(span / 30.44, 1), n_fixes=len(g), med_interval_h=round(gap.median(), 2),
        max_gap_h=round(gap.max(), 1), gaps_gt24h=int((gap > 24).sum()),
        pct_wet=round(100 * (g.season == "wet").mean(), 1), pct_dry=round(100 * (g.season == "dry").mean(), 1),
        centroid_shift_pooled_km=round(np.hypot(*(cw - cd)) / 1000, 2),
        centroid_shift_sy_mean_km=round(np.mean(shifts), 2) if shifts else np.nan, n_sy_pairs=len(shifts)))
t = pd.DataFrame(rows)
t.to_csv("final/data_quality/tables/Table1_centroid_v6.csv", index=False)
print(t.to_string(index=False))
print("total fixes:", t.n_fixes.sum())
