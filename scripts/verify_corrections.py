# -*- coding: utf-8 -*-
"""Verify & archive the numbers corrected during manuscript v2 review (2026-10-18).

Regenerates four traceability artifacts:
  results/transfer_dist_summary.csv     SeSAMe raw-transfer centroid distances (overall)
  results/transfer_dist_by_pred.csv     SeSAMe distances by predicted class
  results/train_dist_quantiles.csv      training-set centroid-distance quantiles (20k z-space)
  results/meta_heterogeneity.json       Cochran Q / I2 / P_het recomputed from meta_high_risk_HR.csv
Also prints the GSE304094 complete-case grade strata (149 = 69+75+5).
"""
import os, json
import pandas as pd
import numpy as np
from scipy.stats import chi2

OUT = r"E:\脑膜瘤npj\results"

# 1) SeSAMe raw-transfer distance summary (corrects the "9,324-12,675" range)
t = pd.read_csv(os.path.join(OUT, "natgen201_transfer_diagnosis.csv"))
q = t["dist"]
pd.DataFrame({"stat": ["min", "P5", "median", "P95", "max"],
              "euclidean": [q.min(), q.quantile(.05), q.median(), q.quantile(.95), q.max()]}
             ).round(0).to_csv(os.path.join(OUT, "transfer_dist_summary.csv"), index=False)
t.groupby("pred")["dist"].agg(["min", "median", "max"]).round(0)\
 .to_csv(os.path.join(OUT, "transfer_dist_by_pred.csv"))
print("SeSAMe dist: min %.0f P5 %.0f median %.0f P95 %.0f max %.0f (n=%d)" %
      (q.min(), q.quantile(.05), q.median(), q.quantile(.95), q.max(), len(t)))

# 2) training-set distance quantiles (confirms THR P99 = 178.03 euclidean / 31,699.97 squared)
m = pd.read_csv(os.path.join(OUT, "cache_m325_sel.csv.gz"), index_col=0)
if not str(m.index[0]).startswith("cg"):
    m = m.T
scaler = pd.read_csv(os.path.join(OUT, "methyl_scaler.csv"), index_col=0)
cent = pd.read_csv(os.path.join(OUT, "methyl_centroids.csv"), index_col=0)
cc = [c for c in m.index if c in scaler.index and c in cent.index]
Z = m.loc[cc].sub(scaler["mean"].reindex(cc), axis=0).div(scaler["sd"].reindex(cc), axis=0).T
M = cent.loc[cc].values.T
d2 = ((Z[cc].values[:, None, :] - M[None, :, :]) ** 2).sum(axis=2)
dmin2, dmin = d2.min(axis=1), np.sqrt(d2.min(axis=1))
rows = [{"space": s, "n": len(v),
         "P50": round(float(np.percentile(v, 50)), 1),
         "P95": round(float(np.percentile(v, 95)), 1),
         "P99": round(float(np.percentile(v, 99)), 2),
         "max": round(float(v.max()), 2)}
        for s, v in [("euclidean", dmin), ("squared", dmin2)]]
pd.DataFrame(rows).to_csv(os.path.join(OUT, "train_dist_quantiles.csv"), index=False)
print(pd.DataFrame(rows).to_string(index=False))

# 3) meta heterogeneity (confirms Q=0.835, I2=0%, P_het=0.66)
mh = pd.read_csv(os.path.join(OUT, "meta_high_risk_HR.csv"))
k = mh[mh["cohort"] != "Pooled (fixed effect)"]
w = 1 / k["se"] ** 2
theta = (w * k["logHR"]).sum() / w.sum()
Q = float((w * (k["logHR"] - theta) ** 2).sum())
df = len(k) - 1
het = {"Q": round(Q, 3), "df": df,
       "I2_percent": round(max(0.0, (Q - df) / Q) * 100, 1),
       "p_heterogeneity": float(chi2.sf(Q, df)),
       "method": "fixed-effect Q on logHR, weights=1/se^2, recomputed from meta_high_risk_HR.csv"}
with open(os.path.join(OUT, "meta_heterogeneity.json"), "w") as fh:
    json.dump(het, fh, indent=1)
print(json.dumps(het, indent=1))

# 4) GSE304094 complete-case grade strata (confirms 149 = 69+75+5)
d = pd.read_csv(r"E:\脑膜瘤npj\results\gse304094_primary_survival.csv")
ccg = d.dropna(subset=["T", "E", "grade", "gtr", "age", "male"])
print("complete-case n=%d, grade x (n, events):" % len(ccg))
print(ccg.groupby("grade")["E"].agg(["size", "sum"]).to_string())
print("VERIFY_CORRECTIONS_DONE")
