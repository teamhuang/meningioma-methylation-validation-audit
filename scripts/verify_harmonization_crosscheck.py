# -*- coding: utf-8 -*-
"""Cross-checks for the UCSF harmonization claims (manuscript v2 review round 2).

1) Same-reference comparison: v11 paths vs the ORIGINAL diagnose k-means
   (9,483-CpG listwise-deleted in-space reference). Shows the 0.545 -> 0.84/0.86
   gap is reference-driven, not path-driven (v11 paths score ~0.53 vs that
   reference too) -> cross-reference comparisons are not like-for-like.
2) Biological gradient: v11 predicted classes vs MOESM20 clinical (MIB1-LI,
   Recurrent) attached positionally to MOESM21 M-IDs (same validation logic as
   analysis_report 9.1). Shows the Hypermitotic-like class is recovered by all
   three paths.

Outputs: results/harmonization_crosscheck.csv (both blocks, long format).
"""
import os
import pandas as pd
import numpy as np
from sklearn.metrics import adjusted_rand_score
from scipy.optimize import linear_sum_assignment

OUT = r"E:\脑膜瘤npj\results"
NG = r"E:\脑膜瘤npj\data\suppl\natgen2022"

m21 = pd.read_excel(os.path.join(NG, "41588_2022_1061_MOESM21_ESM.xlsx"), index_col=0)
m21.columns = [str(c).strip() for c in m21.columns]
mids = list(m21.columns)
dg = pd.read_csv(os.path.join(OUT, "natgen201_transfer_diagnosis.csv"))
assert len(dg) == len(mids), "diagnosis rows must align to MOESM21 columns"
dg["MID"] = mids
pu = pd.read_csv(os.path.join(OUT, "harmonization_pred_ucsfsesame.csv"),
                 dtype={"pred_P1": str, "pred_P2": str, "pred_P3": str})
m = dg.merge(pu, on="MID")
assert len(m) == 200

def bestmap(t, p):
    ct = pd.crosstab(pd.Series(p).astype(str), pd.Series(t).astype(str))
    ri, ci = linear_sum_assignment(-ct.values)
    return ct.values[ri, ci].sum() / ct.values.sum(), str(dict(zip(ct.index[ri], ct.columns[ci])))

rows = []
t_diag = dg["km"].astype(str)
for c in ["pred_P1", "pred_P2", "pred_P3"]:
    ari = adjusted_rand_score(t_diag, m[c].astype(str))
    b, mp = bestmap(t_diag, m[c].astype(str))
    rows.append({"block": "same_reference_vs_diagnose_km", "path": c, "n": len(m),
                 "ARI": round(float(ari), 3), "bestmap_acc": round(float(b), 3), "map": mp})
b, mp = bestmap(t_diag, dg["pred"].astype(str))
rows.append({"block": "same_reference_vs_diagnose_km", "path": "old_direct", "n": len(m),
             "ARI": round(float(adjusted_rand_score(t_diag, dg["pred"].astype(str))), 3),
             "bestmap_acc": round(float(b), 3), "map": mp})

m20 = pd.read_excel(os.path.join(NG, "41588_2022_1061_MOESM20_ESM.xlsx"))
n = min(len(m20), len(m21.columns))
c20 = m20.iloc[:n].copy(); c20["MID"] = mids[:n]
mm = c20.merge(pu, on="MID")
mm["mib1"] = pd.to_numeric(mm["MIB1-LI"], errors="coerce")
mm["rec"] = pd.to_numeric(mm["Recurrent"], errors="coerce")
for c in ["km_ref", "pred_P1", "pred_P2", "pred_P3"]:
    g = mm.groupby(c).agg(n=("mib1", "size"), MIB1=("mib1", "mean"), rec_rate=("rec", "mean")).round(2)
    for cls, r in g.iterrows():
        rows.append({"block": "biological_gradient", "path": f"{c}={cls}", "n": int(r["n"]),
                     "ARI": np.nan, "bestmap_acc": np.nan,
                     "map": f"MIB1={r['MIB1']}, rec_rate={r['rec_rate']}"})

pd.DataFrame(rows).to_csv(os.path.join(OUT, "harmonization_crosscheck.csv"), index=False)
for r in rows:
    print(r)
print("HARMONIZATION_CROSSCHECK_DONE")
