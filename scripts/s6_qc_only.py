# -*- coding: utf-8 -*-
"""S6 only: QC-threshold sensitivity on GSE270371 (994 samples, distribution-level).

Streaming read of intensity columns; keeps the 20k selected CpGs; frozen scaler
normalization; per-CpG missingness thresholds (0.2/0.4/0.6); class distribution,
distance stats and reject rate at the frozen training P99 threshold.
"""
import os, gzip, warnings
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")

OUT = r"E:\脑膜瘤npj\results"
DATA = r"E:\脑膜瘤npj\data"
sel = pd.read_csv(os.path.join(OUT, "methyl_selected_cpgs.csv"))["cpg"].tolist()
scaler = pd.read_csv(os.path.join(OUT, "methyl_scaler.csv"), index_col=0)
cent = pd.read_csv(os.path.join(OUT, "methyl_centroids.csv"), index_col=0)
sel_set = set(sel)
THR = float(np.load(os.path.join(OUT, "methyl_train_thr.npy"))[0])
print("reject THR (training P99):", round(THR, 2))

p = os.path.join(DATA, "GSE270371_Matrix_processed.txt.gz")
with gzip.open(p, "rt") as f:
    hdr = f.readline().rstrip("\n").split("\t")
keep_idx = [i for i, h in enumerate(hdr) if h != "ID_REF" and not h.startswith("Detection Pval")]
keep_names = [hdr[i] for i in keep_idx]

names, rows = [], []
with gzip.open(p, "rt") as f:
    f.readline()
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if parts[0] in sel_set:
            names.append(parts[0])
            rows.append([parts[i] for i in keep_idx])

B = pd.DataFrame(np.array(rows, dtype=np.float32), index=names, columns=keep_names)
print("kept CpGs:", B.shape)

mu = scaler["mean"].reindex(B.index)
sd = scaler["sd"].reindex(B.index)
Z = B.sub(mu, axis=0).div(sd, axis=0).T  # samples x CpGs

# centroid predict on shared CpGs
cc = [c for c in Z.columns if c in cent.index]
M = cent.loc[cc].values.T  # classes x CpGs
d = np.sqrt(((Z[cc].values[:, None, :] - M[None, :, :]) ** 2).sum(axis=2))
pred = cent.columns[np.argmin(d, axis=1)]
dmin = d.min(axis=1)

# per-CpG missingness thresholds
res = []
for thr in (0.2, 0.4, 0.6):
    na_cpg = B.isna().mean(axis=1)
    kept = na_cpg[na_cpg < thr].index
    cc2 = [c for c in kept if c in cent.index]
    M2 = cent.loc[cc2].values.T
    d2 = np.sqrt(((Z[cc2].values[:, None, :] - M2[None, :, :]) ** 2).sum(axis=2))
    pred2 = cent.columns[np.argmin(d2, axis=1)]
    dmin2 = d2.min(axis=1)
    res.append({
        "check": f"S6_na_thr_{thr}",
        "n_cpg_kept": len(cc2),
        "class_dist": pd.Series(pred2).value_counts().to_dict(),
        "dist_P50": round(float(np.nanmedian(dmin2)), 1),
        "dist_P95": round(float(np.nanpercentile(dmin2, 95)), 1),
        "reject_rate": round(float(np.nanmean(dmin2 > THR)), 3),
    })
    print(res[-1])

base = {
    "check": "S6_na_thr_baseline_full_panel",
    "n_cpg_kept": len(cc),
    "class_dist": pd.Series(pred).value_counts().to_dict(),
    "dist_P50": round(float(np.nanmedian(dmin)), 1),
    "dist_P95": round(float(np.nanpercentile(dmin, 95)), 1),
    "reject_rate": round(float(np.nanmean(dmin > THR)), 3),
}
print(base)
rows_out = [base] + res
pd.DataFrame(rows_out).to_csv(os.path.join(OUT, "s6_qc_threshold.csv"), index=False)
print("SAVED s6_qc_threshold.csv")
