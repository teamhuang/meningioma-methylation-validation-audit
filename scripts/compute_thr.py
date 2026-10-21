import pandas as pd, numpy as np, os, json
OUT = r"E:\脑膜瘤npj\results"

# --- THR from cached reduced training matrix (same pipeline as harmonize_transfer.py) ---
m = pd.read_csv(os.path.join(OUT, "cache_m325_sel.csv.gz"), index_col=0)
print("cache shape:", m.shape, "| index head:", [str(x) for x in m.index[:3]],
      "| col head:", [str(c) for c in m.columns[:3]])

scaler = pd.read_csv(os.path.join(OUT, "methyl_scaler.csv"), index_col=0)
cent = pd.read_csv(os.path.join(OUT, "methyl_centroids.csv"), index_col=0)
print("scaler:", scaler.shape, "centroids:", cent.shape)

# orient: rows should be CpGs (20000). If not, transpose.
if not m.index[0].startswith("cg"):
    m = m.T
    print("transposed ->", m.shape)

common = [c for c in m.index if c in scaler.index and c in cent.index]
print("common CpGs:", len(common))

mu = scaler["mean"].reindex(common)
sd = scaler["sd"].reindex(common)
Z = m.loc[common].sub(mu, axis=0).div(sd, axis=0).T  # samples x CpGs
print("Z:", Z.shape, "NaN frac:", float(Z.isna().mean().mean()))

cc = [c for c in Z.columns if c in cent.index]
M = cent.loc[cc].values.T  # classes x CpGs
d = np.sqrt(((Z[cc].values[:, None, :] - M[None, :, :]) ** 2).sum(axis=2))
dmin = d.min(axis=1)
THR = float(np.percentile(dmin, 99))
print(f"TRAIN THR (P99 of min-distance) = {THR:.2f}")
print("train dist P50/P95:", np.percentile(dmin, 50).round(2), np.percentile(dmin, 95).round(2))
np.save(os.path.join(OUT, "methyl_train_thr.npy"), np.array([THR]))
print("saved methyl_train_thr.npy")
