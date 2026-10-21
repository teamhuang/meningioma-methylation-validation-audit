import os, gzip, warnings, time
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")
t0 = time.time()
DATA = r"E:\脑膜瘤npj\data"
OUT = r"E:\脑膜瘤npj\results"
def log(msg): print(f"[{time.time()-t0:6.0f}s] {msg}", flush=True)
scaler = pd.read_csv(os.path.join(OUT, "methyl_scaler.csv"), index_col=0)
centroids = pd.read_csv(os.path.join(OUT, "methyl_centroids.csv"), index_col=0)
sel = pd.read_csv(os.path.join(OUT, "methyl_selected_cpgs.csv"))["cpg"].tolist()
def classify_centroid(Xz_new):
    labs = list(centroids.columns)
    common = [c for c in Xz_new.columns if c in centroids.index]
    M = centroids.loc[common, labs].values.T
    Xv = Xz_new[common].values
    d = ((Xv[:, None, :] - M[None, :, :]) ** 2).sum(axis=2)
    return np.array([labs[i] for i in d.argmin(axis=1)])
def zscore_with_artifacts(bet):
    idx = scaler.index.intersection(bet.index)
    bet = bet.loc[idx]
    mu = scaler["mean"].reindex(idx); sd = scaler["sd"].reindex(idx)
    z = (bet.T - mu).div(sd)
    return z, idx

# ================= GSE270371 (994, streaming 8.3GB gz) =================
log("applying to GSE270371 (streaming) ...")
sel_base = set(str(x).split("_")[0] for x in sel)
path270 = os.path.join(DATA, "GSE270371_Matrix_processed.txt.gz")
rows = {}
with gzip.open(path270, "rt") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    log(f"270371 header cols: {len(hdr)}")
    beta_cols = [i for i, h in enumerate(hdr) if h != "ID_REF" and h != "Detection Pval"]
    beta_names = [hdr[i] for i in beta_cols]
    log(f"beta sample columns: {len(beta_cols)} (first {beta_names[:4]})")
    for line in f:
        parts = line.rstrip("\n").split("\t")
        cid = parts[0].split("_")[0]
        if cid in sel_base:
            vals = []
            ok = True
            for i in beta_cols:
                try:
                    vals.append(float(parts[i]))
                except Exception:
                    ok = False
                    break
            if ok:
                rows[cid] = vals
log(f"stream done: collected {len(rows)} of {len(sel)} selected CpGs")
B270 = pd.DataFrame.from_dict(rows, orient="index", columns=beta_names)
B270 = B270.loc[[c for c in B270.index if c in set(sel)]]
B270 = B270.T
log(f"B270 matrix {B270.shape}")
B270z = (B270 - scaler["mean"].reindex(B270.columns)) / scaler["sd"].reindex(B270.columns)
pred270 = classify_centroid(B270z)
log(f"GSE270371 (n={len(pred270)}) predicted class distribution: {pd.Series(pred270).value_counts().to_dict()}")
pd.DataFrame({"sample": B270z.index, "predicted_class": pred270}).to_csv(
    os.path.join(OUT, "gse270371_classification.csv"), index=False)
log("ANALYSIS_B2_DONE")



