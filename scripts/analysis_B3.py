import os, gzip, warnings, time
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")
t0 = time.time()
DATA = r"E:\脑膜瘤npj\data"
OUT = r"E:\脑膜瘤npj\results"
META = os.path.join(DATA, "meta")
os.makedirs(META, exist_ok=True)
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

# ================= GSE189521 (110 patients, endpoints) =================
p189 = os.path.join(DATA, "GSE189521_series_matrix.txt.gz")
log("parsing GSE189521 series matrix header ...")
geo_acc, titles, chars = [], [], []
with gzip.open(p189, "rt", encoding="utf-8", errors="replace") as f:
    start = None
    for i, line in enumerate(f):
        if line.startswith("!Sample_geo_accession"):
            geo_acc = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_title"):
            titles = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1"):
            chars.append([x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]])
        elif line.startswith("!series_matrix_table_begin"):
            start = i + 1
            break
tab = pd.DataFrame({"GSM": geo_acc, "title": titles})
for i, c in enumerate(chars):
    tab[f"char{i}"] = c
tab.to_csv(os.path.join(META, "GSE189521_sample_map.csv"), index=False)
log(f"sample map: n={len(tab)}; chars={len(chars)}; title example={titles[:2]}")

log("loading beta matrix ...")
B189 = pd.read_csv(p189, sep="\t", index_col=0, skiprows=start, compression="gzip")
B189 = B189[~B189.index.astype(str).str.startswith("!")]
B189.columns = [str(c).strip().strip('"') for c in B189.columns]
B189.index = [str(x).split("_")[0] for x in B189.index]
B189 = B189[~B189.index.duplicated(keep="first")]
log(f"189521 beta matrix {B189.shape}")

B189z, common189 = zscore_with_artifacts(B189)
log(f"189521: {len(common189)}/{len(sel)} classifier CpGs available")
pred189 = classify_centroid(B189z)
cls189 = pd.DataFrame({"GSM": B189z.index,
                       "title": [tab.set_index("GSM").loc[g, "title"] if g in tab.set_index("GSM").index else "" for g in B189z.index],
                       "predicted_methyl_class": pred189})
cls189.to_csv(os.path.join(OUT, "gse189521_classification.csv"), index=False)
log("189521 predicted class distribution: " + pd.Series(pred189).value_counts().to_dict().__str__())
log("ANALYSIS_B3_DONE")




