import os, gzip, warnings, time
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.model_selection import StratifiedKFold

t0 = time.time()
DATA = r"E:\脑膜瘤npj\data"
OUT = r"E:\脑膜瘤npj\results"
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20260317)
TOP_K = 20000

def log(msg):
    print(f"[{time.time()-t0:6.0f}s] {msg}", flush=True)

# ================= helpers =================
def load_processed(path, sample_names=None):
    df = pd.read_csv(path, index_col=0, sep="\t")
    df.columns = [str(c).strip() for c in df.columns]
    return df.astype(np.float32)

def read_sample_map(gse):
    p = os.path.join(DATA, "meta", gse + "_sample_map.csv")
    return pd.read_csv(p, dtype=str)

# ================= 1. GSE292325: labels + matrix =================
log("loading GSE292325 sample map ...")
sm325 = read_sample_map("GSE292325")
def parse_char(col, prefix):
    return col.str.replace(prefix, "", regex=False).str.strip()
sm325["grade"] = parse_char(sm325["char2"], "who grade:")
sm325["choudhury"] = parse_char(sm325["char3"], "choudhury class:")
sm325["kmeans"] = parse_char(sm325["char4"], "kmeans cluster:")
sm325["sex"] = parse_char(sm325["char5"], "Sex:")
log("GSE292325 labels: " + sm325["kmeans"].value_counts().to_dict().__str__())

log("loading GSE292325 processed matrix (1.25GB gz) ...")
m325 = load_processed(os.path.join(DATA, "GSE292325_matrix_processed.txt.gz"))
m325 = m325[~m325.index.duplicated(keep="first")]
log(f"matrix {m325.shape}")
assert m325.shape[1] == len(sm325), "column/sample count mismatch"

# ================= 2. feature selection + QC =================
mad = (m325 - m325.median(axis=1).values[:, None]).abs().median(axis=1)
na_frac = m325.isna().mean(axis=1)
m325 = m325[na_frac < 0.2]
mad = mad[mad.index.isin(m325.index)]
sel = mad.sort_values(ascending=False).index[:TOP_K]
log(f"selected {len(sel)} top-MAD CpGs (from {m325.shape[0]} CpGs)")

X = m325.loc[sel, :].T.copy()
Xz = (X - X.mean()) / X.std()
y = sm325["kmeans"].values

# ================= 3. k-means reproduction =================
km = KMeans(n_clusters=3, n_init=20, random_state=20260317)
pred_km = km.fit_predict(Xz.values)
ari_pub = adjusted_rand_score(y, pred_km)
sil = silhouette_score(Xz.values, pred_km, sample_size=5000, random_state=0)
log(f"k-means ARI vs published labels = {ari_pub:.3f}; silhouette = {sil:.3f}")
ct = pd.crosstab(pd.Series(y, name="published"), pd.Series(pred_km, name="our_km"))
log("crosstab:\n" + ct.to_string())

# ================= 4. nearest-centroid classifier + CV =================
centroids = {}
for lab in np.unique(y):
    centroids[lab] = Xz[y == lab].mean(axis=0)

def classify_centroid(Xz_new):
    labs = list(centroids.keys())
    common = [c for c in Xz_new.columns if c in centroids[labs[0]].index]
    M = np.vstack([centroids[l].reindex(common).values for l in labs])
    Xv = Xz_new[common].values
    d = ((Xv[:, None, :] - M[None, :, :]) ** 2).sum(axis=2)
    return np.array([labs[i] for i in d.argmin(axis=1)])

from sklearn.model_selection import cross_val_predict
def cv_centroid_acc(Xz_df, labels, k=5, seed=1):
    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=seed)
    preds = np.empty(len(labels), dtype=object)
    for tr, te in skf.split(Xz_df, labels):
        cents = {}
        for lab in np.unique(labels[tr]):
            cents[lab] = Xz_df.iloc[tr][labels[tr] == lab].mean(axis=0)
        M = np.vstack([cents[l] for l in sorted(cents)])
        d = ((Xz_df.iloc[te].values[:, None, :] - M[None, :, :]) ** 2).sum(axis=2)
        for i, ti in enumerate(te):
            preds[ti] = sorted(cents)[d[i].argmin()]
    return preds

preds_cv = cv_centroid_acc(Xz, y)
acc_cv = (preds_cv == y).mean()
log(f"nearest-centroid 5-fold CV accuracy vs published labels = {acc_cv:.3f}")
ari_ch = adjusted_rand_score(y, sm325["choudhury"].values)
log(f"published kmeans vs choudhury class ARI = {ari_ch:.3f}")

# save classifier artifacts
pd.Series(sel).to_csv(os.path.join(OUT, "methyl_selected_cpgs.csv"), index=False, header=["cpg"])
pd.DataFrame(centroids).to_csv(os.path.join(OUT, "methyl_centroids.csv"))
scaler = pd.DataFrame({"mean": X.mean(), "sd": X.std().replace(0, 1)})
scaler.to_csv(os.path.join(OUT, "methyl_scaler.csv"))

def zscore_with_artifacts(bet):
    idx = scaler.index.intersection(bet.index)
    bet = bet.loc[idx]
    mu = scaler["mean"].reindex(idx); sd = scaler["sd"].reindex(idx)
    z = (bet.T - mu).div(sd)
    return z, idx

# ================= 5. apply to GSE292326 =================
log("applying to GSE292326 ...")
sm326 = read_sample_map("GSE292326")
# labels are split across char4/char5
def extract_label(row):
    for c in ["char4", "char5", "char3"]:
        v = str(row[c])
        if v.startswith("kmeans cluster:"):
            return v.replace("kmeans cluster:", "").strip()
    return np.nan
sm326["kmeans"] = sm326.apply(extract_label, axis=1)
m326 = load_processed(os.path.join(DATA, "GSE292326_matrix_processed.txt.gz"))
log(f"292326 matrix {m326.shape}")
m326.index = [str(x).split("_")[0] for x in m326.index]
m326 = m326[~m326.index.duplicated(keep="first")]  # cg IDs may carry suffix
bet26, common26 = zscore_with_artifacts(m326)
log(f"292326: {len(common26)}/{len(sel)} classifier CpGs available")
lab26 = sm326["kmeans"].values
pred26 = classify_centroid(bet26)
agree26 = (pred26 == lab26).mean()
log(f"GSE292326 agreement with published labels = {agree26:.3f} (n={len(pred26)})")
pd.DataFrame({"published": lab26, "predicted": pred26}).to_csv(os.path.join(OUT, "gse292326_classification.csv"), index=False)

# ================= 6. apply to GSE304096 + computational pairing =================
log("applying to GSE304096 ...")
sm304 = read_sample_map("GSE304096")
m304 = load_processed(os.path.join(DATA, "GSE304096_longitudinal_processed.txt.gz"))
log(f"304096 matrix {m304.shape}; columns are sentrix barcodes")
m304.index = [str(x).split("_")[0] for x in m304.index]
m304 = m304[~m304.index.duplicated(keep="first")]
bet304, common304 = zscore_with_artifacts(m304)
log(f"304096: {len(common304)}/{len(sel)} classifier CpGs available")
pred304 = classify_centroid(bet304)
# class confidence = margin between best and second centroid distance
labs = sorted(centroids.keys())
common_c = [c for c in bet304.columns if c in centroids[labs[0]].index]
M = np.vstack([centroids[l].reindex(common_c).values for l in labs])
Xv = bet304[common_c].values
d = ((Xv[:, None, :] - M[None, :, :]) ** 2).sum(axis=2)
srt = np.sort(d, axis=1)
conf304 = srt[:, 1] - srt[:, 0]
cls304 = pd.DataFrame({"sample": bet304.index, "predicted_class": pred304, "confidence": conf304})

# computational pairing by correlation on classifier CpGs
C = np.corrcoef(bet304.values)
np.fill_diagonal(C, 0)
n = C.shape[0]
pairs = []
used = set()
flat = [(C[i, j], i, j) for i in range(n) for j in range(i + 1, n)]
flat.sort(reverse=True)
parent = list(range(n))
def find(a):
    while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
    return a
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb: parent[rb] = ra
for corr, i, j in flat:
    ri, rj = find(i), find(j)
    if ri != rj and len(np.where(np.array(parent) == ri)[0] if False else []) >= 0:
        # limit cluster size to 3
        ci = sum(1 for k2 in range(n) if find(k2) == ri)
        cj = sum(1 for k2 in range(n) if find(k2) == rj)
        if ci + cj <= 3:
            union(i, j)
groups = {}
for k2 in range(n):
    groups.setdefault(find(k2), []).append(k2)
pair_rows = []
gid = 0
for r, members in groups.items():
    gid += 1
    names = [bet304.index[m] for m in members]
    for nm in names:
        pair_rows.append({"patient_cluster": f"P{gid:02d}", "sample": nm,
                          "cluster_size": len(members), "predicted_class": pred304[bet304.index.get_loc(nm)]})
pairs_df = pd.DataFrame(pair_rows)
size_counts = pairs_df.groupby("patient_cluster").size().value_counts().to_dict()
log(f"computational pairing: {len(groups)} clusters; size distribution={size_counts}")
cls304 = cls304.merge(pairs_df[["sample", "patient_cluster", "cluster_size"]], on="sample")
# stability within clusters of size>=2
stab_rows = []
for pc, sub in pairs_df.groupby("patient_cluster"):
    if len(sub) < 2: continue
    labs_in = sub["predicted_class"].tolist()
    stab_rows.append({"patient_cluster": pc, "n_samples": len(sub),
                      "classes": "|".join(labs_in), "all_same": int(len(set(labs_in)) == 1)})
stab = pd.DataFrame(stab_rows)
log(f"within-patient class stability: {int(stab['all_same'].sum())}/{len(stab)} clusters fully concordant")
cls304.to_csv(os.path.join(OUT, "gse304096_classification.csv"), index=False)
stab.to_csv(os.path.join(OUT, "gse304096_stability.csv"), index=False)

# paired CpG shifts (primary vs first recurrence within clusters of size>=2)
deltas = []
cluster_ids = pairs_df["patient_cluster"].unique()
for pc in cluster_ids:
    mem = pairs_df[pairs_df["patient_cluster"] == pc]
    if len(mem) < 2: continue
    cols = [np.where(bet304.index == nm)[0][0] for nm in mem["sample"]]
    deltas.append(bet304.values[cols[0], :] - bet304.values[cols[1], :])
if len(deltas) >= 5:
    Dm = np.vstack(deltas)  # rows = patients, cols = CpGs (common304 order)
    from scipy import stats
    tstat, pval = stats.ttest_1samp(Dm, 0, axis=0)
    padj = np.minimum(1, pval * len(pval))
    sig = pd.DataFrame({"cpg": list(map(str, common304)), "mean_delta": Dm.mean(axis=0), "t": tstat, "p": pval, "p_adj": padj})
    sig.sort_values("p_adj").to_csv(os.path.join(OUT, "gse304096_paired_cpg_shifts.csv"), index=False)
    log(f"paired CpG shifts: {(sig['p_adj'] < 0.05).sum()} CpGs BH-significant (n pairs={len(Dm)})")

log("ANALYSIS_B1_DONE")




