import os, gzip, warnings, time
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")
from sklearn.metrics import adjusted_rand_score, balanced_accuracy_score, f1_score

t0 = time.time()
DATA = r"E:\脑膜瘤npj\data"
OUT = r"E:\脑膜瘤npj\results"
def log(m): print(f"[{time.time()-t0:6.0f}s] {m}", flush=True)

# ---------- re-parse GSE292326 series matrix per-sample ----------
p = os.path.join(DATA, "GSE292326_series_matrix.txt.gz")
acc, tit, chars = [], [], []
with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith("!Sample_geo_accession"):
            acc = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_title"):
            tit = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1"):
            chars.append([x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]])
log(f"samples: {len(acc)}; title example: {tit[:3]}")

# per-sample label extraction from all char fields
n = len(acc)
km = [None] * n
ch = [None] * n
gr = [None] * n
for c in chars:
    for i, v in enumerate(c):
        v = str(v).strip()
        if v.startswith("kmeans cluster:"):
            km[i] = v.split(":")[1].strip()
        elif v.startswith("choudhury class:"):
            ch[i] = v.split(":")[1].strip()
        elif v.startswith("who grade:"):
            gr[i] = v.split(":")[1].strip()
log(f"kmeans labels: {sum(x is not None for x in km)}/{n}; choudhury: {sum(x is not None for x in ch)}/{n}")

# ---------- load our classifier predictions for 292326 ----------
pred = pd.read_csv(os.path.join(OUT, "gse292326_classification.csv"))
# saved columns: published (kmeans), predicted — built from the same sample order assumption
pred_v = pred["predicted"].astype(str).values
pub_v = pred["published"].astype(str).values

# rebuild from current parse (to make order explicit)
# our predictions in sample order = rows of the saved file
log(f"saved predictions n={len(pred_v)}")

km_arr = np.array([x if x is not None else "NA" for x in km])
ch_arr = np.array([x if x is not None else "NA" for x in ch])

# mapping correspondence learned from TRAINING (292325): pred class -> choudhury class
sm325 = pd.read_csv(os.path.join(DATA, "meta", "GSE292325_sample_map.csv"), dtype=str)
km325 = sm325["char4"].str.replace("kmeans cluster:", "", regex=False).str.strip()
ch325 = sm325["char3"].str.replace("choudhury class:", "", regex=False).str.strip()
ct325 = pd.crosstab(km325, ch325)
log("292325 kmeans x choudhury (training correspondence):\n" + ct325.to_string())
# derive mapping: for each pred class (same label space as 292325 kmeans), the majority choudhury class
mapping = {}
for k in ct325.index:
    mapping[k] = ct325.loc[k].idxmax()
log("derived mapping pred->choudhury: " + str(mapping))

# TEST 1 (order as-is): agreement of pred vs choudhury through mapping
pred_ch = np.array([mapping.get(p, "?") for p in pred_v])
mask = ch_arr != "NA"
agree = (pred_ch[mask] == ch_arr[mask]).mean()
ari1 = adjusted_rand_score(ch_arr[mask], pred_ch[mask])
ba = balanced_accuracy_score(ch_arr[mask], pred_ch[mask])
f1m = f1_score(ch_arr[mask], pred_ch[mask], average="macro")
log(f"ORDER AS-IS: agreement vs choudhury = {agree:.3f}; ARI = {ari1:.3f}; balanced acc = {ba:.3f}; macro-F1 = {f1m:.3f} (n={mask.sum()})")

# TEST 2: order sensitivity — try a set of plausible alternative permutations
def score_order(perm):
    ch_p = ch_arr[perm]
    m = ch_p != "NA"
    if m.sum() < 20:
        return -1.0
    return adjusted_rand_score(ch_p[m], pred_ch[m])

perms = {
    "as-is": np.arange(n),
    "reversed": np.arange(n)[::-1],
    "sorted-title": np.argsort(tit),
    "sorted-rev": np.argsort(tit)[::-1],
}
rng = np.random.default_rng(0)
rand_scores = []
for _ in range(200):
    rand_scores.append(score_order(rng.permutation(n)))
log("order tests:")
for name, perm in perms.items():
    log(f"  {name}: ARI(pred, choudhury) = {score_order(perm):.3f}")
log(f"  random permutations: mean ARI = {np.mean([s for s in rand_scores if s >= 0]):.3f} (n=200)")

# Test 3: same check against kmeans labels (published) — as-is
m2 = km_arr != "NA"
ari_km = adjusted_rand_score(km_arr[m2], pred_v[m2])
agree_km = (km_arr[m2] == pred_v[m2]).mean()
log(f"kmeans agreement as-is: {agree_km:.3f}; ARI = {ari_km:.3f} (n={m2.sum()})")

# Save the per-sample table for the audit
audit = pd.DataFrame({"GSM": acc, "title": tit, "kmeans": km_arr, "choudhury": ch_arr,
                      "who_grade": gr, "our_pred": pred_v, "pred_as_choudhury": pred_ch})
audit.to_csv(os.path.join(OUT, "gse292326_identity_audit.csv"), index=False)
log("AUDIT_292326_DONE")
