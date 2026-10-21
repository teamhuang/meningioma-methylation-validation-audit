"""Cross-platform transfer harmonization experiment (pre-registered, 3 paths).

Paths (identical on every platform; chosen BEFORE seeing any test/outcome data):
  P1_frozen_raw      : raw frozen centroids, full 20k selected CpG panel
  P2_refnorm_shared  : frozen reference mean/sd scaler + frozen centroids on the
                       shared-CpG subset that is actually present in the target
                       (the diagnosed failure mode: panel mismatch under SeSAMe)
  P3_robust_panel    : a pre-registered ROBUST CpG subset selected from the
                       TRAINING set only (top MAD rank + low-missing rank, k=2000),
                       transferred with frozen centroids resampled to that panel
Plus a reject mechanism: samples whose min distance to any centroid exceeds the
training-set 99th-percentile distance are flagged "low confidence / reject".

Per-platform outputs: confusion matrix, ARI, macro-F1, accuracy, centroid distance
distribution, and reject rate. No platform uses outcome/survival data to pick a
path or threshold — the threshold is derived from the TRAINING set only.
"""
import os, gzip, re, time, warnings
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, f1_score
from scipy.stats import rankdata, percentileofscore

t0 = time.time()
DATA = r"E:\脑膜瘤npj\data"
OUT = r"E:\脑膜瘤npj\results"
NG = os.path.join(DATA, "suppl", "natgen2022")
SEED = 20260317
def log(m): print(f"[{time.time()-t0:6.0f}s] {m}", flush=True)

scaler = pd.read_csv(os.path.join(OUT, "methyl_scaler.csv"), index_col=0)
centroids = pd.read_csv(os.path.join(OUT, "methyl_centroids.csv"), index_col=0)
sel = pd.read_csv(os.path.join(OUT, "methyl_selected_cpgs.csv"))["cpg"].tolist()
labs = list(centroids.columns)
log(f"artifacts: {len(sel)} selected CpGs; labels {labs}; scaler {len(scaler)}")

# ---------- training-set-only robust panel (pre-registered; no test/outcome data) ----------
# cache the reduced training matrix so re-runs are fast
cache325 = os.path.join(OUT, "cache_m325_sel.csv.gz")
if os.path.exists(cache325):
    m325_sel = pd.read_csv(cache325, index_col=0).astype(np.float32)
    log(f"loaded cached m325_sel {m325_sel.shape}")
else:
    m325 = pd.read_csv(os.path.join(DATA, "GSE292325_matrix_processed.txt.gz"),
                       sep="\t", index_col=0).astype(np.float32)
    m325 = m325[~m325.index.duplicated(keep="first")]
    sel325 = [c for c in scaler.index if c in m325.index]
    log(f"selected panel available in 292325: {len(sel325)}/{len(scaler)}")
    m325_sel = m325.loc[sel325, :]
    m325_sel.to_csv(cache325)
    del m325
    log(f"cached m325_sel {m325_sel.shape}")
sel325 = list(m325_sel.index)
Xtr = m325_sel.T
na_frac_tr = Xtr.isna().mean()
mad_tr = (Xtr - Xtr.median()).abs().median()
sp = pd.DataFrame({"mad": mad_tr, "na_frac": na_frac_tr})
sp["rank_mad"] = sp["mad"].rank(ascending=False, method="first")
sp["rank_na"] = sp["na_frac"].rank(ascending=True, method="first")
sp["score"] = sp["rank_mad"] + sp["rank_na"]
ROBUST_K = 2000
robust = list(sp.sort_values("score", ascending=False).index[:ROBUST_K])
log(f"robust training-set panel: {len(robust)} CpGs")
pd.Series(robust).to_csv(os.path.join(OUT, "methyl_robust_cpgs.csv"), header=["cpg"], index=False)

# training-set reference distance threshold (reject cutoff) — derived from TRAINING set
def z_with_scaler(beta_df):
    """beta_df: rows=CpGs, cols=samples. Z: rows=samples, cols=the scaler CpGs
    present in beta_df (frozen reference normalization)."""
    idx = [c for c in scaler.index if c in beta_df.index]
    mu = scaler["mean"].reindex(idx); sd = scaler["sd"].reindex(idx)
    Z = (beta_df.loc[idx].T - mu).div(sd)
    return Z

def predict_centroids(Z, panel):
    cc = [c for c in panel if c in centroids.index and c in Z.columns]
    if len(cc) < 50:
        return (np.full(len(Z), "NA", dtype=object), np.full(len(Z), np.nan),
                np.full(len(Z), np.nan))
    M = centroids.loc[cc, labs].values.T
    Xv = Z[cc].values
    d = ((np.nan_to_num(Xv[:, None, :], nan=0.0) - M[None, :, :]) ** 2).sum(axis=2)
    # NaN CpGs contribute 0 distance; require >=50% finite CpGs per sample
    n_finite = np.isfinite(Xv).sum(axis=1)
    pred = np.array([labs[i] if f >= 0.5 * len(cc) else "NA"
                     for i, f in zip(d.argmin(axis=1), n_finite)], dtype=object)
    dist_min = np.where(n_finite >= 0.5 * len(cc), d.min(axis=1), np.nan)
    dist2nd = np.sort(d, axis=1)[:, 1]
    conf = dist2nd - dist_min
    return pred, dist_min, conf

# training-set P99 min-distance -> global reject threshold
Ztr = z_with_scaler(m325_sel)
ptr, dtr, ctr = predict_centroids(Ztr, sel)
pos = dtr[np.isfinite(dtr)]
if len(pos) == 0:
    raise RuntimeError("training-set distances all NaN — scaler/CpG mismatch")
THR = float(np.percentile(pos, 99))
log(f"training-set min-distance P99 (reject threshold) = {THR:.2f}; "
    f"P50={np.percentile(pos,50):.2f} P95={np.percentile(pos,95):.2f} "
    f"finite={len(pos)}/{len(dtr)}")

# ---- panel-matched training-set P99 reject thresholds (training data only) ----
import json as _json
# P1: raw-beta-space training distances to descaled centroids (deployment metric)
cc_tr = [c for c in sel if c in m325_sel.index]
Cz_tr = centroids.loc[cc_tr, labs].values.T                       # (3, K) z-space
mu_tr = scaler["mean"].reindex(cc_tr).values
sd_tr = scaler["sd"].reindex(cc_tr).values
Cb_tr = Cz_tr * sd_tr[None, :] + mu_tr[None, :]                   # (3, K) raw beta
Xtr_raw = m325_sel.T[cc_tr].values                                # (Ntr, K) raw beta
Xtr_raw = np.where(np.isnan(Xtr_raw), mu_tr[None, :], Xtr_raw)
d1_tr = ((Xtr_raw[:, None, :] - Cb_tr[None, :, :]) ** 2).sum(axis=2)
THR_P1 = float(np.percentile(d1_tr.min(axis=1), 99))
# P3: z-space, robust 2k panel (nan->0 = z-space training mean, matches predict_centroids)
Ztr_full = z_with_scaler(m325_sel)
cc_r = [c for c in robust if c in Ztr_full.columns]
X3 = np.nan_to_num(Ztr_full[cc_r].values, nan=0.0)
d3_tr = ((X3[:, None, :] - centroids.loc[cc_r, labs].values.T[None, :, :]) ** 2).sum(axis=2)
THR_P3 = float(np.percentile(d3_tr.min(axis=1), 99))
THR_P2 = THR  # full 20k z-space P99 (computed above via predict_centroids on Ztr)
log(f"panel-matched training P99 thresholds: P1_raw={THR_P1:.2f} "
    f"P2_z20k={THR_P2:.2f} P3_z2k={THR_P3:.2f} "
    f"(panels: P1={len(cc_tr)}, P2={len(sel)}, P3={len(cc_r)})")
with open(os.path.join(OUT, "methyl_train_thr.json"), "w") as fh:
    _json.dump({"THR_P1_raw": THR_P1, "THR_P2_z20k": THR_P2, "THR_P3_z2k": THR_P3,
                "n_train": int(d1_tr.shape[0]),
                "panels": {"P1": len(cc_tr), "P2": len(sel), "P3": len(cc_r)}}, fh, indent=1)

def report(tag, pred, dist, truth, thr=None, log_ct=True, relabel=None):
    if thr is None:
        thr = THR_P2
    if relabel is not None:
        pred = np.array([relabel.get(x, x) for x in pred], dtype=object)
    truth = np.array([None if (t is None or (isinstance(t, float) and np.isnan(t))) else str(t)
                      for t in truth], dtype=object)
    pred = np.array([str(x) for x in pred], dtype=object)
    lab = np.array([t is not None for t in truth])
    sub = lab & (pred != "NA") & ~np.isnan(dist)
    rej = float(np.nanmean(dist > thr)) if sub.any() else float("nan")
    out = {"path": tag, "n_lab": int(sub.sum()), "reject_rate": round(rej, 3)}
    if sub.sum() >= 2 and len(set(truth[sub])) > 1:
        out["ARI"] = round(adjusted_rand_score(truth[sub], pred[sub]), 3)
        out["macro_F1"] = round(f1_score(truth[sub], pred[sub], average="macro",
                                        zero_division=0), 3)
        out["acc"] = round(float((truth[sub] == pred[sub]).mean()), 3)
        out["dist_P50"] = round(float(np.nanmedian(dist[sub])), 1)
        out["dist_P95"] = round(float(np.nanpercentile(dist[sub], 95)), 1)
    if log_ct:
        out["_ct"] = pd.crosstab(pd.Series(pred[sub]),
                                 pd.Series(truth[sub])).to_dict() if sub.sum() else {}
    return out

PATHS = {"P1_frozen_raw": sel}
P3_ROBUST = robust
def note(): return "P2 uses shared-CpG subset of the frozen 20k panel (panel-mismatch control)"

# =====================================================================
# PLATFORM 1: UCSF/SeSAMe (GSE183647, 200 M-IDs) — reference = unsupervised
# 2000-CpG k-means (same space as analysis_F; label space differs, so report
# crosstab + distance + reject, not raw accuracy vs a single class).
# =====================================================================
log("=== PLATFORM 1: UCSF/SeSAMe (GSE183647) ===")
fl = pd.read_csv(os.path.join(DATA, "suppl", "GSE183647_filelist.txt"), sep="\t",
                 comment="#", header=None, names=["kind", "name", "time", "size", "type"])
fl = fl[fl["kind"] == "File"]
gsm2bar = {}
for nm in fl["name"]:
    parts = str(nm).split("_")
    if parts[0].startswith("GSM") and len(parts) >= 3:
        gsm2bar[parts[0]] = parts[1] + "_" + parts[2]
gsm2mid = {}
with gzip.open(os.path.join(DATA, "GSE183647_series_matrix.txt.gz"), "rt",
               encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith("!Sample_geo_accession"):
            acc = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_title"):
            tit = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
for g, t in zip(acc, tit):
    m2 = re.match(r"UCSF Meningioma (\d+) \(methylation\)", t)
    if m2:
        gsm2mid[g] = "M" + m2.group(1)
bar2mid = {gsm2bar[g]: gsm2mid[g] for g in gsm2mid if g in gsm2bar}
bet = pd.read_csv(os.path.join(DATA, "GSE183647_meningioma_betas.csv.gz"), index_col=0)
bet.columns = [bar2mid.get(str(c).strip(), None) for c in bet.columns]
bet = bet.iloc[:, [i for i, c in enumerate(bet.columns) if c is not None]]
log(f"UCSF betas: {bet.shape}")

m21 = pd.read_excel(os.path.join(NG, "41588_2022_1061_MOESM21_ESM.xlsx"), index_col=0)
m21.columns = [str(c).strip() for c in m21.columns]
m21.index = [str(c).split("_")[0].strip() for c in m21.index]
# reference clustering in the published 2000-CpG fingerprint space (samples = columns)
km21 = KMeans(n_clusters=3, n_init=20, random_state=SEED).fit_predict(
    ((m21.T - m21.T.mean()) / m21.T.std()).fillna(0).values)
ref21 = pd.Series(km21, index=m21.columns)
# UCSF reference for the 200 M-ID samples: k-means on the fingerprint CpGs available
panel200 = [c for c in m21.index if c in bet.index]
X200 = bet.loc[panel200, :].T
km_full = KMeans(n_clusters=3, n_init=20, random_state=SEED).fit_predict(
    ((X200 - X200.mean()) / X200.std()).fillna(0).values)
truth_u = pd.Series(km_full, index=list(bet.columns))
log(f"UCSF reference kmeans: fingerprint CpGs shared={len(panel200)}/{len(m21)}; "
    f"200-sample clusters {pd.Series(truth_u).value_counts().to_dict()}")

def predict_panel(Z, panel):
    """P1: frozen raw-space transfer — descale frozen centroids back to beta space,
    predict on the frozen panel's CpGs present in raw beta_df (no re-scaling)."""
    cc = [c for c in panel if c in centroids.index and c in Z.columns]
    if len(cc) < 50:
        return (np.full(len(Z), "NA", dtype=object), np.full(len(Z), np.nan),
                np.full(len(Z), np.nan))
    # raw-space centroids from scaler mean/sd of the 292325 training set
    Cz = centroids.loc[cc, labs].values.T                      # (3, K) z-space
    mu = scaler["mean"].reindex(cc).values                     # (K,)
    sd = scaler["sd"].reindex(cc).values
    Cb = Cz * sd[None, :] + mu[None, :]                        # (3, K) beta space
    Xv = Z[cc].values                                          # (N, K) raw beta
    n_finite = np.isfinite(Xv).sum(axis=1)
    # NaN imputation: per-CpG training mean (neutral anchor in raw beta space)
    Xv = np.where(np.isnan(Xv), mu[None, :], Xv)
    d = ((Xv[:, None, :] - Cb[None, :, :]) ** 2).sum(axis=2)
    pred = np.array([labs[i] if f >= 0.5 * len(cc) else "NA"
                     for i, f in zip(d.argmin(axis=1), n_finite)], dtype=object)
    dist_min = np.where(n_finite >= 0.5 * len(cc), d.min(axis=1), np.nan)
    dist2nd = np.sort(d, axis=1)[:, 1]
    conf = dist2nd - dist_min
    return pred, dist_min, conf

def predict_shared(Z, sel_full):
    """P2: frozen-reference (scaler) normalization restricted to the shared-CpG
    subset actually present in Z, then z-space centroid prediction."""
    shared = [c for c in sel_full if c in set(Z.columns)]
    if len(shared) < 100:
        return np.full(len(Z), "NA", dtype=object), np.full(len(Z), np.nan), np.full(len(Z), np.nan)
    return predict_centroids(Z, shared)

Zb = z_with_scaler(bet)
avail_u = [c for c in Zb.columns if c in centroids.index]
log(f"UCSF available shared classifier CpGs: {len(avail_u)}/{len(sel)}")
# samples are Zb rows (== bet.columns); truth must align to samples
truth_u_samples = [truth_u.get(x, None) for x in bet.columns]
res_u = []
# P1: frozen raw-space transfer (RAW betas, not z-space)
p, d, c = predict_panel(bet.T, sel)
r = report("P1_frozen_raw", p, d, truth_u_samples, thr=THR_P1)
r["pred_counts"] = pd.Series(p).value_counts().to_dict()
res_u.append(r)
log(f"P1_frozen_raw: pred={pd.Series(p).value_counts().to_dict()} "
    f"ARI={r.get('ARI')} macroF1={r.get('macro_F1')} reject={r['reject_rate']}")
# P2: shared-CpG reference normalization on UCSF
p2, d2, c2 = predict_shared(Zb, sel)
r2 = report("P2_refnorm_shared", p2, d2, truth_u_samples, thr=THR_P2)
r2["pred_counts"] = pd.Series(p2).value_counts().to_dict()
res_u.append(r2)
log(f"P2_refnorm_shared: pred={pd.Series(p2).value_counts().to_dict()} "
    f"ARI={r2.get('ARI')} macroF1={r2.get('macro_F1')} reject={r2['reject_rate']}")
# P3: robust panel on UCSF
p3, d3, c3 = predict_centroids(Zb, [x for x in robust if x in set(Zb.columns)])
r3 = report("P3_robust_panel", p3, d3, truth_u_samples, thr=THR_P3)
r3["pred_counts"] = pd.Series(p3).value_counts().to_dict()
res_u.append(r3)
log(f"P3_robust_panel: pred={pd.Series(p3).value_counts().to_dict()} "
    f"ARI={r3.get('ARI')} macroF1={r3.get('macro_F1')} reject={r3['reject_rate']}")

# =====================================================================
# PLATFORM 2: GSE189521 (EPIC, 110) — ground truth = clinical MenG A/B/C
# =====================================================================
log("=== PLATFORM 2: GSE189521/EPIC ===")
p189 = os.path.join(DATA, "GSE189521_series_matrix.txt.gz")
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
B189 = pd.read_csv(p189, sep="\t", index_col=0, skiprows=start, compression="gzip")
B189 = B189[~B189.index.astype(str).str.startswith("!")]
B189.columns = [str(c).strip().strip('"') for c in B189.columns]
B189.index = [str(x).split("_")[0] for x in B189.index]
B189 = B189[~B189.index.duplicated(keep="first")]
log(f"189521 matrix {B189.shape}")
clin = pd.read_csv(os.path.join(OUT, "clinical_parsed.csv"), dtype={"Sample": str})
clin["Sample"] = clin["Sample"].astype(str).str.strip()
tab = pd.DataFrame({"GSM": geo_acc, "title": titles})
tab["Sample"] = [re.sub(r"\s*\(methylation\)\s*", "", t).strip() for t in titles]
tab["Sample"] = tab["Sample"].astype(str).str.strip()
tab = tab.merge(clin[["Sample", "menG"]], on="Sample", how="left")
log(f"189521 title↔menG matched: {tab['menG'].notna().sum()}/{len(tab)}")
# classifier emits 1/2/3 (training-space labels); menG is A/B/C. Before metrics,
# apply a Hungarian-optimal label mapping fitted ONCE on P2 predictions (the
# pre-registered primary path) — explicitly a data-driven relabeling for
# cross-label-space reporting, not a model change.
tab = tab[tab["menG"].notna()].set_index("Sample")
# map GSM -> Sample via paired geo_acc/titles arrays, then label the matrix columns
_tmp = pd.DataFrame({"GSM": geo_acc, "Sample": [re.sub(r"\s*\(methylation\)\s*", "", t).strip() for t in titles]})
gsm2sample = dict(zip(_tmp["GSM"], _tmp["Sample"]))
truth_189 = [tab.loc[gsm2sample[c], "menG"] if gsm2sample.get(c) in tab.index else None
             for c in B189.columns]
n_lab = sum(1 for c in B189.columns if gsm2sample.get(c) in tab.index)
log(f"189521 labelled cols: {n_lab}")
Z189 = z_with_scaler(B189)
avail_189 = [c for c in Z189.columns if c in centroids.index]
log(f"189521 available shared classifier CpGs: {len(avail_189)}/{len(sel)}")
res_189 = []
truth_189_s = truth_189  # GSM-aligned labels (None where no clinical menG)
# fit the 1/2/3 -> A/B/C label map once, on P2 predictions only (pre-registered)
_p2, _d2, _c2 = predict_shared(Z189, sel)
_m = pd.crosstab(pd.Series([str(x) for x in _p2]), pd.Series([str(t) for t in truth_189_s if t is not None]))
from scipy.optimize import linear_sum_assignment
ri, ci = linear_sum_assignment(-_m.values)
relabel_189 = {str(_m.index[i]): str(_m.columns[j]) for i, j in zip(ri, ci)}
log(f"189521 label map (fit on P2, Hungarian): {relabel_189}; crosstab={_m.to_dict()}")
path_fns = [
    ("P1_frozen_raw", lambda: predict_panel(B189.T, sel), THR_P1),
    ("P2_refnorm_shared", lambda: predict_shared(Z189, sel), THR_P2),
    ("P3_robust_panel", lambda: predict_centroids(Z189, [x for x in robust if x in set(Z189.columns)]), THR_P3),
]
for tag, fn, thr_p in path_fns:
    p, d, c = fn()
    r = report(tag, p, d, truth_189_s, thr=thr_p, relabel=relabel_189)
    r["pred_counts"] = pd.Series(p).value_counts().to_dict()
    res_189.append(r)
    log(f"{tag}: pred={pd.Series(p).value_counts().to_dict()} ARI={r.get('ARI')} "
        f"macroF1={r.get('macro_F1')} reject={r['reject_rate']}")

# =====================================================================
# PLATFORM 3: GSE270371 (Nat Med 994) — no outcome; distribution + reject only
# =====================================================================
log("=== PLATFORM 3: GSE270371 (distribution only) ===")
# GSE270371 file: rows = CpGs, columns alternate intensity (I1..I994) and
# "Detection Pval" (de-duplicated by pandas with .N suffixes). Keep only intensity cols.
B270 = pd.read_csv(os.path.join(DATA, "GSE270371_Matrix_processed.txt.gz"),
                   sep="\t", index_col=0).astype(np.float32)
B270 = B270[[c for c in B270.columns if not str(c).startswith("Detection Pval")]]
B270 = B270[~B270.index.duplicated(keep="first")]
log(f"270371 matrix {B270.shape}")
Z270 = z_with_scaler(B270)
avail_270 = [c for c in Z270.columns if c in centroids.index]
log(f"270371 available shared classifier CpGs: {len(avail_270)}/{len(sel)}")
res_270 = []
path_fns_270 = [
    ("P1_frozen_raw", lambda: predict_panel(B270.T, sel), THR_P1),
    ("P2_refnorm_shared", lambda: predict_shared(Z270, sel), THR_P2),
    ("P3_robust_panel", lambda: predict_centroids(Z270, [x for x in robust if x in set(Z270.columns)]), THR_P3),
]
for tag, fn, thr_p in path_fns_270:
    p, d, c = fn()
    rej = float(np.nanmean(d > thr_p))
    res_270.append({"path": tag, "pred_dist": pd.Series(p).value_counts().to_dict(),
                    "n_reject": int(np.nansum(d > thr_p)), "reject_rate": round(rej, 3),
                    "dist_P50": round(float(np.nanmedian(d)), 1),
                    "dist_P95": round(float(np.nanpercentile(d, 95)), 1)})
    log(f"{tag}: dist={pd.Series(p).value_counts().to_dict()} reject={rej:.3f}")

pd.DataFrame(res_u).drop(columns=[c for c in res_u[0] if c.startswith("_")],
                         errors="ignore").to_csv(os.path.join(OUT, "harmonization_ucsfsesame.csv"), index=False)
pd.DataFrame(res_189).drop(columns=[c for c in res_189[0] if c.startswith("_")],
                           errors="ignore").to_csv(os.path.join(OUT, "harmonization_gse189521.csv"), index=False)
pd.DataFrame(res_270).to_csv(os.path.join(OUT, "harmonization_gse270371.csv"), index=False)
log("HARMONIZATION_DONE")
