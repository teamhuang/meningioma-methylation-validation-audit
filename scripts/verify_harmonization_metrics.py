# -*- coding: utf-8 -*-
"""Recompute harmonization agreement with CONSISTENT, label-space-safe metrics.

Motivated by manuscript review (2026-10-18):
  * harmonization_ucsfsesame.csv acc/macro_F1 were raw string equality between
    frozen labels {1,2,3} and sklearn kmeans labels {0,1,2} -> artifacts (e.g.
    raw-transfer string equality = 0.04 vs Hungarian best-map 0.545).
  * the v11 Hungarian map for GSE189521 was fitted on a crosstab whose truth
    side was index-compressed (pred len 110 vs filtered truth len 100) -> the
    printed crosstab (diagonal 43/100) is misaligned; the fitted map happened
    to be correct but must be re-verified on correctly aligned data.

For every path this script reports, on CORRECTLY aligned pairs:
  ARI            label-invariant headline metric
  bestmap_acc    Hungarian-optimal 1:1 label mapping accuracy (comparable to
                 the 0.545 raw-transfer diagnosis number)
  mapped_acc / mapped_macroF1  under the per-path Hungarian map
  reject_rate    frozen panel-matched P99 threshold
Outputs:
  results/harmonization_agreement_fixed.csv
  results/harmonization_pred_ucsfsesame.csv   (per-sample)
  results/harmonization_pred_gse189521.csv    (per-sample)
"""
import os, gzip, re, time, json, warnings
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, f1_score
from scipy.optimize import linear_sum_assignment

t0 = time.time()
DATA = r"E:\脑膜瘤npj\data"
OUT = r"E:\脑膜瘤npj\results"
NG = os.path.join(DATA, "suppl", "natgen2022")
SEED = 20260317


def log(m): print(f"[{time.time()-t0:6.0f}s] {m}", flush=True)


scaler = pd.read_csv(os.path.join(OUT, "methyl_scaler.csv"), index_col=0)
centroids = pd.read_csv(os.path.join(OUT, "methyl_centroids.csv"), index_col=0)
sel = pd.read_csv(os.path.join(OUT, "methyl_selected_cpgs.csv"))["cpg"].tolist()
robust = pd.read_csv(os.path.join(OUT, "methyl_robust_cpgs.csv"))["cpg"].tolist()
labs = list(centroids.columns)
thr = json.load(open(os.path.join(OUT, "methyl_train_thr.json")))
THR_P1, THR_P2, THR_P3 = thr["THR_P1_raw"], thr["THR_P2_z20k"], thr["THR_P3_z2k"]
log(f"artifacts: {len(sel)} sel / {len(robust)} robust; thr P1={THR_P1:.1f} P2={THR_P2:.1f} P3={THR_P3:.1f}")


def z_with_scaler(beta_df):
    """beta_df rows=CpGs; returns samples x scaler-CpGs z matrix (frozen reference)."""
    idx = [c for c in scaler.index if c in beta_df.index]
    mu = scaler["mean"].reindex(idx); sd = scaler["sd"].reindex(idx)
    return (beta_df.loc[idx].T - mu).div(sd)


def _centroid_predict_on(Xv, cc):
    """Xv: (n_samples, len(cc)) z-space; frozen centroids; NaN->0 (z-space train mean)."""
    M = centroids.loc[cc, labs].values.T
    d = ((np.nan_to_num(Xv[:, None, :], nan=0.0) - M[None, :, :]) ** 2).sum(axis=2)
    n_finite = np.isfinite(Xv).sum(axis=1)
    pred = np.array([labs[i] if f >= 0.5 * len(cc) else "NA"
                     for i, f in zip(d.argmin(axis=1), n_finite)], dtype=object)
    dist_min = np.where(n_finite >= 0.5 * len(cc), d.min(axis=1), np.nan)
    return pred, dist_min


def predict_panel(Bs, panel):
    """P1: Bs = samples x CpGs RAW beta matrix; descaled frozen centroids;
    NaN -> per-CpG training mean (raw-space neutral anchor)."""
    cc = [c for c in panel if c in centroids.index and c in Bs.columns]
    Xv = Bs[cc].values.astype(np.float64)
    n_finite = np.isfinite(Xv).sum(axis=1)
    mu = scaler["mean"].reindex(cc).values
    Cb = centroids.loc[cc, labs].values.T * scaler["sd"].reindex(cc).values[None, :] + mu[None, :]
    Xv = np.where(np.isnan(Xv), mu[None, :], Xv)
    d = ((Xv[:, None, :] - Cb[None, :, :]) ** 2).sum(axis=2)
    pred = np.array([labs[i] if f >= 0.5 * len(cc) else "NA"
                     for i, f in zip(d.argmin(axis=1), n_finite)], dtype=object)
    dist_min = np.where(n_finite >= 0.5 * len(cc), d.min(axis=1), np.nan)
    return pred, dist_min


def predict_shared(Z, panel):
    shared = [c for c in panel if c in set(Z.columns)]
    if len(shared) < 100:
        return np.full(len(Z), "NA", dtype=object), np.full(len(Z), np.nan)
    return _centroid_predict_on(Z[shared].values, shared)


def best_map_acc(truth, pred):
    m = pd.crosstab(pd.Series([str(x) for x in pred]), pd.Series([str(t) for t in truth]))
    ri, ci = linear_sum_assignment(-m.values)
    return m.values[ri, ci].sum() / m.values.sum(), {str(m.index[i]): str(m.columns[j]) for i, j in zip(ri, ci)}


def metrics(tag, platform, truth, pred, dist, thr_p):
    t = np.array([None if x is None or (isinstance(x, float) and np.isnan(x)) else str(x)
                  for x in truth], dtype=object)
    p = np.array([str(x) for x in pred], dtype=object)
    sub = np.array([x is not None for x in t]) & (p != "NA") & np.isfinite(dist)
    out = {"platform": platform, "path": tag, "n": int(sub.sum())}
    fin = np.isfinite(dist) & (p != "NA")
    if sub.sum() >= 2:
        tt, pp = t[sub], p[sub]
        out["ARI"] = round(float(adjusted_rand_score(tt, pp)), 3)
        bm, mp = best_map_acc(tt, pp)
        out["bestmap_acc"] = round(float(bm), 3)
        out["map"] = str(mp)
        pm = np.array([mp.get(x, x) for x in pp], dtype=object)
        out["mapped_acc"] = round(float((tt == pm).mean()), 3)
        out["mapped_macroF1"] = round(float(f1_score(tt, pm, average="macro", zero_division=0)), 3)
    # reject with EXPLICIT denominators (convention note archived with the table)
    rn_all = int(np.nansum((dist > thr_p) & fin))
    den_all = int(fin.sum())
    rn_lab = int(np.nansum(((dist > thr_p) & fin)[sub]))
    out["reject_n_profiled"] = rn_all
    out["reject_denom_profiled"] = den_all
    out["reject_rate_profiled"] = round(rn_all / den_all, 4) if den_all else float("nan")
    out["reject_n_labelled"] = rn_lab
    out["reject_denom_labelled"] = int(sub.sum())
    out["reject_rate_labelled"] = round(rn_lab / int(sub.sum()), 4) if sub.any() else float("nan")
    out["reject_convention"] = ("profiled = all methylation profiles with finite centroid distance; "
                                "labelled = reference-labelled evaluable subset; same convention in both columns")
    log(f"{platform}/{tag}: n={out['n']} ARI={out.get('ARI')} bestmap={out.get('bestmap_acc')} "
        f"acc={out.get('mapped_acc')} F1={out.get('mapped_macroF1')} "
        f"rej(profiled)={rn_all}/{den_all} rej(labelled)={rn_lab}/{int(sub.sum())} map={out.get('map')}")
    return out


rows, preds_u = [], None

# ================= PLATFORM 1: UCSF/SeSAMe (reference = unsupervised kmeans) =================
log("=== UCSF/SeSAMe ===")
fl = pd.read_csv(os.path.join(DATA, "suppl", "GSE183647_filelist.txt"), sep="\t",
                 comment="#", header=None, names=["kind", "name", "time", "size", "type"])
fl = fl[fl["kind"] == "File"]
gsm2bar = {}
for nm in fl["name"]:
    parts = str(nm).split("_")
    if parts[0].startswith("GSM") and len(parts) >= 3:
        gsm2bar[parts[0]] = parts[1] + "_" + parts[2]
gsm2mid = {}
acc_l, tit_l = [], []
with gzip.open(os.path.join(DATA, "GSE183647_series_matrix.txt.gz"), "rt",
               encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith("!Sample_geo_accession"):
            acc_l = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_title"):
            tit_l = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
for g, t in zip(acc_l, tit_l):
    m2 = re.match(r"UCSF Meningioma (\d+) \(methylation\)", t)
    if m2:
        gsm2mid[g] = "M" + m2.group(1)
bar2mid = {gsm2bar[g]: gsm2mid[g] for g in gsm2mid if g in gsm2bar}
bet = pd.read_csv(os.path.join(DATA, "GSE183647_meningioma_betas.csv.gz"), index_col=0)
bet.columns = [bar2mid.get(str(c).strip(), None) for c in bet.columns]
bet = bet.iloc[:, [i for i, c in enumerate(bet.columns) if c is not None]]
log(f"UCSF betas {bet.shape}")
m21 = pd.read_excel(os.path.join(NG, "41588_2022_1061_MOESM21_ESM.xlsx"), index_col=0)
m21.columns = [str(c).strip() for c in m21.columns]
m21.index = [str(c).split("_")[0].strip() for c in m21.index]
panel200 = [c for c in m21.index if c in bet.index]
X200 = bet.loc[panel200, :].T
km_full = KMeans(n_clusters=3, n_init=20, random_state=SEED).fit_predict(
    ((X200 - X200.mean()) / X200.std()).fillna(0).values)
truth_u = pd.Series(km_full, index=list(bet.columns))
log(f"reference kmeans dist: {pd.Series(truth_u).value_counts().sort_index().to_dict()}")
Zb = z_with_scaler(bet)
truth_list = [truth_u.get(x, None) for x in bet.columns]
pu, du = predict_panel(bet.T, sel)
res1 = metrics("P1_frozen_raw", "UCSF-SeSAMe", truth_list, pu, du, THR_P1)
rows.append(res1)
p2, d2 = predict_shared(Zb, sel)
res2 = metrics("P2_refnorm_shared", "UCSF-SeSAMe", truth_list, p2, d2, THR_P2)
rows.append(res2)
cc_u3 = [c for c in robust if c in Zb.columns]
p3, d3 = _centroid_predict_on(Zb[cc_u3].values, cc_u3)
res3 = metrics("P3_robust_panel", "UCSF-SeSAMe", truth_list, p3, d3, THR_P3)
rows.append(res3)
preds_u = pd.DataFrame({"MID": list(bet.columns), "km_ref": truth_list,
                        "pred_P1": pu, "pred_P2": p2, "pred_P3": p3,
                        "dist_P1": du, "dist_P2": d2, "dist_P3": d3})

# ---- POST HOC: supervised evaluation against PUBLISHED Group labels ----
# Label chain (all published data, no labels generated by us):
#   M-ID -> fingerprint column j -> MOESM20 row j (positional bridge, validated
#   unsupervisedly) -> unique MOESM19 record via exact 8-field clinical match
#   (201/201) -> published Group label. Not part of the original analysis plan;
#   labelled as post hoc wherever reported.
m19 = pd.read_excel(os.path.join(NG, "41588_2022_1061_MOESM19_ESM.xlsx"))
m20 = pd.read_excel(os.path.join(NG, "41588_2022_1061_MOESM20_ESM.xlsx"))
shared = [c for c in m20.columns if c in m19.columns]
k19 = m19[shared].astype(str).agg("|".join, axis=1)
idx19 = {k: i for i, k in enumerate(k19)}
mid2j = {str(c).strip(): j for j, c in enumerate(m21.columns)}
pub = {}
for mid, j in mid2j.items():
    if j < len(m20):
        kk = "|".join(str(x) for x in m20.iloc[j][shared])
        if kk in idx19:
            pub[mid] = str(m19["Group"].values[idx19[kk]])
preds_u["pub_group"] = preds_u["MID"].map(pub)
lab_mask = preds_u["pub_group"].notna()
log(f"POST HOC supervised eval: published Group labels available for {int(lab_mask.sum())}/200 profiles; "
    f"group dist {preds_u.loc[lab_mask, 'pub_group'].value_counts().to_dict()}")
sup_rows = []
for path in ["pred_P1", "pred_P2", "pred_P3"]:
    t = preds_u.loc[lab_mask, "pub_group"].astype(str)
    pp = preds_u.loc[lab_mask, path].astype(str)
    keep = pp != "NA"
    t, pp = t[keep], pp[keep]
    ari = round(float(adjusted_rand_score(t, pp)), 3)
    bm, mp = best_map_acc(t, pp)
    pm = np.array([mp.get(x, x) for x in pp], dtype=object)
    ct = pd.crosstab(pp, t)
    # hypermitotic recovery: fraction of Hypermitotic cases in the best-matching predicted class
    sup_rows.append({"path": path, "n": int(len(t)), "ARI": ari,
                     "bestmap_acc": round(float(bm), 3), "map": str(mp),
                     "mapped_acc": round(float((t == pm).mean()), 3),
                     "mapped_macroF1": round(float(f1_score(t, pm, average="macro", zero_division=0)), 3),
                     "crosstab": str(ct.to_dict())})
    log(f"POST HOC {path} vs published Group: n={len(t)} ARI={ari} bestmap={bm:.3f} "
        f"acc={(t == pm).mean():.3f} map={mp}")
    log("crosstab pred x Group:\n" + ct.to_string())
pd.DataFrame(sup_rows).to_csv(os.path.join(OUT, "ucsfsesame_supervised_vs_published.csv"), index=False)
preds_u.to_csv(os.path.join(OUT, "harmonization_pred_ucsfsesame.csv"), index=False)
log("saved supervised-evaluation table + per-sample predictions with published labels")

# raw direct transfer (diagnose_transfer.py output): ARI + bestmap on its pred x km
dg = pd.read_csv(os.path.join(OUT, "natgen201_transfer_diagnosis.csv"))
bm, mp = best_map_acc(dg["km"].astype(str), dg["pred"].astype(str))
raw_ari = round(float(adjusted_rand_score(dg["km"].astype(str), dg["pred"].astype(str))), 3)
rows.append({"platform": "UCSF-SeSAMe", "path": "RAW_direct(diagnose)", "n": int(len(dg)),
             "ARI": raw_ari, "bestmap_acc": round(float(bm), 3), "map": str(mp),
             "mapped_acc": round(float(bm), 3), "mapped_macroF1": np.nan,
             "reject_n_profiled": np.nan, "reject_denom_profiled": int(len(dg)),
             "reject_rate_profiled": np.nan, "reject_n_labelled": np.nan,
             "reject_denom_labelled": np.nan, "reject_rate_labelled": np.nan,
             "reject_convention": "diagnostic run; no reject mechanism applied"})
log(f"UCSF raw direct: ARI={raw_ari} bestmap={bm:.3f} map={mp}")

# ================= PLATFORM 2: GSE189521 (truth = clinical MenG A/B/C) =================
log("=== GSE189521/EPIC ===")
p189 = os.path.join(DATA, "GSE189521_series_matrix.txt.gz")
geo_acc, titles = [], []
start = None
with gzip.open(p189, "rt", encoding="utf-8", errors="replace") as f:
    for i, line in enumerate(f):
        if line.startswith("!Sample_geo_accession"):
            geo_acc = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_title"):
            titles = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
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
tab = clin.set_index("Sample")["menG"]
gsm2sample = {g: re.sub(r"\s*\(methylation\)\s*", "", t).strip() for g, t in zip(geo_acc, titles)}
truth_189 = [tab.get(gsm2sample.get(c), None) for c in B189.columns]
truth_189 = [None if (isinstance(x, float) and np.isnan(x)) else x for x in truth_189]
log(f"189521 labelled cols: {sum(1 for x in truth_189 if x is not None)}/{len(truth_189)}")
Z189 = z_with_scaler(B189)
pj, dj = predict_panel(B189.T, sel)
rows.append(metrics("P1_frozen_raw", "GSE189521", truth_189, pj, dj, THR_P1))
q2, e2 = predict_shared(Z189, sel)
rows.append(metrics("P2_refnorm_shared", "GSE189521", truth_189, q2, e2, THR_P2))
cc_j3 = [c for c in robust if c in Z189.columns]
q3, e3 = _centroid_predict_on(Z189[cc_j3].values, cc_j3)
rows.append(metrics("P3_robust_panel", "GSE189521", truth_189, q3, e3, THR_P3))
preds_j = pd.DataFrame({"GSM": list(B189.columns),
                        "Sample": [gsm2sample.get(g) for g in B189.columns],
                        "menG": truth_189, "pred_P1": pj, "pred_P2": q2, "pred_P3": q3,
                        "dist_P1": dj, "dist_P2": e2, "dist_P3": e3})
# reject counts: denominator = ALL profiled samples (per original report() convention, not the labelled subset)
import json as _json
thr_all = _json.load(open(os.path.join(OUT, "methyl_train_thr.json")))
for tag, dist_col in [("P1", "dist_P1"), ("P2", "dist_P2"), ("P3", "dist_P3")]:
    d = preds_j[dist_col].values
    thr_p = {"P1": THR_P1, "P2": THR_P2, "P3": THR_P3}[tag]
    n_all = int(np.nansum(d > thr_p))
    lab = np.array([x is not None for x in truth_189])
    n_lab = int(np.nansum((d > thr_p)[lab]))
    log(f"189521 {tag}: reject all-profiles {n_all}/110 = {n_all/110:.4f}; within labelled {n_lab}/100 = {n_lab/100:.4f}")

pd.DataFrame(rows).to_csv(os.path.join(OUT, "harmonization_agreement_fixed.csv"), index=False)
preds_u.to_csv(os.path.join(OUT, "harmonization_pred_ucsfsesame.csv"), index=False)
preds_j.to_csv(os.path.join(OUT, "harmonization_pred_gse189521.csv"), index=False)
log("saved harmonization_agreement_fixed.csv + per-sample prediction tables")

# ---- reject-mechanism audit: agreement WITH vs WITHOUT rejected samples ----
log("=== reject-mechanism audit (agreement excluding low-confidence samples) ===")
aud = []

def agree(t, p):
    a = round(float(adjusted_rand_score(t, p)), 3)
    b, mp = best_map_acc(t, p)
    return a, round(float(b), 3)

# UCSF (kmeans reference)
t_all = preds_u["km_ref"].astype(str).values
for path, dcol in [("pred_P1", "dist_P1"), ("pred_P2", "dist_P2"), ("pred_P3", "dist_P3")]:
    p_all = preds_u[path].astype(str).values
    d_all = preds_u[dcol].values
    thr_p = {"P1": THR_P1, "P2": THR_P2, "P3": THR_P3}[path.split("_")[-1]]
    ev = (p_all != "NA") & np.isfinite(d_all)
    keep = ev & ~(d_all > thr_p)
    a_all, b_all = agree(t_all[ev], p_all[ev])
    a_kept, b_kept = agree(t_all[keep], p_all[keep])
    aud.append({"platform": "UCSF-SeSAMe(kmeans ref)", "path": path,
                "n_evaluable": int(ev.sum()), "n_rejected": int((ev & (d_all > thr_p)).sum()),
                "ARI_all": a_all, "bestmap_all": b_all,
                "n_retained": int(keep.sum()), "ARI_retained": a_kept, "bestmap_retained": b_kept,
                "delta_bestmap": round(b_kept - b_all, 3)})
    log(f"UCSF {path}: all bestmap={b_all}; retained bestmap={b_kept} (excluded {int((ev & (d_all > thr_p)).sum())})")
# GSE189521 (MenG reference)
tj = np.array([None if x is None or (isinstance(x, float) and np.isnan(x)) else str(x) for x in preds_j["menG"]], dtype=object)
for path, dcol in [("pred_P1", "dist_P1"), ("pred_P2", "dist_P2"), ("pred_P3", "dist_P3")]:
    p_all = preds_j[path].astype(str).values
    d_all = preds_j[dcol].values
    thr_p = {"P1": THR_P1, "P2": THR_P2, "P3": THR_P3}[path.split("_")[-1]]
    lab = np.array([x is not None for x in tj])
    ev = lab & (p_all != "NA") & np.isfinite(d_all)
    keep = ev & ~(d_all > thr_p)
    a_all, b_all = agree(tj[ev], p_all[ev])
    a_kept, b_kept = agree(tj[keep], p_all[keep])
    aud.append({"platform": "GSE189521(MenG ref)", "path": path,
                "n_evaluable": int(ev.sum()), "n_rejected": int((ev & (d_all > thr_p)).sum()),
                "ARI_all": a_all, "bestmap_all": b_all,
                "n_retained": int(keep.sum()), "ARI_retained": a_kept, "bestmap_retained": b_kept,
                "delta_bestmap": round(b_kept - b_all, 3)})
    log(f"189521 {path}: all bestmap={b_all}; retained bestmap={b_kept} (excluded {int((ev & (d_all > thr_p)).sum())})")
# UCSF vs published Group labels (post hoc chain), with vs without rejects
for path, dcol in [("pred_P1", "dist_P1"), ("pred_P2", "dist_P2"), ("pred_P3", "dist_P3")]:
    m = preds_u["pub_group"].notna() & (preds_u[path].astype(str) != "NA") & np.isfinite(preds_u[dcol])
    t_all = preds_u.loc[m, "pub_group"].astype(str).values
    p_all = preds_u.loc[m, path].astype(str).values
    d_all = preds_u.loc[m, dcol].values
    thr_p = {"P1": THR_P1, "P2": THR_P2, "P3": THR_P3}[path.split("_")[-1]]
    keep = ~(d_all > thr_p)
    a_all, b_all = agree(t_all, p_all)
    a_kept, b_kept = agree(t_all[keep], p_all[keep])
    aud.append({"platform": "UCSF-SeSAMe(published Group, post hoc)", "path": path,
                "n_evaluable": int(m.sum()), "n_rejected": int((~keep).sum()),
                "ARI_all": a_all, "bestmap_all": b_all,
                "n_retained": int(keep.sum()), "ARI_retained": a_kept, "bestmap_retained": b_kept,
                "delta_bestmap": round(b_kept - b_all, 3)})
    log(f"UCSF-published {path}: all bestmap={b_all}; retained bestmap={b_kept} (excluded {int((~keep).sum())})")
pd.DataFrame(aud).to_csv(os.path.join(OUT, "reject_mechanism_audit.csv"), index=False)
log("saved reject_mechanism_audit.csv")
print("HARMONIZATION_METRICS_FIXED_DONE")
