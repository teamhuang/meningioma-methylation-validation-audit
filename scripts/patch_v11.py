# -*- coding: utf-8 -*-
"""Patch harmonize_transfer.py (current half-edited state -> v11).

Current verified state:
- predict_panel docstring says raw-space but callers pass z-space (BUG), NaN->0.0 anchor wrong for raw space
- report() uses global THR (20k z-space P99) for ALL paths (panel mismatch for P3)
- UCSF P1 call already bet.T + thr=THR_P1 (THR_P1 undefined yet)
- 189521/270371 path_fns: P1 gets z-space matrix (BUG), no per-path thr

v11 changes (all thresholds training-only):
1. predict_panel: n_finite guard BEFORE imputation; NaN -> per-CpG training mean (raw-space neutral)
2. Panel-matched training P99 thresholds: THR_P1 (raw, 20k) / THR_P2 (=THR, z-20k) / THR_P3 (z, 2k robust)
   -> results/methyl_train_thr.json
3. report(thr=...); all call sites pass per-path thr
4. 189521/270371: P1 receives B189.T / B270.T (raw betas)
"""
import io

P = r"E:\脑膜瘤npj\scripts\harmonize_transfer.py"
s = io.open(P, encoding="utf-8").read()

# ---- 1. predict_panel: guard before imputation, training-mean anchor ----
old = """    Xv = Z[cc].values                                          # (N, K) raw beta
    d = ((np.nan_to_num(Xv[:, None, :], nan=0.0) - Cb[None, :, :]) ** 2).sum(axis=2)
    n_finite = np.isfinite(Xv).sum(axis=1)"""
new = """    Xv = Z[cc].values                                          # (N, K) raw beta
    n_finite = np.isfinite(Xv).sum(axis=1)
    # NaN imputation: per-CpG training mean (neutral anchor in raw beta space)
    Xv = np.where(np.isnan(Xv), mu[None, :], Xv)
    d = ((Xv[:, None, :] - Cb[None, :, :]) ** 2).sum(axis=2)"""
assert old in s, "block1"; s = s.replace(old, new)

# ---- 2. panel-matched thresholds + report(thr=...) ----
old = """def report(tag, pred, dist, truth, log_ct=True, relabel=None):"""
new = """# ---- panel-matched training-set P99 reject thresholds (training data only) ----
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
        thr = THR_P2"""
assert old in s, "block2"; s = s.replace(old, new)

old = """    rej = float(np.nanmean(dist > THR)) if sub.any() else float("nan")"""
new = """    rej = float(np.nanmean(dist > thr)) if sub.any() else float("nan")"""
assert old in s, "block3"; s = s.replace(old, new)

# ---- 4. UCSF P2/P3 per-path thr (P1 already bet.T + thr=THR_P1) ----
old = """r2 = report("P2_refnorm_shared", p2, d2, truth_u_samples)"""
new = """r2 = report("P2_refnorm_shared", p2, d2, truth_u_samples, thr=THR_P2)"""
assert old in s, "block4"; s = s.replace(old, new)

old = """r3 = report("P3_robust_panel", p3, d3, truth_u_samples)"""
new = """r3 = report("P3_robust_panel", p3, d3, truth_u_samples, thr=THR_P3)"""
assert old in s, "block5"; s = s.replace(old, new)

# ---- 5. 189521: P1 raw matrix; per-path thr ----
old = """path_fns = [
    ("P1_frozen_raw", lambda z: predict_panel(z, sel)),
    ("P2_refnorm_shared", lambda z: predict_shared(z, sel)),
    ("P3_robust_panel", lambda z: predict_centroids(z, [x for x in robust if x in set(z.columns)])),
]
for tag, fn in path_fns:
    p, d, c = fn(Z189)
    r = report(tag, p, d, truth_189_s, relabel=relabel_189)"""
new = """path_fns = [
    ("P1_frozen_raw", lambda: predict_panel(B189.T, sel), THR_P1),
    ("P2_refnorm_shared", lambda: predict_shared(Z189, sel), THR_P2),
    ("P3_robust_panel", lambda: predict_centroids(Z189, [x for x in robust if x in set(Z189.columns)]), THR_P3),
]
for tag, fn, thr_p in path_fns:
    p, d, c = fn()
    r = report(tag, p, d, truth_189_s, thr=thr_p, relabel=relabel_189)"""
assert old in s, "block6"; s = s.replace(old, new)

# ---- 6. 270371: P1 raw matrix; per-path thr ----
old = """path_fns_270 = [
    ("P1_frozen_raw", lambda z: predict_panel(z, sel)),
    ("P2_refnorm_shared", lambda z: predict_shared(z, sel)),
    ("P3_robust_panel", lambda z: predict_centroids(z, [x for x in robust if x in set(z.columns)])),
]
for tag, fn in path_fns_270:
    p, d, c = fn(Z270)
    rej = float(np.nanmean(d > THR))"""
new = """path_fns_270 = [
    ("P1_frozen_raw", lambda: predict_panel(B270.T, sel), THR_P1),
    ("P2_refnorm_shared", lambda: predict_shared(Z270, sel), THR_P2),
    ("P3_robust_panel", lambda: predict_centroids(Z270, [x for x in robust if x in set(Z270.columns)]), THR_P3),
]
for tag, fn, thr_p in path_fns_270:
    p, d, c = fn()
    rej = float(np.nanmean(d > thr_p))"""
assert old in s, "block7"; s = s.replace(old, new)

old = """                    "n_reject": int(np.nansum(d > THR)), "reject_rate": round(rej, 3),"""
new = """                    "n_reject": int(np.nansum(d > thr_p)), "reject_rate": round(rej, 3),"""
assert old in s, "block8"; s = s.replace(old, new)

io.open(P, "w", encoding="utf-8").write(s)
print("PATCH OK: all 8 blocks applied")
