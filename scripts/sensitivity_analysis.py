"""Sensitivity analysis C (P1): does the prognostic increment hold under
alt endpoint definitions, resection, RT, WHO grade version, and QC thresholds?

Pre-registered: the main model is unchanged; these are pre-specified
exploratory checks reporting the same paired delta-C framework, plus
cohort restriction checks. No post-hoc model selection.
"""
import os, time, warnings
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

t0 = time.time()
OUT = r"E:\脑膜瘤npj\results"
DATA = r"E:\脑膜瘤npj\data"
def log(m): print(f"[{time.time()-t0:6.0f}s] {m}", flush=True)
rng = np.random.default_rng(20260317)

def fit(df, cols, pen=0.05):
    d = df[["T", "E"] + cols].dropna().reset_index(drop=True)
    # drop covariates that are constant within this subset (e.g. gtr==1 in the
    # GTR-restricted stratum) — a zero-variance column makes Cox singular
    dropped = [c for c in cols if d[c].nunique() < 2]
    if dropped:
        d = d.drop(columns=dropped)
        cols = [c for c in cols if c not in dropped]
    last = None
    for pen_, step in [(pen, 1.0), (pen, 0.25), (pen*5, 0.25)]:
        try:
            c = CoxPHFitter(penalizer=pen_); c.fit(d, "T", "E", fit_options={"step_size": step})
            return c, d, dropped
        except Exception as e:
            last = e
    raise RuntimeError(f"fit failed {cols} :: {str(last)[:60]}")

def paired_deltaC(df, clin, mol, n_boot=500):
    """Paired bootstrap delta-C; returns dict with 'note' if data insufficient
    or model fails.  Never raises on small/degenerate subsets."""
    d = df[["T","E"] + clin + mol].dropna().reset_index(drop=True)
    if len(d) < 15 or d["E"].sum() < 5:
        return {"n": len(d), "events": int(d["E"].sum()),
                "C_clinical": np.nan, "C_plus_molecular": np.nan,
                "deltaC": np.nan, "deltaC_CI_low": np.nan,
                "deltaC_CI_high": np.nan, "P_delta<=0": np.nan,
                "note": "insufficient data (n<15 or events<5)"}
    try:
        c0, d0, drop0 = fit(d, clin)
        c1, d1, drop1 = fit(d, clin + mol)
        # both models must be scored on the SAME rows; fit() may have dropped
        # constant covariates (and their rows after dropna on other cols)
        common = d0.index.intersection(d1.index)
        if len(common) < len(d0) or len(common) < len(d1):
            d_all = d.loc[common]
            c0, d0, _ = fit(d_all, clin); c1, d1, _ = fit(d_all, clin + mol)
        # d_boot keeps ALL original columns (mol included) so bootstrap fits see
        # the full design; d0 lacks columns that fit() dropped as constants
        d_boot = d.loc[d0.index].copy()
        a0 = concordance_index(d0["T"], -c0.predict_partial_hazard(d0), d0["E"])
        a1 = concordance_index(d1["T"], -c1.predict_partial_hazard(d1), d1["E"])
        delta_point = a1 - a0
    except Exception as e:
        return {"n": len(d), "events": int(d["E"].sum()),
                "C_clinical": np.nan, "C_plus_molecular": np.nan,
                "deltaC": np.nan, "deltaC_CI_low": np.nan,
                "deltaC_CI_high": np.nan, "P_delta<=0": np.nan,
                "note": f"model fit failed: {str(e)[:60]}"}
    deltas = []
    for _ in range(n_boot):
        db = d_boot.iloc[rng.choice(len(d_boot), len(d_boot), replace=True)].reset_index(drop=True)
        try:
            m0, db0, _ = fit(db, clin)
            m1, db1, _ = fit(db, clin + mol)
            deltas.append(concordance_index(db1["T"],-m1.predict_partial_hazard(db1),db1["E"])
                          - concordance_index(db0["T"],-m0.predict_partial_hazard(db0),db0["E"]))
        except Exception:
            continue
    if len(deltas) == 0:
        return {"n": len(d), "events": int(d["E"].sum()),
                "C_clinical": round(a0,3), "C_plus_molecular": round(a1,3),
                "deltaC": round(delta_point,3),
                "deltaC_CI_low": np.nan, "deltaC_CI_high": np.nan,
                "P_delta<=0": np.nan, "note": "all bootstrap fits failed"}
    dl = np.array(deltas)
    lo, hi = np.percentile(dl, [2.5, 97.5])
    p_boot = float((dl <= 0).mean())
    return {"n": len(d), "events": int(d["E"].sum()),
            "C_clinical": round(a0,3), "C_plus_molecular": round(a1,3),
            "deltaC": round(delta_point,3), "deltaC_CI_low": round(lo,3),
            "deltaC_CI_high": round(hi,3), "P_delta<=0": round(p_boot,4)}

rows = []

# =============== COHORT 1: Tübingen GSE304094-primary (n=158, 55 events) ===============
s2 = pd.read_csv(os.path.join(OUT, "gse304094_primary_survival.csv"))
clin2 = ["grade","gtr","age","male"]
mol2 = ["pmc_2","pmc_3"]
log(f"Tübingen: n={len(s2)} events={int(s2['E'].sum())} "
    f"adjuvant_RT={s2['adjuvant_RT'].value_counts().to_dict()} "
    f"WHO16={s2['CNS_WHO_Grade_2016'].value_counts().to_dict()} "
    f"WHO21={s2['CNS_WHO_Grade_2021'].value_counts().to_dict()}")

# S1. RT restriction
for tag, sub in [
    ("S1a_no_adjuvant_RT", s2[s2["adjuvant_RT"]=="no"]),
    ("S1b_adjuvant_RT_only", s2[s2["adjuvant_RT"]=="yes"]),
]:
    r = paired_deltaC(sub, clin2, mol2)
    r["check"] = tag; r["cohort"] = "Tübingen"
    rows.append(r)
    log(f"{tag}: n={len(sub)} ev={int(sub['E'].sum())}")

# S2. WHO grade version (2021 vs 2016)
log(f"S2: WHO2016 vs 2021 concordance: "
    f"{(s2['CNS_WHO_Grade_2016']==s2['CNS_WHO_Grade_2021']).mean():.3f} "
    f"mismatches={int((s2['CNS_WHO_Grade_2016']!=s2['CNS_WHO_Grade_2021']).sum())}")
s2b = s2.copy()
s2b["grade"] = s2["CNS_WHO_Grade_2016"]
r = paired_deltaC(s2b, clin2, mol2)
r["check"] = "S2_WHO2016"; r["cohort"] = "Tübingen"
rows.append(r)
log(f"S2_WHO2016: C0={r['C_clinical']} C1={r['C_plus_molecular']} "
    f"dC={r['deltaC']} CI=({r['deltaC_CI_low']},{r['deltaC_CI_high']})")

# S3. Simpson/resection
for tag, sub in [
    ("S3a_Simpson_GTR", s2[s2["simpson"]==1]),
    ("S3b_Simpson_NON_GTR", s2[s2["simpson"]>1]),
]:
    r = paired_deltaC(sub, clin2, mol2)
    r["check"] = tag; r["cohort"] = "Tübingen"
    rows.append(r)
    log(f"{tag}: n={len(sub)} ev={int(sub['E'].sum())}")

# S4. Grade-restricted (WHO1 only)
sub = s2[s2["CNS_WHO_Grade_2021"]==1]
r = paired_deltaC(sub, ["gtr","age","male"], mol2)
r["check"] = "S4_WHO1_only"; r["cohort"] = "Tübingen"
rows.append(r)
log(f"S4_WHO1_only: n={len(sub)} ev={int(sub['E'].sum())}")

# S5. Endpoint definition
log("S5: Pattern_of_failure distribution:")
print(s2["Pattern_of_failure"].value_counts(dropna=False).to_dict())
if s2["Pattern_of_failure"].notna().any():
    log(f"S5: {int(s2['Pattern_of_failure'].notna().sum())} events have pattern")
    # S5b: endpoint-definition sensitivity — events are all progressions
    # (local / distant / "local and distant"); build (a) any-progression and
    # (b) local-only-recurrence endpoints, both with censored patients kept.
    s2_temp = s2.copy()
    patt = s2_temp["Pattern_of_failure"].astype(str).str.strip().str.lower()
    is_prog = patt.isin(["local", "distant", "local and distant"])
    s2_temp["_localonly"] = patt.eq("local")
    if is_prog.any():
        sub_any = s2_temp.copy()  # main endpoint already = any progression
        sub_loc = s2_temp.copy()
        sub_loc["E"] = np.where(sub_loc["_localonly"], 1, 0)  # local-only events
        for tag, sub in [("S5b_any_progression", sub_any), ("S5b_local_only", sub_loc)]:
            r = paired_deltaC(sub, clin2, mol2)
            r["check"] = tag; r["cohort"] = "Tübingen"
            rows.append(r)
            log(f"{tag}: n={len(sub)} ev={int(sub['E'].sum())}")
    else:
        log("S5b: no recognizable progression patterns")
else:
    log("S5: no Pattern_of_failure recorded")

# S6. CpG QC threshold sensitivity (270371, distribution only, no endpoint)
log("S6: CpG QC threshold sensitivity (GSE270371, distribution only)")
sel = pd.read_csv(os.path.join(OUT, "methyl_selected_cpgs.csv"))["cpg"].tolist()
scaler = pd.read_csv(os.path.join(OUT, "methyl_scaler.csv"), index_col=0)
centroids = pd.read_csv(os.path.join(OUT, "methyl_centroids.csv"), index_col=0)
# GSE270371: rows=CpGs, cols alternate intensity + "Detection Pval.N"; keep intensity only.
B270 = pd.read_csv(os.path.join(DATA, "GSE270371_Matrix_processed.txt.gz"),
                   sep="\t", index_col=0).astype(np.float32)
B270 = B270[[c for c in B270.columns if not str(c).startswith("Detection Pval")]]
B270 = B270[~B270.index.duplicated(keep="first")]
m270 = B270
sel270 = [c for c in sel if c in m270.index]
log(f"270371 available CpGs: {len(sel270)}/{len(sel)}")
for na_thr in [0.10, 0.20, 0.30]:
    # per-CpG missingness across the 994 samples
    sub = m270.loc[sel270]
    na_frac = sub.isna().mean(axis=1)  # axis=1 -> per CpG (row), over samples
    kept = na_frac[na_frac < na_thr]  # index = CpG ids
    if len(kept) < 50:
        log(f"QC<{na_thr}: kept={len(kept)} (too few, skip)")
        continue
    Z270 = ((sub.loc[kept.index] - scaler["mean"].reindex(kept.index))
            / scaler["sd"].reindex(kept.index)).T   # rows=994 samples, cols=CpGs
    common = [c for c in Z270.columns if c in centroids.index]
    M = centroids.loc[common, list(centroids.columns)].values.T   # (3, K)
    Xv = Z270[common].values.astype(np.float32)                   # (994, K)
    # NaN guard: frozen transfer, per-sample distance over finite CpGs
    fin = np.isfinite(Xv)
    Xf = np.nan_to_num(Xv, nan=0.0)
    D = ((Xf[:, None, :] - M[None, :, :])**2).sum(axis=2)         # (994, 3)
    n_fin = fin.sum(axis=1)
    ok = n_fin >= 0.5 * len(common)
    pred = np.where(ok, [list(centroids.columns)[i] for i in D.argmin(axis=1)], "NA")
    d_min = np.where(ok, D.min(axis=1), np.nan)
    log(f"QC<{na_thr}: kept={len(kept)}/{len(sel270)} "
        f"class_dist={pd.Series(pred).value_counts().to_dict()} "
        f"dist_P50={np.nanpercentile(d_min,50):.1f} P95={np.nanpercentile(d_min,95):.1f}")

# =============== COHORT 2: UCSF205 (n=201, 76 events) ===============
log("=== UCSF205 sensitivity ===")
s3 = pd.read_csv(os.path.join(OUT, "natgen205_survival.csv"))
clin3 = ["grade","gtr","male"]
mol3 = ["hm","ie"]
log(f"UCSF205: n={len(s3)} events={int(s3['E'].sum())}")
for tag, sub in [
    ("S1a_Simpson_GTR", s3[s3["gtr"]==1]),
    ("S1b_Simpson_NON_GTR", s3[s3["gtr"]==0]),
]:
    r = paired_deltaC(sub, clin3, mol3)
    r["check"] = tag; r["cohort"] = "UCSF205"
    rows.append(r)
    log(f"{tag}: n={len(sub)} ev={int(sub['E'].sum())}")
sub = s3[s3["grade"]==1]
r = paired_deltaC(sub, ["gtr","male"], mol3)
r["check"] = "S4_WHO1_only"; r["cohort"] = "UCSF205"
rows.append(r)
log(f"UCSF205 WHO1 only: n={len(sub)} ev={int(sub['E'].sum())}")

pd.DataFrame(rows).to_csv(os.path.join(OUT, "sensitivity_analysis.csv"), index=False)
log(f"SAVED {len(rows)} sensitivity checks")
log("SENSITIVITY_DONE")
