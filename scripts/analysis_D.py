import os, re, gzip, warnings, time
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index
from scipy.stats import chi2

t0 = time.time()
DATA = r"E:\脑膜瘤npj\data"
OUT = r"E:\脑膜瘤npj\results"
def log(m): print(f"[{time.time()-t0:6.0f}s] {m}", flush=True)

# ---------- 1. GSM -> barcode (from RAW filelist) ----------
fl = pd.read_csv(os.path.join(DATA, "suppl", "GSE304094_filelist.txt"), sep="\t", comment="#", header=None,
                 names=["kind", "name", "time", "size", "type"])
fl = fl[fl["kind"] == "File"]
gsm2bar = {}
for nm in fl["name"]:
    parts = str(nm).split("_")
    if len(parts) >= 3 and parts[0].startswith("GSM"):
        gsm2bar[parts[0]] = parts[1] + "_" + parts[2]
log(f"filelist: {len(gsm2bar)} GSM->barcode entries")

# ---------- 2. GSM -> T26 ID (series matrix) ----------
gsm2t26 = {}
with gzip.open(os.path.join(DATA, "GSE304094_series_matrix.txt.gz"), "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!Sample_geo_accession"):
            acc = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_title"):
            tit = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!series_matrix_table_begin"):
            break
for g, t in zip(acc, tit):
    gsm2t26[g] = t
log(f"series matrix: {len(gsm2t26)} GSM->title")

# ---------- 3. load clinical (Supp Data 1) ----------
clin = pd.read_excel(os.path.join(DATA, "suppl", "protocadherin_suppldata.xlsx"),
                     sheet_name="Supplementary Data 1", header=None, skiprows=3)
clin.columns = pd.read_excel(os.path.join(DATA, "suppl", "protocadherin_suppldata.xlsx"),
                             sheet_name="Supplementary Data 1", header=None, skiprows=2).iloc[0]
clin = clin.dropna(subset=["Internal_ID"])
clin["T26"] = clin["Internal_ID"].astype(str).str.strip()
log(f"clinical rows: {len(clin)}; primary={int((clin['tumor_setting']=='Primary').sum())}")

# ---------- 4. load matrix, map barcodes -> T26 ----------
mat = pd.read_csv(os.path.join(DATA, "GSE304094_discovery_processed.txt.gz"), sep="\t", index_col=0)
mat.index = [str(x).split("_")[0] for x in mat.index]
mat = mat[~mat.index.duplicated(keep="first")]
mat.columns = [str(c).strip() for c in mat.columns]
log(f"matrix {mat.shape}")
bar2t26 = {b: gsm2t26.get(g) for g, b in gsm2bar.items()}
col_t26 = [bar2t26.get(c) for c in mat.columns]
keep = [i for i, t in enumerate(col_t26) if t is not None]
mat = mat.iloc[:, keep]
mat.columns = [col_t26[i] for i in keep]
log(f"matrix after T26 mapping: {mat.shape}")

# ---------- 5. predict classes with our classifier ----------
scaler = pd.read_csv(os.path.join(OUT, "methyl_scaler.csv"), index_col=0)
centroids = pd.read_csv(os.path.join(OUT, "methyl_centroids.csv"), index_col=0)
sel = pd.read_csv(os.path.join(OUT, "methyl_selected_cpgs.csv"))["cpg"].tolist()
idx = scaler.index.intersection(mat.index)
sub = mat.loc[idx]
Z = ((sub.T - scaler["mean"].reindex(idx)) / scaler["sd"].reindex(idx))
labs = list(centroids.columns)
common = [c for c in Z.columns if c in centroids.index]
M = centroids.loc[common, labs].values.T
Xv = Z[common].values
D = ((Xv[:, None, :] - M[None, :, :]) ** 2).sum(axis=2)
pred = np.array([labs[i] for i in D.argmin(axis=1)])
pred_s = pd.Series(pred, index=Z.index)
log(f"predicted classes: {pd.Series(pred).value_counts().to_dict()} (n={len(pred)})")

# ---------- 6. merge & survival ----------
mg = clin.merge(pred_s.rename("pmc"), left_on="T26", right_index=True, how="inner")
log(f"merged clinical+matrix: {len(mg)}")
mg["T"] = pd.to_numeric(mg["Follow-Up [Months]"], errors="coerce")
mg["E"] = pd.to_numeric(mg["Status (0=censored, 1=progression)"], errors="coerce")
mg["grade"] = pd.to_numeric(mg["CNS_WHO_Grade_2021"], errors="coerce")
mg["simpson"] = pd.to_numeric(mg["Simpson_Grade"], errors="coerce")
mg["gtr"] = np.where(mg["simpson"] <= 3, 1, np.where(mg["simpson"] > 3, 0, np.nan))
mg["age"] = (pd.to_numeric(mg["Age [Years]"], errors="coerce") - pd.to_numeric(mg["Age [Years]"], errors="coerce").mean())
mg["age"] = mg["age"] / mg["age"].std()
mg["male"] = (mg["Sex"].astype(str).str.upper().str[0] == "M").astype(int)
mg["prior_rt"] = mg["prior_RT"].astype(str).str.lower().eq("yes").astype(int)
mg["methhigh"] = mg["Methylation_Cluster"].astype(str).str.lower().eq("methhigh").astype(float)
mg["methhigh"] = mg["methhigh"].where(mg["Methylation_Cluster"].notna())
mg["tor_prolif"] = np.where(mg["Molecular_Group (by Toronto nomenclature)"].astype(str).str.contains("Prolif", case=False, na=False), 1.0,
                    np.where(mg["Molecular_Group (by Toronto nomenclature)"].notna(), 0.0, np.nan))
mg["tor_hypermet"] = np.where(mg["Molecular_Group (by Toronto nomenclature)"].astype(str).str.contains("Hypermet", case=False, na=False), 1.0,
                      np.where(mg["Molecular_Group (by Toronto nomenclature)"].notna(), 0.0, np.nan))
mg["pmc_2"] = (mg["pmc"] == "2").astype(float)
mg["pmc_3"] = (mg["pmc"] == "3").astype(float)

# primary-tumour cohort (first surgery), exclude prior RT
prim = mg[(mg["tumor_setting"] == "Primary") & (mg["prior_rt"] == 0)].copy()
log(f"PRIMARY cohort (no prior RT): n={len(prim)}, events={int(prim['E'].sum())}, "
    f"grades={prim['grade'].value_counts().sort_index().to_dict()}")

CLIN = ["grade", "gtr", "age", "male"]
def fit(df, cols, pen=0.05):
    d = df[["T", "E"] + cols].dropna().reset_index(drop=True)
    for pen_, step in [(pen, 1.0), (pen, 0.25), (pen * 5, 0.25)]:
        try:
            c = CoxPHFitter(penalizer=pen_)
            c.fit(d, "T", "E", fit_options={"step_size": step})
            return c, d
        except Exception:
            continue
    raise RuntimeError("fit failed " + str(cols))

def boot_ci(df, cols, pen=0.05, n=1000):
    cph, d = fit(df, cols, pen)
    c_in = concordance_index(d["T"], -cph.predict_partial_hazard(d), d["E"])
    rng = np.random.default_rng(7); bs = []
    for _ in range(n):
        db = d.iloc[rng.choice(len(d), len(d), replace=True)]
        try:
            m = CoxPHFitter(penalizer=pen); m.fit(db, "T", "E")
            bs.append(concordance_index(db["T"], -m.predict_partial_hazard(db), db["E"]))
        except Exception:
            continue
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return c_in, lo, hi, len(d), int(d["E"].sum()), cph

def cv_c(df, cols, pen=0.05, n_rep=20, k=5):
    rng = np.random.default_rng(11)
    base = df[["T", "E"] + cols].dropna().reset_index(drop=True)
    ss = []
    for _ in range(n_rep):
        idxp = rng.permutation(len(base))
        for f in np.array_split(idxp, k):
            tr, te = base.iloc[np.setdiff1d(idxp, f)], base.iloc[f]
            if te["E"].sum() == 0 or tr["E"].sum() < 3: continue
            try:
                m = CoxPHFitter(penalizer=pen); m.fit(tr, "T", "E")
                ss.append(concordance_index(te["T"], -m.predict_partial_hazard(te), te["E"]))
            except Exception:
                continue
    return float(np.mean(ss)) if ss else np.nan

rows = []
MODELS = {
    "M0_clinical": CLIN,
    "M1_clinical+predMethylClass": CLIN + ["pmc_2", "pmc_3"],
    "M1b_predMethylClass_only": ["pmc_2", "pmc_3"],
    "M2_clinical+METHhigh": CLIN + ["methhigh"],
    "M3_clinical+TorontoGroup": CLIN + ["tor_prolif", "tor_hypermet"],
}
for name, cols in MODELS.items():
    try:
        c_in, lo, hi, n, ev, _ = boot_ci(prim, cols)
        cvv = cv_c(prim, cols)
        rows.append({"model": name, "n": n, "events": ev, "C_insample": round(c_in, 3),
                     "CI": f"({lo:.3f}-{hi:.3f})", "C_cv": round(cvv, 3)})
        log(f"{name}: n={n} ev={ev} C={c_in:.3f} ({lo:.3f}-{hi:.3f}) CV={cvv:.3f}")
    except Exception as e:
        log(f"{name}: FAILED {str(e)[:80]}")
res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "cindex_304094_primary.csv"), index=False)

# LR test M0 vs M1
common_df = prim[["T", "E"] + CLIN + ["pmc_2", "pmc_3"]].dropna()
ca, _ = fit(common_df, CLIN)
cb, _ = fit(common_df, CLIN + ["pmc_2", "pmc_3"])
stat = 2 * (cb.log_likelihood_ - ca.log_likelihood_)
log(f"LR add predicted methyl class: chi2={stat:.2f} df=2 p={chi2.sf(stat, 2):.4g} n={len(common_df)}")

# final cox summary
cph, _ = fit(prim, CLIN + ["pmc_2", "pmc_3"])
cph.summary.to_csv(os.path.join(OUT, "cox_304094_predmethyl_summary.csv"))
log("final model (clinical + predicted class):")
print(cph.summary[["coef", "exp(coef)", "p"]].round(4).to_string())

# grade-stratified validation (drop constant 'grade' within strata)
strat_rows = []
STRAT_CLIN = ["gtr", "age", "male"]
for g, sub in prim.groupby("grade"):
    if len(sub) < 20 or sub["E"].sum() < 5:
        log(f"grade {g}: skipped (n={len(sub)}, ev={int(sub['E'].sum())})")
        continue
    for name, cols in [("M0", STRAT_CLIN), ("M1", STRAT_CLIN + ["pmc_2", "pmc_3"])]:
        try:
            c_in, lo, hi, n, ev, _ = boot_ci(sub, cols, n=500)
            strat_rows.append({"grade": int(g), "model": name, "n": n, "events": ev,
                               "C": round(c_in, 3), "CI": f"({lo:.3f}-{hi:.3f})"})
            log(f"grade {g} {name}: n={n} ev={ev} C={c_in:.3f} ({lo:.3f}-{hi:.3f})")
        except Exception as e:
            log(f"grade {g} {name}: FAILED {str(e)[:60]}")
pd.DataFrame(strat_rows).to_csv(os.path.join(OUT, "cindex_304094_gradestrat.csv"), index=False)

# KM by predicted class
kmf = KaplanMeierFitter()
km_out = []
for k, sub in prim.dropna(subset=["T", "E"]).groupby("pmc"):
    kmf.fit(sub["T"], sub["E"])
    km_out.append({"class": k, "n": len(sub), "events": int(sub["E"].sum()),
                   "median_RFS_months": round(float(kmf.median_survival_time_), 1)})
kmdf = pd.DataFrame(km_out)
kmdf.to_csv(os.path.join(OUT, "km_304094_by_class.csv"), index=False)
log("KM by predicted class: " + kmdf.to_string(index=False))

prim.to_csv(os.path.join(OUT, "gse304094_primary_survival.csv"), index=False)
log("ANALYSIS_D_DONE")
