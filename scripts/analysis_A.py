import os, re, gzip, warnings
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from scipy.stats import chi2

DATA = r"E:\脑膜瘤npj\data"
OUT = r"E:\脑膜瘤npj\results"
os.makedirs(OUT, exist_ok=True)

# ================= 1. Clinical =================
clin = pd.read_csv(os.path.join(OUT, "clinical_parsed.csv"), dtype={"Sample": str})
clin["Sample"] = clin["Sample"].str.strip()

def xlsx_to_counts(s):
    m = re.match(r"^(\d{2})-(\d{2})$", s)
    if m: return f"{m.group(1)}-01-{int(m.group(2)):03d}"
    return s
clin["counts_id"] = clin["Sample"].map(xlsx_to_counts)

# ================= 2. RNA counts + symbol mapping (version-stripped) =================
counts = pd.read_csv(os.path.join(DATA, "GSE189672_RNAseq_raw_counts.csv.gz"), index_col=0)
counts.columns = [c.strip() for c in counts.columns]
counts = counts.loc[:, [c for c in counts.columns if c in set(clin["counts_id"])]]
clin = clin[clin["counts_id"].isin(counts.columns)].copy()
counts = counts.loc[:, clin["counts_id"].tolist()]
counts.columns = clin["Sample"].tolist()

gmap = pd.read_csv(os.path.join(DATA, "gencode_v38_ensg_symbol.csv"), dtype=str)
gmap["ensg_base"] = gmap["ensg_id"].str.split(".").str[0]
gmap = gmap.drop_duplicates("ensg_base").set_index("ensg_base")["symbol"]

base_ids = pd.Index([str(x).split(".")[0] for x in counts.index])
sym = base_ids.map(gmap)
mask = pd.notna(sym)
counts = counts[np.asarray(mask)]
counts.index = np.array(sym)[np.asarray(mask)]
tot = counts.sum(axis=1)
counts = counts.iloc[np.argsort(-tot.values, kind="stable")]
counts = counts[~counts.index.duplicated(keep="first")]
print("counts:", counts.shape)

lib = counts.sum(axis=0)
cpm = counts.div(lib, axis=1) * 1e6
logcpm = np.log2(cpm + 1)
logcpm = logcpm.loc[(cpm > 1).sum(axis=1) >= 10]
print("genes kept:", logcpm.shape)

# ================= 3. Features =================
GENESETS = {
    "myeloid": ["CD14","LST1","FCGR3A","C1QA","C1QB","C1QC","LYZ","S100A8","S100A9","CD68","MS4A7","CST3","TYROBP"],
    "tcell_nk": ["CD3D","CD3E","CD8A","CD2","GZMB","NKG7","GNLY","CCL5","IL7R"],
    "bcell": ["MS4A1","CD79A","MZB1","CD19"],
    "endothelial": ["PECAM1","VWF","CLDN5","CDH5"],
    "ecm_fibro": ["COL1A1","COL1A2","COL3A1","COL6A1","COL6A2","COL6A3","FN1","DCN","LUM","BGN"],
    "proliferation": ["MKI67","TOP2A","PCNA","MCM2","MCM4","MCM6","CCNB1","UBE2C"],
    "hypoxia": ["CA9","VEGFA","SLC2A1","HK2","LDHA","NDRG1"],
}
def zscore(v):
    v = np.asarray(v, dtype=float)
    sd = np.nanstd(v)
    return (v - np.nanmean(v)) / (sd if sd > 1e-12 else 1.0)

def singscore(mat, genes):
    g = [x for x in genes if x in mat.index]
    if len(g) < 3:
        return pd.Series(np.nan, index=mat.columns), 0
    ranks_full = mat.rank(axis=0)
    med_all = ranks_full.median(axis=0)
    score = (ranks_full.loc[g].mean(axis=0) - med_all) / mat.shape[0] * 1000
    return score, len(g)

feat = pd.DataFrame(index=clin.index)
nmap = {}
for k, gs in GENESETS.items():
    s, n = singscore(logcpm, gs)
    feat["score_" + k] = zscore(s.values)
    nmap[k] = n
print("geneset sizes mapped:", nmap)

hv = logcpm.var(axis=1).sort_values(ascending=False).index[:5000]
X_hv = logcpm.loc[hv].T
X_hvz = (X_hv - X_hv.mean()) / X_hv.std()
pcs = PCA(n_components=5, random_state=0).fit_transform(X_hvz.values)
for i in range(5):
    feat[f"PC{i+1}"] = zscore(pcs[:, i])
print("PVE PCs1-5:", np.round(PCA(n_components=5, random_state=0).fit(X_hvz.values).explained_variance_ratio_, 3))

hv2 = logcpm.var(axis=1).sort_values(ascending=False).index[:2000]
X2 = logcpm.loc[hv2].T
X2z = (X2 - X2.mean()) / X2.std()
L = linkage(pdist(X2z.values, metric="correlation"), method="ward")
feat["rna_cluster"] = fcluster(L, t=3, criterion="maxclust")
print("RNA cluster x MenG:")
print(pd.crosstab(feat["rna_cluster"], clin["menG"].values).to_string())

# ================= 4. Survival =================
surv = pd.DataFrame(index=clin.index)
surv["T"] = clin["fu_months"].values
surv["E"] = clin["rec"].values
surv["grade2"] = (clin["grade"] == 2).astype(int).values
surv["gtr"] = clin["gtr"].values
surv["age"] = zscore(clin["age"].values)
surv["male"] = (clin["Gender"].astype(str).str.upper() == "M").astype(int).values
known = clin["menG"].isin(["A", "B", "C"]).values
surv["menG_B"] = np.where(known, (clin["menG"] == "B").astype(float), np.nan)
surv["menG_C"] = np.where(known, (clin["menG"] == "C").astype(float), np.nan)
surv = surv.join(feat)
surv.to_csv(os.path.join(OUT, "survival_input.csv"), index=False)
print("N =", len(surv), "events =", int(surv["E"].sum()),
      "| menG known:", int(known.sum()))

CLIN = ["grade2", "gtr", "age", "male"]
MODELS = {
    "M0_clinical": (CLIN, 0.05),
    "MA_clinical+MenG": (CLIN + ["menG_B", "menG_C"], 0.05),
    "M1_clinical+RNA_PC": (CLIN + ["PC1", "PC2", "PC3"], 0.05),
    "M1b_RNA_PC_only": (["PC1", "PC2", "PC3"], 0.05),
    "M2_clinical+ecosystem": (CLIN + ["score_myeloid", "score_ecm_fibro", "score_proliferation", "score_tcell_nk"], 0.05),
    "M2b_ecosystem_only": (["score_myeloid", "score_ecm_fibro", "score_proliferation", "score_tcell_nk"], 0.05),
}

def fit_cox(df, cols, penalizer=0.05):
    d = df[["T", "E"] + cols].dropna().reset_index(drop=True)
    last_err = None
    for pen, step in [(penalizer, 1.0), (penalizer, 0.25), (penalizer * 5, 0.25)]:
        try:
            cph = CoxPHFitter(penalizer=pen)
            cph.fit(d, duration_col="T", event_col="E", fit_options={"step_size": step})
            return cph, d
        except Exception as e:
            last_err = e
    raise RuntimeError(f"cox fit failed ({last_err}) cols={cols}")

def cv_cindex(df, cols, penalizer, n_rep=20, k=5):
    rng = np.random.default_rng(11)
    base = df[["T", "E"] + cols].dropna().reset_index(drop=True)
    scores = []
    for _ in range(n_rep):
        idx = rng.permutation(len(base))
        for f in np.array_split(idx, k):
            tr, te = base.iloc[np.setdiff1d(idx, f)], base.iloc[f]
            if te["E"].sum() == 0 or tr["E"].sum() < 3:
                continue
            try:
                m = CoxPHFitter(penalizer=penalizer)
                m.fit(tr, "T", "E")
                scores.append(concordance_index(te["T"], -m.predict_partial_hazard(te), te["E"]))
            except Exception:
                continue
    return float(np.mean(scores)), float(np.std(scores)), len(base)

rows = []
for name, (cols, pen) in MODELS.items():
    try:
        cph, d = fit_cox(surv, cols, pen)
        c_in = concordance_index(d["T"], -cph.predict_partial_hazard(d), d["E"])
        rng = np.random.default_rng(7)
        bs = []
        for _ in range(1000):
            db = d.iloc[rng.choice(len(d), len(d), replace=True)]
            try:
                m = CoxPHFitter(penalizer=pen)
                m.fit(db, "T", "E")
                bs.append(concordance_index(db["T"], -m.predict_partial_hazard(db), db["E"]))
            except Exception:
                continue
        lo, hi = (np.nan, np.nan) if len(bs) < 50 else np.percentile(bs, [2.5, 97.5])
        cvm, cvs, n_used = cv_cindex(surv, cols, pen)
        rows.append({"model": name, "n": n_used, "events": int(d["E"].sum()),
                     "C_insample": round(c_in, 3), "C_boot_lo": round(lo, 3), "C_boot_hi": round(hi, 3),
                     "C_cv": round(cvm, 3), "C_cv_sd": round(cvs, 3)})
        print(f"{name}: n={n_used} ev={int(d['E'].sum())} C_in={c_in:.3f} ({lo:.3f}-{hi:.3f}) C_cv={cvm:.3f}+-{cvs:.3f}")
    except Exception as e:
        print(f"{name}: FAILED {type(e).__name__}: {str(e)[:100]}")
        rows.append({"model": name, "n": np.nan, "events": np.nan, "C_insample": np.nan,
                     "C_boot_lo": np.nan, "C_boot_hi": np.nan, "C_cv": np.nan, "C_cv_sd": np.nan})
pd.DataFrame(rows).to_csv(os.path.join(OUT, "cindex_models.csv"), index=False)

def lr(df, cols_a, cols_b):
    sub = df[["T", "E"] + cols_a + cols_b].dropna()
    ca, _ = fit_cox(sub, cols_a)
    cb, _ = fit_cox(sub, cols_b)
    stat = 2 * (cb.log_likelihood_ - ca.log_likelihood_)
    return stat, chi2.sf(stat, len(cols_b) - len(cols_a)), len(sub)

print()
for name, extra in [("MA add MenG", ["menG_B", "menG_C"]),
                    ("M1 add PCs", ["PC1", "PC2", "PC3"]),
                    ("M2 add eco", ["score_myeloid", "score_ecm_fibro", "score_proliferation", "score_tcell_nk"])]:
    try:
        s, p, n = lr(surv, CLIN, CLIN + extra)
        print(f"LR {name}: chi2={s:.2f} df={len(extra)} p={p:.4g} n={n}")
    except Exception as e:
        print(f"LR {name}: FAILED {str(e)[:80]}")

cph, _ = fit_cox(surv, CLIN + ["score_myeloid", "score_ecm_fibro", "score_proliferation", "score_tcell_nk"], 0.05)
cph.summary.to_csv(os.path.join(OUT, "cox_final_summary.csv"))
print("\nFinal clinical+ecosystem model (penalized):")
print(cph.summary[["coef", "exp(coef)", "p"]].round(4).to_string())

cph2, _ = fit_cox(surv, CLIN + ["menG_B", "menG_C"], 0.05)
cph2.summary.to_csv(os.path.join(OUT, "cox_meng_summary.csv"))
print("\nClinical+MenG model (penalized):")
print(cph2.summary[["coef", "exp(coef)", "p"]].round(4).to_string())

# ================= 5. KM =================
kmf = KaplanMeierFitter()
def med(sub):
    kmf.fit(sub["T"], sub["E"]); return kmf.median_survival_time_
print("\nKM median RFS (months): overall =", round(med(surv), 1))
for g, sub in surv.groupby("grade2"):
    print(f"  WHO II={g}: median={med(sub):.1f} n={len(sub)} ev={int(sub['E'].sum())}")
menG_lab = np.where(known, clin["menG"].values, np.nan)
for g in ["A", "B", "C"]:
    sub = surv[menG_lab == g]
    print(f"  MenG {g}: median={med(sub):.1f} n={len(sub)} ev={int(sub['E'].sum())}")
np.save(os.path.join(OUT, "menG_labels.npy"), menG_lab)

sc = surv[[c for c in surv.columns if c.startswith("score_")]].copy()
sc["menG"] = menG_lab
sc["rec"] = surv["E"].values
score_cols = [c for c in surv.columns if c.startswith("score_")]
print("\nScores by MenG (mean):")
print(sc.groupby("menG")[score_cols].mean().round(2).to_string())
print("\nScores by recurrence (0=censored,1=recurred):")
print(sc.groupby("rec")[score_cols].mean().round(2).to_string())
print("\nANALYSIS_A_DONE")
