import os, gzip, warnings, time, re
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index
from scipy.stats import chi2

t0 = time.time()
DATA = r"E:\脑膜瘤npj\data"
OUT = r"E:\脑膜瘤npj\results"
NG = os.path.join(DATA, "suppl", "natgen2022")
def log(m): print(f"[{time.time()-t0:6.0f}s] {m}", flush=True)
rng = np.random.default_rng(20260317)

# ---------- 1. clinical tables ----------
c19 = pd.read_excel(os.path.join(NG, "41588_2022_1061_MOESM19_ESM.xlsx"))
c20 = pd.read_excel(os.path.join(NG, "41588_2022_1061_MOESM20_ESM.xlsx"))
s2e = pd.read_excel(os.path.join(NG, "41588_2022_1061_MOESM6_ESM.xlsx"), sheet_name="Sfigure_2e")
m21 = pd.read_excel(os.path.join(NG, "41588_2022_1061_MOESM21_ESM.xlsx"), index_col=0)
m21.columns = [str(c).strip() for c in m21.columns]
log(f"MOESM19 {c19.shape}; MOESM20 {c20.shape}; Sfig2e {s2e.shape}; MOESM21 {m21.shape}")

# ---------- 2. bridge: barcode -> UCSF M-ID ----------
fl = pd.read_csv(os.path.join(DATA, "suppl", "GSE183647_filelist.txt"), sep="\t", comment="#", header=None,
                 names=["kind", "name", "time", "size", "type"])
fl = fl[fl["kind"] == "File"]
gsm2bar = {}
for nm in fl["name"]:
    parts = str(nm).split("_")
    if parts[0].startswith("GSM") and len(parts) >= 3:
        gsm2bar[parts[0]] = parts[1] + "_" + parts[2]

gsm2mid = {}
with gzip.open(os.path.join(DATA, "GSE183647_series_matrix.txt.gz"), "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith("!Sample_geo_accession"):
            acc = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_title"):
            tit = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
for g, t in zip(acc, tit):
    m = re.match(r"UCSF Meningioma (\d+) \(methylation\)", t)
    if m:
        gsm2mid[g] = "M" + m.group(1)
bar2mid = {gsm2bar[g]: gsm2mid[g] for g in gsm2mid if g in gsm2bar}
log(f"bridge entries (UCSF only): {len(bar2mid)}")

# ---------- 3. load betas, rename columns to M-ID ----------
log("loading 565 betas (3.2GB gz) ...")
bet = pd.read_csv(os.path.join(DATA, "GSE183647_meningioma_betas.csv.gz"), index_col=0)
bet.columns = [bar2mid.get(str(c).strip(), None) for c in bet.columns]
keep = [i for i, c in enumerate(bet.columns) if c is not None]
bet = bet.iloc[:, keep]
log(f"betas (UCSF M-ID columns): {bet.shape}")

# ---------- 4. fingerprint check MOESM21 vs full betas ----------
shared_ids = [c for c in m21.columns if c in bet.columns]
log(f"MOESM21 M-IDs in betas: {len(shared_ids)}/201")
cpgs = [c for c in m21.index[:500] if c in bet.index]
A = m21.loc[cpgs, shared_ids].values.astype(float)
B = bet.loc[cpgs, shared_ids].values.astype(float)
corr = np.corrcoef(A.ravel(), B.ravel())[0, 1]
log(f"fingerprint corr = {corr:.4f}; mean|diff|={np.abs(A - B).mean():.4f}")

# ---------- 5. positional validation of MOESM20 rows vs MOESM21 columns ----------
X21 = m21.T.copy()
X21.index = [str(c).strip() for c in X21.index]
km21 = KMeans(n_clusters=3, n_init=20, random_state=20260317).fit_predict((X21 - X21.mean()) / X21.std())
chk = pd.DataFrame({"km": km21,
                    "mib1": pd.to_numeric(c20["MIB1-LI"], errors="coerce").values[:X21.shape[0]],
                    "rec": pd.to_numeric(c20["Recurrent"], errors="coerce").values[:X21.shape[0]],
                    "grade": pd.to_numeric(c20["WHO grade"], errors="coerce").values[:X21.shape[0]]}).dropna()
tab = chk.groupby("km").agg(n=("mib1", "size"), mean_MIB1=("mib1", "mean"),
                            recurrent_rate=("rec", "mean"), mean_grade=("grade", "mean"))
log("positional check (MOESM21 k-means vs MOESM20 clinical, positionally):\n" + tab.round(2).to_string())

# also positional check on full betas (200 samples, frozen k-means in full space)
sel_full = [c for c in X21.index if c in bet.columns]
Xf = bet.loc[:, sel_full].T
mad_f = (Xf - Xf.median()).abs().median()
naf = Xf.isna().mean(axis=0)
selcpg = mad_f[naf < 0.2].sort_values(ascending=False).index[:20000]
kmf565 = KMeans(n_clusters=3, n_init=20, random_state=20260317).fit_predict(
    ((Xf[selcpg] - Xf[selcpg].mean()) / Xf[selcpg].std()).fillna(0).values)
chk2 = pd.DataFrame({"km": kmf565,
                     "mib1": chk.set_index(chk.index).loc[[X21.index.get_loc(i) for i in Xf.index], "mib1"].values,
                     "rec": chk.set_index(chk.index).loc[[X21.index.get_loc(i) for i in Xf.index], "rec"].values}).dropna()
log("positional check (full-betas k-means vs MOESM20):\n" + chk2.groupby("km").agg(
    n=("mib1", "size"), mean_MIB1=("mib1", "mean"), recurrent_rate=("rec", "mean")).round(2).to_string())

# ---------- 6. frozen-classifier models on the 200 UCSF samples ----------
scaler = pd.read_csv(os.path.join(OUT, "methyl_scaler.csv"), index_col=0)
centroids = pd.read_csv(os.path.join(OUT, "methyl_centroids.csv"), index_col=0)
labs = list(centroids.columns)
selset = [c for c in scaler.index if c in bet.index]
B201 = bet.loc[selset, sel_full]
Zb = ((B201.T - scaler["mean"].reindex(B201.index)) / scaler["sd"].reindex(B201.index)).dropna(axis=1, how='any')
cc2 = [c for c in Zb.columns if c in centroids.index]
M2 = centroids.loc[cc2, labs].values.T
D2 = ((Zb[cc2].values[:, None, :] - M2[None, :, :]) ** 2).sum(axis=2)
pred201 = np.array([labs[i] for i in D2.argmin(axis=1)])
log(f"frozen-classifier prediction (n={len(pred201)}): {pd.Series(pred201).value_counts().to_dict()}")

# ---------- 7. survival (T=LFFR years, E=LF) ----------
surv = pd.DataFrame({
    "T": pd.to_numeric(c20["LFFR"], errors="coerce").values[:len(pred201)],
    "E": pd.to_numeric(c20["LF"], errors="coerce").values[:len(pred201)],
    "grade": pd.to_numeric(c20["WHO grade"], errors="coerce").values[:len(pred201)],
    "gtr": np.where(c20["EOR"].astype(str).str.strip().str.upper().str.startswith("G"), 1,
            np.where(c20["EOR"].astype(str).str.strip().str.upper().str.startswith("S"), 0, np.nan))[:len(pred201)],
    "male": (pd.to_numeric(c20["Sex"], errors="coerce") == 2).astype(int).values[:len(pred201)],
    "mib1": pd.to_numeric(c20["MIB1-LI"], errors="coerce").values[:len(pred201)],
    "pmc_2": (pred201 == "2").astype(float),
    "pmc_3": (pred201 == "3").astype(float),
    "MID": sel_full,
})
surv = surv.dropna(subset=["T", "E"])
log(f"n survival: n={len(surv)}, events={int(surv['E'].sum())}, grades={surv['grade'].value_counts().sort_index().to_dict()}")

CLIN_F = ["grade", "gtr", "male"]
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

rows = []
for name, cols in [("M0_clinical", CLIN_F), ("M1_clinical+predMethylClass", CLIN_F + ["pmc_2", "pmc_3"]),
                   ("M1b_predMethylClass_only", ["pmc_2", "pmc_3"]),
                   ("M2_clinical+MIB1", CLIN_F + ["mib1"])]:
    cph, d = fit(surv, cols)
    c_in = concordance_index(d["T"], -cph.predict_partial_hazard(d), d["E"])
    bs = []
    for _ in range(1000):
        db = d.iloc[rng.choice(len(d), len(d), replace=True)]
        try:
            m = CoxPHFitter(penalizer=0.05); m.fit(db, "T", "E")
            bs.append(concordance_index(db["T"], -m.predict_partial_hazard(db), db["E"]))
        except Exception:
            continue
    lo, hi = np.percentile(bs, [2.5, 97.5])
    rows.append({"model": name, "n": len(d), "events": int(d["E"].sum()),
                 "C": round(c_in, 3), "CI": f"({lo:.3f}-{hi:.3f})"})
    log(f"{name}: n={len(d)} ev={int(d['E'].sum())} C={c_in:.3f} ({lo:.3f}-{hi:.3f})")
pd.DataFrame(rows).to_csv(os.path.join(OUT, "cindex_natgen201.csv"), index=False)

common_df = surv[["T", "E"] + CLIN_F + ["pmc_2", "pmc_3"]].dropna()
ca, _ = fit(common_df, CLIN_F)
cb, _ = fit(common_df, CLIN_F + ["pmc_2", "pmc_3"])
stat = 2 * (cb.log_likelihood_ - ca.log_likelihood_)
log(f"LR add predicted class: chi2={stat:.2f} p={chi2.sf(stat, 2):.4g}")

cph, _ = fit(surv, CLIN_F + ["pmc_2", "pmc_3"])
cph.summary.to_csv(os.path.join(OUT, "cox_natgen201_summary.csv"))
log("final model:\n" + cph.summary[["coef", "exp(coef)", "p"]].round(4).to_string())

kmf = KaplanMeierFitter()
pc_series = pd.Series(pred201, index=surv.index)
for k in ["1", "2", "3"]:
    sub = surv[pc_series == k]
    if len(sub) == 0:
        continue
    kmf.fit(sub["T"], sub["E"])
    log(f"class {k}: n={len(sub)} ev={int(sub['E'].sum())} median={kmf.median_survival_time_:.1f}y")

strat = []
STRAT_C = ["gtr", "male"]
for g, sub in surv.groupby("grade"):
    if len(sub) < 20 or sub["E"].sum() < 5:
        log(f"grade {g}: skipped (n={len(sub)}, ev={int(sub['E'].sum())})")
        continue
    for name, cols in [("M0", STRAT_C), ("M1", STRAT_C + ["pmc_2", "pmc_3"])]:
        cph, d = fit(sub, cols)
        c_in = concordance_index(d["T"], -cph.predict_partial_hazard(d), d["E"])
        strat.append({"grade": int(g), "model": name, "n": len(d), "events": int(d["E"].sum()), "C": round(c_in, 3)})
        log(f"grade {g} {name}: n={len(d)} ev={int(d['E'].sum())} C={c_in:.3f}")
pd.DataFrame(strat).to_csv(os.path.join(OUT, "cindex_natgen201_gradestrat.csv"), index=False)
surv.to_csv(os.path.join(OUT, "natgen201_survival.csv"), index=False)
log("ANALYSIS_F_DONE")

