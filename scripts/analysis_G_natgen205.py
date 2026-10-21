import os, time, warnings
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index
from scipy.stats import chi2

t0 = time.time()
OUT = r"E:\脑膜瘤npj\results"
NG = r"E:\脑膜瘤npj\data\suppl\natgen2022"
def log(m): print(f"[{time.time()-t0:6.0f}s] {m}", flush=True)
rng = np.random.default_rng(20260317)

c19 = pd.read_excel(os.path.join(NG, "41588_2022_1061_MOESM19_ESM.xlsx"))
log(f"MOESM19: {c19.shape}; Group counts: {c19['Group'].value_counts().to_dict()}")

surv = pd.DataFrame({
    "T": pd.to_numeric(c19["LFFR"], errors="coerce"),
    "E": pd.to_numeric(c19["LF"], errors="coerce"),
    "grade": pd.to_numeric(c19["WHO grade"], errors="coerce"),
    "gtr": np.where(c19["EOR"].astype(str).str.strip().str.upper().str.startswith("G"), 1,
            np.where(c19["EOR"].astype(str).str.strip().str.upper().str.startswith("S"), 0, np.nan)),
    "male": (pd.to_numeric(c19["Sex"], errors="coerce") == 2).astype(int),
    "mib1": pd.to_numeric(c19["MIB1-LI"], errors="coerce"),
    "hm": (c19["Group"] == "Hypermitotic").astype(float),
    "ie": (c19["Group"] == "Immune-enriched").astype(float),
})
surv = surv.dropna(subset=["T", "E", "hm", "ie"])
log(f"n={len(surv)}, events={int(surv['E'].sum())}, grades={surv['grade'].value_counts().sort_index().to_dict()}")

CLIN = ["grade", "gtr", "male"]
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
for name, cols in [("M0_clinical", CLIN), ("M1_clinical+publishedGroup", CLIN + ["hm", "ie"]),
                   ("M2_clinical+MIB1", CLIN + ["mib1"]),
                   ("M1b_publishedGroup_only", ["hm", "ie"])]:
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
pd.DataFrame(rows).to_csv(os.path.join(OUT, "cindex_natgen205.csv"), index=False)

common_df = surv[["T", "E"] + CLIN + ["hm", "ie"]].dropna()
ca, _ = fit(common_df, CLIN)
cb, _ = fit(common_df, CLIN + ["hm", "ie"])
stat = 2 * (cb.log_likelihood_ - ca.log_likelihood_)
log(f"LR add published Group: chi2={stat:.2f} p={chi2.sf(stat, 2):.4g}")

cph, _ = fit(surv, CLIN + ["hm", "ie"])
cph.summary.to_csv(os.path.join(OUT, "cox_natgen205_summary.csv"))
log("final model (published labels):\n" + cph.summary[["coef", "exp(coef)", "p"]].round(4).to_string())

# grade-stratified (published labels)
strat = []
STRAT_C = ["gtr", "male"]
for g, sub in surv.groupby("grade"):
    if len(sub) < 20 or sub["E"].sum() < 5:
        log(f"grade {g}: skipped (n={len(sub)}, ev={int(sub['E'].sum())})")
        continue
    for name, cols in [("M0", STRAT_C), ("M1", STRAT_C + ["hm", "ie"])]:
        cph, d = fit(sub, cols)
        c_in = concordance_index(d["T"], -cph.predict_partial_hazard(d), d["E"])
        strat.append({"grade": int(g), "model": name, "n": len(d), "events": int(d["E"].sum()), "C": round(c_in, 3)})
        log(f"grade {g} {name}: n={len(d)} ev={int(d['E'].sum())} C={c_in:.3f}")
pd.DataFrame(strat).to_csv(os.path.join(OUT, "cindex_natgen205_gradestrat.csv"), index=False)

kmf = KaplanMeierFitter()
for gname, col in [("Merlin-intact", "#2a78d6"), ("Immune-enriched", "#1baf7a"), ("Hypermitotic", "#eda100")]:
    sub = surv[c19.loc[surv.index, "Group"] == gname]
    if len(sub) == 0: continue
    kmf.fit(sub["T"], sub["E"])
    log(f"{gname}: n={len(sub)} ev={int(sub['E'].sum())} median={kmf.median_survival_time_:.1f}y")

surv.to_csv(os.path.join(OUT, "natgen205_survival.csv"), index=False)
log("ANALYSIS_G_DONE")
