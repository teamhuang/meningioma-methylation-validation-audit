import os, warnings
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index
from sklearn.metrics import adjusted_rand_score
from scipy.stats import chi2

OUT = r"E:\脑膜瘤npj\results"
clin = pd.read_csv(os.path.join(OUT, "clinical_parsed.csv"), dtype={"Sample": str})
cls = pd.read_csv(os.path.join(OUT, "gse189521_classification.csv"), dtype={"title": str})
cls["Sample"] = cls["title"].str.replace(r"\s*\(methylation\)\s*", "", regex=True).str.strip()
mg = cls.merge(clin, on="Sample", how="inner")
print("linked:", len(mg), "of", len(cls))

mg["pmc"] = mg["predicted_methyl_class"].astype(int).astype(str)
print("\nPredicted methyl class x clinical MenG:")
ct = pd.crosstab(mg["pmc"], mg["menG"])
print(ct.to_string())
both = mg.dropna(subset=["menG"])
ari = adjusted_rand_score(both["menG"], both["pmc"])
agree = (both["menG"] == both["pmc"]).mean()
print(f"\nARI(predicted class, MenG) = {ari:.3f}; exact agreement = {agree:.3f} (n={len(both)})")

# RFS by predicted methylation class
surv = pd.DataFrame({
    "T": mg["fu_months"], "E": mg["rec"],
    "grade2": (mg["grade"] == 2).astype(int),
    "gtr": mg["gtr"],
    "age": (mg["age"] - mg["age"].mean()) / mg["age"].std(),
    "male": (mg["Gender"].astype(str).str.upper() == "M").astype(int),
    "pmc_2": (mg["pmc"] == "2").astype(float),
    "pmc_3": (mg["pmc"] == "3").astype(float),
})

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

CLIN = ["grade2", "gtr", "age", "male"]
rows = []
for name, cols in [("M0_clinical", CLIN),
                   ("M3_clinical+predMethylClass", CLIN + ["pmc_2", "pmc_3"]),
                   ("M3b_predMethylClass_only", ["pmc_2", "pmc_3"])]:
    cph, d = fit(surv, cols)
    c_in = concordance_index(d["T"], -cph.predict_partial_hazard(d), d["E"])
    rng = np.random.default_rng(7); bs = []
    for _ in range(1000):
        db = d.iloc[rng.choice(len(d), len(d), replace=True)]
        try:
            m = CoxPHFitter(penalizer=0.05); m.fit(db, "T", "E"); bs.append(concordance_index(db["T"], -m.predict_partial_hazard(db), db["E"]))
        except Exception: continue
    lo, hi = np.percentile(bs, [2.5, 97.5])
    rows.append({"model": name, "n": len(d), "events": int(d["E"].sum()),
                 "C_insample": round(c_in, 3), "CI": f"({lo:.3f}-{hi:.3f})"})
    print(f"{name}: n={len(d)} ev={int(d['E'].sum())} C={c_in:.3f} ({lo:.3f}-{hi:.3f})")
pd.DataFrame(rows).to_csv(os.path.join(OUT, "cindex_predmethyl_models.csv"), index=False)

# LR: add predicted class to clinical
common = surv[["T", "E"] + CLIN + ["pmc_2", "pmc_3"]].dropna()
ca, _ = fit(common, CLIN)
cb, _ = fit(common, CLIN + ["pmc_2", "pmc_3"])
stat = 2 * (cb.log_likelihood_ - ca.log_likelihood_)
print(f"LR add predicted methyl class: chi2={stat:.2f} df=2 p={chi2.sf(stat, 2):.4g} n={len(common)}")

cph, _ = fit(surv, CLIN + ["pmc_2", "pmc_3"])
cph.summary.to_csv(os.path.join(OUT, "cox_predmethyl_summary.csv"))
print("\nClinical + predicted methyl class:")
print(cph.summary[["coef", "exp(coef)", "p"]].round(4).to_string())

kmf = KaplanMeierFitter()
def med(sub):
    kmf.fit(sub["T"], sub["E"]); return kmf.median_survival_time_
print("\nKM median RFS by predicted methyl class:")
for k, sub in surv.join(mg[["pmc"]]).groupby("pmc"):
    print(f"  class {k}: median={med(sub):.1f} n={len(sub)} ev={int(sub['E'].sum())}")

mg.to_csv(os.path.join(OUT, "gse189521_linked_endpoints.csv"), index=False)
print("\nANALYSIS_C_DONE")
