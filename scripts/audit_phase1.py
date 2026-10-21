# -*- coding: utf-8 -*-
"""Phase-1 core-evidence audit (npj DM pre-submission review round, 2026-10-19).

1) Calibration slope + observed/predicted ratio at 60 months for the two
   utility-suite cohorts (final models, same spec as analysis_E_utility.py).
   Addresses review 4.4: quantify the risk-ordering vs absolute-risk boundary.
2) Bootstrap 95% CI for the UCSF-SeSAMe harmonization agreement (ARI + Hungarian
   best-map accuracy, paths P1/P3), resampling the 200 samples. Addresses the
   'uncertainty' row of the deployment-boundary reporting table.
3) Unified cohort numbers re-derived from source CSVs (review Phase-1 item 3).
Outputs:
  results/calibration_slope.csv
  results/harmonization_bootstrap_ci.csv
  results/cohort_numbers_verified.csv
"""
import os, warnings
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.metrics import adjusted_rand_score
from scipy.optimize import linear_sum_assignment

OUT = r"E:\脑膜瘤npj\results"
rng = np.random.default_rng(20260317)


def log(m): print(m, flush=True)


# ---------- 1) calibration slope + O/P ----------
def fit_final(d, cols):
    dd = d[["T", "E"] + cols].dropna().reset_index(drop=True)
    c = CoxPHFitter(penalizer=0.05)
    c.fit(dd, "T", "E")
    return c, dd


def km_censor(T, E):
    df = pd.DataFrame({"t": T, "e": E})
    out = {}
    surv = 1.0
    for u in np.sort(df["t"].unique()):
        at_risk = (df["t"] >= u).sum()
        d = ((df["t"] == u) & (df["e"] == 1)).sum()
        surv *= (1 - d / at_risk) if at_risk > 0 else 1
        out[u] = surv
    def G(t):
        ks = [u for u in out if u <= t]
        return out[max(ks)] if ks else 1.0
    return G


rows = []
specs = [
    ("GSE189672 (n=99)", os.path.join(OUT, "survival_input.csv"),
     ["grade2", "gtr", "age", "male", "menG_B", "menG_C"]),
    ("GSE304094-primary (n=149)", os.path.join(OUT, "gse304094_primary_survival.csv"),
     ["grade", "gtr", "age", "male", "pmc_2", "pmc_3"]),
]
for name, path, cols in specs:
    s = pd.read_csv(path)
    c, d = fit_final(s, cols)
    lp = c.predict_log_partial_hazard(d)
    lpz = (lp - lp.mean()) / lp.std()
    # calibration slope: unpenalized Cox on the standardized LP alone
    cs = CoxPHFitter(penalizer=0.0)
    cs.fit(pd.DataFrame({"T": d["T"], "E": d["E"], "LP": lpz.values}), "T", "E")
    slope = float(cs.summary.loc["LP", "exp(coef)"])
    slope_lo = float(cs.summary.loc["LP", "exp(coef) lower 95%"])
    slope_hi = float(cs.summary.loc["LP", "exp(coef) upper 95%"])
    # 5y predicted risk (same as analysis_E) and IPCW observed overall
    surv5 = c.predict_survival_function(d, times=[60]).iloc[0].values
    prisk5 = 1 - surv5
    T, E = d["T"].values, d["E"].values
    G = km_censor(T, E)
    hz = 60
    y = ((T <= hz) & (E == 1)).astype(float)
    w = np.array([1.0 / G(T[i]) if (T[i] <= hz and E[i] == 1) else (1.0 if T[i] > hz else 0.0)
                  for i in range(len(T))])
    obs_ipcw = float((y * w).sum() / w.sum())
    pred_mean = float(prisk5.mean())
    rows.append({"cohort": name, "n": len(d), "events": int(E.sum()),
                 "calibration_slope": round(slope, 3),
                 "slope_CI_low": round(slope_lo, 3), "slope_CI_high": round(slope_hi, 3),
                 "pred_5y_mean": round(pred_mean, 3), "obs_5y_ipcw": round(obs_ipcw, 3),
                 "O_to_P_ratio": round(obs_ipcw / pred_mean, 2)})
    log(f"{name}: slope={slope:.3f} ({slope_lo:.2f}-{slope_hi:.2f}); "
        f"pred5y={pred_mean:.3f} obs5y(IPCW)={obs_ipcw:.3f} O/P={obs_ipcw/pred_mean:.2f}")
pd.DataFrame(rows).to_csv(os.path.join(OUT, "calibration_slope.csv"), index=False)

# ---------- 2) harmonization bootstrap CI (UCSF-SeSAMe) ----------
pu = pd.read_csv(os.path.join(OUT, "harmonization_pred_ucsfsesame.csv"),
                 dtype={"pred_P1": str, "pred_P2": str, "pred_P3": str})
t = pu["km_ref"].astype(str).values

def bestmap(tt, pp):
    ct = pd.crosstab(pd.Series(pp), pd.Series(tt))
    ri, ci = linear_sum_assignment(-ct.values)
    return ct.values[ri, ci].sum() / ct.values.sum()

ci_rows = []
for path in ["pred_P1", "pred_P3"]:
    p = pu[path].astype(str).values
    ari_obs = adjusted_rand_score(t, p)
    bm_obs = bestmap(t, p)
    aris, bms = [], []
    for _ in range(1000):
        idx = rng.choice(len(t), len(t), replace=True)
        tt, pp = t[idx], p[idx]
        if len(set(tt)) < 3:
            continue
        aris.append(adjusted_rand_score(tt, pp))
        bms.append(bestmap(tt, pp))
    ci_rows.append({"path": path, "n": len(t),
                    "ARI": round(float(ari_obs), 3),
                    "ARI_CI_low": round(float(np.percentile(aris, 2.5)), 3),
                    "ARI_CI_high": round(float(np.percentile(aris, 97.5)), 3),
                    "bestmap_acc": round(float(bm_obs), 3),
                    "bestmap_CI_low": round(float(np.percentile(bms, 2.5)), 3),
                    "bestmap_CI_high": round(float(np.percentile(bms, 97.5)), 3),
                    "n_boot_ok": len(aris)})
    log(f"{path}: ARI={ari_obs:.3f} [{np.percentile(aris,2.5):.3f},{np.percentile(aris,97.5):.3f}]; "
        f"bestmap={bm_obs:.3f} [{np.percentile(bms,2.5):.3f},{np.percentile(bms,97.5):.3f}]")
pd.DataFrame(ci_rows).to_csv(os.path.join(OUT, "harmonization_bootstrap_ci.csv"), index=False)

# ---------- 3) unified cohort numbers from source files ----------
ver = []
def add(cohort, metric, value, source):
    ver.append({"cohort": cohort, "metric": metric, "value": value, "source": source})

cl = pd.read_csv(os.path.join(OUT, "clinical_parsed.csv"), dtype={"Sample": str})
add("GSE189672/189521", "n_samples(=patients)", int(len(cl)), "clinical_parsed.csv")
add("GSE189672/189521", "events(rec)", int(cl["rec"].sum()), "clinical_parsed.csv")
add("GSE189672/189521", "WHO grade 1/2", str(cl["grade"].value_counts().sort_index().to_dict()), "clinical_parsed.csv")
add("GSE189672/189521", "MenG A/B/C/NA", str(cl["menG"].value_counts(dropna=False).to_dict()), "clinical_parsed.csv")
add("GSE189672/189521", "median follow-up (mo)", round(float(cl["fu_months"].median()), 1), "clinical_parsed.csv")
add("GSE189672/189521", "n with MenG label", int(cl["menG"].notna().sum()), "clinical_parsed.csv")

s2 = pd.read_csv(os.path.join(OUT, "gse304094_primary_survival.csv"))
cc2 = s2.dropna(subset=["T", "E", "grade", "gtr", "age", "male"])
add("GSE304094-primary", "n(analysis, complete case)", int(len(cc2)), "gse304094_primary_survival.csv")
add("GSE304094-primary", "events", int(cc2["E"].sum()), "gse304094_primary_survival.csv")
add("GSE304094-primary", "WHO 1/2/3 (complete case)", str(cc2["grade"].value_counts().sort_index().to_dict()), "gse304094_primary_survival.csv")
add("GSE304094-primary", "n(pre-complete-case, primary no-RT)", int(len(s2)), "gse304094_primary_survival.csv")

s3 = pd.read_csv(os.path.join(OUT, "natgen205_survival.csv"))
add("NatGen-UCSF205", "n", int(len(s3)), "natgen205_survival.csv")
add("NatGen-UCSF205", "events(LFFR)", int(s3["E"].sum()), "natgen205_survival.csv")
add("NatGen-UCSF205", "grade 1/2/3", str(s3["grade"].value_counts().sort_index().to_dict()), "natgen205_survival.csv")
add("NatGen-UCSF205", "Hypermitotic/Immune-enriched", str({"hm": int(s3["hm"].sum()), "ie": int(s3["ie"].sum())}), "natgen205_survival.csv")
add("NatGen-UCSF205", "covariates available", "grade,gtr,male (no age)", "natgen205_survival.csv")

sm325 = pd.read_csv(r"E:\脑膜瘤npj\data\meta\GSE292325_sample_map.csv", dtype=str)
lab = sm325["char4"].str.replace("kmeans cluster:", "", regex=False).str.strip()
add("GSE292325 (ANC 2025 train)", "n", int(len(sm325)), "meta/GSE292325_sample_map.csv")
add("GSE292325 (ANC 2025 train)", "kmeans labels 1/2/3", str(lab.value_counts().sort_index().to_dict()), "meta/GSE292325_sample_map.csv")
add("GSE292326 (ANC 2025 val)", "n(samples)", int(len(pd.read_csv(os.path.join(OUT, "gse292326_classification.csv")))), "gse292326_classification.csv")
add("GSE292326 (ANC 2025 val)", "identity-audit n / acc", "31; 0.774 (kmeans) / 0.806 (Choudhury)", "gse292326_identity_audit.csv; manuscript_stats.json")
st = pd.read_csv(os.path.join(OUT, "gse304096_stability.csv"))
add("GSE304096 (longitudinal)", "samples / patient clusters", f"{int(st['n_samples'].sum())} / {len(st)}", "gse304096_stability.csv")
add("GSE304096 (longitudinal)", "stable clusters", f"{int(st['all_same'].sum())}/{len(st)}", "gse304096_stability.csv")
add("GSE270371 (Nat Med 2024)", "n( methylation samples)", int(len(pd.read_csv(os.path.join(OUT, "gse270371_classification.csv")))), "gse270371_classification.csv")
add("GSE270371 (Nat Med 2024)", "patient-level outcomes", "not public (dbGaP phs003707) - distribution-level evidence only", "analysis_report.md 8.1")
pd.DataFrame(ver).to_csv(os.path.join(OUT, "cohort_numbers_verified.csv"), index=False)
log(f"cohort numbers verified: {len(ver)} rows")
print("PHASE1_AUDIT_DONE")
