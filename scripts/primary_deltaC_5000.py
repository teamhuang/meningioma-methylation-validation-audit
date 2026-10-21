# -*- coding: utf-8 -*-
"""Primary ΔC comparisons at n_boot = 5,000 with archived bootstrap distributions.

Recomputes the five pre-specified primary comparisons (published-label routes
plus the transferred-classifier secondary routes) with 5,000 paired bootstrap
resamples (both models refitted on each resample of the same patients), and
archives:
  results/primary_deltaC_5000.csv          summary: one-sided P, two-sided P, CI
  results/primary_deltaC_distributions.npz per-comparison ΔC resample vectors

Data preparation replicates paired_deltaC.py exactly (same dropna, same
covariates, same seed family 20260317).
"""
import os, warnings, time
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

t0 = time.time()
OUT = r"E:\脑膜瘤npj\results"
rng = np.random.default_rng(20260317)
N_BOOT = 5000


def log(m): print(f"[{time.time()-t0:7.0f}s] {m}", flush=True)


def fit(df, cols, pen=0.05):
    d = df[["T", "E"] + cols].dropna().reset_index(drop=True)
    last = None
    for pen_, step in [(pen, 1.0), (pen, 0.25), (pen * 5, 0.25)]:
        try:
            c = CoxPHFitter(penalizer=pen_)
            c.fit(d, "T", "E", fit_options={"step_size": step})
            return c, d
        except Exception as e:
            last = e
    raise RuntimeError("fit failed " + str(cols) + " :: " + str(last)[:80])


def paired_deltaC(df, clin, mol, label):
    d = df[["T", "E"] + clin + mol].dropna().reset_index(drop=True)
    c0, _ = fit(d, clin)
    c1, _ = fit(d, clin + mol)
    c0_obs = concordance_index(d["T"], -c0.predict_partial_hazard(d), d["E"])
    c1_obs = concordance_index(d["T"], -c1.predict_partial_hazard(d), d["E"])
    deltas = np.empty(N_BOOT)
    ok = 0
    for b in range(N_BOOT):
        db = d.iloc[rng.choice(len(d), len(d), replace=True)].reset_index(drop=True)
        try:
            m0 = CoxPHFitter(penalizer=0.05); m0.fit(db[["T", "E"] + clin], "T", "E")
            m1 = CoxPHFitter(penalizer=0.05); m1.fit(db[["T", "E"] + clin + mol], "T", "E")
            a = concordance_index(db["T"], -m0.predict_partial_hazard(db), db["E"])
            bb = concordance_index(db["T"], -m1.predict_partial_hazard(db), db["E"])
            deltas[ok] = bb - a
            ok += 1
        except Exception:
            continue
    deltas = deltas[:ok]
    lo, hi = np.percentile(deltas, [2.5, 97.5])
    p_one = float((deltas <= 0).mean())
    p_two = float(min(1.0, 2 * min((deltas <= 0).mean(), (deltas >= 0).mean())))
    res = {"comparison": label, "n": len(d), "events": int(d["E"].sum()),
           "C_clinical": round(c0_obs, 3), "C_plus_molecular": round(c1_obs, 3),
           "deltaC": round(c1_obs - c0_obs, 3),
           "CI_low": round(lo, 3), "CI_high": round(hi, 3),
           "P_one_sided": round(p_one, 4), "P_two_sided": round(p_two, 4),
           "n_boot_ok": ok}
    log(f"{label}: n={len(d)} dC={c1_obs-c0_obs:+.3f} ({lo:+.3f},{hi:+.3f}) "
        f"P1={p_one:.4f} P2={p_two:.4f}")
    return res, deltas


# ---- cohort 1: GSE189672 (replicates paired_deltaC.py data prep) ----
lk = pd.read_csv(os.path.join(OUT, "gse189521_linked_endpoints.csv"))
cl0 = pd.read_csv(os.path.join(OUT, "clinical_parsed.csv"), dtype={"Sample": str})
cl0["grade2"] = (cl0["grade"] == 2).astype(int)
cl0_m = cl0.merge(lk[["Sample", "pmc"]], on="Sample", how="left")
known_mg = cl0_m["menG"].isin(["A", "B", "C"])
df1 = pd.DataFrame({
    "T": cl0_m["fu_months"], "E": cl0_m["rec"],
    "grade2": cl0_m["grade2"], "gtr": cl0_m["gtr"],
    "age": (cl0_m["age"] - cl0_m["age"].mean()) / cl0_m["age"].std(),
    "male": (cl0_m["Gender"].astype(str).str.upper() == "M").astype(int),
    "menG_B": np.where(known_mg, (cl0_m["menG"] == "B").astype(float), np.nan),
    "menG_C": np.where(known_mg, (cl0_m["menG"] == "C").astype(float), np.nan),
    "pmc_2": (pd.to_numeric(cl0_m["pmc"], errors="coerce") == 2).astype(float),
    "pmc_3": (pd.to_numeric(cl0_m["pmc"], errors="coerce") == 3).astype(float),
})
clin1 = ["grade2", "gtr", "age", "male"]

# ---- cohort 2: GSE304094-primary ----
s2 = pd.read_csv(os.path.join(OUT, "gse304094_primary_survival.csv"))
clin2 = ["grade", "gtr", "age", "male"]

# ---- cohort 3: UCSF205 ----
s3 = pd.read_csv(os.path.join(OUT, "natgen205_survival.csv"))
clin3 = ["grade", "gtr", "male"]

res, dists = [], {}
r, dv = paired_deltaC(df1, clin1, ["menG_B", "menG_C"], "GSE189672: +published MenG (primary)")
res.append(r); dists["gse189672_meng"] = dv
r, dv = paired_deltaC(df1, clin1, ["pmc_2", "pmc_3"], "GSE189672: +transferred class (secondary)")
res.append(r); dists["gse189672_transferred"] = dv
r, dv = paired_deltaC(s2, clin2, ["pmc_2", "pmc_3"], "GSE304094: +transferred class (secondary)")
res.append(r); dists["gse304094_transferred"] = dv
r, dv = paired_deltaC(s2, clin2, ["tor_prolif", "tor_hypermet"], "GSE304094: +published Toronto (primary)")
res.append(r); dists["gse304094_toronto"] = dv
r, dv = paired_deltaC(s3, clin3, ["hm", "ie"], "UCSF205: +published Group (primary)")
res.append(r); dists["ucsf205_group"] = dv

pd.DataFrame(res).to_csv(os.path.join(OUT, "primary_deltaC_5000.csv"), index=False)
np.savez_compressed(os.path.join(OUT, "primary_deltaC_distributions.npz"), **dists)
log("saved primary_deltaC_5000.csv + primary_deltaC_distributions.npz")
print("PRIMARY_DELTA_C_5000_DONE")
