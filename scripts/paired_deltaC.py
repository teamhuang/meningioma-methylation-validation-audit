import os, warnings, time
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

t0 = time.time()
OUT = r"E:\脑膜瘤npj\results"
def log(m): print(f"[{time.time()-t0:6.0f}s] {m}", flush=True)
rng = np.random.default_rng(20260317)

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
    raise RuntimeError("fit failed " + str(cols) + " :: " + str(last)[:60])

def paired_deltaC(df, clin, mol, n_boot=1000, label=""):
    d = df[["T", "E"] + clin + mol].dropna().reset_index(drop=True)
    c0, _ = fit(d, clin)
    c1, _ = fit(d, clin + mol)
    c0_obs = concordance_index(d["T"], -c0.predict_partial_hazard(d), d["E"])
    c1_obs = concordance_index(d["T"], -c1.predict_partial_hazard(d), d["E"])
    deltas = []
    for _ in range(n_boot):
        db = d.iloc[rng.choice(len(d), len(d), replace=True)].reset_index(drop=True)
        try:
            m0 = CoxPHFitter(penalizer=0.05); m0.fit(db[["T", "E"] + clin], "T", "E")
            m1 = CoxPHFitter(penalizer=0.05); m1.fit(db[["T", "E"] + clin + mol], "T", "E")
            a = concordance_index(db["T"], -m0.predict_partial_hazard(db), db["E"])
            b = concordance_index(db["T"], -m1.predict_partial_hazard(db), db["E"])
            deltas.append(b - a)
        except Exception:
            continue
    deltas = np.array(deltas)
    lo, hi = np.percentile(deltas, [2.5, 97.5])
    p_boot = float((deltas <= 0).mean())
    res = {"analysis": label, "n": len(d), "events": int(d["E"].sum()),
           "C_clinical": round(c0_obs, 3), "C_plus_molecular": round(c1_obs, 3),
           "deltaC": round(c1_obs - c0_obs, 3),
           "deltaC_CI_low": round(lo, 3), "deltaC_CI_high": round(hi, 3),
           "P_delta<=0": round(p_boot, 4), "n_boot_ok": len(deltas)}
    log(f"{label}: n={len(d)} ev={int(d['E'].sum())} C0={c0_obs:.3f} C1={c1_obs:.3f} "
        f"dC={c1_obs-c0_obs:+.3f} ({lo:+.3f}, {hi:+.3f}) p_boot={p_boot:.3f}")
    return res

rows = []

# ---------- cohort 1: GSE189672 ----------
s1 = pd.read_csv(os.path.join(OUT, "survival_input.csv"))
lk = pd.read_csv(os.path.join(OUT, "gse189521_linked_endpoints.csv"))
lk_map = lk.set_index("Sample")["pmc"]
clin1 = ["grade2", "gtr", "age", "male"]
d1 = s1[["T", "E"] + clin1 + ["menG_B", "menG_C"]].dropna().reset_index(drop=True)
# transfer class for the same 99 patients: need Sample->pmc; s1 lost Sample id; use lk order: replicate via survival_input? rebuild from clinical_parsed
cl0 = pd.read_csv(os.path.join(OUT, "clinical_parsed.csv"), dtype={"Sample": str})
cl0["grade2"] = (cl0["grade"] == 2).astype(int)
# align with s1 rows: s1 came from feat/clin without Sample; instead recompute directly from cl0
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
rows.append(paired_deltaC(df1, clin1, ["menG_B", "menG_C"], label="GSE189672: +published MenG"))
rows.append(paired_deltaC(df1, clin1, ["pmc_2", "pmc_3"], label="GSE189672: +transferred class"))

# ---------- cohort 2: GSE304094-primary ----------
s2 = pd.read_csv(os.path.join(OUT, "gse304094_primary_survival.csv"))
clin2 = ["grade", "gtr", "age", "male"]
rows.append(paired_deltaC(s2, clin2, ["pmc_2", "pmc_3"], label="GSE304094: +transferred class"))
s2b = s2.copy()
s2b["tor_hypermet"] = s2b["tor_hypermet"].fillna(0)  # if available
if "methhigh" in s2b.columns:
    rows.append(paired_deltaC(s2b, clin2, ["methhigh"], label="GSE304094: +published METHhigh"))
rows.append(paired_deltaC(s2b, clin2, ["tor_prolif", "tor_hypermet"], label="GSE304094: +published Toronto"))

# ---------- cohort 3: UCSF205 ----------
s3 = pd.read_csv(os.path.join(OUT, "natgen205_survival.csv"))
clin3 = ["grade", "gtr", "male"]
rows.append(paired_deltaC(s3, clin3, ["hm", "ie"], label="UCSF205: +published Group"))
# MIB1 comparison
if "mib1" in s3.columns:
    rows.append(paired_deltaC(s3, clin3, ["mib1"], label="UCSF205: +MIB1"))

# ---------- grade-stratified paired dC ----------
print()
for g, sub, clin, mol, tag in [
    (1, s2, ["gtr", "age", "male"], ["pmc_2", "pmc_3"], "GSE304094 WHO1: +transferred"),
    (2, s2, ["gtr", "age", "male"], ["pmc_2", "pmc_3"], "GSE304094 WHO2: +transferred"),
    (1, s3, ["gtr", "male"], ["hm", "ie"], "UCSF205 WHO1: +published Group"),
    (2, s3, ["gtr", "male"], ["hm", "ie"], "UCSF205 WHO2: +published Group"),
]:
    ss = sub[sub["grade"] == g]
    if len(ss) >= 20 and ss["E"].sum() >= 8:
        rows.append(paired_deltaC(ss, clin, mol, label=tag))
    else:
        log(f"{tag}: skipped (n={len(ss)}, ev={int(ss['E'].sum())})")

pd.DataFrame(rows).to_csv(os.path.join(OUT, "paired_deltaC_all.csv"), index=False)
log("PAIRED_DELTAC_DONE")
