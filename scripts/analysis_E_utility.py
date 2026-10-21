import os, gzip, warnings, time
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index
from scipy.stats import chi2

t0 = time.time()
OUT = r"E:\脑膜瘤npj\results"
def log(m): print(f"[{time.time()-t0:6.0f}s] {m}", flush=True)
rng = np.random.default_rng(20260317)

def km_censor(times, events):
    """Kaplan-Meier estimate of censoring survival G(t) (for IPCW)."""
    df = pd.DataFrame({"t": times, "e": 1 - np.asarray(events, dtype=int)}).sort_values("t")
    uniq = np.unique(df["t"])
    n = len(df); surv = 1.0; out = {}
    for u in uniq:
        at_risk = (df["t"] >= u).sum()
        d = ((df["t"] == u) & (df["e"] == 1)).sum()
        surv *= (1 - d / at_risk) if at_risk > 0 else 1
        out[u] = surv
    def G(t):
        ks = [u for u in out if u <= t]
        return out[max(ks)] if ks else 1.0
    return G

def time_auc(T, E, risk, horizon):
    """Cumulative/dynamic AUC with IPCW (Heagerty-style, simple)."""
    T = np.asarray(T, float); E = np.asarray(E, int); risk = np.asarray(risk, float)
    G = km_censor(T, E)
    cases = (T <= horizon) & (E == 1)
    controls = T > horizon
    if cases.sum() < 3 or controls.sum() < 3:
        return np.nan, int(cases.sum()), int(controls.sum())
    w = []
    for i in np.where(cases)[0]:
        gi = G(T[i])
        w.append(1.0 / gi if gi > 0 else 0.0)
    w = np.array(w)
    if w.sum() == 0:
        return np.nan, int(cases.sum()), int(controls.sum())
    # weighted concordance: case risk > control risk
    num = 0.0; den = 0.0
    ci = np.where(cases)[0]
    co = np.where(controls)[0]
    for a, wi in zip(ci, w):
        comp = risk[co] < risk[a]
        ties = risk[co] == risk[a]
        num += wi * (comp.sum() + 0.5 * ties.sum())
        den += wi * len(co)
    auc = num / den if den > 0 else np.nan
    return float(auc), int(cases.sum()), int(controls.sum())

def dca_net_benefit(T, E, pred_risk, thresholds, horizon):
    """DCA with IPCW: net benefit for threshold-guided policy vs treat-all/none."""
    T = np.asarray(T, float); E = np.asarray(E, int)
    pr = np.asarray(pred_risk, float)
    G = km_censor(T, E)
    N = len(T)
    out = []
    for thr in thresholds:
        # outcome "event by horizon" with IPCW for censored before horizon
        y = np.zeros(N)
        w = np.zeros(N)
        for i in range(N):
            if T[i] <= horizon and E[i] == 1:
                y[i] = 1
                gi = G(T[i]); w[i] = 1.0 / gi if gi > 0 else 0.0
            elif T[i] > horizon:
                y[i] = 0; w[i] = 1.0
            else:
                y[i] = 0; w[i] = 0.0  # censored before horizon: uninformative
        flagged = pr >= thr
        if flagged.sum() == 0:
            nb = 0.0
        else:
            tp = w[flagged & (y == 1)].sum()
            fp = w[flagged & (y == 0)].sum()
            nb = tp / N - (fp / N) * (thr / (1 - thr))
        # treat all
        ev_all = w[(y == 1)].sum()
        nb_all = ev_all / N - ((N - ev_all) / N) * (thr / (1 - thr))
        out.append({"threshold": thr, "NB_model": nb, "NB_treat_all": nb_all, "NB_treat_none": 0.0})
    return pd.DataFrame(out)

def calib_by_tertile(T, E, pred_risk, horizon):
    T = np.asarray(T, float); E = np.asarray(E, int)
    pr = np.asarray(pred_risk, float)
    G = km_censor(T, E)
    q = pd.qcut(pr, 3, labels=["low", "mid", "high"], duplicates="drop")
    rows = []
    for g in ["low", "mid", "high"]:
        m = (q == g)
        if m.sum() == 0: continue
        kmf = KaplanMeierFitter(); kmf.fit(T[m], E[m])
        obs = float(1 - kmf.predict(horizon))
        # IPCW-adjusted observed event rate by horizon
        y = ((T <= horizon) & (E == 1)).astype(float)
        w = np.array([1.0 / G(T[i]) if (T[i] <= horizon and E[i] == 1) else (1.0 if T[i] > horizon else 0.0) for i in np.where(m)[0]])
        wsum = w.sum()
        obs_ipcw = float((y[np.where(m)[0]] * w).sum() / wsum) if wsum > 0 else np.nan
        rows.append({"tertile": g, "n": int(m.sum()), "pred_5y_risk": float(pr[m].mean()),
                     "obs_5y_event_ipcw": obs_ipcw})
    return pd.DataFrame(rows)

def risk_from_cph(cph, df):
    lp = cph.predict_log_partial_hazard(df)
    return (lp - lp.mean()) / lp.std()

# ================= cohort 1: GSE189672 =================
log("cohort 1: GSE189672")
s1 = pd.read_csv(os.path.join(OUT, "survival_input.csv"))
c1, d1 = None, None
cols1 = ["T", "E", "grade2", "gtr", "age", "male", "menG_B", "menG_C"]
d1 = s1[cols1].dropna().reset_index(drop=True)
c1 = CoxPHFitter(penalizer=0.05); c1.fit(d1, "T", "E")
# full-cohort model with predicted methyl class instead of MenG labels (to use all 109)
s1b = s1[["T", "E", "grade2", "gtr", "age", "male", "menG_B", "menG_C"]].copy()
# use clinical+MenG for cohort1 utility (as in report)
lp1 = c1.predict_log_partial_hazard(d1)
risk1 = ((lp1 - lp1.mean()) / lp1.std()).values
# 5-year predicted risk from baseline survival
bdf1 = d1.copy()
surv5_1 = c1.predict_survival_function(bdf1, times=[60]).iloc[0].values
prisk5_1 = 1 - surv5_1

# ================= cohort 2: GSE304094 primary =================
log("cohort 2: GSE304094 primary")
s2 = pd.read_csv(os.path.join(OUT, "gse304094_primary_survival.csv"))
cols2 = ["T", "E", "grade", "gtr", "age", "male", "pmc_2", "pmc_3"]
d2 = s2[cols2].dropna().reset_index(drop=True)
c2 = CoxPHFitter(penalizer=0.05); c2.fit(d2, "T", "E")
risk2 = risk_from_cph(c2, d2)
surv5_2 = c2.predict_survival_function(d2, times=[60]).iloc[0].values
prisk5_2 = 1 - surv5_2

res_all = []
for cname, d, risk, prisk5, cph in [("GSE189672 (n=99)", d1, risk1, prisk5_1, c1),
                                     ("GSE304094-primary (n=149)", d2, risk2, prisk5_2, c2)]:
    T, E = d["T"].values, d["E"].values
    for hz in [36, 60]:
        a, nc, nctl = time_auc(T, E, risk, hz)
        a0, _, _ = time_auc(T, E, np.random.default_rng(1).uniform(size=len(T)), hz)
        res_all.append({"cohort": cname, "metric": f"timeAUC_{hz}mo", "value": round(a, 3),
                        "n_cases": nc, "n_controls": nctl, "note": f"random-ref AUC={a0:.2f}"})
        log(f"{cname} timeAUC {hz}mo = {a:.3f} (cases={nc})")
    hz = 60
    cal = calib_by_tertile(T, E, prisk5, hz)
    cal.insert(0, "cohort", cname); cal["metric"] = f"calibration_{hz}mo"
    res_all.append(cal.assign(value=np.nan, note=""))
    log(f"{cname} calibration:\n{cal.to_string(index=False)}")
    thresholds = np.arange(0.05, 0.55, 0.05)
    dca = dca_net_benefit(T, E, prisk5, thresholds, hz)
    dca.insert(0, "cohort", cname); dca["metric"] = f"DCA_{hz}mo"; dca["value"] = np.nan; dca["note"] = ""
    res_all.append(dca)
    dca.to_csv(os.path.join(OUT, f"dca_{cname.split()[0]}.csv"), index=False)
    cal.to_csv(os.path.join(OUT, f"calib_{cname.split()[0]}.csv"), index=False)

pd.DataFrame([r for r in res_all if isinstance(r, dict)]).to_csv(
    os.path.join(OUT, "time_auc_results.csv"), index=False)

# ================= meta-pooled HR for high-risk class =================
log("meta-pool high-risk class HR")
def fit_hr(df, cols, ref_col_pairs):
    d = df[["T", "E"] + cols].dropna().reset_index(drop=True)
    c = CoxPHFitter(penalizer=0.05); c.fit(d, "T", "E")
    s = c.summary
    return s

# cohort1: menG_C vs menG_A (from earlier cox_meng_summary): recompute
s1x = s1[["T", "E", "grade2", "gtr", "age", "male", "menG_B", "menG_C"]].dropna().reset_index(drop=True)
c1m = CoxPHFitter(penalizer=0.05); c1m.fit(s1x, "T", "E")
hr1 = float(c1m.summary.loc["menG_C", "exp(coef)"])
se1 = float((np.log(c1m.summary.loc["menG_C", "exp(coef)"])) / c1m.summary.loc["menG_C", "z"])
# cohort2: pmc_3 vs pmc_1
s2x = s2[["T", "E", "grade", "gtr", "age", "male", "pmc_2", "pmc_3"]].dropna().reset_index(drop=True)
c2m = CoxPHFitter(penalizer=0.05); c2m.fit(s2x, "T", "E")
hr2 = float(c2m.summary.loc["pmc_3", "exp(coef)"])
se2 = float(np.log(hr2) / c2m.summary.loc["pmc_3", "z"])
b1, b2 = np.log(hr1), np.log(hr2)
w1, w2 = 1 / se1 ** 2, 1 / se2 ** 2
b_pool = (b1 * w1 + b2 * w2) / (w1 + w2)
se_pool = float(np.sqrt(1 / (w1 + w2)))
hr_pool = float(np.exp(b_pool))
ci_lo = float(np.exp(b_pool - 1.96 * se_pool)); ci_hi = float(np.exp(b_pool + 1.96 * se_pool))
z_pool = b_pool / se_pool
from scipy.stats import norm
p_pool = float(2 * norm.sf(abs(z_pool)))
Q = w1 * (b1 - b_pool) ** 2 + w2 * (b2 - b_pool) ** 2
I2 = max(0.0, (Q - 1) / Q * 100) if Q > 1 else 0.0
meta = pd.DataFrame({
    "cohort": ["GSE189672 (MenG C vs A)", "GSE304094-primary (Class 3 vs 1)", "Pooled (fixed effect)"],
    "HR": [round(hr1, 3), round(hr2, 3), round(hr_pool, 3)],
    "CI_low": [round(float(np.exp(b1 - 1.96 * se1)), 3), round(float(np.exp(b2 - 1.96 * se2)), 3), round(ci_lo, 3)],
    "CI_high": [round(float(np.exp(b1 + 1.96 * se1)), 3), round(float(np.exp(b2 + 1.96 * se2)), 3), round(ci_hi, 3)],
    "p": [float(2 * norm.sf(abs(b1 / se1))), float(2 * norm.sf(abs(b2 / se2))), p_pool],
})
meta.attrs["I2"] = round(I2, 1)
meta.to_csv(os.path.join(OUT, "meta_high_risk_HR.csv"), index=False)
log(f"meta HR: {hr_pool:.2f} ({ci_lo:.2f}-{ci_hi:.2f}) p={p_pool:.2e} I2={I2:.0f}%")
log("ANALYSIS_E_DONE")
