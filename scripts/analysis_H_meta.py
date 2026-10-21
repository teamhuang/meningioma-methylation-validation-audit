import os, warnings
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")
from scipy.stats import norm

OUT = r"E:\脑膜瘤npj\results"
parts = [
    ("GSE189672 (MenG C vs A)", os.path.join(OUT, "cox_meng_summary.csv"), "menG_C"),
    ("GSE304094-primary (Class 3 vs 1)", os.path.join(OUT, "cox_304094_predmethyl_summary.csv"), "pmc_3"),
    ("NatGen-UCSF205 (Hypermitotic vs Merlin-intact)", os.path.join(OUT, "cox_natgen205_summary.csv"), "hm"),
]
rows = []
for name, path, cov in parts:
    s = pd.read_csv(path, index_col=0)
    hr = float(s.loc[cov, "exp(coef)"])
    z = abs(float(s.loc[cov, "z"]))
    se = float(np.log(hr) / z)
    rows.append({"cohort": name, "logHR": np.log(hr), "se": se,
                 "HR": round(hr, 3), "CI_low": round(float(np.exp(np.log(hr) - 1.96 * se)), 3),
                 "CI_high": round(float(np.exp(np.log(hr) + 1.96 * se)), 3),
                 "p": float(2 * norm.sf(np.log(hr) / se))})
b = np.array([r["logHR"] for r in rows]); w = np.array([1 / r["se"] ** 2 for r in rows])
b_pool = float((b * w).sum() / w.sum()); se_pool = float(np.sqrt(1 / w.sum()))
Q = float((w * (b - b_pool) ** 2).sum())
I2 = max(0.0, (Q - (len(b) - 1)) / Q * 100) if Q > 0 else 0.0
rows.append({"cohort": "Pooled (fixed effect)", "logHR": b_pool, "se": se_pool,
             "HR": round(float(np.exp(b_pool)), 3),
             "CI_low": round(float(np.exp(b_pool - 1.96 * se_pool)), 3),
             "CI_high": round(float(np.exp(b_pool + 1.96 * se_pool)), 3),
             "p": float(2 * norm.sf(b_pool / se_pool))})
meta = pd.DataFrame(rows)
meta.attrs["I2"] = round(I2, 1)
meta.to_csv(os.path.join(OUT, "meta_high_risk_HR.csv"), index=False)
print(meta.drop(columns=["logHR", "se"]).to_string(index=False))
print(f"I2 = {I2:.0f}%  Q = {Q:.2f} (df={len(b)-1})")
