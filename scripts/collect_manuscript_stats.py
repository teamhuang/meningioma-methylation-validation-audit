# -*- coding: utf-8 -*-
"""Collect precise manuscript statistics into results/manuscript_stats.json."""
import os, json
import pandas as pd
import numpy as np

OUT = r"E:\脑膜瘤npj\results"
stats = {}

# a) UCSF-SeSAMe transfer diagnosis: crosstab frozen-classifier pred (1/2/3) vs official kmeans (0/1/2)
d = pd.read_csv(os.path.join(OUT, "natgen201_transfer_diagnosis.csv"))
ct = pd.crosstab(d["pred"], d["km"])
stats["sesame_n"] = int(len(d))
stats["sesame_crosstab_pred_x_km"] = {str(k): {str(c): int(v) for c, v in row.items()} for k, row in ct.to_dict("index").items()}
stats["sesame_km_dist"] = d["km"].value_counts().sort_index().to_dict()
stats["sesame_pred_dist"] = d["pred"].value_counts().sort_index().to_dict()
# best one-to-one mapping accuracy (greedy by row)
ctn = ct.values
# try mapping pred->km for accuracy via hungarian (scipy)
from scipy.optimize import linear_sum_assignment
r, c = linear_sum_assignment(-ctn)
acc_map = ctn[r, c].sum() / ctn.sum()
stats["sesame_bestmap_acc"] = round(float(acc_map), 3)
stats["sesame_bestmap"] = {str(int(pr)): str(int(km)) for pr, km in zip(r, c)}

# b) GSE292326 identity audit agreement rates
a = pd.read_csv(os.path.join(OUT, "gse292326_identity_audit.csv"))
a = a[~a["kmeans"].astype(str).isin(["NA", "nan", ""])]
same_kmeans = (a["our_pred"].astype(int).astype(str) == a["kmeans"].astype(int).astype(str)).mean()
same_chou = (a["pred_as_choudhury"].str.lower() == a["choudhury"].str.lower()).mean()
stats["gse292326_n"] = int(len(a))
stats["gse292326_ourpred_vs_kmeans_acc"] = round(float(same_kmeans), 3)
stats["gse292326_predvschoudhury_acc"] = round(float(same_chou), 3)
stats["gse292326_choudhury_dist"] = a["choudhury"].value_counts().to_dict()

# c) longitudinal stability
s = pd.read_csv(os.path.join(OUT, "gse304096_stability.csv"))
stats["long_n_patients"] = int(len(s))
stats["long_n_samples"] = int(s["n_samples"].sum())
stats["long_stable_computed"] = int(s["all_same"].sum())
stats["long_unstable_patients"] = s.loc[s["all_same"] == 0, "patient_cluster"].tolist()

# d) paired CpG shifts significance
sh = pd.read_csv(os.path.join(OUT, "gse304096_paired_cpg_shifts.csv"))
stats["cpg_shifts_n"] = int(len(sh))
stats["cpg_shifts_min_padj"] = float(sh["p_adj"].min())
stats["cpg_shifts_sig_005"] = int((sh["p_adj"] < 0.05).sum())

# e) GSE189521 linked endpoints cohort description
lk = pd.read_csv(os.path.join(OUT, "gse189521_linked_endpoints.csv"))
stats["gse189521_n"] = int(len(lk))
stats["gse189521_menG_dist"] = lk["menG"].value_counts(dropna=False).to_dict()
stats["gse189521_recurrence_rate"] = round(float(lk["rec"].mean()), 3)
stats["gse189521_median_fu_mo"] = round(float(lk["fu_months"].median()), 1)
stats["gse189521_grade_dist"] = lk["grade"].value_counts(dropna=False).to_dict()

# f) S6 QC threshold results
s6 = pd.read_csv(os.path.join(OUT, "s6_qc_threshold.csv"))
stats["s6"] = s6.to_dict("records")

# g) frozen reject threshold
stats["train_thr_p99"] = round(float(np.load(os.path.join(OUT, "methyl_train_thr.npy"))[0]), 2)

# h) sensitivity rows
sens = pd.read_csv(os.path.join(OUT, "sensitivity_analysis.csv"))
stats["sens_n_checks"] = int(len(sens))
stats["sens_sig"] = sens[sens["P_delta<=0"] < 0.05][["check", "cohort", "n", "deltaC", "deltaC_CI_low", "deltaC_CI_high", "P_delta<=0"]].to_dict("records")

with open(os.path.join(OUT, "manuscript_stats.json"), "w", encoding="utf-8") as f:
    json.dump(stats, f, ensure_ascii=False, indent=1, default=str)
print(json.dumps(stats, ensure_ascii=False, indent=1, default=str))
