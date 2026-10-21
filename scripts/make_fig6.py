import os, warnings
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
STYLE = Path(r"C:\Users\Admin\AppData\Roaming\com.ai4s.workbench\runtime\xdg-config\opencode\skills\publication-figures\openscience.mplstyle")
plt.style.use(str(STYLE)) if STYLE.exists() else None
from lifelines import KaplanMeierFitter

OUT = r"E:\脑膜瘤npj\results"
FIG = r"E:\脑膜瘤npj\figures"
prim = pd.read_csv(os.path.join(OUT, "gse304094_primary_survival.csv"))
res = pd.read_csv(os.path.join(OUT, "cindex_304094_primary.csv"))
strat = pd.read_csv(os.path.join(OUT, "cindex_304094_gradestrat.csv"))

fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.7), constrained_layout=True)
# A: KM by transferred class
ax = axes[0]
kmf = KaplanMeierFitter()
for k, col in [("1", "#2a78d6"), ("2", "#1baf7a"), ("3", "#eda100")]:
    sub = prim[(prim["pmc"] == int(k))].dropna(subset=["T", "E"])
    if len(sub) == 0: continue
    kmf.fit(sub["T"], sub["E"], label=f"Class {k} (n={len(sub)}, ev={int(sub['E'].sum())})")
    kmf.plot_survival_function(ax=ax, color=col, linewidth=2, ci_show=False)
ax.set_xlabel("Months from surgery"); ax.set_ylabel("Progression-free probability")
ax.legend(fontsize=8, frameon=False, loc="lower left")
ax.set_title("A  Tübingen primary cohort (n=149, 55 events)", loc="left", fontsize=10)

# B: model C comparison
ax = axes[1]
NAME = {"M0_clinical": "Clinical only", "M1_clinical+predMethylClass": "+ transferred methyl class",
        "M1b_predMethylClass_only": "Transferred class only", "M2_clinical+METHhigh": "+ published METHhigh/low",
        "M3_clinical+TorontoGroup": "+ published Toronto 3-group"}
d = res.copy()
d["label"] = d["model"].map(NAME)
def cilo(s):
    s = str(s).strip("()"); a, b = s.split("-"); return float(a), float(b)
d["lo"] = d["CI"].map(lambda s: cilo(s)[0]); d["hi"] = d["CI"].map(lambda s: cilo(s)[1])
d = d.sort_values("C_insample")
y = np.arange(len(d))
ax.hlines(y, d["lo"], d["hi"], color="#9ec5f4", linewidth=3)
ax.scatter(d["C_insample"], y, s=42, color="#2a78d6", zorder=3, linewidths=0)
ax.scatter(d["C_cv"], y, s=26, color="#e34948", marker="D", zorder=3)
ax.axvline(0.5, color="#bbbbbb", linewidth=1, linestyle="--")
ax.set_yticks(y); ax.set_yticklabels(d["label"], fontsize=9)
ax.set_xlabel("Harrell's C-index (progression-free survival)")
ax.scatter([], [], s=42, color="#2a78d6", label="In-sample C (95% bootstrap CI)")
ax.scatter([], [], s=26, color="#e34948", marker="D", label="Cross-validated C")
ax.legend(fontsize=8, frameon=False, loc="lower right")
ax.set_title("B  Model comparison (n=149)", loc="left", fontsize=10)

# C: grade-stratified
ax = axes[2]
gg = strat.pivot(index="grade", columns="model", values="C")
x = np.arange(len(gg)); w = 0.34
ax.bar(x - w / 2, gg["M0"], width=w, color="#9ec5f4", label="Clinical only")
ax.bar(x + w / 2, gg["M1"], width=w, color="#2a78d6", label="+ transferred methyl class")
for i, g in enumerate(gg.index):
    m0 = gg.loc[g, "M0"]; m1 = gg.loc[g, "M1"]
    ax.text(i - w / 2, m0 + 0.012, f"{m0:.3f}", ha="center", fontsize=8)
    ax.text(i + w / 2, m1 + 0.012, f"{m1:.3f}", ha="center", fontsize=8)
ax.set_xticks(x); ax.set_xticklabels([f"WHO grade {int(g)}" for g in gg.index])
ax.set_ylim(0.45, 0.8); ax.set_ylabel("C-index (in-sample)")
ax.legend(fontsize=8, frameon=False, loc="upper left")
ax.set_title("C  Grade-stratified incremental value", loc="left", fontsize=10)
fig.savefig(os.path.join(FIG, "fig6_external_validation.png"), bbox_inches="tight")
fig.savefig(os.path.join(FIG, "fig6_external_validation.pdf"), bbox_inches="tight")
plt.close(fig)
print("fig6 done")
