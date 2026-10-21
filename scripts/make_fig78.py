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

OUT = r"E:\脑膜瘤npj\results"
FIG = r"E:\脑膜瘤npj\figures"

# ---------- fig7: clinical utility ----------
calA = pd.read_csv(os.path.join(OUT, "calib_GSE189672.csv"))
calB = pd.read_csv(os.path.join(OUT, "calib_GSE304094-primary.csv"))
dcaA = pd.read_csv(os.path.join(OUT, "dca_GSE189672.csv"))
dcaB = pd.read_csv(os.path.join(OUT, "dca_GSE304094-primary.csv"))
tauc = pd.read_csv(os.path.join(OUT, "time_auc_results.csv"))
tauc = tauc[tauc["metric"].str.startswith("timeAUC")].copy()

fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.7), constrained_layout=True)
ax = axes[0]
for cal, name, col, mk in [(calA, "GSE189672 (n=99)", "#2a78d6", "o"), (calB, "GSE304094-primary (n=149)", "#1baf7a", "s")]:
    ax.plot(cal["pred_5y_risk"], cal["obs_5y_event_ipcw"], marker=mk, color=col, linewidth=0,
            markersize=9, label=name)
xs = np.linspace(0, 0.9, 10)
ax.plot(xs, xs, color="#bbbbbb", linewidth=1, linestyle="--", label="Ideal")
ax.set_xlabel("Predicted 5-year event probability (tertile mean)")
ax.set_ylabel("Observed 5-year event rate (IPCW)")
ax.legend(fontsize=8, frameon=False, loc="upper left")
ax.set_title("A  Calibration by risk tertile (5-year)", loc="left", fontsize=10)

ax = axes[1]
for dca, name, col in [(dcaA, "GSE189672", "#2a78d6"), (dcaB, "GSE304094-primary", "#1baf7a")]:
    ax.plot(dca["threshold"], dca["NB_model"], color=col, linewidth=2, label=name + " — model")
    ax.plot(dca["threshold"], dca["NB_treat_all"], color=col, linewidth=1.2, linestyle=":", label=name + " — treat all")
ax.axhline(0, color="#888888", linewidth=1)
ax.set_xlabel("Threshold probability (5-year)")
ax.set_ylabel("Net benefit")
ax.set_ylim(-0.08, 0.32)
ax.legend(fontsize=7.5, frameon=False)
ax.set_title("B  Decision curve analysis (5-year)", loc="left", fontsize=10)

ax = axes[2]
lab = []
for _, r in tauc.iterrows():
    lab.append(r["cohort"].split(" (")[0] + "\n" + r["metric"].replace("timeAUC_", "").replace("mo", "") + " months")
vals = tauc["value"].values
cases = tauc["n_cases"].values
cols = ["#2a78d6" if "189672" in c else "#1baf7a" for c in tauc["cohort"]]
y = np.arange(len(tauc))
ax.barh(y, vals, color=cols, height=0.55)
for yi, v, n in zip(y, vals, cases):
    ax.text(v + 0.01, yi, f"{v:.3f} ({n} ev)", va="center", fontsize=8)
ax.axvline(0.5, color="#bbbbbb", linewidth=1, linestyle="--")
ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=8)
ax.set_xlim(0.4, 1.0); ax.set_xlabel("Time-dependent AUC (IPCW)")
ax.set_title("C  Time-dependent discrimination", loc="left", fontsize=10)
fig.savefig(os.path.join(FIG, "fig7_clinical_utility.png"), bbox_inches="tight")
fig.savefig(os.path.join(FIG, "fig7_clinical_utility.pdf"), bbox_inches="tight")
plt.close(fig)
print("fig7 done")

# ---------- fig8: meta forest ----------
meta = pd.read_csv(os.path.join(OUT, "meta_high_risk_HR.csv"))
fig, ax = plt.subplots(figsize=(8.2, 3.0), constrained_layout=True)
y = np.arange(len(meta))[::-1]
for yi, (_, r) in zip(y, meta.iterrows()):
    pooled = "Pooled" in r["cohort"]
    ax.plot([r["CI_low"], r["CI_high"]], [yi, yi], color="#184f95" if pooled else "#2a78d6",
            linewidth=2.4 if pooled else 1.8)
    ax.scatter(r["HR"], yi, s=150 if pooled else 60, color="#184f95" if pooled else "#2a78d6",
               marker="D" if pooled else "o", zorder=3, linewidths=0)
    ax.text(r["CI_high"] + 0.06, yi, f"HR={r['HR']:.2f} ({r['CI_low']:.2f}-{r['CI_high']:.2f})", va="center", fontsize=8.5)
ax.axvline(1, color="#bbbbbb", linewidth=1, linestyle="--")
ax.set_yticks(y); ax.set_yticklabels(meta["cohort"], fontsize=9)
ax.set_xlim(0.4, 11.0)
ax.set_xlabel("Progression/recurrence-free survival, high-risk class vs reference (adjusted Cox)")
ax.set_title(f"High-risk molecular class across cohorts (pooled p=2.1e-05)", loc="left", fontsize=10)
fig.savefig(os.path.join(FIG, "fig8_meta_forest.png"), bbox_inches="tight")
fig.savefig(os.path.join(FIG, "fig8_meta_forest.pdf"), bbox_inches="tight")
plt.close(fig)
print("fig8 done")


