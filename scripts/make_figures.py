import os, json, warnings
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
os.makedirs(FIG, exist_ok=True)
C = {"1": "#2a78d6", "2": "#1baf7a", "3": "#eda100"}
LAB3 = {"1": "Class 1 (Merlin-intact-like / A)", "2": "Class 2 (Immune-enriched-like / B)", "3": "Class 3 (Hypermitotic-like / C)"}
SHORT = {"1": "Class 1", "2": "Class 2", "3": "Class 3"}

# ---------- data ----------
surv = pd.read_csv(os.path.join(OUT, "survival_input.csv"))
clin = pd.read_csv(os.path.join(OUT, "clinical_parsed.csv"))
linked = pd.read_csv(os.path.join(OUT, "gse189521_linked_endpoints.csv"))
cls304 = pd.read_csv(os.path.join(OUT, "gse304096_classification.csv"))
stab = pd.read_csv(os.path.join(OUT, "gse304096_stability.csv"))
dist994 = pd.read_csv(os.path.join(OUT, "gse270371_classification.csv"))
sm325 = pd.read_csv(os.path.join(OUT, "..", "data", "meta", "GSE292325_sample_map.csv"), dtype=str)
y325 = sm325["char4"].str.replace("kmeans cluster:", "", regex=False).str.strip()
ci_a = pd.read_csv(os.path.join(OUT, "cindex_models.csv"))
ci_b = pd.read_csv(os.path.join(OUT, "cindex_predmethyl_models.csv"))
menG_lab = np.load(os.path.join(OUT, "menG_labels.npy"), allow_pickle=True)
menG_lab = np.array(["nan" if (x is None or (isinstance(x, float) and np.isnan(x))) else str(x) for x in menG_lab], dtype=object)
linked["pmc"] = linked["pmc"].astype(str)

# ================= FIG 1: endpoint audit =================
fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.4), constrained_layout=True)
ax = axes[0]
for val, lab, col in [(0, "No recurrence (censored)", "#2a78d6"), (1, "Recurrence (event)", "#1baf7a")]:
    sub = surv[surv["E"] == val]["T"].dropna()
    ax.scatter(np.full(len(sub), val) + np.random.default_rng(1).uniform(-0.13, 0.13, len(sub)),
               sub, s=14, alpha=0.55, color=col, linewidths=0, label=f"{lab} (n={len(sub)})")
ax.set_xticks([0, 1]); ax.set_xticklabels(["Censored", "Recurrence"])
ax.set_ylabel("Follow-up / event time (months)")
ax.legend(loc="upper left", fontsize=8, frameon=False)
ax.set_title("A  Endpoint audit: time by status", loc="left", fontsize=10)

ax = axes[1]
ctg = pd.crosstab(clin["menG"], clin["grade"]).reindex(["A", "B", "C"]).fillna(0)
bottom = np.zeros(3)
for grade, col in [(1, "#2a78d6"), (2, "#eda100")]:
    vals = ctg[grade].values if grade in ctg else np.zeros(3)
    ax.bar(ctg.index, vals, bottom=bottom, color=col, label=f"WHO grade {'I' if grade==1 else 'II'}", width=0.62)
    bottom += vals
for i, tot in enumerate(bottom):
    ev = clin[(clin["menG"] == ctg.index[i]) & (clin["rec"] == 1)].shape[0]
    ax.text(i, tot + 1.2, f"{int(ev)}/{int(tot)} recurred", ha="center", fontsize=8)
ax.set_ylim(0, 68); ax.set_xlabel("Published MenG subtype"); ax.set_ylabel("Patients")
ax.legend(fontsize=8, frameon=False)
ax.set_title("B  MenG subtype vs grade (events/total)", loc="left", fontsize=10)

ax = axes[2]
kmf = KaplanMeierFitter()
kmf.fit(surv["T"], surv["E"], label="All patients")
kmf.plot_survival_function(ax=ax, color="#2a78d6", linewidth=2, ci_show=False)
ax.set_xlabel("Months from surgery"); ax.set_ylabel("Recurrence-free survival")
ax.legend(fontsize=8, frameon=False, loc="lower left")
ax.set_title("C  Overall RFS, GSE189672 (n=110, 19 events)", loc="left", fontsize=10)
fig.savefig(os.path.join(FIG, "fig1_endpoint_audit.png"), bbox_inches="tight")
fig.savefig(os.path.join(FIG, "fig1_endpoint_audit.pdf"), bbox_inches="tight")
plt.close(fig)
print("fig1 done")

# ================= FIG 2: KM by subtype =================
fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6), constrained_layout=True)
ax = axes[0]
for g, col in [("A", "#2a78d6"), ("B", "#1baf7a"), ("C", "#eda100")]:
    sub = surv[menG_lab == g]
    if len(sub) == 0: continue
    kmf.fit(sub["T"], sub["E"], label=f"MenG {g} (n={len(sub)}, ev={int(sub['E'].sum())})")
    kmf.plot_survival_function(ax=ax, color=col, linewidth=2, ci_show=False)
ax.set_xlabel("Months from surgery"); ax.set_ylabel("Recurrence-free probability")
ax.legend(fontsize=8, frameon=False, loc="lower left")
ax.set_title("A  Published MenG classification", loc="left", fontsize=10)

ax = axes[1]
for k, col in [("1", "#2a78d6"), ("2", "#1baf7a"), ("3", "#eda100")]:
    sub = linked[linked["pmc"] == k]
    kmf.fit(sub["fu_months"], sub["rec"], label=f"{SHORT[k]} (n={len(sub)}, ev={int(sub['rec'].sum())})")
    kmf.plot_survival_function(ax=ax, color=col, linewidth=2, ci_show=False)
ax.set_xlabel("Months from surgery"); ax.set_ylabel("Recurrence-free probability")
ax.legend(fontsize=8, frameon=False, loc="lower left")
ax.set_title("B  Transferred methylation classifier (ANC-2025-trained)", loc="left", fontsize=10)
fig.savefig(os.path.join(FIG, "fig2_km_subtypes.png"), bbox_inches="tight")
fig.savefig(os.path.join(FIG, "fig2_km_subtypes.pdf"), bbox_inches="tight")
plt.close(fig)
print("fig2 done")

# ================= FIG 3: C-index comparison =================
rows = []
for _, r in ci_a.iterrows():
    rows.append(r)
for _, r in ci_b.iterrows():
    if r["model"] != "M0_clinical":
        rows.append(r)
ddf = pd.DataFrame(rows)
def lo_hi(s):
    s = str(s).strip("()")
    a, b = s.split("-")
    return float(a), float(b)
def get_lo(r):
    if pd.notna(r.get("C_boot_lo", np.nan)): return r["C_boot_lo"]
    if pd.notna(r.get("CI", np.nan)): return lo_hi(r["CI"])[0]
    return np.nan
def get_hi(r):
    if pd.notna(r.get("C_boot_hi", np.nan)): return r["C_boot_hi"]
    if pd.notna(r.get("CI", np.nan)): return lo_hi(r["CI"])[1]
    return np.nan
ddf["lo"] = ddf.apply(get_lo, axis=1)
ddf["hi"] = ddf.apply(get_hi, axis=1)
NAME = {
    "M0_clinical": "Clinical only (grade, resection, age, sex)",
    "MA_clinical+MenG": "Clinical + published MenG",
    "M1_clinical+RNA_PC": "Clinical + RNA PCs (top 5k genes)",
    "M2_clinical+ecosystem": "Clinical + cell-ecosystem scores",
    "M1b_RNA_PC_only": "RNA PCs only",
    "M2b_ecosystem_only": "Cell-ecosystem scores only",
    "M3_clinical+predMethylClass": "Clinical + transferred methyl class",
    "M3b_predMethylClass_only": "Transferred methyl class only",
}
ddf["label"] = ddf["model"].map(NAME)
ddf = ddf.sort_values("C_insample")
fig, ax = plt.subplots(figsize=(8.6, 4.2), constrained_layout=True)
ypos = np.arange(len(ddf))
ax.hlines(ypos, ddf["lo"], ddf["hi"], color="#9ec5f4", linewidth=3)
ax.scatter(ddf["C_insample"], ypos, s=42, color="#2a78d6", zorder=3, linewidths=0)
for y, (cv, n, ev) in enumerate(zip(ddf["C_cv"], ddf["n"], ddf["events"])):
    if pd.notna(cv):
        ax.scatter(cv, y, s=26, color="#e34948", zorder=3, linewidths=0, marker="D")
ax.axvline(0.5, color="#bbbbbb", linewidth=1, linestyle="--")
ax.set_yticks(ypos); ax.set_yticklabels(ddf["label"], fontsize=9)
ax.set_xlabel("Harrell's C-index for recurrence-free survival")
ax.scatter([], [], s=42, color="#2a78d6", label="In-sample C (bootstrap 95% CI)")
ax.scatter([], [], s=26, color="#e34948", marker="D", label="Cross-validated C (20x 5-fold)")
ax.legend(fontsize=8, frameon=False, loc="lower right")
ax.set_title("Model comparison, GSE189672 (n=110, 19 events; MenG models n=99, 17 events)", loc="left", fontsize=10)
fig.savefig(os.path.join(FIG, "fig3_cindex_models.png"), bbox_inches="tight")
fig.savefig(os.path.join(FIG, "fig3_cindex_models.pdf"), bbox_inches="tight")
plt.close(fig)
print("fig3 done")

# ================= FIG 4: classifier transfer =================
fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.6), constrained_layout=True)
ax = axes[0]
ctm = pd.crosstab(linked["pmc"], linked["menG"]).reindex(index=["1", "2", "3"], columns=["A", "B", "C"]).fillna(0)
im = ax.imshow(ctm.values, cmap="Blues", vmin=0)
for i in range(ctm.shape[0]):
    for j in range(ctm.shape[1]):
        ax.text(j, i, int(ctm.values[i, j]), ha="center", va="center", fontsize=10,
                color="white" if ctm.values[i, j] > ctm.values.max() * 0.6 else "#184f95")
ax.set_xticks(range(3)); ax.set_xticklabels(["MenG A", "MenG B", "MenG C"])
ax.set_yticks(range(3)); ax.set_yticklabels(["Class 1", "Class 2", "Class 3"])
ax.set_xlabel("Published subtype (Sci Adv 2022)"); ax.set_ylabel("Transferred classifier (ANC 2025-trained)")
ax.set_title("A  Cross-study agreement (ARI = 0.780, n=110)", loc="left", fontsize=10)

ax = axes[1]
dist = pd.DataFrame({
    "Dev cohort (n=186)": pd.Series(y325).value_counts().reindex(["1", "2", "3"]),
    "Nat Med 2024 cohort (n=994)": dist994["predicted_class"].value_counts().reindex(["1", "2", "3"]),
    "Endpoint cohort (n=110)": linked["pmc"].value_counts().reindex(["1", "2", "3"]),
})
dist = dist.div(dist.sum(axis=0), axis=1) * 100
x = np.arange(3); w = 0.26
for i, (cname, col) in enumerate(zip(dist.columns, ["#9ec5f4", "#3987e5", "#184f95"])):
    ax.bar(x + (i - 1) * w, dist[cname].values, width=w, color=col, label=cname)
ax.set_xticks(x); ax.set_xticklabels(["Class 1", "Class 2", "Class 3"])
ax.set_ylabel("% of cohort")
ax.legend(fontsize=8, frameon=False)
ax.set_title("B  Class proportions across cohorts", loc="left", fontsize=10)
fig.savefig(os.path.join(FIG, "fig4_transfer.png"), bbox_inches="tight")
fig.savefig(os.path.join(FIG, "fig4_transfer.pdf"), bbox_inches="tight")
plt.close(fig)
print("fig4 done")

# ================= FIG 5: longitudinal stability =================
fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.4), constrained_layout=True)
ax = axes[0]
cls_map = {"1": ("Class 1", "#2a78d6"), "2": ("Class 2", "#1baf7a"), "3": ("Class 3", "#eda100")}
order = sorted(cls304["patient_cluster"].unique())
yidx = {p: i for i, p in enumerate(order)}
rng5 = np.random.default_rng(3)
discordant = set(stab.loc[stab["all_same"] == 0, "patient_cluster"])
for _, r in cls304.iterrows():
    lab, col = cls_map[str(r["predicted_class"])]
    j = yidx[r["patient_cluster"]]
    ax.scatter(int(r["predicted_class"]) + rng5.uniform(-0.16, 0.16), j, s=46,
               color=col, linewidths=0.8, edgecolors="#184f95" if r["patient_cluster"] in discordant else "none",
               zorder=3)
ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=7)
ax.set_xticks([1, 2, 3]); ax.set_xticklabels(["Class 1", "Class 2", "Class 3"], fontsize=9)
ax.set_xlim(0.5, 3.5)
ax.set_xlabel("Predicted methylation class (GSE304096, 42 paired samples)")
ax.set_title("A  Within-patient class concordance (16/17 stable)", loc="left", fontsize=10)

ax = axes[1]
sig = pd.read_csv(os.path.join(OUT, "gse304096_paired_cpg_shifts.csv"))
ax.hist(sig["mean_delta"].dropna(), bins=60, color="#2a78d6", linewidth=0)
ax.axvline(0, color="#e34948", linewidth=1.2, linestyle="--")
ax.set_xlabel("Mean paired $\\Delta\\beta$ (primary $-$ recurrence)")
ax.set_ylabel("Number of CpGs")
ax.set_title("B  Paired methylation shift (0 CpGs BH-significant)", loc="left", fontsize=10)
fig.savefig(os.path.join(FIG, "fig5_longitudinal.png"), bbox_inches="tight")
fig.savefig(os.path.join(FIG, "fig5_longitudinal.pdf"), bbox_inches="tight")
plt.close(fig)
print("fig5 done")
print("ALL_FIGURES_DONE")

