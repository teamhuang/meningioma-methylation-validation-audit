import os, re, json
import pandas as pd
import numpy as np

DATA = r"E:\脑膜瘤npj\data"
OUT = r"E:\脑膜瘤npj\results"
os.makedirs(OUT, exist_ok=True)

# ---------- load clinical ----------
raw = pd.read_excel(os.path.join(DATA, "GSE189672_Clinical_data.xlsx"), header=None)
hdr_idx = raw.index[raw[0] == "Sample"][0]
clin = raw.iloc[hdr_idx + 1:].copy()
clin.columns = [str(x).strip() for x in raw.iloc[hdr_idx]]
clin = clin.dropna(subset=["Sample"]).reset_index(drop=True)
clin["Sample"] = clin["Sample"].astype(str).str.strip()
print("n patients in xlsx:", len(clin))
print("columns:", list(clin.columns))

# ---------- parse fields ----------
def yn(x):
    if pd.isna(x): return np.nan
    s = str(x).strip().upper()
    if s in ("Y", "YES"): return 1
    if s in ("N", "NO"): return 0
    return np.nan

clin["rec"] = clin["Recurrence"].map(yn)
clin["bi"] = clin["Brain Invasion"].map(yn)
clin["nec"] = clin["Necrosis"].map(yn)
clin["age"] = pd.to_numeric(clin["Age"], errors="coerce")
clin["mib"] = pd.to_numeric(clin["MIB"], errors="coerce")
clin["fu_months"] = pd.to_numeric(clin["MonthsFollowUp"], errors="coerce")
clin["grade"] = clin["WHO grade"].astype(str).str.extract(r"(I+)")[0].map({"I": 1, "II": 2, "III": 3})
clin["simpson"] = pd.to_numeric(clin["Simpson_grade"], errors="coerce")
clin["gtr"] = np.where(clin["simpson"] <= 3, 1, np.where(clin["simpson"] > 3, 0, np.nan))
clin["tertp"] = clin["TERT Promoter Mutation"].map(yn)
clin["menG"] = clin["MenG"].astype(str).str.extract(r"(A|B|C)")[0]
for g in ["NF2", "SMARCB1", "TRAF7", "KLF4", "AKT1", "PIK3CA", "SMO", "POLR2A", "ARID1A"]:
    clin["mut_" + g] = clin[g].map(yn)
clin["n_mut_any"] = clin[["mut_NF2", "mut_SMARCB1", "mut_TRAF7", "mut_KLF4", "mut_AKT1", "mut_PIK3CA", "mut_SMO", "mut_POLR2A", "mut_ARID1A"]].sum(axis=1)

print("\n--- Endpoint audit ---")
print("Recurrence value counts (raw):", clin["Recurrence"].value_counts(dropna=False).to_dict())
print("rec: 1 = %d, 0 = %d, NA = %d" % ((clin["rec"] == 1).sum(), (clin["rec"] == 0).sum(), clin["rec"].isna().sum()))
print("MonthsFollowUp: min=%.2f  median=%.2f  max=%.2f  NA=%d" % (
    clin["fu_months"].min(), clin["fu_months"].median(), clin["fu_months"].max(), clin["fu_months"].isna().sum()))
print("\nBy recurrence status:")
print(clin.groupby("rec")["fu_months"].describe().to_string())
print("\nGrade counts:", clin["grade"].value_counts(dropna=False).sort_index().to_dict())
print("Simpson counts:", clin["simpson"].value_counts(dropna=False).sort_index().to_dict())
print("MenG counts:", clin["menG"].value_counts(dropna=False).to_dict())
print("TERTp mutant:", clin["tertp"].sum(), " NA:", clin["tertp"].isna().sum())
print("\nEvents among WHO I:", int(((clin["grade"] == 1) & (clin["rec"] == 1)).sum()), "/", int((clin["grade"] == 1).sum()))
print("Events among WHO II:", int(((clin["grade"] == 2) & (clin["rec"] == 1)).sum()), "/", int((clin["grade"] == 2).sum()))
print("Events among WHO III:", int(((clin["grade"] == 3) & (clin["rec"] == 1)).sum()), "/", int((clin["grade"] == 3).sum()))

# Recurrence=Yes with tiny follow-up (suspicious encoding)?
sus = clin[(clin["rec"] == 1) & (clin["fu_months"] < 1)]
print("\nRecurrence=Yes with fu_months<1: n=", len(sus))
print(sus[["Sample", "grade", "fu_months", "Simpson_grade"]].to_string() if len(sus) else " none")

# zero follow-up patients
print("\nPatients with fu_months==0:", clin.loc[clin["fu_months"] == 0, ["Sample", "rec", "grade"]].to_dict("records"))

# ---------- RNA sample ID mapping ----------
counts_path = os.path.join(DATA, "GSE189672_RNAseq_raw_counts.csv.gz")
import gzip
with gzip.open(counts_path, "rt") as f:
    cols = f.readline().rstrip("\n").split(",")[1:]
print("\nn RNA samples:", len(cols))

def xlsx_to_counts(s):
    m = re.match(r"^(\d{2})-(\d{2})$", s)
    if m: return f"{m.group(1)}-01-{m.group(2)}"
    return s

clin["counts_id"] = clin["Sample"].map(xlsx_to_counts)
mapped = clin["counts_id"].isin(cols)
print("xlsx patients mapped to counts columns:", int(mapped.sum()), "/", len(clin))
print("Unmapped:", clin.loc[~mapped, "Sample"].tolist())
extra = [c for c in cols if c not in set(clin["counts_id"])]
print("counts cols not in xlsx:", extra)

clin.loc[mapped, "in_rna"] = True
clin.to_csv(os.path.join(OUT, "clinical_parsed.csv"), index=False)

# ---------- MenG x endpoint cross-tab ----------
print("\nMenG x Recurrence:")
print(pd.crosstab(clin["menG"], clin["rec"]).to_string())
print("\nMenG x Grade:")
print(pd.crosstab(clin["menG"], clin["grade"]).to_string())
print("\nAUDIT_A_DONE")
