# -*- coding: utf-8 -*-
"""UCSF/Toronto per-sample identity chain (review P0: auditable mapping).

Builds, for every fingerprint column, the full published-data identity chain:
  M-ID (fingerprint col) -> GSM (series-matrix title 'UCSF Meningioma N') ->
  raw barcode (GEO filelist) -> methylation-profile column in GSE183647 betas ->
  MOESM20 row (positional bridge, disclosed assumption) -> MOESM19 record
  (exact 8-field clinical match) -> published Group label.
Flags: in-201 complete-covariate models; in-200 transfer analysis; the one
complete-covariate patient without a fingerprint column is listed separately.

Output: results/ucsftoronto_sample_mapping.csv (overwrites with the full chain)
"""
import os, gzip, re
import pandas as pd
import numpy as np

DATA = r"E:\脑膜瘤npj\data"
OUT = r"E:\脑膜瘤npj\results"
NG = os.path.join(DATA, "suppl", "natgen2022")

# 1) GSM -> title (series matrix)
acc_l, tit_l = [], []
with gzip.open(os.path.join(DATA, "GSE183647_series_matrix.txt.gz"), "rt",
               encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith("!Sample_geo_accession"):
            acc_l = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_title"):
            tit_l = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
gsm2mid = {}
for g, t in zip(acc_l, tit_l):
    m = re.match(r"UCSF Meningioma (\d+) \(methylation\)", t)
    if m:
        gsm2mid[g] = "M" + m.group(1)
mid2gsm = {v: k for k, v in gsm2mid.items()}

# 2) GSM -> raw barcode (filelist)
fl = pd.read_csv(os.path.join(DATA, "suppl", "GSE183647_filelist.txt"), sep="\t",
                 comment="#", header=None, names=["kind", "name", "time", "size", "type"])
fl = fl[fl["kind"] == "File"]
gsm2bar = {}
for nm in fl["name"]:
    parts = str(nm).split("_")
    if parts[0].startswith("GSM") and len(parts) >= 3:
        gsm2bar.setdefault(parts[0], parts[1] + "_" + parts[2])

# 3) fingerprint columns (MOESM21)
m21 = pd.read_excel(os.path.join(NG, "41588_2022_1061_MOESM21_ESM.xlsx"), index_col=0)
m21.columns = [str(c).strip() for c in m21.columns]
mids = list(m21.columns)

# 4) clinical chain: M20 row j <-> fingerprint col j (positional, disclosed) -> M19 record
m19 = pd.read_excel(os.path.join(NG, "41588_2022_1061_MOESM19_ESM.xlsx"))
m20 = pd.read_excel(os.path.join(NG, "41588_2022_1061_MOESM20_ESM.xlsx"))
shared = [c for c in m20.columns if c in m19.columns]
k19 = m19[shared].astype(str).agg("|".join, axis=1)
idx19 = {k: i for i, k in enumerate(k19)}

# 5) methylation beta columns (confirm profile availability per barcode)
bet = pd.read_csv(os.path.join(DATA, "GSE183647_meningioma_betas.csv.gz"), index_col=0, nrows=5)
beta_cols = set(str(c).strip() for c in bet.columns)

rows = []
for j, mid in enumerate(mids):
    gsm = mid2gsm.get(mid)
    bar = gsm2bar.get(gsm)
    kk = "|".join(str(x) for x in m20.iloc[j][shared])
    r19 = idx19.get(kk)
    rows.append({
        "MID": mid,
        "fingerprint_col_index": j,
        "GSM": gsm,
        "raw_barcode": bar,
        "barcode_in_beta_matrix": bar in beta_cols if bar is not None else False,
        "M20_row_index": j,
        "M19_record_index": r19 if r19 is not None else "NO_MATCH",
        "published_Group": m19["Group"].values[r19] if r19 is not None else "",
        "in_201_complete_covariate_models": True,
        "in_200_transfer_analysis": True,
        "linkage": "MID<->GSM<->barcode: ID-based (series title + filelist); M20row<->fingerprint col: positional (disclosed); M20row->M19 record: exact 8-field match",
    })
mp = pd.DataFrame(rows)
# the 201st complete-covariate patient (M20 row 200) has no fingerprint column
kk200 = "|".join(str(x) for x in m20.iloc[200][shared])
r200 = idx19.get(kk200)
extra = pd.DataFrame([{
    "MID": None, "fingerprint_col_index": None, "GSM": None, "raw_barcode": None,
    "barcode_in_beta_matrix": False, "M20_row_index": 200,
    "M19_record_index": r200, "published_Group": m19["Group"].values[r200],
    "in_201_complete_covariate_models": True, "in_200_transfer_analysis": False,
    "linkage": "complete-covariate patient without fingerprint column (MOESM21 provides 200 columns)",
}])
mp = pd.concat([mp, extra], ignore_index=True)
mp.to_csv(os.path.join(OUT, "ucsftoronto_sample_mapping.csv"), index=False)
print("rows:", len(mp))
print("barcode_in_beta_matrix:", mp.loc[mp.MID.notna(), "barcode_in_beta_matrix"].value_counts().to_dict())
print("Group dist (200):", mp.loc[mp.MID.notna(), "published_Group"].value_counts().to_dict())
print("extra row (M20 row 200):", extra.to_dict("records"))
print("UCSFTORONTO_IDENTITY_CHAIN_DONE")
