import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")

p = r"E:\脑膜瘤npj\data\suppl\protocadherin_suppldata.xlsx"
df = pd.read_excel(p, sheet_name="Supplementary Data 4", header=None, skiprows=3)
df.columns = ["Patient", "Sex", "Setting", "Setting_detail", "WHO_2021", "Heidelberg", "Toronto_Group", "MethCluster"]
df = df.dropna(subset=["Patient"])
df["Patient"] = df["Patient"].astype(int)
print("rows:", len(df), " patients:", df["Patient"].nunique())
print("setting detail counts:", df["Setting_detail"].value_counts().to_dict())
conc = df.groupby("Patient")["MethCluster"].nunique()
print("\nMethCluster per patient: consistent =", int((conc == 1).sum()), "/", len(conc),
      "; mixed =", int((conc > 1).sum()))
sizes = df.groupby("Patient").size()
print("timepoints per patient:", sizes.value_counts().to_dict())
# who switched
for pid, sub in df.groupby("Patient"):
    if sub["MethCluster"].nunique() > 1:
        print("switch patient", pid, ":", sub[["Setting_detail", "WHO_2021", "Toronto_Group", "MethCluster"]].to_string(index=False))
# save tidy table
df.to_csv(r"E:\脑膜瘤npj\results\gse304096_official_longitudinal_table.csv", index=False)
print("\nSD4_OUT_DONE")
