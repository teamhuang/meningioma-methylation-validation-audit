import os, warnings, time
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")

t0 = time.time()
NG = r"E:\脑膜瘤npj\data\suppl\natgen2022"
def log(m): print(f"[{time.time()-t0:6.0f}s] {m}", flush=True)

c19 = pd.read_excel(os.path.join(NG, "41588_2022_1061_MOESM19_ESM.xlsx"))
c20 = pd.read_excel(os.path.join(NG, "41588_2022_1061_MOESM20_ESM.xlsx"))
log(f"c19 {c19.shape}; c20 {c20.shape}")

# match rows by MIB1 + LF + LFFR + WHO + EOR fingerprint
def key(df):
    return (pd.to_numeric(df["MIB1-LI"], errors="coerce").round(3).astype(str) + "|" +
            pd.to_numeric(df["LF"], errors="coerce").astype(str) + "|" +
            pd.to_numeric(df["LFFR"], errors="coerce").round(3).astype(str) + "|" +
            pd.to_numeric(df["WHO grade"], errors="coerce").astype(str) + "|" +
            df["EOR"].astype(str).str.strip().str.upper())

k19 = key(c19)
k20 = key(c20)
# index each key in c19
from collections import defaultdict
idx19 = defaultdict(list)
for i, k in enumerate(k19):
    idx19[k].append(i)
log(f"c19 unique keys: {len(idx19)}/{len(k19)}")
matches = []
missing = 0
used = set()
for j, k in enumerate(k20):
    cand = [i for i in idx19.get(k, []) if i not in used]
    if cand:
        i = cand[0]
        used.add(i)
        matches.append((j, i))
    else:
        missing += 1
log(f"matched {len(matches)}/{len(k20)}; unmatched c20 rows: {missing}")
mi = np.array([m[1] for m in matches])
mono = np.all(np.diff(mi) > 0)
log(f"embedding monotone (order-preserving): {mono}")
if not mono:
    # allow order check anyway
    inversions = int((np.diff(mi) < 0).sum())
    log(f"inversions: {inversions}")

# Map Group from c19 to c20 rows
grp = c19["Group"].astype(str).values
map_df = pd.DataFrame({"c20_row": [m[0] for m in matches], "c19_row": mi,
                       "Group": grp[mi]})
log("Group counts among matched c20 rows: " + str(map_df["Group"].value_counts().to_dict()))
map_df.to_csv(r"E:\脑膜瘤npj\results\natgen_c19_c20_rowmap.csv", index=False)
log("ROWMATCH_DONE")
