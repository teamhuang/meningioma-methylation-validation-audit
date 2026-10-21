import os, gzip, re, time, warnings
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")
from sklearn.cluster import KMeans

t0 = time.time()
DATA = r"E:\脑膜瘤npj\data"
OUT = r"E:\脑膜瘤npj\results"
NG = os.path.join(DATA, "suppl", "natgen2022")
def log(m): print(f"[{time.time()-t0:6.0f}s] {m}", flush=True)

c20 = pd.read_excel(os.path.join(NG, "41588_2022_1061_MOESM20_ESM.xlsx"))
m21 = pd.read_excel(os.path.join(NG, "41588_2022_1061_MOESM21_ESM.xlsx"), index_col=0)
m21.columns = [str(c).strip() for c in m21.columns]
X21 = m21.T
X21.index = [str(c).strip() for c in X21.index]

# rebuild bridge + betas (fast enough at ~60s)
fl = pd.read_csv(os.path.join(DATA, "suppl", "GSE183647_filelist.txt"), sep="\t", comment="#", header=None,
                 names=["kind", "name", "time", "size", "type"])
fl = fl[fl["kind"] == "File"]
gsm2bar = {}
for nm in fl["name"]:
    parts = str(nm).split("_")
    if parts[0].startswith("GSM") and len(parts) >= 3:
        gsm2bar[parts[0]] = parts[1] + "_" + parts[2]
gsm2mid = {}
with gzip.open(os.path.join(DATA, "GSE183647_series_matrix.txt.gz"), "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith("!Sample_geo_accession"):
            acc = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_title"):
            tit = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
for g, t in zip(acc, tit):
    m = re.match(r"UCSF Meningioma (\d+) \(methylation\)", t)
    if m:
        gsm2mid[g] = "M" + m.group(1)
bar2mid = {gsm2bar[g]: gsm2mid[g] for g in gsm2mid if g in gsm2bar}
bet = pd.read_csv(os.path.join(DATA, "GSE183647_meningioma_betas.csv.gz"), index_col=0)
bet.columns = [bar2mid.get(str(c).strip(), None) for c in bet.columns]
bet = bet.iloc[:, [i for i, c in enumerate(bet.columns) if c is not None]]
log(f"betas: {bet.shape}")

# our classifier prediction (same as analysis_F)
scaler = pd.read_csv(os.path.join(OUT, "methyl_scaler.csv"), index_col=0)
centroids = pd.read_csv(os.path.join(OUT, "methyl_centroids.csv"), index_col=0)
labs = list(centroids.columns)
sel_full = [c for c in X21.index if c in bet.columns]
selset = [c for c in scaler.index if c in bet.index]
B = bet.loc[selset, sel_full]
Zb = ((B.T - scaler["mean"].reindex(B.index)) / scaler["sd"].reindex(B.index)).dropna(axis=1, how="any")
cc = [c for c in Zb.columns if c in centroids.index]
log(f"Zb: {Zb.shape}; shared centroid CpGs: {len(cc)}")
M2 = centroids.loc[cc, labs].values.T
D2 = ((Zb[cc].values[:, None, :] - M2[None, :, :]) ** 2).sum(axis=2)
pred = np.array([labs[i] for i in D2.argmin(axis=1)])
dist_min = D2.min(axis=1)

# k-means in the same feature space (unsupervised reference)
km = KMeans(n_clusters=3, n_init=20, random_state=20260317).fit_predict(Zb[cc].values)

mib1 = pd.to_numeric(c20["MIB1-LI"], errors="coerce").values[:len(pred)]
rec = pd.to_numeric(c20["Recurrent"], errors="coerce").values[:len(pred)]
lf = pd.to_numeric(c20["LF"], errors="coerce").values[:len(pred)]
lffr = pd.to_numeric(c20["LFFR"], errors="coerce").values[:len(pred)]

diag = pd.DataFrame({"pred": pred, "km": km, "mib1": mib1, "rec": rec, "lf": lf, "dist": dist_min})
log("MIB1 / LF-rate by OUR classifier class:\n" + diag.groupby("pred").agg(
    n=("mib1", "size"), mean_MIB1=("mib1", "mean"), LF_rate=("lf", "mean"),
    mean_dist=("dist", "mean")).round(2).to_string())
log("crosstab our-class x kmeans:\n" + pd.crosstab(diag["pred"], diag["km"]).to_string())
log("MIB1 / LF-rate by kmeans cluster:\n" + diag.groupby("km").agg(
    n=("mib1", "size"), mean_MIB1=("mib1", "mean"), LF_rate=("lf", "mean")).round(2).to_string())
# distance distributions per predicted class (batch-mismatch indicator)
log("pred class dist ranges:\n" + diag.groupby("pred")["dist"].describe().round(1).to_string())
diag.to_csv(os.path.join(OUT, "natgen201_transfer_diagnosis.csv"), index=False)
log("DIAG_DONE")
