import os, warnings, time
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")

t0 = time.time()
DATA = r"E:\脑膜瘤npj\data"
OUT = r"E:\脑膜瘤npj\results"
def log(m): print(f"[{time.time()-t0:6.0f}s] {m}", flush=True)

m = pd.read_csv(os.path.join(DATA, "GSE292325_matrix_processed.txt.gz"), sep="\t", index_col=0)
m = m[~m.index.duplicated(keep="first")]
sm = pd.read_csv(os.path.join(DATA, "meta", "GSE292325_sample_map.csv"), dtype=str)
y = sm["char4"].str.replace("kmeans cluster:", "", regex=False).str.strip().values
log(f"matrix {m.shape}; labels {pd.Series(y).value_counts().to_dict()}")

def centroid_predict(tr_idx, te_idx, top_k=20000):
    X_tr_raw = m.iloc[:, tr_idx]
    mad = (X_tr_raw - X_tr_raw.median(axis=1).values[:, None]).abs().median(axis=1)
    na_frac = X_tr_raw.isna().mean(axis=1)
    sel = mad[(na_frac < 0.2)].sort_values(ascending=False).index[:top_k]
    X_tr = X_tr_raw.loc[sel].T
    mu = X_tr.mean()
    sd = X_tr.std().replace(0, 1)
    Z_tr = (X_tr - mu) / sd
    cents = {}
    ytr = y[tr_idx]
    for lab in np.unique(ytr):
        cents[lab] = Z_tr[ytr == lab].mean(axis=0)
    labs = sorted(cents)
    M = np.vstack([cents[l] for l in labs])
    X_te = m.iloc[:, te_idx].loc[sel].T
    Z_te = (X_te - mu) / sd
    Z_te = Z_te.dropna(axis=1, how="any")
    common = [c for c in Z_te.columns if c in cents[labs[0]].index]
    M2 = np.vstack([cents[l].reindex(common).values for l in labs])
    D = ((Z_te[common].values[:, None, :] - M2[None, :, :]) ** 2).sum(axis=2)
    pred = np.array([labs[i] for i in D.argmin(axis=1)])
    return pred

from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, adjusted_rand_score, balanced_accuracy_score, f1_score

all_true, all_pred = [], []
n_rep = 10
for seed in range(n_rep):
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    for tr, te in skf.split(np.zeros(len(y)), y):
        pred = centroid_predict(tr, te)
        all_true.extend(y[te])
        all_pred.extend(pred)
all_true = np.array(all_true); all_pred = np.array(all_pred)
acc = accuracy_score(all_true, all_pred)
ari = adjusted_rand_score(all_true, all_pred)
ba = balanced_accuracy_score(all_true, all_pred)
f1m = f1_score(all_true, all_pred, average="macro")
log(f"NESTED 10x5-fold CV (feature selection + scaler + centroids per fold):")
log(f"  accuracy={acc:.3f}; ARI={ari:.3f}; balanced acc={ba:.3f}; macro-F1={f1m:.3f} (n={len(all_true)})")
per_class = pd.crosstab(pd.Series(all_true, name="true"), pd.Series(all_pred, name="pred"), normalize="index").round(3)
log("per-class recall:\n" + per_class.to_string())
pd.DataFrame({"metric": ["accuracy", "ARI", "balanced_acc", "macro_F1"],
              "value": [round(acc, 3), round(ari, 3), round(ba, 3), round(f1m, 3)],
              "n": [len(all_true)] * 4}).to_csv(os.path.join(OUT, "nested_cv_metrics.csv"), index=False)
log("NESTED_CV_DONE")
