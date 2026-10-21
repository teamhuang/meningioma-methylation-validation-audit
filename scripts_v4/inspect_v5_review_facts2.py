import os, glob, pandas as pd, numpy as np
from sklearn.metrics import cohen_kappa_score
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def P(*a): return os.path.join(ROOT, *a)

# ---- grade/sex/label positional mapping: kappa + chance agreement ----
fp = P("identity_preprocessing_20261019/results/ucsf_old_positional_mapping_rejected.csv")
d = pd.read_csv(fp)
print("[positional] cols:", list(d.columns), "n=", len(d))
for c in d.columns:
    print("  ", c, "->", d[c].head(4).tolist())

def chance_and_kappa(a, b, name):
    a=np.asarray(a); b=np.asarray(b)
    m=~(pd.isna(a)|pd.isna(b))
    a,b=a[m],b[m]
    obs=np.mean(a==b)
    cats=np.unique(np.concatenate([a,b]))
    pa=np.array([(a==k).mean() for k in cats]); pb=np.array([(b==k).mean() for k in cats])
    exp=np.sum(pa*pb)
    kappa=(obs-exp)/(1-exp) if exp<1 else np.nan
    print(f"  {name}: n={len(a)} observed={obs:.3f} expected_chance={exp:.3f} kappa={kappa:.3f} 'worse-than-chance'={obs<exp}")

# find grade-like and group-like paired columns (pred vs published)
for c in d.columns:
    lc=c.lower()
    if lc in ('grade_pred','grade_published','grade','true_grade','positional_grade','recorded_grade'):
        pass
print("---- attempt kappa on plausible label/grade/sex pairs ----")
# print all unique pairs by looking for _pred/_published or agree cols
print(d.head(3).T.to_string())

# ---- Table 3 denominators: HK hypermitotic total vs uniquely assigned ----
print("\n[HKG confusion / linked predictions]")
for name in ["hku_raw_author_confusion_all_assigned.csv","hku_raw_author_ID_linked_predictions.csv"]:
    fp=P("identity_preprocessing_20261019/results",name)
    if os.path.exists(fp):
        x=pd.read_csv(fp)
        print("\n",name,"cols=",list(x.columns),"n=",len(x))
        print(x.head(6).to_string(index=False))
