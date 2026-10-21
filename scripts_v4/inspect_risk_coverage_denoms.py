import os, pandas as pd, numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def P(*a): return os.path.join(ROOT, *a)
files = {
 "Tübingen":"identity_preprocessing_20261019/results/tubingen_raw_author_ID_linked_predictions.csv",
 "UCSF":"identity_preprocessing_20261019/results/ucsf_raw_author_ID_linked_predictions.csv",
 "Hong Kong":"identity_preprocessing_20261019/results/hku_raw_author_ID_linked_predictions.csv",
}
for coh,rel in files.items():
    fp=P(rel)
    d=pd.read_csv(fp)
    gcol='Group'; conf='confidence'; tied='tied_top'
    if gcol not in d.columns:
        print(coh,"cols",list(d.columns)); continue
    print(f"\n===== {coh}: n_all={len(d)} =====")
    order = [c for c in d[gcol].dropna().unique()]
    for thr in [0.75,0.95]:
        print(f"  -- threshold {thr} --")
        for g in sorted(order):
            sub=d[d[gcol]==g]
            N_all=len(sub)
            # assigned = not tied
            assigned = sub[~sub[tied].astype(bool)]
            # retained at threshold = among all true-g, pass confidence>=thr and not tied and assigned
            retained = sub[(~sub[tied].astype(bool)) & (sub[conf]>=thr)]
            print(f"    {g:16s} N_all={N_all:3d} assigned={len(assigned):3d} retained@thr={len(retained):3d} "
                  f"retention_of_ASSIGNED={len(retained)/max(len(assigned),1):.3f} retention_of_ALL={len(retained)/N_all:.3f} "
                  f"lost_of_ALL={1-len(retained)/N_all:.3f}")
