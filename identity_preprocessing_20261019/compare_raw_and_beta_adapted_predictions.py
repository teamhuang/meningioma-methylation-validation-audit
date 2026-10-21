from pathlib import Path
import pandas as pd,numpy as np
from sklearn.metrics import adjusted_rand_score
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
raw=pd.concat([pd.read_csv(f) for f in (OUT/'raw_preprocessed').glob('*_summary.csv')],ignore_index=True)
rows=[]
for cohort,n in [('baylor',110),('ucsf',200)]:
    r=raw[raw.cohort==cohort]
    if len(r)!=n:print('Deferred',cohort,len(r),n);continue
    old=pd.read_csv(ROOT.parent/('strengthening_20261019/results/official_knn_predictions_'+cohort+'.csv'))
    key='GSM' if cohort=='baylor' else 'barcode';old=old.rename(columns={'id':key})
    d=r.merge(old,on=key,validate='one_to_one',suffixes=('_raw','_beta_adapted'))
    assert len(d)==n
    d.to_csv(OUT/(cohort+'_paired_raw_vs_beta_adapted_predictions.csv'),index=False)
    mask=d.predicted_class_raw.notna()&d.predicted_class_beta_adapted.notna()
    rows.append({'cohort':cohort,'n_profiles':len(d),'n_both_assigned':int(mask.sum()),'n_changed_assigned_classes':int(d.loc[mask,'predicted_class_raw'].ne(d.loc[mask,'predicted_class_beta_adapted']).sum()),
       'ARI_raw_vs_processed_among_both_assigned':adjusted_rand_score(d.loc[mask,'predicted_class_raw'],d.loc[mask,'predicted_class_beta_adapted']),
       'raw_ties':int(d.tied_top_raw.sum()),'beta_adapted_ties':int(d.tied_top_beta_adapted.sum()),'raw_max_imputed_probes':0,
       'beta_adapted_max_imputed_probes':int(d.imputed_features.max()),'inference':'paired descriptive preprocessing comparison; not evidence of improved clinical performance'})
if rows:pd.DataFrame(rows).to_csv(OUT/'paired_raw_vs_beta_adapted_summary.csv',index=False);print(pd.DataFrame(rows).to_string(index=False))
