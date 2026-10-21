from pathlib import Path
import pandas as pd,json
from scipy.optimize import linear_sum_assignment
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
ref=pd.read_csv(ROOT.parent/'strengthening_20261019/results/official_model/official_knn_training.csv')
meta=pd.read_csv(OUT/'raw_idat_manifest.csv');meta=meta[meta.cohort.isin(['author186','author37'])]
d=ref[['sample','group']].merge(meta,on=None,left_on='sample',right_on='barcode',validate='one_to_one')
assert len(d)==141
ct=pd.crosstab(d.group,d['choudhury class']);a,b=linear_sum_assignment(-ct.to_numpy())
mapping={str(ct.index[i]):str(ct.columns[j]) for i,j in zip(a,b)}
ct.to_csv(OUT/'training_only_semantic_contingency.csv');d.to_csv(OUT/'training_reference_ID_proven_labels.csv',index=False)
(OUT/'training_only_semantic_map.json').write_text(json.dumps({'n_reference':len(d),'mapping':mapping,'source':'141 author reference barcodes joined literally to GEO training metadata; target UCSF labels never used to select mapping','agreement':int(ct.to_numpy()[a,b].sum())/len(d)},indent=2),encoding='utf-8')
print('TRAINING SEMANTICS',mapping,'contingency\n',ct.to_string())
