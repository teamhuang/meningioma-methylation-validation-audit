from pathlib import Path
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
ref=pd.read_csv(ROOT.parent/'strengthening_20261019/results/official_model/official_knn_training.csv').set_index('sample')
meta=pd.read_csv(OUT/'requested_raw_author186.csv');rows=[]
for _,r in meta.iterrows():
    path=OUT/'raw_preprocessed'/(r.GSM+'_403beta.csv')
    if not path.exists():raise ValueError('anchor not processed '+r.GSM)
    b=pd.read_csv(path).set_index('cpg').beta;v=ref.loc[r.barcode,b.index].astype(float)
    diff=np.abs(b-v)
    rows.append({'GSM':r.GSM,'barcode':r.barcode,'reference_class':r.group,'n_cpgs':len(b),'max_abs_beta_difference':float(diff.max()),
                'mean_abs_beta_difference':float(diff.mean()),'pearson_correlation':float(np.corrcoef(b,v)[0,1]),'pass_tolerance_1e12':bool(diff.max()<=1e-12)})
d=pd.DataFrame(rows);d.to_csv(OUT/'reference_anchor_raw_vs_original_403beta.csv',index=False);print(d.to_string(index=False));assert d.pass_tolerance_1e12.all()
