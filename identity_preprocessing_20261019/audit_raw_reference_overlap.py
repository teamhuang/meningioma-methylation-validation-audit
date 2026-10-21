from pathlib import Path
import pandas as pd,numpy as np,json
from scipy.spatial.distance import cdist
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
ref=pd.read_csv(ROOT.parent/'strengthening_20261019/results/official_model/official_knn_training.csv').set_index('sample')
cpgs=ref.columns.drop('group').tolist();R=ref[cpgs].to_numpy(float)
meta=pd.read_csv(OUT/'raw_idat_manifest.csv').set_index('GSM');items=[];rows=[]
for f in (OUT/'raw_preprocessed').glob('*_403beta.csv'):
    gsm=f.name.split('_')[0];b=pd.read_csv(f).set_index('cpg').beta.reindex(cpgs)
    if not np.isfinite(b).all():continue
    x=b.to_numpy(float);dist=np.sqrt(np.mean((R-x)**2,axis=1));k=int(dist.argmin())
    rows.append({'GSM':gsm,'cohort':meta.loc[gsm,'cohort'],'barcode':meta.loc[gsm,'barcode'],'nearest_author_reference_barcode':ref.index[k],
                 'min_reference_RMSE':dist[k],'max_abs_difference_to_nearest_reference':np.max(abs(R[k]-x)),
                 'barcode_in_author_reference':str(meta.loc[gsm,'barcode']) in set(ref.index.astype(str)),
                 'exact_beta_reference_match_tolerance_1e12':bool(np.max(abs(R[k]-x))<=1e-12)})
d=pd.DataFrame(rows);d.to_csv(OUT/'raw_author_reference_overlap_audit.csv',index=False)
print(d.groupby('cohort').agg(n=('GSM','size'),min_RMSE=('min_reference_RMSE','min'),barcode_overlap=('barcode_in_author_reference','sum'),exact_beta_matches=('exact_beta_reference_match_tolerance_1e12','sum')).to_string())
print('This checks literal barcode and exact classifier-profile duplicates; absence of these does not prove absence of different tumours from the same person in the reference.')
