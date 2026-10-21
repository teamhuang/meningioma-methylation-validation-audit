from pathlib import Path
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from download_idats import fetch
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
ref=pd.read_csv(ROOT.parent/'strengthening_20261019/results/official_model/official_knn_training.csv')[['sample','group']]
meta=pd.read_csv(OUT/'raw_idat_manifest.csv');meta=meta[meta.cohort=='author186']
d=meta.merge(ref,left_on='barcode',right_on='sample',validate='one_to_one').groupby('group',as_index=False).head(1)
assert len(d)==3 and set(d.group)=={1,2,3}
d.to_csv(OUT/'requested_raw_author186.csv',index=False)
tasks=[(r,c) for r in d.to_dict('records') for c in ['Grn','Red']]
with ThreadPoolExecutor(max_workers=1) as pool:rows=list(pool.map(fetch,tasks))
pd.DataFrame(rows).to_csv(OUT/'download_log_reference_anchors.csv',index=False)
print(pd.DataFrame(rows)[['GSM','color','status']].to_string(index=False))
