from pathlib import Path
import json,re
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results';BASE=ROOT.parent
ref=pd.read_csv(BASE/'strengthening_20261019/results/official_model/official_knn_training.csv',index_col=0)
a=pd.read_csv(OUT/'raw_preprocessed/GSM8856328_403beta.csv').set_index('cpg').beta
raw=pd.read_csv(OUT/'requested_raw_pilot.csv');r=raw[raw.GSM=='GSM8856328'].iloc[0]
b=ref.loc[r.barcode,a.index].astype(float)
delta=np.abs(a-b)
result={'sample':r.GSM,'barcode':r.barcode,'features':len(a),'max_abs_beta_diff':delta.max(),'mean_abs_beta_diff':delta.mean(),'correlation':np.corrcoef(a,b)[0,1]}
(OUT/'author_preprocessing_pilot_comparison.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(result)
html=(OUT/'baylor_public_box.html').read_text(encoding='utf-8')
for mat in re.finditer(r'"typedID"',html):print(html[mat.start():mat.start()+650])
for x in ['Box.postStreamData','window.Box','accessToken','requestToken']:
    at=html.find(x);print('SCRIPT',x,at,html[max(0,at-80):at+80])
