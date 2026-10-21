import json,hashlib,zipfile,py_compile
from pathlib import Path
import numpy as np,pandas as pd
import run_strengthening as m
from run_strengthening import DEST,OUT,BASE

checks=[]
def check(label,ok):
    checks.append({'check':label,'passed':bool(ok)});assert ok,label
raw=m.datasets();m.USE_SNAPSHOT=True;snap=m.datasets()
for name,(df,clin,mol,h) in raw.items():
    expected=df[['id','T','E']+clin+mol].dropna().reset_index(drop=True)
    actual=snap[name][0]
    check(name+' snapshot ids identical',expected.id.astype(str).tolist()==actual.id.astype(str).tolist())
    cols=['T','E']+clin+mol
    check(name+' snapshot numbers identical within decimal serialization tolerance',np.allclose(expected[cols],actual[cols],rtol=1e-12,atol=1e-12))
for file in DEST.glob('*.py'):py_compile.compile(str(file),doraise=True)
check('no new figure files',not any(DEST.rglob('*.png')) and not any(DEST.rglob('*.pdf')) and not any(DEST.rglob('*.svg')))
check('all fourteen comparisons retained',len(snap)==14)
ph=pd.read_csv(OUT/'PH_diagnostics_all.csv')
check('PH summary has no duplicated model variable records',not ph.duplicated(['comparison','index']).any())
rank=pd.read_csv(OUT/'robust_panel_rank_audit.csv')
check('corrected panel includes exactly 2000 lowest rank sums',rank.selected_corrected.sum()==2000 and rank.loc[rank.selected_corrected,'score'].max()<=rank.loc[~rank.selected_corrected,'score'].min())
(OUT/'snapshot_validation.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
print('Snapshot and release checks passed:',len(checks))

files=[f for f in DEST.rglob('*') if f.is_file() and '__pycache__' not in f.parts and f.name!='SHA256SUMS.csv']
manifest=[]
for f in sorted(files):manifest.append({'path':f.relative_to(DEST).as_posix(),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
pd.DataFrame(manifest).to_csv(DEST/'SHA256SUMS.csv',index=False,encoding='utf-8-sig')
archive=BASE/'npjDM_data_strengthening_20261019.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for f in files+[DEST/'SHA256SUMS.csv']:z.write(f,arcname='strengthening_20261019/'+f.relative_to(DEST).as_posix())
with zipfile.ZipFile(archive) as z:
    check('zip CRC test',z.testzip() is None)
    check('zip count agrees with hash manifest',len(z.namelist())==len(manifest)+1)
print('Archive:',archive,'files',len(manifest)+1,'bytes',archive.stat().st_size)
