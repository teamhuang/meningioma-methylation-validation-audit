from pathlib import Path
import json,requests,hashlib
import pandas as pd
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results';DEST=ROOT/'source_data/zenodo_11501567';DEST.mkdir(exist_ok=True)
z=json.loads((OUT/'zenodo_11501566_record.json').read_text(encoding='utf-8'));rows=[]
for f in z['files']:
    if f['key'] not in ['MeningiomaClinical.xlsx','MeningiomaSegs.zip','hg38.14_chrom_arms_manual.txt']:continue
    path=DEST/f['key'];url=f['links']['self']
    if not path.exists():
        r=requests.get(url,timeout=120);r.raise_for_status();path.write_bytes(r.content)
    assert path.stat().st_size==f['size']
    expected=f.get('checksum','')
    if expected.startswith('md5:'):assert hashlib.md5(path.read_bytes()).hexdigest()==expected.split(':',1)[1]
    rows.append({'file':str(path.relative_to(ROOT)),'url':url,'bytes':f['size'],'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'source_checksum':expected,'checksum_pass':True})
(OUT/'zenodo_identity_download_audit.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
for sh,d in pd.read_excel(DEST/'MeningiomaClinical.xlsx',sheet_name=None).items():
    print(sh,d.shape,d.columns.tolist());print(d.head(5).to_string(index=False));d.to_csv(OUT/('zenodo_clinical_'+sh.replace('/','_')+'.csv'),index=False)
