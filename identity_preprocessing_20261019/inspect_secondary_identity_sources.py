from pathlib import Path
import pandas as pd,json,requests,re
from recover_longitudinal_centroid_links import sheet
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
for file in (ROOT/'source_data').glob('*.xlsx'):
    if file.name.startswith(('401_','41467_2025')):
        book=pd.ExcelFile(file)
        for sh in book.sheet_names:
            d=sheet(file,sh)
            print('SHEET',file.name,sh,d.shape,repr(d.iloc[:4,:25].fillna('').values.tolist())[:3200],flush=True)
for rec in ['11501566','11501381']:
    r=requests.get('https://zenodo.org/api/records/'+rec,timeout=45);r.raise_for_status();z=r.json()
    (OUT/('zenodo_'+rec+'_record.json')).write_text(json.dumps(z,indent=2),encoding='utf-8')
    print('ZENODO',rec,z.get('metadata',{}).get('title'),[(f['key'],f.get('size'),f.get('links')) for f in z.get('files',[])],flush=True)
