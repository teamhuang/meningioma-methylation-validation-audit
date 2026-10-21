from pathlib import Path
import hashlib,json,re,requests
import pandas as pd
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent;OUT=ROOT/'results'
files=[*list((BASE/'data/suppl/natgen2022').glob('*.xlsx')),
    BASE/'data/suppl/protocadherin_suppldata.xlsx',BASE/'data/suppl/protocadherin_sourcedata.xlsx',
    ROOT/'source_data/zenodo_11501567/MeningiomaClinical.xlsx',ROOT/'source_data/401_2024_2708_MOESM1_ESM.xlsx']
rows=[]
for f in files:
    doi='10.1038/s41588-022-01061-8' if 'natgen2022' in f.parts else '10.1038/s41467-026-77170-3' if f.name.startswith('protocadherin') else '10.1038/s41467-025-61063-y' if 'MeningiomaClinical' in f.name else '10.1007/s00401-024-02708-y'
    rows.append({'local_file':str(f.relative_to(BASE)),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'source_article':'https://doi.org/'+doi})
pd.DataFrame(rows).to_csv(OUT/'original_source_workbook_hashes.csv',index=False)
url='https://www.nature.com/articles/s41467-026-77170-3'
r=requests.get(url,timeout=40)
urls=sorted(set(requests.compat.urljoin(url,x.replace('&amp;','&')) for x in re.findall(r'href=["\']([^"\']+\.xlsx)["\']',r.text)))
(OUT/'protocadherin_original_supplement_urls.json').write_text(json.dumps({'article':url,'status':r.status_code,'xlsx_urls':urls},indent=2),encoding='utf-8')
print('Hashed original workbooks',len(rows),'original supplement URLs',urls)
