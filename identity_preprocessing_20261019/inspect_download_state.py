from pathlib import Path
import pandas as pd,requests,time,json
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
for cohort in ['baylor','tubingen','ucsf']:
    folder=ROOT/'raw_idat'/cohort
    raw=list(folder.glob('*.idat'));parts=list(folder.glob('*.part'))
    print(cohort,'complete files',len(raw),'partial',len(parts),'partial sizes',[(p.name,p.stat().st_size) for p in parts[:3]],flush=True)
    for p in parts[:1]:
        d=pd.read_csv(OUT/'raw_idat_manifest.csv');row=d[(d.cohort==cohort)&d.barcode.eq(p.name.split('_Grn')[0].split('_Red')[0])].iloc[0]
        color='Grn' if '_Grn' in p.name else 'Red'
        offset=p.stat().st_size
        try:
            r=requests.get(row[color+'_URL'],headers={'Range':'bytes='+str(offset)+'-','Accept-Encoding':'identity'},stream=True,timeout=15)
            print('range diagnostic',r.status_code,dict((k,v) for k,v in r.headers.items() if k.lower() in ['content-length','content-range','accept-ranges','content-encoding']),flush=True);r.close()
        except Exception as e:print(str(e))
