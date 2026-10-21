from pathlib import Path
import gzip,re
import pandas as pd
BASE=Path(__file__).resolve().parent.parent
OUT=Path(__file__).resolve().parent/'results';OUT.mkdir(exist_ok=True)
rows=[]
for cohort,acc in [('baylor','GSE189521'),('tubingen','GSE304094'),('ucsf','GSE183647'),('longitudinal','GSE304096'),('author186','GSE292325'),('author37','GSE292326')]:
    fields=[]
    with gzip.open(BASE/'data'/(acc+'_series_matrix.txt.gz'),'rt',encoding='utf-8',errors='replace') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):break
            bits=[x.strip().strip('"') for x in line.rstrip().split('\t')]
            if bits[0].startswith('!Sample_'):fields.append(bits)
    gsm=next(x[1:] for x in fields if x[0]=='!Sample_geo_accession')
    titles=next(x[1:] for x in fields if x[0]=='!Sample_title')
    chars=[x[1:] for x in fields if x[0]=='!Sample_characteristics_ch1']
    urls=[x[1:] for x in fields if x[0]=='!Sample_supplementary_file']
    for i,(g,t) in enumerate(zip(gsm,titles)):
        d={'cohort':cohort,'accession':acc,'GSM':g,'title':t}
        for vals in chars:
            if ':' in vals[i]:
                key,value=vals[i].split(':',1);d[key.strip().lower()]=value.strip()
        for vals in urls:
            u=vals[i]
            if '.idat' not in u:continue
            match=re.search(r'_(\d+_R\d+C\d+)_(Grn|Red)\.idat',u)
            if match:d['barcode']=match.group(1);d[match.group(2)+'_URL']=u.replace('ftp://','https://')
        rows.append(d)
df=pd.DataFrame(rows)
df.loc[(df.cohort=='ucsf')&df.institution.eq('HKU'),'cohort']='hku'
df.to_csv(OUT/'raw_idat_manifest.csv',index=False,encoding='utf-8-sig')
print(df.groupby('cohort').agg(samples=('GSM','size'),idat_pairs=('Grn_URL','count')).to_string())
