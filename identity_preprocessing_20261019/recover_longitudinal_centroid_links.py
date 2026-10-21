from pathlib import Path
from xml.etree import ElementTree as ET
from itertools import combinations
import zipfile,json,re
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent;OUT=ROOT/'results'
NS={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
def sheet(book,name):
    with zipfile.ZipFile(book) as z:
        wb=ET.fromstring(z.read('xl/workbook.xml'));rels=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        links={x.attrib['Id']:x.attrib['Target'] for x in rels};shared=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            for _,e in ET.iterparse(z.open('xl/sharedStrings.xml'),events=['end']):
                if e.tag.endswith('}si'):shared.append(''.join(e.itertext()));e.clear()
        sh=next(x for x in wb.find('m:sheets',NS) if x.attrib['name']==name)
        ident=sh.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'];path=links[ident].lstrip('/')
        if not path.startswith('xl/'):path='xl/'+path
        rows=[]
        for _,e in ET.iterparse(z.open(path),events=['end']):
            if e.tag.endswith('}row'):
                row={}
                for c in e:
                    v=c.find('m:v',NS);value=v.text if v is not None else ''.join(c.itertext())
                    if c.attrib.get('t')=='s' and v is not None:value=shared[int(value)]
                    address=re.match('[A-Z]+',c.attrib.get('r','A')).group(0)
                    n=0
                    for letter in address:n=n*26+ord(letter)-64
                    row[n-1]=value
                rows.append(row);e.clear()
        maxcol=max(max(x,default=0) for x in rows)
        return pd.DataFrame([[r.get(c,np.nan) for c in range(maxcol+1)] for r in rows])

if __name__=='__main__':
    book=BASE/'data/suppl/protocadherin_sourcedata.xlsx'
    individual=sheet(book,'Fig. 5f');centroids=sheet(book,'Supplementary Fig. 15d')
    for name,d in [('pca_individual',individual),('pca_centroids',centroids)]:d.to_csv(OUT/('longitudinal_source_'+name+'.csv'),index=False,header=False)
    individual=individual.iloc[2:,1:10].copy();individual.columns=['ID','PC1','PC2','curve_x','curve_y','signature','grade','setting','group']
    centroids=centroids.iloc[2:,1:6].copy();centroids.columns=['centroid_id','PC1','PC2','signature','condition']
    manifest=pd.read_csv(OUT/'raw_idat_manifest.csv');ids=manifest[manifest.cohort=='longitudinal'].title
    d=individual[individual.ID.isin(ids)].copy()
    for col in ['PC1','PC2','signature']:d[col]=pd.to_numeric(d[col],errors='coerce');centroids[col]=pd.to_numeric(centroids[col],errors='coerce')
    d=d.dropna(subset=['PC1','PC2','signature']).reset_index(drop=True)
    c=centroids[centroids.condition.astype(str).str.contains('corresponding',case=False)].dropna(subset=['PC1','PC2','signature']).copy()
    print('Source coordinates, longitudinal profiles',len(d),'patient centroids',len(c));print(c.to_string(index=False))
    X=d[['PC1','PC2','signature']].to_numpy()
    combos=[ix for n in [2,3] for ix in combinations(range(len(d)),n)]
    means=np.array([X[list(ix)].mean(axis=0) for ix in combos]);rows=[]
    for _,r in c.iterrows():
        target=r[['PC1','PC2','signature']].to_numpy(float)
        delta=np.max(np.abs(means-target)/np.array([1.,1.,.01]),axis=1)
        best=np.argsort(delta)[:3]
        for rank,k in enumerate(best,1):
            rows.append({'centroid_id':r.centroid_id,'rank':rank,'n':len(combos[k]),'IDs':'|'.join(d.iloc[list(combos[k])].ID),
                'max_PC_difference':np.max(np.abs(means[k,:2]-target[:2])),
                'signature_difference':abs(means[k,2]-target[2]),'scaled_max_difference':delta[k],
                'exact_numeric_candidate':np.max(np.abs(means[k,:2]-target[:2]))<1e-6 and abs(means[k,2]-target[2])<1e-8})
    pd.DataFrame(rows).to_csv(OUT/'longitudinal_centroid_identity_candidates.csv',index=False)
    print(pd.DataFrame(rows).query('rank==1').to_string(index=False))
