from pathlib import Path
import json,re,requests,gzip
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results';BASE=ROOT.parent
html=(OUT/'baylor_public_box.html').read_text(encoding='utf-8')
start=html.index('Box.postStreamData = ')+len('Box.postStreamData = ')
obj,_=json.JSONDecoder().raw_decode(html[start:])
folder=obj['/app-api/enduserapp/shared-folder']
items=folder['items'];print('BOX files',[(x['id'],x['name'],x.get('itemSize')) for x in items])
(OUT/'baylor_box_file_inventory.json').write_text(json.dumps([{'id':x['id'],'name':x['name'],'bytes':x.get('itemSize')} for x in items],indent=2),encoding='utf-8')
for x in items:
    if any(w in x['name'].lower() for w in ['clinical','metadata','sample','information']):
        u='https://bcm.box.com/s/nsd4pwofiitxvp0zh1sp56way2z8mkti/file/'+str(x['id'])
        r=requests.get(u,timeout=40);(OUT/('baylor_box_file_'+str(x['id'])+'.html')).write_text(r.text,encoding='utf-8')
        print('FILE',x['name'],r.status_code,len(r.content))
        for match in re.finditer(r'(?i)download.{0,20}',r.text):
            at=match.start();snippet=r.text[max(0,at-50):at+350]
            if 'http' in snippet:print('LINK',snippet)

# Look up source-native, explicitly identified group assignments.
book=BASE/'data/suppl/natgen2022/41588_2022_1061_MOESM4_ESM.xlsx'
native=[]
for sheet in ['3','4','5']:
    d=pd.read_excel(book,sheet_name=sheet,header=None)
    for i,row in d.iterrows():
        for j,value in row.items():
            if re.fullmatch(r'M\d+',str(value).strip()):
                neighbors=row.iloc[max(0,j-1):min(len(row),j+4)].astype(str).tolist()
                group=next((q for q in neighbors if q in ['Merlin-intact','Immune-enriched','Hypermitotic']),None)
                if group:native.append({'MID':str(value).strip(),'Group':group,'source_sheet':sheet,'excel_row':i+1})
nd=pd.DataFrame(native).drop_duplicates()
nd.to_csv(OUT/'ucsf_ID_proven_published_groups.csv',index=False)
print('NATIVE UCSF group IDs',nd.MID.nunique(),nd.groupby('MID').Group.nunique().value_counts().to_dict())

# Test the previously assumed UCSF clinical-row order against independently labelled GEO sex and grade.
manifest=pd.read_csv(OUT/'raw_idat_manifest.csv');u=manifest[manifest.cohort=='ucsf'].copy()
u['MID']=u.title.str.extract(r'UCSF Meningioma (\d+)')[0].map(lambda v:'M'+v if pd.notna(v) else None)
u=u[u.MID.notna()].copy()
chain=pd.read_csv(BASE/'results/ucsftoronto_sample_mapping.csv').dropna(subset=['MID'])
chain=chain.merge(u[['MID','sex','grade','age']],on='MID',validate='one_to_one')
clin=pd.read_excel(BASE/'data/suppl/natgen2022/41588_2022_1061_MOESM20_ESM.xlsx')
chain['positional_sex']=chain.M20_row_index.map(clin.Sex)
chain['positional_grade']=chain.M20_row_index.map(clin['WHO grade'])
chain['sex_agrees']=pd.to_numeric(chain.sex,errors='coerce').eq(chain.positional_sex)
chain['grade_agrees']=pd.to_numeric(chain.grade,errors='coerce').eq(chain.positional_grade)
chain.to_csv(OUT/'ucsf_positional_bridge_independent_check.csv',index=False)
print('POSITIONAL UCSF independent sex/grade concordance',chain.sex_agrees.mean(),chain.grade_agrees.mean(),(chain.sex_agrees&chain.grade_agrees).mean())

# Age metadata may support unambiguous deterministic linkage, never forced nearest-neighbour matching.
os=pd.read_excel(BASE/'data/suppl/natgen2022/41588_2022_1061_MOESM6_ESM.xlsx',sheet_name='Sfigure_2e')
links=[]
for _,r in u.dropna(subset=['MID']).iterrows():
    age=float(r.age);grade=float(r.grade)
    cand=os.index[(np.abs(pd.to_numeric(os.Age,errors='coerce')-age)<=.050000001)&(pd.to_numeric(os['WHO grade'],errors='coerce')==grade)].tolist()
    links.append({'MID':r.MID,'GSM':r.GSM,'GEO_age':age,'GEO_grade':grade,'candidate_rows':json.dumps(cand),'n_candidates':len(cand),
      'candidate_groups':'|'.join(os.iloc[cand].Group.astype(str))})
ll=pd.DataFrame(links);ll.to_csv(OUT/'ucsf_age_grade_linkage_candidates.csv',index=False)
print('UCSF age-grade candidate counts',ll.n_candidates.value_counts().to_dict())
