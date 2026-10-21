from pathlib import Path
import pandas as pd,numpy as np,json
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent;OUT=ROOT/'results'
manifest=pd.read_csv(OUT/'raw_idat_manifest.csv')
u=manifest[manifest.cohort=='ucsf'].copy()
u['MID']=u.title.str.extract(r'UCSF Meningioma (\d+)')[0].map(lambda v:'M'+v if pd.notna(v) else None)
u=u[u.MID.notna()].copy()
native=pd.read_csv(OUT/'ucsf_ID_proven_published_groups.csv')
assert native.groupby('MID').Group.nunique().eq(1).all()
native=native.groupby('MID',as_index=False).agg(Group=('Group','first'),source_rows=('excel_row',lambda x:'|'.join(str(q) for q in x)),source_sheets=('source_sheet',lambda x:'|'.join(str(q) for q in x)))
u=u.merge(native,on='MID',validate='one_to_one')
fingerprint=pd.read_excel(BASE/'data/suppl/natgen2022/41588_2022_1061_MOESM21_ESM.xlsx',index_col=0)
u['in_fingerprint']=u.MID.isin(fingerprint.columns)
assert len(u)==200 and u.in_fingerprint.all()
u.to_csv(OUT/'ucsf_ID_proven_molecular_linkage_200.csv',index=False)
old=pd.read_csv(BASE/'results/ucsftoronto_sample_mapping.csv').dropna(subset=['MID'])
comp=old.merge(u[['MID','Group','sex','grade','age']],on='MID',validate='one_to_one')
cl=pd.read_excel(BASE/'data/suppl/natgen2022/41588_2022_1061_MOESM20_ESM.xlsx')
comp['sex_agrees']=pd.to_numeric(comp.sex).eq(comp.M20_row_index.map(cl.Sex))
comp['grade_agrees']=pd.to_numeric(comp.grade).eq(comp.M20_row_index.map(cl['WHO grade']))
comp['published_group_agrees']=comp.published_Group.eq(comp.Group)
comp.to_csv(OUT/'ucsf_old_positional_mapping_rejected.csv',index=False)
print('UCSF direct unique ID links',len(u),'group distribution',u.Group.value_counts().to_dict())
print('Old positional agreement sex/grade/group',comp[['sex_agrees','grade_agrees','published_group_agrees']].mean().to_dict())
stats=[]
for source,idcol,cols in [(BASE/'results/harmonization_pred_ucsfsesame.csv','MID',['pred_P1','pred_P2','pred_P3'])]:
    p=pd.read_csv(source).merge(u[['MID','Group']],on='MID',validate='one_to_one')
    for col in cols:
        ct=pd.crosstab(p[col],p.Group);ri,ci=linear_sum_assignment(-ct.to_numpy())
        stats.append({'path':col,'n':len(p),'ARI':adjusted_rand_score(p.Group,p[col]),'bestmap_agreement':ct.to_numpy()[ri,ci].sum()/len(p),'reference':'published_Group_direct_MID_source_Table3'})
        ct.to_csv(OUT/('ucsf_supervised_confusion_'+col+'.csv'))
pd.DataFrame(stats).to_csv(OUT/'ucsf_ID_proven_supervised_legacy_transfer.csv',index=False)
print(pd.DataFrame(stats).to_string(index=False))
os=pd.read_excel(BASE/'data/suppl/natgen2022/41588_2022_1061_MOESM6_ESM.xlsx',sheet_name='Sfigure_2e')
links=[]
for _,r in u.iterrows():
    candidates=os.index[(np.abs(pd.to_numeric(os.Age,errors='coerce')-float(r.age))<=.050000001)&(pd.to_numeric(os['WHO grade'],errors='coerce')==float(r.grade))&os.Group.eq(r.Group)].tolist()
    links.append({'MID':r.MID,'Group':r.Group,'GSM':r.GSM,'GEO_age':r.age,'GEO_grade':r.grade,'candidate_rows':json.dumps(candidates),'n_candidates':len(candidates)})
ll=pd.DataFrame(links);ll.to_csv(OUT/'ucsf_age_grade_group_linkage_candidates.csv',index=False)
print('Age-grade-group matches',ll.n_candidates.value_counts().to_dict())
