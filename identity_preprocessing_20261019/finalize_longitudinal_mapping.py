from pathlib import Path
import pandas as pd,numpy as np,json
from recover_longitudinal_centroid_links import sheet
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent;OUT=ROOT/'results'
links=pd.read_csv(OUT/'longitudinal_centroid_identity_candidates.csv')
best=links[links['rank']==1].copy();assert len(best)==18 and best.exact_numeric_candidate.all()
second=links[links['rank']==2].set_index('centroid_id')
assert second.scaled_max_difference.min()>1e-4
members=[]
for _,r in best.iterrows():
    for sample in r.IDs.split('|'):
        members.append({'Internal_ID':sample,'source_patient_centroid':int(r.centroid_id),
           'PC_max_abs_error':r.max_PC_difference,'signature_abs_error':r.signature_difference,
           'linkage':'unique exact source centroid decomposition of PC1 PC2 and methylation signature; not correlation clustering'})
mp=pd.DataFrame(members);assert len(mp)==42 and mp.Internal_ID.is_unique
manifest=pd.read_csv(OUT/'raw_idat_manifest.csv');lm=manifest[manifest.cohort=='longitudinal']
mp=mp.merge(lm[['GSM','title','barcode']],left_on='Internal_ID',right_on='title',validate='one_to_one')
clinical=pd.read_excel(BASE/'data/suppl/protocadherin_suppldata.xlsx',sheet_name='Supplementary Data 4',header=2)
source=sheet(BASE/'data/suppl/protocadherin_sourcedata.xlsx','Fig. 5f').iloc[2:,1:10].copy()
source.columns=['Internal_ID','PC1','PC2','curve_x','curve_y','signature','grade','setting','group']
mp=mp.merge(source[['Internal_ID','grade','group']],on='Internal_ID',validate='one_to_one')
checks=[]
for pid,sub in mp.groupby('source_patient_centroid'):
    cl=clinical[pd.to_numeric(clinical.Patient,errors='coerce')==pid]
    observed=sorted(zip(pd.to_numeric(sub.grade),sub.group))
    expected=sorted(zip(pd.to_numeric(cl.WHO_Grade_2021),cl['Molecular_Group (by Toronto nomenclature)']))
    checks.append({'patient':pid,'n_samples':len(sub),'count_agrees':len(sub)==len(cl),'grade_group_multiset_agrees':observed==expected})
check=pd.DataFrame(checks);check.to_csv(OUT/'longitudinal_source_patient_verification.csv',index=False)
pred=pd.read_csv(BASE/'results/gse304096_classification.csv')
print('Predicted table columns',pred.columns.tolist())
barcode_col=next(x for x in pred if 'sample' in x.lower() or 'barcode' in x.lower())
mp=mp.merge(pred,left_on='barcode',right_on=barcode_col,validate='one_to_one')
mp.to_csv(OUT/'longitudinal_verified_patient_mapping_42.csv',index=False)
cohort=manifest[manifest.cohort=='tubingen']
mp2=mp.merge(cohort[['GSM','title','barcode']],left_on='Internal_ID',right_on='title',how='left',suffixes=('','_discovery'))
mp2['same_physical_barcode_in_discovery']=mp2.barcode.eq(mp2.barcode_discovery)
mp2.to_csv(OUT/'longitudinal_discovery_sample_overlap.csv',index=False)
print('Source patient checks',check.to_string(index=False))
print('Overlap discovery series',mp2.GSM_discovery.notna().sum(),'same physical barcode',mp2.same_physical_barcode_in_discovery.sum())
print('42 profiles partitioned as',mp.groupby('source_patient_centroid').size().value_counts().to_dict())
