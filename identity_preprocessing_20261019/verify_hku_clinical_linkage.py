from pathlib import Path
import pandas as pd,numpy as np,json,zipfile
from recover_longitudinal_centroid_links import sheet
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent;OUT=ROOT/'results'
clinical=pd.read_excel(ROOT/'source_data/zenodo_11501567/MeningiomaClinical.xlsx')
meta=pd.read_csv(OUT/'raw_idat_manifest.csv');meta=meta[meta.institution.eq('HKU')].copy()
meta['HKU_case_number']=meta.title.str.extract(r'^HKU Meningioma (\d+) \(methylation\)$')[0].astype(int)
clinical['HKU_case_number']=clinical.ID.str.extract(r'^QM(\d+)$')[0].astype(int)
assert not clinical.HKU_case_number.duplicated().any() and not meta.HKU_case_number.duplicated().any()
d=meta.merge(clinical,on='HKU_case_number',validate='one_to_one')
assert len(d)==365
d['independent_sex_agrees']=d.sex.astype(int).map({1:'M',2:'F'}).eq(d.Sex)
fusion=sheet(ROOT/'source_data/401_2024_2708_MOESM1_ESM.xlsx','SupplementaryTable4')
fusion.columns=fusion.iloc[1];fusion=fusion.iloc[2:].copy()
f=fusion[['Sample','Age','Sex','WHO_2021','DNA_methylation_group','Local_failure','LFFR','Death','OS']].copy()
d=d.merge(f,left_on='ID',right_on='Sample',how='left',validate='one_to_one',suffixes=('','_fusion2024'))
d['independent_age_abs_difference']=abs(pd.to_numeric(d.age)-pd.to_numeric(d.Age))
d['independent_grade_agrees']=pd.to_numeric(d.grade).eq(pd.to_numeric(d.WHO_2021))
assert d.Sample.notna().sum()==302
for col,other in [('LF','Local_failure'),('LFFR','LFFR_fusion2024'),('VS','Death'),('OS','OS_fusion2024')]:
    mask=d[other].notna();d['source_agrees_'+col]=False;d.loc[mask,'source_agrees_'+col]=np.isclose(pd.to_numeric(d.loc[mask,col]),pd.to_numeric(d.loc[mask,other]))
d['source_group_agrees']=d.DNA_methylation_groups.eq(d.DNA_methylation_group)
d['identity_demographic_conflict']=~d.independent_sex_agrees | (d.Sample.notna()&(d.independent_age_abs_difference>.0500001))
d['identity_use_primary']=~d.identity_demographic_conflict
d['grade_version_concordant_sensitivity']=d.identity_use_primary & (~d.Sample.notna()|d.independent_grade_agrees)
d['identity_method']='source_cohort_case_number_QMnn_HKU_Meningioma_nn; exact case-number normalization; independently checked demographics'
d['clinical_source']='Zenodo 11501567 MeningiomaClinical.xlsx; referenced by Nguyen et al Nat Commun 2025 doi:10.1038/s41467-025-61063-y'
d['GEO_male']=d.sex.astype(int).map({1:1,2:0})
d.to_csv(OUT/'hku_ID_linked_clinical_365.csv',index=False)
r={'n_ID_linked':len(d),'n_identity_concordant_primary':int(d.identity_use_primary.sum()),'n_demographic_conflicts':int(d.identity_demographic_conflict.sum()),'independent_sex_matches':int(d.independent_sex_agrees.sum()),'fusion2024_overlapping_cases':int(d.Sample.notna().sum()),
   'independent_age_max_abs_difference':d.independent_age_abs_difference.max(),'independent_grade_matches':int(d.loc[d.Sample.notna(),'independent_grade_agrees'].sum()),
   'event_count':int(d.LF.sum()),'time_unit':'years','endpoint':'local freedom from recurrence','stage':'posthoc identity repair external validation; no claim of prospective prespecification',
   'caveat':'QM to HKU title uses a cohort-prefix normalization of literal case numbers; all independent demographic discrepancies are flagged and excluded from primary analysis; person IDs are not supplied'}
(OUT/'hku_clinical_linkage_audit.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(r)
