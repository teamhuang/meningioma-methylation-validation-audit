from pathlib import Path
import pandas as pd,json
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
d=pd.read_csv(OUT/'hku_ID_linked_clinical_365.csv');snp=pd.read_csv(OUT/'SNP_donor_clusters_current.csv')[['GSM','callable_SNPs','donor_overlap_guard_key']]
p=pd.concat([pd.read_csv(f) for f in (OUT/'raw_preprocessed').glob('*_summary.csv')],ignore_index=True);p=p[p.cohort=='hku']
d=d.merge(snp,on='GSM',validate='one_to_one').merge(p[['GSM','predicted_class','tied_top']],on='GSM',validate='one_to_one')
d['SNP_guard_infeasible']=d.callable_SNPs.lt(35)
d['eligible_before_classification']=d.identity_use_primary&~d.SNP_guard_infeasible
d['clinical_covariates_complete']=d[['LFFR','LF','age','grade','GEO_male']].notna().all(axis=1)&d.EOR.isin(['GTR','STR'])
d['eligible_author_KNN_primary']=d.eligible_before_classification&~d.tied_top&d.clinical_covariates_complete
d['eligible_published_primary']=d.eligible_before_classification&d.clinical_covariates_complete
d.to_csv(OUT/'hku_final_primary_exclusion_ledger.csv',index=False)
z={'source_records':len(d),'demographic_conflicts':int(d.identity_demographic_conflict.sum()),'SNP_below_35':int(d.SNP_guard_infeasible.sum()),
    'identity_and_SNP_eligible':int(d.eligible_before_classification.sum()),'author_primary_complete_assigned':int(d.eligible_author_KNN_primary.sum()),
    'published_primary_complete':int(d.eligible_published_primary.sum()),'event_count_author_primary':int(d.loc[d.eligible_author_KNN_primary,'LF'].sum())}
(OUT/'hku_final_primary_selection_summary.json').write_text(json.dumps(z,indent=2),encoding='utf-8');print(z)
