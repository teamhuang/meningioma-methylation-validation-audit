from pathlib import Path
import pandas as pd,numpy as np,json
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
mp=pd.read_csv(OUT/'longitudinal_verified_patient_mapping_42.csv')
summ=pd.concat([pd.read_csv(f) for f in (OUT/'raw_preprocessed').glob('*_summary.csv')],ignore_index=True)
sub=summ[summ.cohort=='longitudinal']
if len(sub)!=42:raise SystemExit('Deferred: only '+str(len(sub))+'/42 longitudinal profiles processed')
d=mp.merge(sub[['GSM','predicted_class','confidence','tied_top']],on='GSM',validate='one_to_one',suffixes=('_legacy','_author'))
d.to_csv(OUT/'longitudinal_source_patient_and_raw_author_predictions.csv',index=False)
pairs=pd.read_csv(OUT/'SNP_pairwise_identity_metrics.csv');pairmap={tuple(sorted([r.GSM1,r.GSM2])):r for _,r in pairs.iterrows()}
rows=[]
from itertools import combinations
for patient,g in d.groupby('source_patient_centroid'):
    for a,b in combinations(g.GSM,2):
        p=pairmap[tuple(sorted([a,b]))]
        rows.append({'source_patient':patient,'GSM1':a,'GSM2':b,'common_highconfidence_SNPs':p.common_highconfidence_SNPs,'discordant_SNPs':p.discordant_SNPs,
          'discordance_fraction':p.discordance_fraction,'SNP_beta_correlation':p.SNP_beta_correlation,'primary_SNP_threshold_agrees':p.same_donor_candidate,
          'strict_SNP_threshold_agrees':p.strict_candidate,'relaxed_SNP_threshold_agrees':p.relaxed_candidate})
r=pd.DataFrame(rows);r.to_csv(OUT/'longitudinal_source_groups_independent_SNP_validation.csv',index=False)
source_group=d.set_index('GSM').source_patient_centroid.to_dict()
between=pairs[pairs.GSM1.isin(source_group)&pairs.GSM2.isin(source_group)&pairs.relaxed_candidate].copy()
between['source_group1']=between.GSM1.map(source_group);between['source_group2']=between.GSM2.map(source_group)
between=between[between.source_group1.ne(between.source_group2)]
between.to_csv(OUT/'longitudinal_cross_source_group_SNP_candidates.csv',index=False)
patient=[]
for pid,g in d.groupby('source_patient_centroid'):
    pr=r[r.source_patient==pid]
    assigned=g.predicted_class_author.notna().all()
    patient.append({'source_patient':pid,'n_samples':len(g),'all_assigned_author':assigned,'author_class_concordant_if_all_assigned':g.predicted_class_author.nunique()==1 if assigned else np.nan,
        'author_classes':'|'.join(g.predicted_class_author.fillna(-1).astype(int).astype(str)),
        'legacy_centroid_concordant':g.predicted_class_legacy.nunique()==1,'all_pairs_SNP_primary_agree':pr.primary_SNP_threshold_agrees.all(),
        'all_pairs_SNP_relaxed_agree':pr.relaxed_SNP_threshold_agrees.all(),'max_pair_SNP_discordance':pr.discordance_fraction.max(),
        'source_grade_group_conflict':pid==18})
pd.DataFrame(patient).to_csv(OUT/'longitudinal_raw_patient_stability_summary.csv',index=False)
print(pd.DataFrame(patient).to_string(index=False))
status={'source_labelled_groups_recovered':18,'profiles':42,'n_source_groups_with_all_pairs_relaxed_SNP_support':int(pd.DataFrame(patient).all_pairs_SNP_relaxed_agree.sum()),
    'cross_source_group_candidates':sorted(set(tuple(sorted([int(a),int(b)])) for a,b in zip(between.source_group1,between.source_group2))),
    'independent_patient_count_confirmed':False,'visit_order_confirmed':False,
    'interpretation':'Literal/numeric source-record grouping recovered; independent donor identity remains discordant. Do not infer 18 independent patients or temporal CpG changes.'}
(OUT/'longitudinal_identity_resolution_status.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
print(status)
print('Every-visit majority-class assignment gives 18/18 concordant patient groups; raw class stability is descriptive and cannot outperform that trivial baseline.')
