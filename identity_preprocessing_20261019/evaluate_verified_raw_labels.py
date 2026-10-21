from pathlib import Path
import json
import pandas as pd,numpy as np
from sklearn.metrics import adjusted_rand_score,accuracy_score,balanced_accuracy_score,f1_score
from scipy.optimize import linear_sum_assignment
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
summaries=pd.concat([pd.read_csv(f) for f in (OUT/'raw_preprocessed').glob('*_summary.csv')],ignore_index=True)
mapping=json.loads((OUT/'training_only_semantic_map.json').read_text(encoding='utf-8'))['mapping']
canonical={'Merlin Intact':'Merlin-intact','Immune Enriched':'Immune-enriched','Hypermitotic':'Hypermitotic'}
mapping={float(k):canonical[v] for k,v in mapping.items()}
rows=[]
for cohort,expected in [('ucsf',200),('hku',365),('tubingen',149)]:
    if cohort=='ucsf':truth=pd.read_csv(OUT/'ucsf_ID_proven_molecular_linkage_200.csv')[['MID','GSM','Group']]
    elif cohort=='hku':truth=pd.read_csv(OUT/'hku_ID_linked_clinical_365.csv')[['ID','GSM','DNA_methylation_groups','identity_use_primary']].rename(columns={'DNA_methylation_groups':'Group'})
    else:
        cl=pd.read_csv(ROOT.parent/'results/gse304094_primary_survival.csv')
        selected=pd.read_csv(OUT/'requested_raw_tubingen.csv')[['GSM','title']]
        truth=selected.merge(cl[['Internal_ID','UCSF_Methylation_Group']],left_on='title',right_on='Internal_ID',validate='one_to_one')
        truth['Group']=truth.UCSF_Methylation_Group.map({'Merlinintact':'Merlin-intact','Immuneenriched':'Immune-enriched','Hypermitotic':'Hypermitotic'})
    p=summaries[summaries.cohort==cohort]
    if len(p)!=expected:
        print(cohort,'deferred: raw processed',len(p),'/',expected);continue
    d=truth.merge(p,on='GSM',validate='one_to_one');assert len(d)==expected and d.missing_classifier_cpgs.eq(0).all()
    d['semantic_prediction']=d.predicted_class.map(mapping)
    d.to_csv(OUT/(cohort+'_raw_author_ID_linked_predictions.csv'),index=False)
    for subset in ['all_assigned','highconfidence_075']+(['identity_concordant_assigned','identity_concordant_highconfidence_075'] if cohort=='hku' else []):
        a=d[d.semantic_prediction.notna()].copy()
        if 'highconfidence_075' in subset:a=a[a.confidence>=.75]
        if subset.startswith('identity_concordant'):a=a[a.identity_use_primary]
        ct=pd.crosstab(a.Group,a.semantic_prediction);ri,ci=linear_sum_assignment(-ct.to_numpy())
        ct.to_csv(OUT/(cohort+'_raw_author_confusion_'+subset+'.csv'))
        rows.append({'cohort':cohort,'subset':subset,'total_profiles':len(d),'assigned_profiles':len(a),'unassigned_top_vote_ties_all':int(d.tied_top.sum()),
            'accuracy_fixed_training_semantics':accuracy_score(a.Group,a.semantic_prediction),'balanced_accuracy_fixed_training_semantics':balanced_accuracy_score(a.Group,a.semantic_prediction),
            'macroF1_fixed_training_semantics':f1_score(a.Group,a.semantic_prediction,average='macro'),'ARI':adjusted_rand_score(a.Group,a.semantic_prediction),
            'bestmap_agreement_descriptive_only':ct.to_numpy()[ri,ci].sum()/len(a),'confidence_075_is_secondary':True})
if rows:
    pd.DataFrame(rows).to_csv(OUT/'raw_author_supervised_label_metrics.csv',index=False);print(pd.DataFrame(rows).to_string(index=False))
