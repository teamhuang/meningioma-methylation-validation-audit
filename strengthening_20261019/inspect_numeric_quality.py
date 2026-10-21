import json
from pathlib import Path
import numpy as np,pandas as pd
from run_strengthening import OUT,BASE,CORE,save

rows=[]
for name in CORE+['baylor_author_knn_beta_adapted']:
    d=pd.read_csv(OUT/('unpenalized_coefficients_'+name+'.csv'))
    key={'baylor_published':'menG_C','baylor_transferred':'pmc_3','tubingen_transferred':'pmc_3',
      'tubingen_published':'tor_prolif','ucsf_published':'hm','baylor_author_knn_beta_adapted':'official_3'}[name]
    r=d[d['covariate']==key].iloc[0]
    rows.append({'comparison':name,'contrast_column':key,'HR':r['exp(coef)'],'CI_low':r['exp(coef) lower 95%'],
        'CI_high':r['exp(coef) upper 95%'],'P_Wald_exploratory':r.p})
save(pd.DataFrame(rows),'cohort_specific_unpenalized_HR.csv')
checks=[]
def check(label,ok):
    checks.append({'check':label,'passed':bool(ok)});assert ok,label
s=pd.read_csv(OUT/'validation_summary.csv')
check('all planned fourteen comparisons saved',len(s)==14)
check('every comparison completed all twenty repeated CV runs',(s.cv_complete_repeats==20).all())
check('no CV fits failed',(s.cv_failed_folds==0).all())
check('six paired optimism analyses each have 1000 valid bootstrap samples',(s.loc[s.bootstrap_attempted>0,'bootstrap_success']==1000).all())
for name in s.comparison:
    pred=pd.read_csv(OUT/('oof_predictions_'+name+'.csv'))
    check(name+' each record predicted exactly once per repeat',pred.groupby(['repeat','id']).size().eq(1).all() and len(pred)==int(s[s.comparison==name].n.iloc[0])*20)
    check(name+' predictions finite and valid probabilities',np.isfinite(pred[['LP_clinical','LP_molecular','risk_clinical','risk_molecular']]).all().all() and pred[['risk_clinical','risk_molecular']].ge(0).all().all() and pred[['risk_clinical','risk_molecular']].le(1).all().all())
    inp=pd.read_csv(OUT/('analysis_input_'+name+'.csv'))
    check(name+' preserved complete-case denominator',len(inp)==int(s[s.comparison==name].n.iloc[0]))
    ass=pd.read_csv(OUT/('fold_assignments_'+name+'.csv'))
    check(name+' assignment and predicted test rows identical',set(zip(ass['repeat'],ass['fold'],ass['id']))==set(zip(pred['repeat'],pred['fold'],pred['id'])))
for _,r in pd.read_csv(OUT/'original_hashes.csv').iterrows():
    import hashlib
    check('original preserved '+r.path,hashlib.sha256((BASE/r.path).read_bytes()).hexdigest()==r.sha256)
training=pd.read_csv(OUT/'official_training_release_identity.csv')
check('121 training fingerprints unique in release186',training.in_official_training.sum()==121 and training.loc[training.in_official_training,'reference_nearest_barcode'].nunique()==121)
held=pd.read_csv(OUT/'official_knn_predictions_release37.csv')
check('20 training fingerprints unique in release37',held.in_official_training.sum()==20 and held.loc[held.in_official_training,'reference_nearest_barcode'].nunique()==20)
check('combined source releases cover all 141 unique model-training references',len(set(training.loc[training.in_official_training,'reference_nearest_barcode'])|set(held.loc[held.in_official_training,'reference_nearest_barcode']))==141)
tt=training.loc[training.in_official_training]
hh=held.loc[held.in_official_training]
check('all 141 fingerprint anchors agree with claimed series-order barcodes',tt.reference_nearest_barcode.eq(tt.barcode_positional).all() and hh.reference_nearest_barcode.eq(hh.barcode_positional).all())
(OUT/'package_validation.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
print(pd.DataFrame(rows).to_string(index=False));print('PASSED',len(checks),'checks')
