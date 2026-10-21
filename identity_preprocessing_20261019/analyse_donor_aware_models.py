from pathlib import Path
import sys,json,time,warnings,argparse
import pandas as pd,numpy as np
from sklearn.model_selection import StratifiedGroupKFold
from lifelines import CoxPHFitter
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent
sys.path.insert(0,str(BASE/'strengthening_20261019'))
import donor_model_utilities as legacy
OUT=ROOT/'results/donor_aware';OUT.mkdir(exist_ok=True)
SEED=2026101902

def raw_summary(cohort):
    files=list((ROOT/'results/raw_preprocessed').glob('*_summary.csv'))
    d=pd.concat([pd.read_csv(f) for f in files],ignore_index=True)
    return d[d.cohort==cohort]

def prepare(cohort):
    snapshot=ROOT/'results/model_snapshot'/('prepared_'+cohort+'.csv')
    if snapshot.exists():
        specs=json.loads((snapshot.parent/'specifications.json').read_text(encoding='utf-8'))[cohort]
        return pd.read_csv(snapshot,float_precision='round_trip'),specs['clinical'],specs['layers'],specs['horizon']
    import run_strengthening as previous
    ds=previous.datasets()
    if cohort=='baylor':
        d,clin,mol,h=ds['baylor_published'];_,_,surrogate,_=ds['baylor_transferred']
        mapping=pd.read_csv(BASE/'results/gse189521_linked_endpoints.csv')[['Sample','GSM']]
        d=d.merge(mapping,left_on='id',right_on='Sample',validate='one_to_one');expect=110
    elif cohort=='tubingen':
        d,clin,mol,h=ds['tubingen_published'];_,_,surrogate,_=ds['tubingen_transferred']
        mapping=pd.read_csv(ROOT/'results/raw_idat_manifest.csv');mapping=mapping[mapping.cohort=='tubingen'][['title','GSM']]
        d=d.merge(mapping,left_on='id',right_on='title',validate='one_to_one');expect=149
        d=d.dropna(subset=['T','E']+clin+mol)
    elif cohort=='hku':
        d=pd.read_csv(ROOT/'results/hku_ID_linked_clinical_365.csv')
        d['id']=d.ID;d['T']=d.LFFR;d['E']=d.LF;d['gtr']=d.EOR.map({'GTR':1.,'STR':0.})
        d['male']=d.GEO_male;d['age']=pd.to_numeric(d.age);d['grade']=pd.to_numeric(d.grade)
        d['published_immune']=d.DNA_methylation_groups.eq('Immune-enriched').astype(float)
        d['published_hyper']=d.DNA_methylation_groups.eq('Hypermitotic').astype(float)
        clin=['grade','gtr','age','male'];mol=['published_immune','published_hyper'];surrogate=[];h=5.;expect=365
    else:raise ValueError(cohort)
    snp=pd.read_csv(ROOT/'results/SNP_donor_clusters_current.csv')[['GSM','donor_strict_key','donor_overlap_guard_key','callable_SNPs']]
    d=d.merge(snp,on='GSM',validate='one_to_one')
    if len(d)!=expect:raise ValueError(f'{cohort}: only {len(d)}/{expect} SNP samples; full-cohort analysis deferred')
    preds=raw_summary(cohort)
    d=d.merge(preds[['GSM','predicted_class','missing_classifier_cpgs','confidence','tied_top']],on='GSM',validate='one_to_one')
    if len(d)!=expect or d.missing_classifier_cpgs.ne(0).any():raise ValueError('raw preprocessing incomplete')
    for k in [2,3]:d['author_'+str(k)]=np.where(d.predicted_class.notna(),d.predicted_class.eq(k).astype(float),np.nan)
    if cohort=='hku':d=d[d.identity_use_primary&d.callable_SNPs.ge(35)].copy()
    return d,clin,{'published':mol,'centroid_surrogate':surrogate,'author_KNN_raw':['author_2','author_3']},h

def analyse(cohort,layer,group_rule,one_record=False,common_set=False,grade_concordant=False):
    started=time.time();d,clin,layers,h=prepare(cohort);mol=layers[layer]
    donor='donor_overlap_guard_key' if group_rule=='overlap_guard' else 'donor_strict_key'
    if grade_concordant:
        if cohort!='hku':raise ValueError('grade concordance sensitivity is defined only for HKU')
        d=d[d.grade_version_concordant_sensitivity].copy()
    if common_set:d=d.dropna(subset=sorted(set(clin+sum(layers.values(),[]))))
    d=d[['id','GSM','T','E',donor]+clin+mol].dropna().reset_index(drop=True)
    d=d.rename(columns={donor:'donor'})
    if one_record:d=d.sort_values('GSM').drop_duplicates('donor',keep='first').reset_index(drop=True)
    name=cohort+'_'+layer+'_'+group_rule+('_one_record' if one_record else '')+('_common_set' if common_set else '')+('_grade_concordant' if grade_concordant else '')
    d.to_csv(OUT/('input_'+name+'.csv'),index=False)
    c0=legacy.fit(d,clin);c1=legacy.fit(d,clin+mol)
    app0=legacy.concord(d,c0.predict_log_partial_hazard(d))[0];app1=legacy.concord(d,c1.predict_log_partial_hazard(d))[0]
    repeats=[];predictions=[];failures=[]
    for rep in range(20):
        cv=StratifiedGroupKFold(5,shuffle=True,random_state=SEED+rep);tot=np.zeros(3)
        for fold,(tr,te) in enumerate(cv.split(d,d.E,d.donor)):
            a,b=d.iloc[tr],d.iloc[te]
            assert not set(a.donor)&set(b.donor)
            try:
                m0=legacy.fit(a,clin);m1=legacy.fit(a,clin+mol)
                x=m0.predict_log_partial_hazard(b).to_numpy();y=m1.predict_log_partial_hazard(b).to_numpy()
                _,c,n=legacy.concord(b,x);_,cc,nn=legacy.concord(b,y);assert n==nn;tot+=np.array([c,cc,n])
                for (_,r),s0,s1 in zip(b.iterrows(),x,y):predictions.append({'repeat':rep,'fold':fold,'id':r.id,'donor':r.donor,'LP_clinical':s0,'LP_molecular':s1})
            except Exception as exc:failures.append({'stage':'CV','repeat':rep,'fold':fold,'error':str(exc)})
        repeats.append({'repeat':rep,'C_clinical':tot[0]/tot[2],'C_molecular':tot[1]/tot[2],'deltaC':(tot[1]-tot[0])/tot[2]})
    boot=[];rng=np.random.default_rng(SEED);units=d.donor.unique();by={x:d[d.donor==x] for x in units}
    nboot=1000 if group_rule=='overlap_guard' and not one_record and not common_set and not grade_concordant else 0
    for k in range(nboot):
        drawn=rng.choice(units,len(units),replace=True)
        db=pd.concat([by[x] for x in drawn],ignore_index=True)
        try:
            a=legacy.fit(db,clin);b=legacy.fit(db,clin+mol)
            train0=legacy.concord(db,a.predict_log_partial_hazard(db))[0];train1=legacy.concord(db,b.predict_log_partial_hazard(db))[0]
            test0=legacy.concord(d,a.predict_log_partial_hazard(d))[0];test1=legacy.concord(d,b.predict_log_partial_hazard(d))[0]
            boot.append({'bootstrap':k,'n_records':len(db),'n_events':int(db.E.sum()),'delta_app':train1-train0,
                'optimism_clinical':train0-test0,'optimism_molecular':train1-test1,'optimism_delta':(train1-train0)-(test1-test0)})
        except Exception as exc:failures.append({'stage':'bootstrap','bootstrap':k,'error':str(exc)})
    rp=pd.DataFrame(repeats);bp=pd.DataFrame(boot)
    for obj,suffix in [(rp,'cv'),(pd.DataFrame(predictions),'oof'),(bp,'bootstrap'),(pd.DataFrame(failures),'failures')]:obj.to_csv(OUT/(suffix+'_'+name+'.csv'),index=False)
    r={'comparison':name,'records':len(d),'events':int(d.E.sum()),'donor_clusters':d.donor.nunique(),'C_app_clinical':app0,'C_app_molecular':app1,'delta_app':app1-app0,
       'C_cv_clinical':rp.C_clinical.mean(),'C_cv_molecular':rp.C_molecular.mean(),'delta_cv':rp.deltaC.mean(),
       'delta_cv_split_p025':rp.deltaC.quantile(.025),'delta_cv_split_p975':rp.deltaC.quantile(.975),'positive_repeats':int((rp.deltaC>0).sum()),
       'bootstrap_success':len(bp),'bootstrap_attempted':nboot,'n_failures':len(failures)}
    if nboot and len(bp):r.update(delta_app_CI_low=bp.delta_app.quantile(.025),delta_app_CI_high=bp.delta_app.quantile(.975),
       C_corrected_clinical=app0-bp.optimism_clinical.mean(),C_corrected_molecular=app1-bp.optimism_molecular.mean(),
       delta_corrected=(app1-app0)-bp.optimism_delta.mean())
    # Report robust individual coefficients; do not substitute an independence-based global LR test.
    try:
        m=CoxPHFitter().fit(d[['T','E','donor']+clin+mol],'T','E',cluster_col='donor',robust=True)
        m.summary.to_csv(OUT/('cluster_robust_cox_'+name+'.csv'))
    except Exception as exc:r['cluster_robust_cox_error']=str(exc)
    pd.DataFrame([r]).to_csv(OUT/('summary_'+name+'.csv'),index=False)
    print('DONE',json.dumps(r),'seconds',time.time()-started,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cohort',required=True);p.add_argument('--layer',required=True);p.add_argument('--group-rule',default='overlap_guard',choices=['overlap_guard','strict']);p.add_argument('--one-record-per-donor',action='store_true');p.add_argument('--common-complete-set',action='store_true');p.add_argument('--grade-concordant-only',action='store_true');a=p.parse_args()
    analyse(a.cohort,a.layer,a.group_rule,a.one_record_per_donor,a.common_complete_set,a.grade_concordant_only)
