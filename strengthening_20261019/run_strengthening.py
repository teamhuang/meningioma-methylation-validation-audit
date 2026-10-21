"""Auditable retrospective statistical repair. Does not alter source results or draw figures."""
from pathlib import Path
import argparse, json, hashlib, platform, sys, time, warnings, traceback
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from scipy.stats import chi2
from sklearn.model_selection import StratifiedKFold
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.exceptions import ConvergenceWarning
from lifelines.statistics import proportional_hazard_test
from lifelines.utils.concordance import _concordance_summary_statistics
from statsmodels.stats.multitest import multipletests

BASE = Path(__file__).resolve().parent.parent
DEST = Path(__file__).resolve().parent
OUT = DEST / 'results'
OUT.mkdir(exist_ok=True)
SEED, REPEATS, FOLDS, NBOOT = 20261019, 20, 5, 1000
USE_SNAPSHOT=False

def init_snapshot(enabled):
    global USE_SNAPSHOT
    USE_SNAPSHOT=enabled

def save(df, name):
    df.to_csv(OUT / name, index=False, encoding='utf-8-sig')

def datasets():
    if USE_SNAPSHOT or not (BASE/'results/clinical_parsed.csv').exists():
        snapshot=json.loads((OUT/'specifications.json').read_text(encoding='utf-8'))
        return {name:(pd.read_csv(OUT/('analysis_input_'+name+'.csv')),x['clinical'],x['molecular'],x['horizon']) for name,x in snapshot.items()}
    c = pd.read_csv(BASE/'results/clinical_parsed.csv', dtype={'Sample':str})
    p = pd.read_csv(BASE/'results/gse189521_linked_endpoints.csv', dtype={'Sample':str})
    assert c.Sample.is_unique and p.Sample.is_unique
    c = c.merge(p[['Sample','pmc']], on='Sample', how='left', validate='one_to_one')
    d = pd.DataFrame({'id':c.Sample, 'T':c.fu_months, 'E':c.rec,
        'grade2':(c.grade==2).astype(float),'gtr':c.gtr,'age':c.age,
        'male':c.Gender.str.upper().map({'M':1.,'F':0.})})
    for label in ['B','C']:
        d['menG_'+label] = np.where(c.menG.isin(['A','B','C']), (c.menG==label).astype(float), np.nan)
    for k in [2,3]:
        d['pmc_'+str(k)] = np.where(c.pmc.notna(),(pd.to_numeric(c.pmc)==k).astype(float),np.nan)
    t = pd.read_csv(BASE/'results/gse304094_primary_survival.csv')
    t['id'] = t.Internal_ID
    t['age'] = pd.to_numeric(t['Age [Years]'],errors='coerce')
    group = t['Molecular_Group (by Toronto nomenclature)'].astype('string')
    for col, name in [('tor_immunogenic','Immunogenic'),('tor_prolif','Proliferative'),('tor_hypermet','Hypermetabolic')]:
        t[col] = np.where(group.notna(),group.str.lower().eq(name.lower()).fillna(False).astype(float),np.nan)
    for k in [2,3]:
        t['grade_'+str(k)] = np.where(t.grade.notna(),(t.grade==k).astype(float),np.nan)
    u = pd.read_csv(BASE/'results/natgen205_survival.csv')
    u['id'] = ['MOESM19_row_'+str(i) for i in range(len(u))]
    for k in [2,3]:
        u['grade_'+str(k)] = np.where(u.grade.notna(),(u.grade==k).astype(float),np.nan)
    specs = {
      'baylor_published':(d,['grade2','gtr','age','male'],['menG_B','menG_C'],60.),
      'baylor_transferred':(d,['grade2','gtr','age','male'],['pmc_2','pmc_3'],60.),
      'tubingen_transferred':(t,['grade','gtr','age','male'],['pmc_2','pmc_3'],60.),
      'tubingen_published':(t,['grade','gtr','age','male'],['tor_prolif','tor_hypermet'],60.),
      'ucsf_published':(u,['grade','gtr','male'],['hm','ie'],5.),
      'tubingen_four_groups':(t,['grade','gtr','age','male'],['tor_immunogenic','tor_hypermet','tor_prolif'],60.),
      'tubingen_categorical_grade':(t,['grade_2','grade_3','gtr','age','male'],['tor_prolif','tor_hypermet'],60.),
      'ucsf_categorical_grade':(u,['grade_2','grade_3','gtr','male'],['hm','ie'],5.),
      'ucsf_mib1':(u,['grade','gtr','male'],['mib1'],5.),
      'ucsf_groups_beyond_mib1':(u,['grade','gtr','male','mib1'],['hm','ie'],5.),
      'ucsf_positive_followup':(u[u['T']>0].copy(),['grade','gtr','male'],['hm','ie'],5.),
      'baylor_published_positive_followup':(d[d['T']>0].copy(),['grade2','gtr','age','male'],['menG_B','menG_C'],60.),
    }
    fp=OUT/'official_knn_predictions_baylor.csv'
    if fp.exists():
        pp=pd.read_csv(fp)
        lk=pd.read_csv(BASE/'results/gse189521_linked_endpoints.csv',dtype={'Sample':str})
        pp=lk[['GSM','Sample']].merge(pp,left_on='GSM',right_on='id',validate='one_to_one')
        dd=d.merge(pp[['Sample','predicted_class','confidence','tied_top']],left_on='id',right_on='Sample',validate='one_to_one')
        for k in [2,3]: dd['official_'+str(k)]=np.where(dd.predicted_class.notna(),(dd.predicted_class==k).astype(float),np.nan)
        specs['baylor_author_knn_beta_adapted']=(dd,['grade2','gtr','age','male'],['official_2','official_3'],60.)
        specs['baylor_author_knn_beta_highconfidence']=(dd[dd.confidence>=.75].copy(),['grade2','gtr','age','male'],['official_2','official_3'],60.)
    return specs

CORE=['baylor_published','baylor_transferred','tubingen_transferred','tubingen_published','ucsf_published']

def fit(d,cols,pen=.05):
    active=[x for x in cols if d[x].nunique()>1]
    with warnings.catch_warnings():
        warnings.simplefilter('error', ConvergenceWarning)
        m=CoxPHFitter(penalizer=pen)
        m.fit(d[['T','E']+active], 'T','E')
    if not np.all(np.isfinite(m.params_)): raise ValueError('nonfinite coefficients')
    return m

def concord(d,score):
    a,b,n=_concordance_summary_statistics(d['T'].to_numpy(),-np.asarray(score),d.E.to_numpy())
    if n==0: raise ValueError('no comparable pairs')
    return float((a+.5*b)/n), float(a+.5*b),int(n)

def G_at(train, times, left=False):
    """KM censor survival with event-first handling of event/censor ties."""
    km=KaplanMeierFitter().fit(train['T'],1-train.E)
    x=km.survival_function_.index.to_numpy(); y=km.survival_function_.iloc[:,0].to_numpy()
    idx=np.searchsorted(x,np.asarray(times),side='left' if left else 'right')-1
    return np.where(idx<0,1.,y[np.maximum(idx,0)])

def weights(train,test,h):
    ev=(test['T'].to_numpy()<=h)&(test.E.to_numpy()==1)
    ct=test['T'].to_numpy()>h
    w=np.zeros(len(test)); ge=G_at(train,test.loc[ev,'T'].to_numpy(),True)
    gh=float(G_at(train,[h])[0])
    if h>train['T'].max() or gh<=0 or np.any(ge<=0):
        return None,ev,ct,gh
    w[ev]=1/ge; w[ct]=1/gh
    return w,ev,ct,gh

def horizon_metrics(train,test,r0,r1,h,tag):
    w,y,ct,gh=weights(train,test,h)
    base={'horizon':h,'n':len(test),'cases':int(y.sum()),'controls':int(ct.sum()),'G_horizon':gh,
          'at_risk':int((test['T']>=h).sum()),'supported':w is not None,'evaluation':tag}
    if w is None: return base
    base['observed_ipcw']=float(np.sum(w*y)/len(test))
    for name,r in [('clinical',r0),('molecular',r1)]:
        base['predicted_'+name]=float(np.mean(r))
        base['brier_'+name]=float(np.sum(w*(y-r)**2)/len(test))
        base['OP_'+name]=base['observed_ipcw']/np.mean(r) if np.mean(r)>0 else np.nan
        # Pairwise cumulative/dynamic AUC; control weights constant within a test fold.
        a=np.asarray(r)[y,None]; b=np.asarray(r)[ct][None,:]
        den=w[y].sum()*w[ct].sum()
        base['auc_'+name]=float(np.sum(((a>b)+.5*(a==b))*w[y,None]*w[ct][None,:])/den) if den>0 else np.nan
    base['delta_brier']=base['brier_molecular']-base['brier_clinical']
    return base

def audit():
    specs=datasets(); rows=[]
    (OUT/'specifications.json').write_text(json.dumps({name:{'clinical':c,'molecular':m,'horizon':h} for name,(d,c,m,h) in specs.items()},indent=2),encoding='utf-8')
    for name,(d,c,m,h) in specs.items():
        dd=d[['id','T','E']+c+m].dropna()
        rows.append({'comparison':name,'input_records':len(d),'complete_records':len(dd),
          'events':int(dd.E.sum()),'excluded_missing':len(d)-len(dd),'zero_time':int((dd['T']==0).sum()),
          'negative_time':int((dd['T']<0).sum()),'duplicate_ids':int(dd.id.duplicated().sum()),
          'invalid_events':int((~dd.E.isin([0,1])).sum()),'predictor_parameters':len(c+m),
          'events_per_parameter':dd.E.sum()/len(c+m),'clinical_columns':','.join(c),'molecular_columns':','.join(m)})
        miss=d[['T','E']+c+m].isna().sum().rename_axis('field').reset_index(name='missing')
        save(miss,f'missing_{name}.csv')
        save(dd,f'analysis_input_{name}.csv')
    save(pd.DataFrame(rows),'denominators.csv')
    ng=BASE/'data/suppl/natgen2022'
    a=pd.read_excel(ng/'41588_2022_1061_MOESM19_ESM.xlsx')
    b=pd.read_excel(ng/'41588_2022_1061_MOESM20_ESM.xlsx')
    fields=[x for x in b.columns if x in a.columns]
    ka=a[fields].astype(str).agg('|'.join,axis=1); kb=b[fields].astype(str).agg('|'.join,axis=1)
    links=[]
    for i,key in enumerate(kb):
        candidates=[int(x) for x in np.flatnonzero(ka==key)]
        links.append({'M20_row':i,'n_candidate_M19':len(candidates),'candidate_M19_rows':json.dumps(candidates),
          'candidate_groups':'|'.join(a.iloc[candidates].Group.astype(str)),
          'identity_proven':False,'unique_clinical_key':len(candidates)==1})
    save(pd.DataFrame(links),'ucsf_clinical_key_cardinality.csv')
    same=kb[kb.duplicated(keep=False)]
    save(b.loc[same.index].assign(M20_row=same.index),'ucsf_duplicate_clinical_keys.csv')
    clin=pd.read_csv(BASE/'results/clinical_parsed.csv')
    # Candidate duplicates are never collapsed without a patient identifier.
    save(clin[clin.duplicated(['Gender','age','Location','fu_months','rec'],keep=False)],'baylor_candidate_duplicate_records.csv')
    long=pd.read_csv(BASE/'results/gse304096_stability.csv')
    save(pd.DataFrame([{'computed_clusters':len(long),'observed_concordant':int(long.all_same.sum()),
          'all_class3_concordant_clusters':len(long),'observed_all_class3_clusters':int((long.classes.str.replace('|','',regex=False).str.strip('3')=='').sum()),
          'interpretation':'Constant-label assignment has 100% within-cluster agreement; 15/17 is observed all-class3 clusters, not a constant-class baseline.'}]),'longitudinal_baseline_correction.csv')
    for f in ['harmonization_gse270371.csv','harmonization_gse189521.csv','harmonization_ucsfsesame.csv']:
        d=pd.read_csv(BASE/'results'/f)
        n=994 if '270371' in f else 110 if '189521' in f else 200
        if 'n_reject' in d:
            d['profile_denominator']=n; d['reject_percent_recomputed']=100*d.n_reject/n
        save(d,'audited_'+f)
    legacy=pd.read_csv(BASE/'results/primary_deltaC_5000.csv').rename(columns={
       'P_one_sided':'descriptive_bootstrap_fraction_delta_le_zero','P_two_sided':'twice_smaller_descriptive_tail'})
    legacy['valid_null_calibrated_P']=False
    save(legacy,'legacy_bootstrap_reinterpreted.csv')
    files=[BASE/'manuscript_npjDM.md',BASE/'results/primary_deltaC_5000.csv',BASE/'results/calibration_slope.csv']
    save(pd.DataFrame([{'path':str(f.relative_to(BASE)),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in files]),'original_hashes.csv')
    env={'python':sys.version,'platform':platform.platform(),'seed':SEED,'repeats':REPEATS,'folds':FOLDS,'bootstrap':NBOOT}
    import lifelines, sklearn, scipy, statsmodels
    for mod in [np,pd,lifelines,sklearn,scipy,statsmodels]: env[mod.__name__]=mod.__version__
    (OUT/'environment.json').write_text(json.dumps(env,indent=2),encoding='utf-8')
    print(pd.DataFrame(rows)[['comparison','complete_records','events','zero_time']].to_string(index=False),flush=True)
    print('UCSF clinical keys:',pd.Series([r['n_candidate_M19'] for r in links]).value_counts().to_dict(),flush=True)

def run(name):
    started=time.time(); df,clin,mol,hz=datasets()[name]
    d=df[['id','T','E']+clin+mol].dropna().reset_index(drop=True)
    assert d.id.is_unique and d.E.isin([0,1]).all() and d['T'].ge(0).all()
    c0=fit(d,clin); c1=fit(d,clin+mol)
    app0=concord(d,c0.predict_log_partial_hazard(d))[0]; app1=concord(d,c1.predict_log_partial_hazard(d))[0]
    coef=c1.summary.reset_index(); coef['comparison']=name; save(coef,f'penalized_coefficients_{name}.csv')
    preds=[]; splits=[]; reps=[]; metrics=[]; dc=[]; failures=[]; slopes=[]
    # Same split seed across comparable specifications, enabling paired checks.
    for rep in range(REPEATS):
        cv=StratifiedKFold(FOLDS,shuffle=True,random_state=SEED+rep)
        sums=np.zeros(3); risk0=np.full(len(d),np.nan); risk1=risk0.copy(); folds=np.zeros(len(d),int)
        lp0=risk0.copy(); lp1=risk0.copy(); nfailed=0
        for fold,(tr,te) in enumerate(cv.split(d,d.E)):
            train,test=d.iloc[tr],d.iloc[te]
            for ix in te: splits.append({'comparison':name,'repeat':rep,'fold':fold,'id':d.id.iloc[ix],'event':d.E.iloc[ix]})
            try:
                a,b=fit(train,clin),fit(train,clin+mol)
                x=a.predict_log_partial_hazard(test).to_numpy(); y=b.predict_log_partial_hazard(test).to_numpy()
                ac,an,pairs=concord(test,x); bc,bn,pairs2=concord(test,y); assert pairs==pairs2
                sums+=np.array([an,bn,pairs]); lp0[te]=x; lp1[te]=y; folds[te]=fold
                p0=1-a.predict_survival_function(test,times=[hz]).iloc[0].to_numpy()
                p1=1-b.predict_survival_function(test,times=[hz]).iloc[0].to_numpy()
                risk0[te]=p0; risk1[te]=p1
                for ix,x0,x1,r0,r1 in zip(te,x,y,p0,p1):
                    preds.append({'comparison':name,'repeat':rep,'fold':fold,'id':d.id.iloc[ix],
                        'T':d['T'].iloc[ix],'E':d.E.iloc[ix],'LP_clinical':x0,'LP_molecular':x1,
                        'risk_clinical':r0,'risk_molecular':r1,'horizon':hz})
                for h in [hz*.6,hz]:
                    rr0=1-a.predict_survival_function(test,times=[h]).iloc[0].to_numpy()
                    rr1=1-b.predict_survival_function(test,times=[h]).iloc[0].to_numpy()
                    mm=horizon_metrics(train,test,rr0,rr1,h,'out_of_fold')
                    metrics.append(dict(comparison=name,repeat=rep,fold=fold,**mm))
                w,ev,ct,gh=weights(train,test,hz)
                if w is not None:
                    for threshold in np.arange(.05,.501,.05):
                        odds=threshold/(1-threshold)
                        dc.append({'comparison':name,'repeat':rep,'fold':fold,'threshold':threshold,'n':len(test),
                          'NB_clinical':float(np.sum(w*(p0>=threshold)*(ev-ct*odds))/len(test)),
                          'NB_molecular':float(np.sum(w*(p1>=threshold)*(ev-ct*odds))/len(test)),
                          'NB_all':float(np.sum(w*(ev-ct*odds))/len(test)),'NB_none':0.})
            except Exception as exc:
                failures.append({'stage':'cv','repeat':rep,'fold':fold,'error':str(exc)}); nfailed+=1
        reps.append({'comparison':name,'repeat':rep,'C_clinical':sums[0]/sums[2],
           'C_molecular':sums[1]/sums[2],'deltaC':(sums[1]-sums[0])/sums[2],
           'comparable_pairs':int(sums[2]),'failed_folds':nfailed,'complete_repeat':nfailed==0})
        if nfailed==0:
            for label,lp in [('clinical',lp0),('molecular',lp1)]:
                try:
                    sd=pd.DataFrame({'T':d['T'],'E':d.E,'LP':lp,'fold':folds})
                    cm=CoxPHFitter().fit(sd,'T','E',strata='fold')
                    ss=cm.summary.loc['LP']
                    slopes.append({'comparison':name,'repeat':rep,'model':label,'slope':ss['coef'],
                        'CI_low':ss['coef lower 95%'],'CI_high':ss['coef upper 95%']})
                except Exception as exc: failures.append({'stage':'cv_slope','repeat':rep,'error':str(exc)})
    save(pd.DataFrame(preds),f'oof_predictions_{name}.csv')
    save(pd.DataFrame(splits),f'fold_assignments_{name}.csv')
    save(pd.DataFrame(reps),f'cv_repeats_{name}.csv')
    save(pd.DataFrame(metrics),f'oof_horizon_metrics_{name}.csv')
    save(pd.DataFrame(dc),f'oof_decision_metrics_{name}.csv')
    save(pd.DataFrame(slopes),f'oof_calibration_slopes_{name}.csv')
    # Correct apparent calibration on the original, unstandardized LP; slope is coef, not HR.
    apparent=[]
    for h in [hz*.6,hz]:
        p0=1-c0.predict_survival_function(d,times=[h]).iloc[0].to_numpy()
        p1=1-c1.predict_survival_function(d,times=[h]).iloc[0].to_numpy()
        mm=horizon_metrics(d,d,p0,p1,h,'apparent')
        km=KaplanMeierFitter().fit(d['T'],d.E)
        mm['observed_KM']=float(1-km.predict(h))
        mm['KM_risk_CI_low']=float(1-km.confidence_interval_survival_function_.iloc[np.searchsorted(km.timeline,h,side='right')-1,1])
        mm['KM_risk_CI_high']=float(1-km.confidence_interval_survival_function_.iloc[np.searchsorted(km.timeline,h,side='right')-1,0])
        apparent.append(dict(comparison=name,**mm))
    save(pd.DataFrame(apparent),f'appar ent_horizon_{name}.csv'.replace(' ',''))
    appsl=[]
    for label,m in [('clinical',c0),('molecular',c1)]:
        lp=m.predict_log_partial_hazard(d).to_numpy()
        cm=CoxPHFitter().fit(pd.DataFrame({'T':d['T'],'E':d.E,'LP':lp}),'T','E')
        ss=cm.summary.loc['LP']
        appsl.append({'comparison':name,'model':label,'slope':ss['coef'],'CI_low':ss['coef lower 95%'],
            'CI_high':ss['coef upper 95%'],'LP_SD':np.std(lp,ddof=1)})
    save(pd.DataFrame(appsl),f'apparent_calibration_slopes_{name}.csv')
    # Exploratory unpenalized association testing: never use penalized likelihood differences.
    lr={'comparison':name,'n':len(d),'events':int(d.E.sum()),'status':'failed'}
    try:
        a,b=fit(d,clin,0),fit(d,clin+mol,0)
        stat=2*(b.log_likelihood_-a.log_likelihood_); dfree=len(b.params_)-len(a.params_)
        lr.update(status='ok',LR=stat,df=dfree,P=chi2.sf(stat,dfree))
        ss=b.summary.reset_index(); ss['comparison']=name; save(ss,f'unpenalized_coefficients_{name}.csv')
        ph=proportional_hazard_test(b,d[['T','E']+list(b.params_.index)],time_transform='rank').summary.reset_index()
        ph['comparison']=name; save(ph,f'PH_diagnostics_{name}.csv')
    except Exception as exc: lr['error']=str(exc)
    save(pd.DataFrame([lr]),f'unpenalized_LRT_{name}.csv')
    # More expensive paired optimism correction only for the five core comparisons.
    boot=[]; rng=np.random.default_rng(SEED)
    if name in CORE+['baylor_author_knn_beta_adapted']:
        for k in range(NBOOT):
            ix=rng.integers(0,len(d),len(d)); db=d.iloc[ix].reset_index(drop=True)
            try:
                a,b=fit(db,clin),fit(db,clin+mol)
                train0=concord(db,a.predict_log_partial_hazard(db))[0]; train1=concord(db,b.predict_log_partial_hazard(db))[0]
                test0=concord(d,a.predict_log_partial_hazard(d))[0]; test1=concord(d,b.predict_log_partial_hazard(d))[0]
                boot.append({'bootstrap':k,'C_train_clinical':train0,'C_train_molecular':train1,
                    'deltaC_apparent':train1-train0,'optimism_clinical':train0-test0,
                    'optimism_molecular':train1-test1,'optimism_delta':(train1-train0)-(test1-test0)})
            except Exception as exc: failures.append({'stage':'bootstrap','bootstrap':k,'error':str(exc)})
            if (k+1)%250==0: print(name,'bootstrap',k+1,'elapsed',round(time.time()-started),flush=True)
    save(pd.DataFrame(boot),f'bootstrap_optimism_{name}.csv')
    save(pd.DataFrame(failures,columns=['stage','repeat','fold','bootstrap','error']),f'failures_{name}.csv')
    reps=pd.DataFrame(reps); valid=reps[reps.complete_repeat]
    result={'comparison':name,'n':len(d),'events':int(d.E.sum()),'C_app_clinical':app0,'C_app_molecular':app1,
       'delta_app':app1-app0,'cv_complete_repeats':len(valid),'cv_C_clinical':valid.C_clinical.mean(),
       'cv_C_molecular':valid.C_molecular.mean(),'cv_delta_mean':valid.deltaC.mean(),
       'cv_delta_split_p025':valid.deltaC.quantile(.025),'cv_delta_split_p975':valid.deltaC.quantile(.975),
       'cv_delta_split_SD':valid.deltaC.std(),'cv_positive_repeats':int((valid.deltaC>0).sum()),
       'cv_failed_folds':int(reps.failed_folds.sum()),'bootstrap_attempted':NBOOT if name in CORE+['baylor_author_knn_beta_adapted'] else 0,'bootstrap_success':len(boot)}
    if boot:
        bs=pd.DataFrame(boot)
        result.update(C_corrected_clinical=app0-bs.optimism_clinical.mean(),
            C_corrected_molecular=app1-bs.optimism_molecular.mean(),delta_corrected=(app1-app0)-bs.optimism_delta.mean(),
            optimism_delta_MCSE=bs.optimism_delta.std()/np.sqrt(len(bs)),
            delta_app_CI_low=bs.deltaC_apparent.quantile(.025),delta_app_CI_high=bs.deltaC_apparent.quantile(.975),
            descriptive_fraction_delta_le_zero=(bs.deltaC_apparent<=0).mean())
    save(pd.DataFrame([result]),f'summary_{name}.csv')
    print(name,'DONE',json.dumps(result),flush=True)
    return result

def combine():
    sums=pd.concat([pd.read_csv(f) for f in sorted(OUT.glob('summary_*.csv'))],ignore_index=True)
    save(sums,'validation_summary.csv')
    lr=pd.concat([pd.read_csv(f) for f in sorted(OUT.glob('unpenalized_LRT_*.csv'))],ignore_index=True)
    lr['P_Holm_core_family']=np.nan; mask=lr.comparison.isin(CORE)&lr.P.notna()
    lr.loc[mask,'P_Holm_core_family']=multipletests(lr.loc[mask,'P'],method='holm')[1]
    lr['P_Holm_all_tested']=np.nan; mask=lr.P.notna()
    lr.loc[mask,'P_Holm_all_tested']=multipletests(lr.loc[mask,'P'],method='holm')[1]
    save(lr,'unpenalized_association_tests.csv')
    allmet=[]; allnb=[]; allsl=[]
    for name in sums.comparison:
        met=pd.read_csv(OUT/f'oof_horizon_metrics_{name}.csv')
        # Only complete, supported five-fold repeats. Fold AUC not averaged into cohort AUC.
        for (rep,h),dd in met.groupby(['repeat','horizon']):
            if len(dd)!=FOLDS or not dd.supported.all(): continue
            row={'comparison':name,'repeat':rep,'horizon':h,'n':int(dd.n.sum()),'min_G_horizon':dd.G_horizon.min()}
            for col in ['observed_ipcw','predicted_clinical','predicted_molecular','brier_clinical','brier_molecular','delta_brier']:
                row[col]=np.average(dd[col],weights=dd.n)
            row['OP_clinical']=row['observed_ipcw']/row['predicted_clinical']
            row['OP_molecular']=row['observed_ipcw']/row['predicted_molecular']
            allmet.append(row)
        nb=pd.read_csv(OUT/f'oof_decision_metrics_{name}.csv')
        for (rep,t),dd in nb.groupby(['repeat','threshold']):
            if len(dd)!=FOLDS: continue
            row={'comparison':name,'repeat':rep,'threshold':t}
            for col in ['NB_clinical','NB_molecular','NB_all','NB_none']:row[col]=np.average(dd[col],weights=dd.n)
            allnb.append(row)
        sl=pd.read_csv(OUT/f'apparent_calibration_slopes_{name}.csv');allsl.append(sl)
    save(pd.DataFrame(allmet),'oof_risk_metrics_by_repeat.csv')
    save(pd.DataFrame(allnb),'oof_decision_metrics_by_repeat.csv')
    save(pd.concat(allsl,ignore_index=True),'apparent_calibration_slopes.csv')
    print(sums[['comparison','n','events','cv_delta_mean','delta_corrected']].to_string(index=False),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--audit',action='store_true');parser.add_argument('--route');parser.add_argument('--all',action='store_true');parser.add_argument('--combine',action='store_true');parser.add_argument('--snapshot',action='store_true')
    args=parser.parse_args()
    USE_SNAPSHOT=args.snapshot
    if args.audit:audit()
    if args.route:run(args.route)
    if args.all:
        with ProcessPoolExecutor(max_workers=3,initializer=init_snapshot,initargs=(USE_SNAPSHOT,)) as pool: list(pool.map(run,datasets()))
        combine()
    if args.combine:combine()
