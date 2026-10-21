from pathlib import Path
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results';rows=[]
snp=pd.read_csv(OUT/'SNP_donor_clusters_current.csv')[['GSM','donor_overlap_guard_key']]
labels=['Merlin-intact','Immune-enriched','Hypermitotic'];index={x:i for i,x in enumerate(labels)}
for cohort in ['tubingen','ucsf','hku']:
    full=pd.read_csv(OUT/(cohort+'_raw_author_ID_linked_predictions.csv')).merge(snp,on='GSM',validate='one_to_one')
    for subset in ['all_assigned','highconfidence_075']+(['identity_concordant_assigned','identity_concordant_highconfidence_075'] if cohort=='hku' else []):
        d=full[full.semantic_prediction.notna()].copy()
        if 'highconfidence_075' in subset:d=d[d.confidence>=.75]
        if subset.startswith('identity_concordant'):d=d[d.identity_use_primary]
        units=d.donor_overlap_guard_key.unique();ix={x:i for i,x in enumerate(units)}
        unit_counts=np.zeros((len(units),9),dtype=int)
        for r in d.itertuples():unit_counts[ix[r.donor_overlap_guard_key],3*index[r.Group]+index[r.semantic_prediction]]+=1
        rng=np.random.default_rng(2026101905);weights=rng.multinomial(len(units),np.repeat(1/len(units),len(units)),size=1000)
        ct=(weights@unit_counts).reshape(-1,3,3);n=ct.sum(axis=(1,2));correct=ct[:,np.arange(3),np.arange(3)].sum(axis=1)
        accuracy=correct/n;row=ct.sum(axis=2);col=ct.sum(axis=1);diagonal=ct[:,np.arange(3),np.arange(3)]
        with np.errstate(divide='ignore',invalid='ignore'):
            balanced=np.mean(diagonal/row,axis=1);macro=np.mean(2*diagonal/(row+col),axis=1)
            comb=lambda a:a*(a-1)/2
            sums=comb(ct).sum(axis=(1,2));rs=comb(row).sum(axis=1);cs=comb(col).sum(axis=1);expect=rs*cs/comb(n)
            ari=(sums-expect)/((rs+cs)/2-expect)
        metric={'accuracy':accuracy,'balanced_accuracy':balanced,'macroF1':macro,'ARI':ari}
        eligible=int(full.identity_use_primary.sum()) if subset.startswith('identity_concordant') else len(full)
        z={'cohort':cohort,'subset':subset,'profiles':len(d),'donor_clusters':len(units),'coverage_all_profiles':len(d)/len(full),'eligible_profiles':eligible,'coverage_eligible_profiles':len(d)/eligible,'bootstrap_attempted':1000}
        output=pd.DataFrame(metric);output.to_csv(OUT/(cohort+'_label_cluster_bootstrap_'+subset+'.csv'),index=False)
        for key,value in metric.items():z[key+'_CI_low'],z[key+'_CI_high']=np.nanquantile(value,[.025,.975]);z[key+'_finite_bootstraps']=int(np.isfinite(value).sum())
        rows.append(z)
result=pd.DataFrame(rows);result.to_csv(OUT/'raw_author_label_cluster_bootstrap_intervals.csv',index=False)
point=pd.read_csv(OUT/'raw_author_supervised_label_metrics.csv');point.merge(result,on=['cohort','subset'],validate='one_to_one').to_csv(OUT/'raw_author_supervised_label_metrics_with_CI.csv',index=False)
print(result[['cohort','subset','profiles','donor_clusters','coverage_all_profiles','accuracy_CI_low','accuracy_CI_high']].to_string(index=False))
