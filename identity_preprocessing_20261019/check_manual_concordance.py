from pathlib import Path
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results/donor_aware'
names=[c+'_'+m+'_overlap_guard' for c,ms in [('baylor',['published','centroid_surrogate','author_KNN_raw']),('tubingen',['published','centroid_surrogate','author_KNN_raw']),('hku',['published','author_KNN_raw'])] for m in ms]
rows=[]
def pair_counts(d,score):
    t=d['T'].to_numpy();e=d.E.to_numpy();s=np.asarray(score)
    i,j=np.triu_indices(len(d),1)
    early_i=((t[i]<t[j])&(e[i]==1))|((t[i]==t[j])&(e[i]==1)&(e[j]==0))
    early_j=((t[j]<t[i])&(e[j]==1))|((t[i]==t[j])&(e[j]==1)&(e[i]==0))
    valid=early_i|early_j;diff=np.where(early_i,s[i]-s[j],s[j]-s[i])[valid]
    return np.sum(diff>0)+.5*np.sum(diff==0),int(valid.sum())
for name in names:
    d=pd.read_csv(OUT/('input_'+name+'.csv'),float_precision='round_trip');p=pd.read_csv(OUT/('oof_'+name+'.csv'),float_precision='round_trip')
    cv=pd.read_csv(OUT/('cv_'+name+'.csv'),float_precision='round_trip').set_index('repeat')
    z=p.merge(d[['id','T','E']],on='id',validate='many_to_one')
    for rep,g in z.groupby('repeat'):
        counts=np.zeros(3)
        for fold,b in g.groupby('fold'):
            a,n=pair_counts(b,b.LP_clinical);aa,nn=pair_counts(b,b.LP_molecular);assert n==nn;counts+=np.array([a,aa,n])
        c0,c1=counts[:2]/counts[2];ref=cv.loc[rep]
        err=max(abs(c0-ref.C_clinical),abs(c1-ref.C_molecular),abs(c1-c0-ref.deltaC))
        rows.append({'model':name,'repeat':rep,'n_comparable_pairs':int(counts[2]),'C_manual_clinical':c0,'C_manual_molecular':c1,'max_absolute_error':err,'passed':err<=1e-12})
result=pd.DataFrame(rows);result.to_csv(ROOT/'results/manual_pairwise_concordance_verification.csv',index=False)
print('Independent pair-count checks',len(result),'max error',result.max_absolute_error.max());assert result.passed.all()
