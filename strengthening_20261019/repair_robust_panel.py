"""Correct ascending rank selection. Report alongside, never replace legacy P3 silently."""
import json,time,ast,re,gzip
from concurrent.futures import ProcessPoolExecutor
import numpy as np,pandas as pd
from sklearn.metrics import adjusted_rand_score
from scipy.optimize import linear_sum_assignment
from run_strengthening import BASE,OUT,save
import extract_official_panels as extractor

def choose():
    train=pd.read_csv(BASE/'results/cache_m325_sel.csv.gz',index_col=0).astype(np.float32).T
    st=pd.DataFrame({'MAD':(train-train.median()).abs().median(),'missing_fraction':train.isna().mean()})
    st['rank_mad']=st.MAD.rank(ascending=False,method='first')
    st['rank_missing']=st.missing_fraction.rank(ascending=True,method='first')
    st['score']=st.rank_mad+st.rank_missing
    st=st.sort_values('score',ascending=True)
    st['selected_corrected']=np.arange(len(st))<2000
    old=set(pd.read_csv(BASE/'results/methyl_robust_cpgs.csv').cpg)
    st['selected_legacy']=st.index.isin(old)
    save(st.reset_index(names='cpg'),'robust_panel_rank_audit.csv')
    panel=st.head(2000).index.tolist()
    (OUT/'corrected_robust_cpgs.txt').write_text('\n'.join(panel),encoding='utf-8')
    print('Corrected / legacy panel overlap:',len(old.intersection(panel)),flush=True)
    return panel

def extract_one(name):
    extractor.targets=set((OUT/'corrected_robust_cpgs.txt').read_text().split())
    spec={'baylor_correctedP3':('GSE189521_series_matrix.txt.gz','\t',True),
          'ucsf_correctedP3':('GSE183647_meningioma_betas.csv.gz',',',False),
          'multicentre_correctedP3':('GSE270371_Matrix_processed.txt.gz','\t',False)}
    extractor.extract(name,*spec[name])

def stats(pred,truth):
    ct=pd.crosstab(pd.Series(pred),pd.Series(truth));ri,ci=linear_sum_assignment(-ct.to_numpy())
    return adjusted_rand_score(truth,pred),ct.to_numpy()[ri,ci].sum()/len(pred)

def evaluate():
    panel=(OUT/'corrected_robust_cpgs.txt').read_text().split()
    train=pd.read_csv(BASE/'results/cache_m325_sel.csv.gz',index_col=0).reindex(panel).T
    scale=pd.read_csv(BASE/'results/methyl_scaler.csv',index_col=0).reindex(panel)
    cent=pd.read_csv(BASE/'results/methyl_centroids.csv',index_col=0).reindex(panel)
    def predict(beta,cc):
        z=(beta.loc[:,cc]-scale.loc[cc,'mean']).div(scale.loc[cc,'sd']).to_numpy()
        finite=np.isfinite(z).sum(axis=1);z=np.nan_to_num(z,nan=0)
        dist=((z[:,None,:]-cent.loc[cc].to_numpy().T[None,:,:])**2).sum(axis=2)
        pred=np.asarray(cent.columns)[dist.argmin(axis=1)]
        return pred,dist.min(axis=1),finite
    rows=[]
    for name in ['baylor_correctedP3','ucsf_correctedP3','multicentre_correctedP3']:
        bet=pd.read_csv(OUT/('official_panel_'+name+'.csv.gz'),index_col=0)
        # GSE270371 contains detection-p columns as well as beta columns.
        if name.startswith('multicentre'):bet=bet.loc[:,~bet.columns.str.startswith('Detection Pval')]
        present=bet.notna().any(axis=1);cc=bet.index[present].tolist()
        ptr,dtr,ftr=predict(train,cc);threshold=float(np.percentile(dtr,99))
        pr,dist,finite=predict(bet.T,cc)
        eligible=finite>=.5*len(cc)
        data=pd.DataFrame({'id':bet.columns,'predicted_class':pr,'distance':dist,'finite_CpGs':finite,
            'evaluable':eligible,'rejected':dist>threshold})
        row={'cohort':name,'profiled':len(data),'shared_cpgs':len(cc),'n_evaluable':int(eligible.sum()),
             'threshold_P99_training':threshold,'n_rejected':int((dist>threshold).sum()),'reject_percent':100*np.mean(dist>threshold),
             'class_counts':json.dumps(pd.Series(pr).value_counts().to_dict())}
        if name.startswith('baylor'):
            lk=pd.read_csv(BASE/'results/gse189521_linked_endpoints.csv',dtype={'GSM':str})
            truth=lk.set_index('GSM').menG.reindex(data.id)
        elif name.startswith('ucsf'):
            mapping=pd.read_csv(BASE/'results/ucsftoronto_sample_mapping.csv').dropna(subset=['MID'])
            bm=dict(zip(mapping.raw_barcode,mapping.MID));data['MID']=data.id.map(bm)
            ref=pd.read_csv(BASE/'results/harmonization_pred_ucsfsesame.csv').set_index('MID').km_ref
            truth=ref.reindex(data.MID)
        else:truth=pd.Series(np.nan,index=data.index)
        truth=np.asarray(truth);known=pd.notna(truth)&eligible
        if known.sum()>1:
            row['reference_n']=int(known.sum());row['ARI'],row['best_map_agreement']=stats(pr[known],truth[known])
            keep=known&~(dist>threshold)
            row['retained_ARI'],row['retained_best_map_agreement']=stats(pr[keep],truth[keep])
        save(data,'corrected_P3_predictions_'+name+'.csv');rows.append(row)
        if name.startswith('ucsf'):
            ss=data[data.MID.notna()].copy()
            row2=dict(row,cohort='ucsf_discovery200_correctedP3',profiled=len(ss),n_evaluable=int(ss.evaluable.sum()),
                n_rejected=int(ss.rejected.sum()),reject_percent=100*ss.rejected.mean(),
                class_counts=json.dumps(ss.predicted_class.value_counts().to_dict()))
            rows.append(row2)
    save(pd.DataFrame(rows),'corrected_P3_transfer_summary.csv');print(pd.DataFrame(rows).to_string(index=False),flush=True)

if __name__=='__main__':
    import sys
    if '--evaluate-only' in sys.argv:
        evaluate();sys.exit(0)
    choose()
    with ProcessPoolExecutor(max_workers=3) as pool:list(pool.map(extract_one,['baylor_correctedP3','ucsf_correctedP3','multicentre_correctedP3']))
    evaluate()
