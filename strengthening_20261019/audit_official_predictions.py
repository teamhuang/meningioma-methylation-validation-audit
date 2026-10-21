import json,gzip
from scipy.spatial.distance import cdist
import numpy as np,pandas as pd
from sklearn.metrics import adjusted_rand_score
from scipy.optimize import linear_sum_assignment
from run_strengthening import BASE,OUT,save
ref=pd.read_csv(OUT/'official_model/official_knn_training.csv',index_col=0)
rows=[]
def sample_map(accession):
    values={}
    with gzip.open(BASE/'data'/('GSE'+accession+'_series_matrix.txt.gz'),'rt',encoding='utf-8',errors='replace') as f:
        for line in f:
            if line.startswith('!Sample_geo_accession'):gs=[x.strip().strip('"') for x in line.rstrip().split('\t')[1:]]
            if line.startswith('!Sample_title'):tt=[x.strip().strip('"') for x in line.rstrip().split('\t')[1:]]
            if line.startswith('!series_matrix_table_begin'):break
    return gs,tt
for name in ['baylor','release186','release37','ucsf']:
    p=pd.read_csv(OUT/('official_knn_predictions_'+name+'.csv'))
    if name.startswith('release'):
        gs,tt=sample_map('292325' if name=='release186' else '292326')
        assert len(gs)==len(p)
        p['GSM_positional']=gs
        p['barcode_positional']=tt
        p['label_identity_status']='SAMPLE-number-to-series-order assumption; not an ID-proven label join'
        panel=pd.read_csv(OUT/('official_panel_'+name+'.csv.gz'),index_col=0).T
        xx=panel.loc[p.id,ref.columns[:-1]].to_numpy(float)
        rr=ref.loc[:,ref.columns[:-1]].to_numpy(float)
        dist=cdist(xx,rr,'sqeuclidean')
        closest=dist.argmin(axis=1)
        mad=np.array([np.max(np.abs(x-rr[j])) for x,j in zip(xx,closest)])
        p['reference_min_squared_distance']=dist.min(axis=1)
        p['reference_nearest_max_abs_difference']=mad
        p['reference_nearest_barcode']=ref.index.to_numpy()[closest]
        p['in_official_training']=mad<1e-8
        p['fingerprint_barcode_agrees_with_positional_title']=p.reference_nearest_barcode.eq(p.barcode_positional).where(p.in_official_training)
        save(p,'official_knn_predictions_'+name+'.csv')
    rows.append({'cohort':name,'n_scored':len(p),'n_ties':int(p.tied_top.sum()),'n_in_model_training':int(p.in_official_training.sum()),
       'n_confidence_ge_075':int(p.confidence_ge_075.sum()),'max_imputed_CpGs':int(p.imputed_features.max())})
    rows[-1]['training_overlap_check']='exact_403_probe_fingerprint' if name.startswith('release') else 'comparable_barcode_IDs' if name=='ucsf' else 'not_established_GSM_identifiers_not_comparable_to_training_barcodes'
    if name=='baylor':rows[-1]['n_in_model_training']=np.nan
    if name=='release186':
        save(p[['id','GSM_positional','barcode_positional','reference_nearest_barcode','reference_min_squared_distance','reference_nearest_max_abs_difference','in_official_training']],'official_training_release_identity.csv')
    if name=='release37':
        a=pd.read_csv(BASE/'results/gse292326_identity_audit.csv',dtype={'title':str})
        z=a.merge(p,left_on='GSM',right_on='GSM_positional',how='left',validate='one_to_one')
        z['official_group_numeric']=pd.to_numeric(z.predicted_class,errors='coerce')
        z['published_kmeans_numeric']=pd.to_numeric(z.kmeans,errors='coerce')
        z['held_out_from_official_training']=~z.in_official_training.fillna(True)
        save(z,'official_release37_label_audit.csv')
        use=z.dropna(subset=['official_group_numeric','published_kmeans_numeric'])
        rows[-1]['n_known_nontied_labels']=len(use)
        rows[-1]['ARI_vs_kmeans']=adjusted_rand_score(use.official_group_numeric,use.published_kmeans_numeric)
        # Fixed author mapping: risk output1=low corresponds kmeans2; output2=medium corresponds kmeans1.
        ct=pd.crosstab(use.official_group_numeric,use.published_kmeans_numeric)
        ri,ci=linear_sum_assignment(-ct.to_numpy())
        rows[-1]['bestmap_agreement_vs_kmeans']=ct.to_numpy()[ri,ci].sum()/len(use) if len(use) else np.nan
        held=use[use.held_out_from_official_training]
        rows[-1]['n_nontraining_known_nontied_labels']=len(held)
        rows[-1]['nontraining_ARI_vs_kmeans']=adjusted_rand_score(held.official_group_numeric,held.published_kmeans_numeric) if len(held)>1 else np.nan
        rows[-1]['label_linkage']='SAMPLE-number/series-order assumption; not fully verified individual identifiers'
        # Learn the cross-nomenclature mapping only on verified training fingerprints.
        inref=z[z.in_official_training.fillna(False)].dropna(subset=['published_kmeans_numeric'])
        inref['author_reference_group']=inref.reference_nearest_barcode.map(ref.group)
        trct=pd.crosstab(inref.author_reference_group,inref.published_kmeans_numeric)
        save(trct.reset_index(),'official_group_vs_kmeans_training_crosstab.csv')
        if len(trct):
            ri,ci=linear_sum_assignment(-trct.to_numpy());mapping={int(trct.index[i]):float(trct.columns[j]) for i,j in zip(ri,ci)}
            rows[-1]['training_only_label_map']=json.dumps(mapping)
            rows[-1]['nontraining_agreement_training_map']=float(held.official_group_numeric.map(mapping).eq(held.published_kmeans_numeric).mean())
save(pd.DataFrame(rows),'official_knn_audit_summary.csv')
print(pd.DataFrame(rows).to_string(index=False))
