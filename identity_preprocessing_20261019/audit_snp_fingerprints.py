from pathlib import Path
import json
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results';RAW=OUT/'raw_preprocessed'
manifest=pd.read_csv(OUT/'raw_idat_manifest.csv')
summaries=[];vectors=[]
for file in RAW.glob('*_snp.csv'):
    gsm=file.name.split('_')[0]
    meta=manifest[manifest.GSM==gsm]
    if len(meta)!=1:continue
    s=pd.read_csv(file).set_index('snp').beta.rename(gsm)
    vectors.append(s);summaries.append(meta.iloc[0].to_dict())
if not vectors:raise SystemExit('No SNP data yet')
mat=pd.concat(vectors,axis=1).T.sort_index();meta=pd.DataFrame(summaries).set_index('GSM').reindex(mat.index)
call=pd.DataFrame(np.nan,index=mat.index,columns=mat.columns)
call[mat<=.2]=0;call[(mat>=.35)&(mat<=.65)]=1;call[mat>=.8]=2
rows=[]
for i in range(len(mat)):
    a=call.iloc[i].to_numpy();ba=mat.iloc[i].to_numpy()
    for j in range(i+1,len(mat)):
        b=call.iloc[j].to_numpy();bb=mat.iloc[j].to_numpy();common=np.isfinite(a)&np.isfinite(b)
        n=int(common.sum());dif=int((a[common]!=b[common]).sum());rate=dif/n if n else np.nan
        finite=np.isfinite(ba)&np.isfinite(bb)
        rows.append({'GSM1':mat.index[i],'GSM2':mat.index[j],'cohort1':meta.cohort.iloc[i],'cohort2':meta.cohort.iloc[j],
           'title1':meta.title.iloc[i],'title2':meta.title.iloc[j],'common_highconfidence_SNPs':n,'discordant_SNPs':dif,'discordance_fraction':rate,
           'SNP_beta_correlation':np.corrcoef(ba[finite],bb[finite])[0,1] if finite.sum()>3 else np.nan,
           'same_donor_candidate':n>=40 and dif<=1 and rate<=.02,
           'strict_candidate':n>=40 and dif==0,'relaxed_candidate':n>=35 and dif<=2 and rate<=.05})
pairs=pd.DataFrame(rows);pairs.to_csv(OUT/'SNP_pairwise_identity_metrics.csv',index=False)
mat.to_csv(OUT/'SNP_beta_all_processed.csv');call.to_csv(OUT/'SNP_highconfidence_genotypes.csv')
parent={g:g for g in mat.index}
def find(g):
    while parent[g]!=g:parent[g]=parent[parent[g]];g=parent[g]
    return g
for _,r in pairs[pairs.same_donor_candidate].iterrows():parent[find(r.GSM1)]=find(r.GSM2)
groups={g:find(g) for g in parent}
labels={g:'D'+str(k+1).zfill(3) for k,g in enumerate(sorted(set(groups.values())))}
meta['SNP_donor_cluster']=[labels[groups[g]] for g in meta.index]
meta['donor_primary_key']=[groups[g] for g in meta.index]
strict_parent={g:g for g in mat.index}
def strict_find(g):
    while strict_parent[g]!=g:strict_parent[g]=strict_parent[strict_parent[g]];g=strict_parent[g]
    return g
for _,r in pairs[pairs.strict_candidate].iterrows():
    a,b=strict_find(r.GSM1),strict_find(r.GSM2);strict_parent[max(a,b)]=min(a,b)
meta['donor_strict_key']=[strict_find(g) for g in meta.index]
parent2={g:g for g in mat.index}
def find2(g):
    while parent2[g]!=g:parent2[g]=parent2[parent2[g]];g=parent2[g]
    return g
for _,r in pairs[pairs.relaxed_candidate].iterrows():
    a,b=find2(r.GSM1),find2(r.GSM2);parent2[max(a,b)]=min(a,b)
meta['donor_overlap_guard_key']=[find2(g) for g in meta.index]
meta['callable_SNPs']=call.notna().sum(axis=1)
meta.index.name='GSM'
meta.to_csv(OUT/'SNP_donor_clusters_current.csv')
print('Processed SNP samples',len(meta),'by cohort',meta.cohort.value_counts().to_dict())
print('Candidate same donor pairs\n',pairs[pairs.same_donor_candidate].to_string(index=False))
