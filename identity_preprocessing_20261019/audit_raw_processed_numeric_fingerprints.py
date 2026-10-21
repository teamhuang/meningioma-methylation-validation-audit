from pathlib import Path
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results';OLD=ROOT.parent/'strengthening_20261019/results'
meta=pd.read_csv(OUT/'raw_idat_manifest.csv').set_index('GSM');rows=[]
for cohort in ['baylor','ucsf']:
    processed=pd.read_csv(OLD/('official_panel_'+cohort+'.csv.gz'),index_col=0)
    p=processed.T;cpgs=processed.index
    for file in (OUT/'raw_preprocessed').glob('*_403beta.csv'):
        gsm=file.name.split('_')[0]
        if meta.loc[gsm,'cohort']!=cohort:continue
        key=gsm if cohort=='baylor' else str(meta.loc[gsm,'barcode'])
        if key not in p.index:raise ValueError('numeric source column missing '+gsm)
        x=pd.read_csv(file).set_index('cpg').beta.reindex(cpgs).to_numpy(float)
        y=p.loc[key].to_numpy(float);mask=np.isfinite(x)&np.isfinite(y)
        # A source-profile check only: normalization differences prevent an equality requirement.
        rows.append({'cohort':cohort,'GSM':gsm,'source_profile_identifier':key,'n_shared_finite_CpGs':int(mask.sum()),
            'raw_vs_published_beta_Pearson':float(np.corrcoef(x[mask],y[mask])[0,1]),'mean_absolute_beta_difference':float(np.mean(abs(x[mask]-y[mask]))),
            'RMSE':float(np.sqrt(np.mean((x[mask]-y[mask])**2))),
            'interpretation':'source-barcode-aligned numeric concordance; not a clinical-outcome crosswalk or an exact-preprocessing equality check'})
d=pd.DataFrame(rows);d.to_csv(OUT/'raw_vs_published_processed_beta_fingerprint_audit.csv',index=False)
print(d.groupby('cohort').agg(n=('GSM','size'),min_shared=('n_shared_finite_CpGs','min'),median_Pearson=('raw_vs_published_beta_Pearson','median'),min_Pearson=('raw_vs_published_beta_Pearson','min'),max_RMSE=('RMSE','max')).to_string())
