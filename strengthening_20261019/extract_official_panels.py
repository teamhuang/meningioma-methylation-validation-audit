"""Stream only the author's 403 CpGs from local processed beta matrices."""
from pathlib import Path
import gzip,io,time,argparse,re
import pandas as pd
from run_strengthening import BASE,OUT,save

targets=set((OUT/'official_model/official_knn_cpgs.txt').read_text().split())
def extract(name,file,sep,series=False):
    start=time.time(); rows=[]; header=None; scanned=0
    with gzip.open(BASE/'data'/file,'rt',encoding='utf-8',errors='replace') as f:
        if series:
            for line in f:
                if line.startswith('!series_matrix_table_begin'):break
        header=f.readline()
        for line in f:
            if line.startswith('!'):continue
            key=line.split(sep,1)[0].strip().strip('"').split('_')[0]
            if key in targets:rows.append(line)
            scanned+=1
    tab=pd.read_csv(io.StringIO(header+''.join(rows)),sep=sep,index_col=0)
    tab.index=[str(x).strip().split('_')[0] for x in tab.index]
    tab.columns=[str(x).strip() for x in tab.columns]
    # EPICv2 duplicate probes are averaged once after collapsing to base CpG identifiers.
    if tab.index.duplicated().any():tab=tab.groupby(level=0,sort=False).mean()
    tab=tab.reindex(sorted(targets))
    tab.to_csv(OUT/('official_panel_'+name+'.csv.gz'))
    record={'cohort':name,'source_file':file,'n_samples':tab.shape[1],'n_official_features':len(targets),
      'n_features_with_any_data':int(tab.notna().any(axis=1).sum()),'min_features_per_sample':int(tab.notna().sum().min()),
      'max_features_per_sample':int(tab.notna().sum().max()),'n_complete_samples':int(tab.notna().all(axis=0).sum()),'scanned_rows':scanned,
      'seconds':time.time()-start}
    save(pd.DataFrame([record]),'official_feature_coverage_'+name+'.csv')
    print(record,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--cohort',required=True);args=parser.parse_args()
    files={'baylor':('GSE189521_series_matrix.txt.gz','\t',True),
        'tubingen':('GSE304094_discovery_processed.txt.gz','\t',False),
        'release186':('GSE292325_matrix_processed.txt.gz','\t',False),
        'release37':('GSE292326_matrix_processed.txt.gz','\t',False),
        'ucsf':('GSE183647_meningioma_betas.csv.gz',',',False)}
    extract(args.cohort,*files[args.cohort])
