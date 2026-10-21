from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results';rows=[]
snp=pd.read_csv(OUT/'SNP_donor_clusters_current.csv').set_index('GSM')
for f in (OUT/'donor_aware').glob('input_*.csv'):
    d=pd.read_csv(f);rule='donor_strict_key' if '_strict' in f.stem else 'donor_overlap_guard_key'
    expected=d.GSM.map(snp[rule]);ok=d.donor.eq(expected).all()
    rows.append({'model':f.stem[6:],'records':len(d),'all_final_donor_keys_unchanged':ok,'n_changed_keys':int(d.donor.ne(expected).sum())})
result=pd.DataFrame(rows);result.to_csv(OUT/'final_model_donor_key_consistency.csv',index=False)
print(result.to_string(index=False));assert result.all_final_donor_keys_unchanged.all()
