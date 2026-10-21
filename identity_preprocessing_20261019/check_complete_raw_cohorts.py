from pathlib import Path
import pandas as pd,json
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results';rows=[]
for cohort,n in [('author186',3),('baylor',110),('tubingen',149),('ucsf',200),('hku',365),('longitudinal',42)]:
    d=pd.read_csv(OUT/('requested_raw_'+cohort+'.csv'));assert len(d)==n
    for _,r in d.iterrows():
        prefix=ROOT/'raw_idat'/cohort/str(r.barcode)
        raw=all(Path(str(prefix)+'_'+color+'.idat').exists() for color in ['Grn','Red'])
        summary=OUT/'raw_preprocessed'/(r.GSM+'_summary.csv');qc=OUT/'raw_qc'/(r.GSM+'_qc.csv')
        rows.append({'cohort':cohort,'GSM':r.GSM,'barcode':r.barcode,'raw_pair_complete':raw,'author_preprocessed':summary.exists(),'raw_qc':qc.exists()})
z=pd.DataFrame(rows);z.to_csv(OUT/'raw_cohort_completion_status.csv',index=False)
print(z.groupby('cohort').agg(expected=('GSM','size'),raw_pair=('raw_pair_complete','sum'),author_preprocessed=('author_preprocessed','sum'),raw_qc=('raw_qc','sum')).to_string())
