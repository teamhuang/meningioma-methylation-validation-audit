import pandas as pd,numpy as np
from run_strengthening import OUT,save
summary=pd.read_csv(OUT/'validation_summary.csv')
metrics=pd.read_csv(OUT/'oof_risk_metrics_by_repeat.csv')
cols=['observed_ipcw','predicted_clinical','predicted_molecular','brier_clinical','brier_molecular','delta_brier','OP_molecular','min_G_horizon']
g=metrics.groupby(['comparison','horizon'])[cols].agg(['mean','min','max'])
g.columns=['_'.join(x) for x in g.columns];save(g.reset_index(),'oof_risk_summary.csv')
ph=[]
for file in OUT.glob('PH_diagnostics_*.csv'):
    if file.name=='PH_diagnostics_all.csv':continue
    d=pd.read_csv(file); d['source']=file.name;ph.append(d)
ph=pd.concat(ph,ignore_index=True)
from statsmodels.stats.multitest import multipletests
ph['P_BH_all_diagnostics']=multipletests(ph.p,method='fdr_bh')[1]
save(ph,'PH_diagnostics_all.csv')
sl=[]
for file in OUT.glob('oof_calibration_slopes_*.csv'):
    d=pd.read_csv(file);sl.append(d)
sl=pd.concat(sl,ignore_index=True)
save(sl.groupby(['comparison','model']).slope.agg(['mean','median','min','max']).reset_index(),'oof_slope_summary.csv')
print('OOF risk metrics (5-year horizon)')
gg=g.reset_index();print(gg[(gg.horizon==60)|(gg.horizon==5)][['comparison','OP_molecular_mean','brier_clinical_mean','brier_molecular_mean','delta_brier_mean']].to_string(index=False))
print('PH nominal P<0.05\n',ph[ph.p<.05][['comparison','index','p','P_BH_all_diagnostics']].to_string(index=False))
