import json
import numpy as np,pandas as pd
import run_strengthening as m
name='ucsf_groups_beyond_mib1'
before=pd.read_csv(m.OUT/('summary_'+name+'.csv')).iloc[0]
m.USE_SNAPSHOT=True
after=m.run(name)
keys=['C_app_clinical','C_app_molecular','cv_C_clinical','cv_C_molecular','cv_delta_mean']
differences={key:abs(float(before[key])-float(after[key])) for key in keys}
assert max(differences.values())<1e-10
(m.OUT/'snapshot_execution_check.json').write_text(json.dumps({'comparison':name,'max_absolute_differences':differences,'tolerance':1e-10,'passed':True},indent=2),encoding='utf-8')
print('Prepared-input end-to-end rerun agrees within 1e-10.')
