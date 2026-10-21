from pathlib import Path
import shutil,subprocess,json
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parent;TEST=ROOT/'snapshot_reproduction_check';TEST.mkdir(exist_ok=True)
for name in ['analyse_donor_aware_models.py','donor_model_utilities.py']:shutil.copy2(ROOT/name,TEST/name)
(TEST/'results/model_snapshot').mkdir(parents=True,exist_ok=True)
for f in (ROOT/'results/model_snapshot').iterdir():shutil.copy2(f,TEST/'results/model_snapshot'/f.name)
cmd=['rtk','proxy','python',str(TEST/'analyse_donor_aware_models.py'),'--cohort','tubingen','--layer','author_KNN_raw','--common-complete-set']
subprocess.run(cmd,cwd=TEST,check=True)
name='summary_tubingen_author_KNN_raw_overlap_guard_common_set.csv'
old=pd.read_csv(ROOT/'results/donor_aware'/name);new=pd.read_csv(TEST/'results/donor_aware'/name)
cols=['records','events','donor_clusters','C_app_clinical','C_app_molecular','delta_app','C_cv_clinical','C_cv_molecular','delta_cv','positive_repeats']
delta=np.abs(old[cols].to_numpy(float)-new[cols].to_numpy(float));assert delta.max()<=1e-12
result={'portable_input_snapshot_pass':True,'old_workspace_inputs_available_in_test':False,'max_absolute_metric_difference':float(delta.max()),'comparison':name,
        'tested_metrics':cols,'note':'Executed from a separate minimal folder with only the current analysis scripts and frozen prepared inputs; older workspace clinical sources were not copied.'}
(ROOT/'results/portable_model_snapshot_verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(result)
