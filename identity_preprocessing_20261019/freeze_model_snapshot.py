from pathlib import Path
import json,hashlib
import pandas as pd
import analyse_donor_aware_models as analysis
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results/model_snapshot';OUT.mkdir(exist_ok=True)
specs={};manifest=[]
for cohort in ['baylor','tubingen','hku']:
    d,clinical,layers,horizon=analysis.prepare(cohort)
    path=OUT/('prepared_'+cohort+'.csv');d.to_csv(path,index=False)
    specs[cohort]={'clinical':clinical,'layers':layers,'horizon':horizon,'records_before_layer_complete_case_filter':len(d)}
    manifest.append({'file':path.name,'rows':len(d),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
(OUT/'specifications.json').write_text(json.dumps(specs,indent=2),encoding='utf-8')
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print(manifest)
