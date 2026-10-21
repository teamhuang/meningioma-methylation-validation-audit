from pathlib import Path
import json,re
import pandas as pd
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent;OUT=ROOT/'results'
data=json.loads((OUT/'workbook_inventory.json').read_text(encoding='utf-8'))
for x in data:
    if any(a in x['file'] for a in ['MOESM28_','MOESM29_','MOESM30_','MOESM31_']):
        print(x['file'],x['sheet'],repr(x['preview'][:4])[:2200])
mp=pd.read_csv(OUT/'longitudinal_verified_patient_mapping_42.csv')
cl=pd.read_excel(BASE/'data/suppl/protocadherin_suppldata.xlsx',sheet_name='Supplementary Data 4',header=2)
print('Patient18 source annotations\n',cl[cl.Patient==18].to_string(index=False))
print('Patient18 source coordinates\n',mp[mp.source_patient_centroid==18][['Internal_ID','grade','group']].to_string(index=False))
print('True cluster counts\n',mp.groupby('source_patient_centroid').predicted_class.agg(['size','nunique']).to_string())
manifest=pd.read_csv(OUT/'raw_idat_manifest.csv')
print('GEO institutions\n',manifest[manifest.cohort=='ucsf'].institution.value_counts().to_dict())
