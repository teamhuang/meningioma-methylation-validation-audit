from pathlib import Path
import pandas as pd,json
from recover_longitudinal_centroid_links import sheet
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent
for name in ['41588_2022_1061_MOESM19_ESM.xlsx','41588_2022_1061_MOESM27_ESM.xlsx']:
    d=pd.read_excel(BASE/'data/suppl/natgen2022'/name)
    print(name,'shape',d.shape,'columns',d.columns.tolist());print(d.head(3).to_string(index=False))
for name in ['41591_2023_2586_MOESM3_ESM.xlsx','41591_2023_2586_MOESM4_ESM.xlsx']:
    path=ROOT/'source_data'/name;book=pd.ExcelFile(path)
    for sh in book.sheet_names:
        if any(k in sh.lower() for k in ['4f','4g','3c','3d','2a']):
            d=sheet(path,sh);print(name,sh,d.shape);print(d.iloc[:7,:15].to_string(index=False,header=False))
with __import__('gzip').open(BASE/'data/GSE183647_meningioma_betas.csv.gz','rt') as f:print('FULL BETA HEADER',f.readline()[:1600])
