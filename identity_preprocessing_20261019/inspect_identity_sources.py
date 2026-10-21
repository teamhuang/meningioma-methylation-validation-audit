from pathlib import Path
import gzip,json,re
import pandas as pd
BASE=Path(__file__).resolve().parent.parent
OUT=Path(__file__).resolve().parent/'results';OUT.mkdir(parents=True,exist_ok=True)
report=[]
for file in sorted((BASE/'data/suppl/natgen2022').glob('*.xlsx')):
    sheets=pd.ExcelFile(file).sheet_names
    for sheet in sheets:
        d=pd.read_excel(file,sheet_name=sheet,header=None,nrows=8)
        report.append({'file':file.name,'sheet':sheet,'preview':d.iloc[:8,:14].fillna('').astype(str).values.tolist()})
for file in [BASE/'data/GSE189672_Clinical_data.xlsx',BASE/'data/GSE189521_Clinical_data.xlsx',BASE/'data/suppl/protocadherin_suppldata.xlsx']:
    for sheet in pd.ExcelFile(file).sheet_names:
        d=pd.read_excel(file,sheet_name=sheet,header=None,nrows=10)
        report.append({'file':file.name,'sheet':sheet,'preview':d.iloc[:10,:18].fillna('').astype(str).values.tolist()})
(OUT/'workbook_inventory.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
for file in ['GSE183647_series_matrix.txt.gz','GSE304096_series_matrix.txt.gz','GSE189672_series_matrix.txt.gz','GSE189521_series_matrix.txt.gz']:
    rows=[]
    with gzip.open(BASE/'data'/file,'rt',encoding='utf-8',errors='replace') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):break
            if line.startswith(('!Sample_','!Series_')):
                vals=[v.strip().strip('"') for v in line.rstrip().split('\t')]
                if vals[0] in ['!Sample_title','!Sample_geo_accession','!Sample_characteristics_ch1','!Sample_source_name_ch1','!Sample_supplementary_file','!Series_supplementary_file','!Series_relation']:
                    rows.append({'field':vals[0],'first_values':vals[1:14]})
    (OUT/(file+'.header.json')).write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
    print(file,json.dumps(rows,ensure_ascii=False)[:14000])
for x in report:
    if 'protocadherin_suppldata' in x['file']:
        print(x['file'],x['sheet'],json.dumps(x['preview'],ensure_ascii=False)[:3500])
