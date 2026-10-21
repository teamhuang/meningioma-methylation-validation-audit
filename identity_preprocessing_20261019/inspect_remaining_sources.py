from pathlib import Path
import pandas as pd,requests,json,gzip
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent;OUT=ROOT/'results'
book=BASE/'data/suppl/protocadherin_sourcedata.xlsx'
report=[]
for sheet in pd.ExcelFile(book).sheet_names:
    df=pd.read_excel(book,sheet_name=sheet,header=None,nrows=9)
    report.append({'sheet':sheet,'preview':df.iloc[:9,:20].fillna('').astype(str).values.tolist()})
(OUT/'protocadherin_source_inventory.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
for r in report:print(r['sheet'],repr(r['preview'])[:1400])
book=BASE/'data/suppl/natgen2022/41588_2022_1061_MOESM4_ESM.xlsx'
for sheet in ['3','4','5']:
    df=pd.read_excel(book,sheet_name=sheet,header=None)
    df.to_csv(OUT/('ucsf_native_identified_table_'+sheet+'.csv'),index=False,header=False)
    print('UCSF identified table',sheet,df.head(12).to_string(index=False,header=False))
for item in [('jbayley','meningiomas','Analysis _script_Bayley_et_al.R'),('william-c-chen','Meningioma_GE_Biomarker','Files/Combined_TPM_CSVs/combined.tpms.UCSF.csv')]:
    user,repo,path=item;url='https://raw.githubusercontent.com/'+user+'/'+repo+'/main/'+requests.utils.quote(path)
    r=requests.get(url,timeout=40);r.raise_for_status();target=OUT/(repo+'_'+Path(path).name);target.write_text(r.text,encoding='utf-8')
    print('Downloaded',target.name,len(r.content))
    if path.endswith('.csv'):print(r.text[:2000])
box='https://bcm.box.com/s/nsd4pwofiitxvp0zh1sp56way2z8mkti'
r=requests.get(box,timeout=40);(OUT/'baylor_public_box.html').write_text(r.text,encoding='utf-8')
print('BOX status',r.status_code,'bytes',len(r.content))
