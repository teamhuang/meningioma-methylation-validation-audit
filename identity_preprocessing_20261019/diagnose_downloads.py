from pathlib import Path
import requests,pandas as pd,json
ROOT=Path(__file__).resolve().parent
d=pd.read_csv(ROOT/'results/requested_raw_tubingen.csv')
for _,r in d.head(2).iterrows():
    u=r.Grn_URL
    try:
        q=requests.get(u,stream=True,timeout=30)
        print(r.GSM,u,q.status_code,q.headers.get('Content-Type'),q.raw.read(200) if q.status_code!=200 else 'OK')
        q.close()
    except Exception as exc:print(type(exc).__name__,str(exc))
box=(ROOT/'results/baylor_public_box.html').read_text(encoding='utf-8')
for text in ['typedID','itemID','Clinical','clinical','download','Patient','patient']:
    at=box.find(text)
    print('BOX',text,at,box[max(0,at-100):at+600] if at>=0 else '')
p=ROOT/'results/Meningioma_GE_Biomarker_combined.tpms.UCSF.csv'
df=pd.read_csv(p)
print('UCSF GE matrix',df.shape,df.columns.tolist()[:30]);print(df.iloc[:5,:12].to_string(index=False))
