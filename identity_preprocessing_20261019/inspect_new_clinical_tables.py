from pathlib import Path
import pandas as pd,requests,json,re
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results';BASE=ROOT.parent
for file in (ROOT/'source_data').glob('*.xlsx'):
    book=pd.ExcelFile(file)
    for sh in book.sheet_names:
        d=pd.read_excel(book,sheet_name=sh,header=None,nrows=8)
        print(file.name,sh,d.shape,repr(d.iloc[:6,:18].fillna('').astype(str).values.tolist())[:4000])
r=requests.get('https://www.nature.com/articles/s41588-022-01061-8',timeout=40)
(OUT/'natgen2022_article.html').write_text(r.text,encoding='utf-8')
for x in ['MOESM19','MOESM20','MOESM3','Supplementary Data 1','Supplementary Data 2']:
    at=r.text.find(x);print('SOURCE ROLE',x,re.sub('<[^>]+>',' ',r.text[max(0,at-250):at+350]) if at>=0 else 'NOT FOUND')
