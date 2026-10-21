from pathlib import Path
import urllib.request,re
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent;OUT=ROOT/'results'
url='https://open-proposals.ucsf.edu/sites/default/files/chv/cw-projects-2025/1709303946/choudhurynatgenet2022.pdf'
file=ROOT/'source_data/choudhurynatgenet2022.pdf'
if not file.exists():
    print('Full paper URL unavailable; reading retained supplementary PDF.')
for pdf in [p for p in [file,BASE/'data/suppl/natgen2022/41588_2022_1061_MOESM1_ESM.pdf'] if p.exists()]:
    rd=PdfReader(pdf);pages=[]
    for i,p in enumerate(rd.pages):
        text=p.extract_text() or '';pages.append('PAGE '+str(i+1)+'\n'+text)
        if any(w in text.lower() for w in ['nomogram','201 meningioma','205 meningioma','supplementary figure 14','supplementary figure 15']):
            for term in ['nomogram','201','205','Supplementary Figure 14','Supplementary Figure 15']:
                for mat in list(re.finditer(re.escape(term),text,re.I))[:3]:print(pdf.name,'page',i+1,text[max(0,mat.start()-250):mat.start()+900].replace('\n',' '))
    (OUT/(pdf.stem+'_text.txt')).write_text('\n\n'.join(pages),encoding='utf-8')
