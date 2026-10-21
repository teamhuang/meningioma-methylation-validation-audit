from pypdf import PdfReader
from pathlib import Path
p=Path(r'E:\脑膜瘤npj\npjDM_submission_candidate_v6')
for i in range(1,5):
    page=PdfReader(str(p/f'Supplementary_Figure_S{i}.pdf')).pages[0]
    print(i,float(page.mediabox.width),float(page.mediabox.height))
