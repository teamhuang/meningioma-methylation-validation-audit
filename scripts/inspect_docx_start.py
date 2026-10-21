from docx import Document
from pathlib import Path
p=Path(r'E:\脑膜瘤npj\npjDM_submission_candidate_v6\manuscript_npjDM_submission_candidate_v6.docx')
d=Document(p)
for i,x in enumerate(d.paragraphs[:12]):
    print(i,repr(x.text))
    if x.runs: print(' runs:',[(r.text,r.italic,r.bold) for r in x.runs])
