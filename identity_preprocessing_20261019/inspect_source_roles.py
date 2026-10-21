from pathlib import Path
import re,html,json
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
raw=(OUT/'natgen2022_article.html').read_text(encoding='utf-8')
text=html.unescape(re.sub('<[^>]+>',' ',raw))
text=re.sub(r'\s+',' ',text)
(OUT/'natgen2022_article_text.txt').write_text(text,encoding='utf-8')
for term in ['nomogram','201 meningiomas','205','Supplementary Fig. 14','Supplementary Fig. 15','63 meningiomas','Toronto']:
    print('\nTERM',term)
    for match in list(re.finditer(re.escape(term),text,re.I))[:4]:print(text[max(0,match.start()-450):match.start()+950])
