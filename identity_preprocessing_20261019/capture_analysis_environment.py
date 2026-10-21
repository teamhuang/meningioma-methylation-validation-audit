from pathlib import Path
import importlib.metadata as md,json,platform,sys,hashlib,requests
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
packages={p:md.version(p) for p in ['numpy','pandas','scipy','scikit-learn','lifelines','openpyxl','requests']}
(OUT/'analysis_environment.json').write_text(json.dumps({'python':sys.version,'executable':sys.executable,'platform':platform.platform(),'packages':packages,
  'R':'4.3.3; complete sessions in author_batch_session_* and raw_qc_session_*','random_seeds':{'donor_model':2026101902,'label_cluster_bootstrap':2026101905}},indent=2),encoding='utf-8')
dois=['10.1038/s41588-022-01061-8','10.1126/sciadv.abm6247','10.1186/s40478-025-02184-9','10.1038/s41467-026-77170-3','10.1038/s41467-025-61063-y','10.1007/s00401-024-02708-y']
items=[]
for doi in dois:
    r=requests.get('https://api.crossref.org/works/'+doi,timeout=40);r.raise_for_status();z=r.json()['message']
    items.append({'DOI':doi,'title':z.get('title',[''])[0],'journal':z.get('container-title',[''])[0],'author':z.get('author',[]),'published':z.get('published'),
      'volume':z.get('volume'),'issue':z.get('issue'),'page':z.get('page'),'article_number':z.get('article-number'),'URL':z.get('URL')})
(OUT/'primary_reference_metadata.json').write_text(json.dumps(items,indent=2),encoding='utf-8');print([(z['DOI'],z['title']) for z in items])
