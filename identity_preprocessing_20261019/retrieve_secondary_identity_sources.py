from pathlib import Path
import requests,re,json
from xml.etree import ElementTree as ET
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results';DEST=ROOT/'source_data'
sources=[('pancancer2025','https://www.nature.com/articles/s41467-025-61063-y','PMC12222647'),('fusions2024','https://link.springer.com/article/10.1007/s00401-024-02708-y','PMC10954959')]
ledger=[]
for name,url,pmc in sources:
    r=requests.get(url,timeout=45);(OUT/(name+'_page.html')).write_text(r.text,encoding='utf-8')
    urls=[requests.compat.urljoin(url,x.replace('&amp;','&')) for x in re.findall(r'href=["\']([^"\']+\.(?:xlsx|csv|zip|rds)(?:\?[^"\']*)?)["\']',r.text,re.I)]
    print(name,'status',r.status_code,'github',sorted(set(re.findall(r'https://github.com/[^\s<>"\']+',r.text))),'downloads',urls,flush=True)
    rr=requests.get('https://www.ebi.ac.uk/europepmc/webservices/rest/'+pmc+'/fullTextXML',timeout=45)
    if rr.ok:
        (OUT/(pmc+'.xml')).write_text(rr.text,encoding='utf-8')
        tree=ET.fromstring(rr.content)
        for el in tree.iter('supplementary-material'):
            print(name,'SUPPLEMENT',''.join(el.itertext())[:350],flush=True)
            for child in el.iter():
                href=child.attrib.get('{http://www.w3.org/1999/xlink}href')
                if href and href.endswith(('.xlsx','.csv','.zip','.rds')):urls.append(href if href.startswith('http') else 'https://pmc.ncbi.nlm.nih.gov/articles/'+pmc+'/bin/'+href)
        for el in tree.iter():
            href=el.attrib.get('{http://www.w3.org/1999/xlink}href','')
            if 'github.com' in href or 'zenodo' in href:print('CODE LINK',href,flush=True)
    for u in sorted(set(urls)):
        file=DEST/u.split('/')[-1].split('?')[0]
        try:
            if not file.exists():
                res=requests.get(u,timeout=90);res.raise_for_status()
                if file.suffix in ['.xlsx','.zip'] and not res.content.startswith(b'PK'):raise ValueError('not an XLSX/ZIP')
                file.write_bytes(res.content)
            ledger.append({'source':name,'url':u,'file':file.name,'bytes':file.stat().st_size,'status':'downloaded'})
        except Exception as exc:ledger.append({'source':name,'url':u,'status':'failed','error':str(exc)})
(OUT/'secondary_identity_source_retrieval.json').write_text(json.dumps(ledger,indent=2),encoding='utf-8')
