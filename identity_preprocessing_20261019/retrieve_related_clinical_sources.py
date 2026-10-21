from pathlib import Path
import requests,re,json
from xml.etree import ElementTree as ET
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results';dest=ROOT/'source_data';dest.mkdir(exist_ok=True)
base='https://www.nature.com/articles/s41591-023-02586-z'
r=requests.get(base,timeout=40);(OUT/'natmed2023_page.html').write_text(r.text,encoding='utf-8')
urls=[requests.compat.urljoin(base,x.replace('&amp;','&')) for x in re.findall(r'href=["\']([^"\']+\.(?:xlsx|zip)(?:\?[^"\']*)?)["\']',r.text)]
print('Article',r.status_code,'supplement URLs',urls,flush=True)
search=requests.get('https://www.ebi.ac.uk/europepmc/webservices/rest/search',params={'query':'EXT_ID:37944590 AND SRC:MED','format':'json','resultType':'core'},timeout=40)
search.raise_for_status();z=search.json();(OUT/'natmed2023_metadata.json').write_text(json.dumps(z,indent=2),encoding='utf-8')
pmc=z['resultList']['result'][0].get('pmcid')
if pmc:
    rr=requests.get('https://www.ebi.ac.uk/europepmc/webservices/rest/'+pmc+'/fullTextXML',timeout=40)
    if rr.ok:
        (OUT/(pmc+'.xml')).write_text(rr.text,encoding='utf-8')
        tree=ET.fromstring(rr.content)
        for el in tree.iter('supplementary-material'):
            print('SUPPLEMENT',''.join(el.itertext())[:1000],el.attrib,flush=True)
            for child in el.iter():
                href=child.attrib.get('{http://www.w3.org/1999/xlink}href')
                if href and href.endswith(('.xlsx','.zip')):
                    urls.append(href if href.startswith('http') else 'https://pmc.ncbi.nlm.nih.gov/articles/'+pmc+'/bin/'+href)
rows=[]
for u in sorted(set(urls)):
    try:
        rr=requests.get(u,timeout=90);rr.raise_for_status()
        if not rr.content.startswith(b'PK'):raise ValueError('not an XLSX/ZIP payload')
        file=dest/u.split('/')[-1].split('?')[0];file.write_bytes(rr.content)
        rows.append({'url':u,'file':file.name,'bytes':len(rr.content),'status':'downloaded'})
    except Exception as exc:rows.append({'url':u,'status':'failed','error':str(exc)})
(OUT/'related_clinical_retrieval.json').write_text(json.dumps(rows,indent=2),encoding='utf-8');print(rows,flush=True)
