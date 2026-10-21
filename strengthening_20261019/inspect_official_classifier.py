import requests,json
from run_strengthening import OUT
url='https://api.github.com/repos/NWM-MGP/Public_Resources/contents/meningioma_risk_classifier?ref=develop'
r=requests.get(url,timeout=40);r.raise_for_status();obj=r.json()
(OUT/'official_classifier_repository_listing.json').write_text(json.dumps(obj,indent=2),encoding='utf-8')
for x in obj:
    print(x['name'],x['type'],x['size'],x.get('download_url'))
    if x['type']=='file' and x['size']<150000 and x['name'].lower().endswith(('.md','.r','.py','.txt','dockerfile')):
        s=requests.get(x['download_url'],timeout=40);s.raise_for_status()
        (OUT/('official_'+x['name'])).write_text(s.text,encoding='utf-8');print(s.text[:16000])
