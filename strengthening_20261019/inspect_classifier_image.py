import requests,json
from run_strengthening import OUT
repo='bioinformaticsnm/public_resources/meningioma-risk-classifier'
s=requests.Session()
u='https://quay.io/v2/'+repo+'/manifests/v0.1.0'
h={'Accept':'application/vnd.docker.distribution.manifest.v2+json, application/vnd.oci.image.manifest.v1+json, application/vnd.docker.distribution.manifest.list.v2+json'}
r=s.get(u,headers=h,timeout=30)
if r.status_code==401:
    auth=r.headers.get('WWW-Authenticate','')
    import re
    parts=dict(re.findall(r'(\w+)="([^"]+)"',auth))
    tr=s.get(parts.pop('realm'),params=parts,timeout=30);tr.raise_for_status()
    h['Authorization']='Bearer '+tr.json()['token']
    r=s.get(u,headers=h,timeout=30)
print('manifest',r.status_code,r.text[:4000])
if r.ok:
    (OUT/'official_classifier_image_manifest.json').write_text(r.text,encoding='utf-8')
    obj=r.json()
    cfg=s.get('https://quay.io/v2/'+repo+'/blobs/'+obj['config']['digest'],headers=h,timeout=30)
    cfg.raise_for_status();(OUT/'official_classifier_image_config.json').write_text(cfg.text,encoding='utf-8')
    hist=cfg.json()['history']; k=0
    for row in hist:
        if row.get('empty_layer'):continue
        layer=obj['layers'][k];k+=1
        print('LAYER',k,layer['size'],row.get('created_by',''))
