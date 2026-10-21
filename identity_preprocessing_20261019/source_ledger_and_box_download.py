from pathlib import Path
import json,re,requests
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
target=OUT/'baylor_box_file_891275442848.html'
html=target.read_text(encoding='utf-8')
start=html.index('Box.postStreamData = ')+len('Box.postStreamData = ')
obj,_=json.JSONDecoder().raw_decode(html[start:])
file=next(v for k,v in obj.items() if k.startswith('/app-api/enduserapp/item/'))
url=file.get('authenticated_download_url')
for u in [url,'https://bcm.box.com/index.php?rm=box_download_shared_file&shared_name=nsd4pwofiitxvp0zh1sp56way2z8mkti&file_id=891275442848']:
    if not u:continue
    r=requests.get(u,headers={'BoxApi':'shared_link=https://bcm.box.com/s/nsd4pwofiitxvp0zh1sp56way2z8mkti'},timeout=40)
    print('Public Box download',r.status_code,r.headers.get('content-type'),len(r.content))
    if r.ok and ('csv' in r.headers.get('content-type','') or r.text.startswith(',')):
        (OUT/'Baylor_Methylation_to_Sample_Name.csv').write_bytes(r.content);print('Baylor sample mapping downloaded');break
ledger=[
 {'issue':'UCSF per-case molecular labels','status':'resolved','evidence':'MOESM4 Supplementary Tables 3/4/5, literal M IDs; 200 IDs with a single consistent published group; GSM and raw barcode from GEO'},
 {'issue':'UCSF old clinical-row to fingerprint-column positional bridge','status':'rejected','evidence':'Independent GEO sex concordance 53.5%, grade 41.0%, native molecular group 26.5%; do not use for outcomes, MIB1 or biological gradient'},
 {'issue':'UCSF clinical outcome to individual methylation sample','status':'under_investigation','evidence':'173 unique age-grade-group candidates are covariate linkage only, not a verified source ID'},
 {'issue':'Baylor tumour to patient clusters','status':'raw_SNP_audit_underway','evidence':'Full raw-IDAT cohort requested for 59 SNP fingerprint assessment'},
 {'issue':'Author preprocessing pilot','status':'passed','evidence':'Original 403-probe reference matched exactly from raw IDAT using extracted author functions and minfi 1.43.1'},
 {'issue':'EPICv2 author preprocessing pilot','status':'passed','evidence':'Author minfi fork and manifest yield all 403 classifier probes; native standard minfi does not support this array correctly'}]
(OUT/'resolution_ledger.json').write_text(json.dumps(ledger,ensure_ascii=False,indent=2),encoding='utf-8')
