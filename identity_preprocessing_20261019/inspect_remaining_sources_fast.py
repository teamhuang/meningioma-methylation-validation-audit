from pathlib import Path
import zipfile,json,requests,shutil
from xml.etree import ElementTree as ET
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent;OUT=ROOT/'results'
NS={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
with zipfile.ZipFile(BASE/'data/suppl/protocadherin_sourcedata.xlsx') as z:
    wb=ET.fromstring(z.read('xl/workbook.xml'));rels=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
    links={x.attrib['Id']:x.attrib['Target'] for x in rels}
    shared=[]
    if 'xl/sharedStrings.xml' in z.namelist():
        for _,el in ET.iterparse(z.open('xl/sharedStrings.xml'),events=['end']):
            if el.tag.endswith('}si'):shared.append(''.join(el.itertext()));el.clear()
    report=[]
    for sh in wb.find('m:sheets',NS):
        ident=sh.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'];path=links[ident].lstrip('/')
        if not path.startswith('xl/'):path='xl/'+path
        rows=[]
        for _,el in ET.iterparse(z.open(path),events=['end']):
            if el.tag.endswith('}row'):
                row=[]
                for c in el:
                    v=c.find('m:v',NS);val=v.text if v is not None else ''.join(c.itertext())
                    if c.attrib.get('t')=='s' and v is not None:val=shared[int(val)]
                    row.append((c.attrib.get('r'),val))
                rows.append(row[:24]);el.clear()
                if len(rows)>=7:break
        report.append({'sheet':sh.attrib['name'],'preview':rows})
    (OUT/'protocadherin_source_inventory_fast.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    for r in report:print(r['sheet'],repr(r['preview'][:4])[:1500])
print('Disk free GB',shutil.disk_usage(BASE).free/1024**3)
for user,repo,path in [('jbayley','meningiomas','Analysis _script_Bayley_et_al.R'),('william-c-chen','Meningioma_GE_Biomarker','Files/Combined_TPM_CSVs/combined.tpms.UCSF.csv')]:
    u='https://raw.githubusercontent.com/'+user+'/'+repo+'/main/'+requests.utils.quote(path)
    r=requests.get(u,timeout=40);r.raise_for_status();f=OUT/(repo+'_'+Path(path).name);f.write_text(r.text,encoding='utf-8');print('Downloaded',f.name,len(r.content))
r=requests.get('https://bcm.box.com/s/nsd4pwofiitxvp0zh1sp56way2z8mkti',timeout=40)
(OUT/'baylor_public_box.html').write_text(r.text,encoding='utf-8');print('Box',r.status_code,len(r.content))
