from pathlib import Path
import zipfile,json,hashlib,importlib.metadata as md
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent;OUT=ROOT/'results'
checks=__import__('pandas').read_csv(OUT/'identity_validation_checks.csv')
assert checks.passed.all() and (ROOT/'数据身份与原始预处理验证报告_20261019.md').exists()
entries={}
def add(path,name=None):
    path=Path(path);entries[name or str(path.relative_to(ROOT)).replace('\\','/')]=path.read_bytes()
for f in ROOT.iterdir():
    if f.is_file() and f.suffix in ['.py','.R','.md']:add(f)
for f in (ROOT/'R_library').rglob('*'):
    if f.is_file():add(f)
for f in OUT.rglob('*'):
    if not f.is_file():continue
    if f.suffix.lower() not in ['.csv','.json','.rds','.log','.txt','.gz']:continue
    if f.name.endswith('_text.txt') or f.name.startswith(('baylor_public_box','baylor_box_file_')):continue
    add(f)
for name in ['knn_reference_fc2.rds','official_knn_training.csv','official_knn_cpgs.txt','run_meningioma_risk_analysis.R']:
    add(BASE/'strengthening_20261019/results/official_model'/name,'author_model/'+name)
for name in ['official_classifier_image_manifest.json','official_classifier_image_config.json','original_hashes.csv']:
    add(BASE/'strengthening_20261019/results'/name,'provenance_previous_release/'+name)
for name in ['41588_2022_1061_MOESM4_ESM.xlsx','41588_2022_1061_MOESM6_ESM.xlsx','41588_2022_1061_MOESM19_ESM.xlsx','41588_2022_1061_MOESM20_ESM.xlsx']:
    add(BASE/'data/suppl/natgen2022'/name,'source_data/natgen2022/'+name)
add(BASE/'data/suppl/protocadherin_suppldata.xlsx','source_data/protocadherin_suppldata.xlsx')
for name in ['401_2024_2708_MOESM1_ESM.xlsx','zenodo_11501567/MeningiomaClinical.xlsx']:
    add(ROOT/'source_data'/name)
# Only filesystem lookup changes in the release copy. Scientific functions and model parameters are untouched.
for name in ['run_author_batch.R','run_author_preprocessing.R']:
    if name in entries:
        source=entries[name].decode('utf-8')
        source=source.replace('author<-file.path(root,"..","strengthening_20261019/results/official_model")',
            'author<-file.path(root,"author_model");if(!dir.exists(author))author<-file.path(root,"..","strengthening_20261019/results/official_model")')
        entries[name]=source.encode('utf-8')
req='\n'.join(k+'=='+md.version(k) for k in ['numpy','pandas','scipy','scikit-learn','lifelines','openpyxl','requests','tabulate'])+'\n'
entries['requirements.txt']=req.encode()
entries['NOTICE.md']=('This package contains derived analyses and copies of public research reference data and R packages. Original package copyright and license files are retained. Minfi and its dependencies retain their original licenses. Author model and source-publication data must be cited using the sources in README.md and the provenance manifests. No large raw IDAT files are redistributed here. No clinical outcome identity is asserted for UCSF by row order. No figures were generated.\n').encode()
manifest=[{'file':k,'bytes':len(v),'sha256':hashlib.sha256(v).hexdigest()} for k,v in sorted(entries.items())]
entries['RELEASE_CONTENTS.json']=json.dumps(manifest,ensure_ascii=False,indent=2).encode('utf-8')
path=BASE/'npjDM_identity_preprocessing_validated_20261019.zip'
with zipfile.ZipFile(path,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for name,data in sorted(entries.items()):z.writestr('identity_preprocessing_20261019/'+name,data)
with zipfile.ZipFile(path) as z:
    assert z.testzip() is None
    for item in manifest:assert hashlib.sha256(z.read('identity_preprocessing_20261019/'+item['file'])).hexdigest()==item['sha256']
sha=hashlib.sha256(path.read_bytes()).hexdigest();path.with_suffix('.zip.sha256.txt').write_text(sha+'  '+path.name+'\n',encoding='utf-8')
result={'archive':str(path),'files':len(entries),'bytes':path.stat().st_size,'sha256':sha,'zip_CRC_pass':True,'all_member_hashes_pass':True}
(OUT/'final_release_validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(result)
