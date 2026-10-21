from pathlib import Path,PurePosixPath
import json,hashlib,requests,tarfile,io
BASE=Path(__file__).resolve().parent.parent
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
lib=ROOT/'R_library';lib.mkdir(exist_ok=True)
manifest=json.loads((BASE/'strengthening_20261019/results/official_classifier_image_manifest.json').read_text())
rows=[]
for ix in [22,23]:
    layer=manifest['layers'][ix];digest=layer['digest']
    r=requests.get('https://quay.io/v2/bioinformaticsnm/public_resources/meningioma-risk-classifier/blobs/'+digest,timeout=120)
    r.raise_for_status();assert hashlib.sha256(r.content).hexdigest()==digest.split(':')[1]
    with tarfile.open(fileobj=io.BytesIO(r.content),mode='r:gz') as tf:
        for member in tf:
            if not member.isfile():continue
            pos=PurePosixPath(member.name)
            parts=pos.parts
            package='IlluminaHumanMethylationEPICv2manifest' if ix==22 else 'minfi'
            if package not in parts:continue
            start=parts.index(package);tail=parts[start:]
            if '..' in tail:raise ValueError('invalid archive path')
            target=lib.joinpath(*tail);target.parent.mkdir(parents=True,exist_ok=True)
            payload=tf.extractfile(member).read();target.write_bytes(payload)
            rows.append({'file':str(target.relative_to(ROOT)),'bytes':len(payload),'sha256':hashlib.sha256(payload).hexdigest(),'source_layer':digest})
    print(package,'copied from author image',flush=True)
(OUT/'author_R_package_extraction.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
print('Native libraries:',[(str(p.relative_to(ROOT)),p.stat().st_size) for p in lib.rglob('*') if p.suffix in ['.so','.dll']])
