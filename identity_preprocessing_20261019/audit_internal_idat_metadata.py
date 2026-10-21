from pathlib import Path
import struct,json,hashlib
import pandas as pd
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'

def read_string(f):
    length=0;shift=0
    while True:
        c=f.read(1)[0];length|=(c&127)<<shift
        if c<128:break
        shift+=7
        if shift>35:raise ValueError('invalid string length')
    return f.read(length).decode('utf-8')

def metadata(path):
    with path.open('rb') as f:
        assert f.read(4)==b'IDAT'
        version=struct.unpack('<Q',f.read(8))[0];n=struct.unpack('<I',f.read(4))[0]
        fields={}
        for _ in range(n):
            code,offset=struct.unpack('<HQ',f.read(10));fields[code]=offset
        out={'idat_version':version,'n_fields':n}
        for code,label in [(402,'internal_slide_barcode'),(403,'internal_chip_type')]:
            f.seek(fields[code]);out[label]=read_string(f)
        f.seek(fields[1000]);out['n_beads']=struct.unpack('<I',f.read(4))[0]
        return out

if __name__=='__main__':
    manifest=pd.read_csv(OUT/'raw_idat_manifest.csv');rows=[]
    for _,r in manifest.iterrows():
        for color in ['Grn','Red']:
            path=ROOT/'raw_idat'/r.cohort/(str(r.barcode)+'_'+color+'.idat')
            if not path.exists():continue
            z=metadata(path);z.update(cohort=r.cohort,GSM=r.GSM,color=color,barcode=r.barcode,file=str(path.relative_to(ROOT)),bytes=path.stat().st_size)
            z['internal_barcode_matches']=str(z['internal_slide_barcode'])==str(r.barcode).split('_')[0]
            z['sha256']=hashlib.sha256(path.read_bytes()).hexdigest();rows.append(z)
    d=pd.DataFrame(rows);d.to_csv(OUT/'raw_internal_barcode_and_hash_audit.csv',index=False)
    assert d.internal_barcode_matches.all()
    print(d.groupby('cohort').agg(channels=('GSM','size'),samples=('GSM','nunique'),barcode_matches=('internal_barcode_matches','sum')).to_string())
    print('Chip types',d.internal_chip_type.value_counts().to_dict())
