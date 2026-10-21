from pathlib import Path
import argparse,gzip,hashlib,json,time,threading
from concurrent.futures import ThreadPoolExecutor,as_completed
import requests,pandas as pd
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results';RAW=ROOT/'raw_idat';RAW.mkdir(exist_ok=True)
_local=threading.local()

def session():
    if not hasattr(_local,'http'):_local.http=requests.Session()
    return _local.http

def fetch(task):
    row,color=task;folder=RAW/row['cohort'];folder.mkdir(exist_ok=True)
    target=folder/(row['barcode']+'_'+color+'.idat')
    url=row[color+'_URL'];cache=target.with_suffix('.idat.gz')
    started=time.time()
    if not target.exists():
        last=None
        for attempt in range(6):
            try:
                if not cache.exists():
                    part=Path(str(cache)+'.part')
                    while not cache.exists():
                        offset=part.stat().st_size if part.exists() else 0
                        headers={'Accept-Encoding':'identity','Range':'bytes='+str(offset)+'-'+str(offset+2*1024*1024-1)}
                        r=session().get(url,stream=True,timeout=(20,30),headers=headers)
                        if r.status_code==416 and r.headers.get('Content-Range','')=='bytes */'+str(offset):part.replace(cache);break
                        r.raise_for_status()
                        if offset and r.status_code!=206:offset=0
                        total=None
                        if r.status_code==206:
                            content_range=r.headers.get('Content-Range','')
                            if not content_range.startswith('bytes '+str(offset)+'-'):raise ValueError('unexpected range')
                            total=int(content_range.split('/')[-1])
                        with open(part,'ab' if offset else 'wb') as f:
                            for chunk in r.iter_content(64*1024):f.write(chunk)
                        r.close()
                        if total is None or part.stat().st_size==total:part.replace(cache)
                        elif part.stat().st_size>total:raise ValueError('download exceeds source length')
                with gzip.open(cache,'rb') as src,open(str(target)+'.part','wb') as dest:
                    head=src.read(4)
                    if head!=b'IDAT':raise ValueError('not an IDAT: '+repr(head))
                    dest.write(head)
                    while True:
                        chunk=src.read(1024*1024)
                        if not chunk:break
                        dest.write(chunk)
                Path(str(target)+'.part').replace(target);last=None;break
            except Exception as e:
                last=e
                if cache.exists():cache.unlink()
                time.sleep(1+attempt*2)
        if last is not None:
            print('DOWNLOAD_ERROR',row['GSM'],color,type(last).__name__,str(last),flush=True)
            return {'cohort':row['cohort'],'GSM':row['GSM'],'barcode':row['barcode'],'color':color,'url':url,'status':'failed','error':str(last)}
    return {'cohort':row['cohort'],'GSM':row['GSM'],'barcode':row['barcode'],'color':color,'url':url,'status':'success',
      'file':str(target.relative_to(ROOT)),'bytes':target.stat().st_size,'sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'seconds':round(time.time()-started,2)}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cohort');p.add_argument('--pilot',action='store_true');p.add_argument('--workers',type=int,default=8);p.add_argument('--limit',type=int);args=p.parse_args()
    d=pd.read_csv(OUT/'raw_idat_manifest.csv').dropna(subset=['barcode','Grn_URL','Red_URL'])
    if args.pilot:
        dd=[]
        ref=pd.read_csv(ROOT.parent/'strengthening_20261019/results/official_model/official_knn_training.csv')
        for cohort in ['author186','baylor','tubingen','ucsf']:
            sub=d[d.cohort==cohort]
            if cohort=='author186':sub=sub[sub.barcode.isin(ref['sample'])]
            dd.append(sub.head(1))
        d=pd.concat(dd)
    elif args.cohort:
        d=d[d.cohort==args.cohort]
        if args.cohort=='tubingen':
            clinical=pd.read_csv(ROOT.parent/'results/gse304094_primary_survival.csv')
            ids=set(clinical.dropna(subset=['T','E','grade','gtr','age','male']).Internal_ID)
            d=d[d.title.isin(ids)]
        if args.cohort=='ucsf':d=d[d.institution=='UCSF']
    if args.limit:d=d.head(args.limit)
    d.to_csv(OUT/('requested_raw_'+('pilot' if args.pilot else args.cohort)+'.csv'),index=False)
    tasks=[(row,color) for row in d.to_dict('records') for color in ['Grn','Red']]
    rows=[]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(fetch,task) for task in tasks]
        for k,f in enumerate(as_completed(futures),1):
            row=f.result();rows.append(row)
            with open(OUT/('download_live_'+('pilot' if args.pilot else args.cohort)+'.jsonl'),'a',encoding='utf-8') as log:log.write(json.dumps(row)+'\n')
            if k%20==0 or row['status']!='success' or len(tasks)<=8:print(k,'/',len(tasks),row['cohort'],row['GSM'],row['color'],row['status'],row.get('bytes'),flush=True)
    file=OUT/('download_log_'+('pilot' if args.pilot else args.cohort)+'.csv')
    pd.DataFrame(rows).to_csv(file,index=False)
    print('DONE',len(rows),'failed',sum(x['status']!='success' for x in rows),flush=True)
