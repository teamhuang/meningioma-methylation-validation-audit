from pathlib import Path
import subprocess,time,sys,datetime
import pandas as pd
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent;cohort=sys.argv[1]
preprocess_only='--preprocess-only' in sys.argv
qc_partition=None
if '--qc-partition' in sys.argv:
    at=sys.argv.index('--qc-partition');qc_partition=tuple(map(int,sys.argv[at+1:at+3]))
d=pd.read_csv(ROOT/'results'/('requested_raw_'+cohort+'.csv'))
if qc_partition:d=d.iloc[qc_partition[0]-1::qc_partition[1]].copy()
R=Path('C:/Program Files/R/R-4.3.3/bin/Rscript.exe')
start=time.time();last=None
log=ROOT/'results'/('pipeline_'+cohort+('_qc_part'+str(qc_partition[0]) if qc_partition else '')+'.log')
while True:
    raw=[];done=[];qc=[]
    for r in d.itertuples():
        if all((ROOT/'raw_idat'/cohort/(r.barcode+'_'+c+'.idat')).exists() for c in ['Grn','Red']):raw.append(r.GSM)
        if (ROOT/'results/raw_preprocessed'/(r.GSM+'_summary.csv')).exists():done.append(r.GSM)
        if (ROOT/'results/raw_qc'/(r.GSM+'_qc.csv')).exists():qc.append(r.GSM)
    state=(len(raw),len(done),len(qc),len(d))
    if state!=last:print(datetime.datetime.now().isoformat(timespec='seconds'),cohort,'raw/preprocessed/QC/expected',state,flush=True);last=state
    if qc_partition and len(qc)==len(d):break
    if not qc_partition and len(done)==len(d) and (preprocess_only or len(qc)==len(d)):break
    if time.time()-start>3*3600:raise TimeoutError('cohort pipeline incomplete; inspect download and R logs')
    for script,needed in [('run_author_batch.R',set(raw)-set(done)),('run_raw_qc.R',set(raw)-set(qc))]:
        if preprocess_only and script=='run_raw_qc.R':continue
        if qc_partition and script=='run_author_batch.R':continue
        if not needed:continue
        cmd=['rtk','proxy',str(R),str(ROOT/script),str(ROOT),cohort]
        if qc_partition:cmd=['rtk','proxy',str(R),str(ROOT/'run_raw_qc_partitioned.R'),str(ROOT),cohort]+list(map(str,qc_partition))
        if script=='run_author_batch.R':cmd+=['8']
        with log.open('a',encoding='utf-8') as f:
            f.write('\nPHASE '+script+' '+datetime.datetime.now().isoformat()+'\n');f.flush()
            result=subprocess.run(cmd,cwd=BASE,stdout=f,stderr=subprocess.STDOUT)
        if result.returncode:raise RuntimeError(script+' failed; see '+str(log))
    time.sleep(10)
print('COMPLETE',cohort,'expected',len(d),flush=True)
