from pathlib import Path
import json, requests, gzip, re
import pandas as pd
from run_strengthening import BASE,OUT,save

ng=BASE/'data/suppl/natgen2022'
a=pd.read_excel(ng/'41588_2022_1061_MOESM19_ESM.xlsx')
b=pd.read_excel(ng/'41588_2022_1061_MOESM20_ESM.xlsx')
c=pd.read_csv(BASE/'results/clinical_parsed.csv')
t=pd.read_csv(BASE/'results/gse304094_primary_survival.csv')
report=[]
report.append('UCSF source columns: '+repr(a.columns.tolist()))
report.append('UCSF source endpoint/sex/resection preview:\n'+a[['LFFR','LF','Sex','EOR','Group']].head(8).to_string(index=False))
for name,df,cols in [('Baylor',c,['Gender','ResectionType','Simpson_grade','Recurrence']),('UCSF',a,['Sex','EOR','LF']),('Tubingen',t,['Sex','Molecular_Group (by Toronto nomenclature)','tumor_setting','prior_RT'])]:
    for col in cols:report.append(name+' '+col+' '+repr(df[col].value_counts(dropna=False).to_dict()))
gtr=c.ResectionType.astype(str).str.strip().str.upper().map({'GTR':1.,'STR':0.})
bad=c[(gtr.notna()) & c.gtr.notna() & (gtr!=c.gtr)]
save(bad,'baylor_resection_disagreements.csv')
report.append('Baylor resection disagreements: '+str(len(bad)))
report.append('Baylor repeated demographics:\n'+pd.read_csv(OUT/'baylor_candidate_duplicate_records.csv').to_string(index=False))
report.append('Zero time UCSF:\n'+a[pd.to_numeric(a.LFFR,errors='coerce')==0][['LFFR','LF','Sex','EOR','Group']].to_string(index=False))
report.append('Zero time Baylor:\n'+c[c.fu_months==0][['Sample','rec','grade','menG']].to_string(index=False))
(OUT/'source_field_audit.txt').write_text('\n\n'.join(report),encoding='utf-8')
print('\n\n'.join(report))

refs=[]
for pmid in ['35108039','41466325','35534562','39169220','42668313','38626948']:
    try:
        rr=requests.get('https://www.ebi.ac.uk/europepmc/webservices/rest/search',params={'query':'EXT_ID:'+pmid+' AND SRC:MED','format':'json','resultType':'core'},timeout=40)
        rr.raise_for_status(); obj=rr.json()
        (OUT/('reference_'+pmid+'.json')).write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')
        result=obj.get('resultList',{}).get('result',[])
        if result:
            z=result[0]; refs.append({'PMID':pmid,'title':z.get('title'),'authors':z.get('authorString'),'journal':z.get('journalInfo',{}).get('journal',{}).get('title'),'year':z.get('pubYear'),'DOI':z.get('doi'),'PMCID':z.get('pmcid')})
    except Exception as exc:refs.append({'PMID':pmid,'error':str(exc)})
save(pd.DataFrame(refs),'verified_reference_metadata.csv')
for pmcid in ['PMC9374001','PMC11313601']+[r.get('PMCID') for r in refs if r.get('PMID') in ['41466325','42668313'] and r.get('PMCID')]:
    try:
        rr=requests.get('https://www.ebi.ac.uk/europepmc/webservices/rest/'+pmcid+'/fullTextXML',timeout=40)
        rr.raise_for_status(); (OUT/(pmcid+'.xml')).write_text(rr.text,encoding='utf-8')
    except Exception as exc:print(pmcid,exc)
