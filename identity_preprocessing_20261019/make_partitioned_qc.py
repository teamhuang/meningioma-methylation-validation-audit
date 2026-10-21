from pathlib import Path
ROOT=Path(__file__).resolve().parent
s=(ROOT/'run_raw_qc.R').read_text(encoding='utf-8')
s=s.replace('root<-args[1];cohort<-args[2]','root<-args[1];cohort<-args[2]\npartition<-as.integer(args[3]);partitions<-as.integer(args[4])')
s=s.replace('out<-file.path(root,"results/raw_qc")','d<-d[((seq_len(nrow(d))-1L)%%partitions)+1L==partition,,drop=FALSE]\nout<-file.path(root,"results/raw_qc")')
s=s.replace('cat("QC",cohort,start','cat("QC",cohort,partition,start')
s=s.replace('write.csv(result,file.path(out,paste0(chunk$GSM[i],"_qc.csv")),row.names=FALSE)',
    'target<-file.path(out,paste0(chunk$GSM[i],"_qc.csv"));temp<-paste0(target,".part")\n    write.csv(result,temp,row.names=FALSE);stopifnot(file.rename(temp,target))')
s=s.replace('paste0("raw_qc_session_",cohort,".txt")','paste0("raw_qc_session_",cohort,"_part",partition,"of",partitions,".txt")')
(ROOT/'run_raw_qc_partitioned.R').write_text(s,encoding='utf-8')
