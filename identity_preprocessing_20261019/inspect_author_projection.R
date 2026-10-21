Sys.setlocale("LC_CTYPE",".UTF-8")
root<-commandArgs(trailingOnly=TRUE)[1]
for(gsm in c("GSM8856328","GSM9143148","GSM5566283")){
 f<-file.path(root,"results/raw_preprocessed",paste0(gsm,"_author.rds"))
 if(!file.exists(f))next
 x<-readRDS(f)
 cat(gsm,"n_projection",nrow(x$probe_projection),"duplicate_base",sum(duplicated(x$probe_projection$base_cpg)),"array",x$summary$array_type,"\n")
 print(head(x$probe_projection));print(x$summary)
}
