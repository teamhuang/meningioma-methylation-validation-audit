# Execute the author's preprocessing functions, with no plotting or target-data fitting.
args<-commandArgs(trailingOnly=TRUE)
Sys.setlocale("LC_CTYPE",".UTF-8")
root<-args[1]
mode<-if(length(args)>=2)args[2] else "pilot"
.libPaths(c(file.path(root,"R_library"),.libPaths()))
suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(caret))
cat("R",R.version.string,"minfi",as.character(packageVersion("minfi")),"\n")
cat("Loaded minfi path",find.package("minfi"),"\n")
author.script<-file.path(root,"..","strengthening_20261019","results","official_model","run_meningioma_risk_analysis.R")
expressions<-parse(author.script)
functions<-c("MNPpreprocessIllumina","MNPnormalize.illumina.control","find_failed_colors","calculate_matrix")
for(expr in expressions){
  if(is.call(expr)&&identical(expr[[1]],as.name("<-"))&&is.symbol(expr[[2]])&&as.character(expr[[2]])%in%functions)eval(expr)
}
stopifnot(all(vapply(functions,exists,logical(1))))
ref<-readRDS(file.path(root,"..","strengthening_20261019","results","official_model","knn_reference_fc2.rds"))
cpgs<-colnames(ref)[colnames(ref)!="group"]
knn<-caret::knn3(x=ref[,cpgs],y=ref$group,k=10)
manifest<-read.csv(file.path(root,"results",paste0("requested_raw_",mode,".csv")),check.names=FALSE)
out<-file.path(root,"results","raw_preprocessed");dir.create(out,showWarnings=FALSE,recursive=TRUE)
rows<-list()
for(i in seq_len(nrow(manifest))){
  r<-manifest[i,];prefix<-file.path(root,"raw_idat",r$cohort,r$barcode)
  target<-file.path(out,paste0(r$GSM,"_author.rds"))
  cat("START",i,nrow(manifest),r$cohort,r$GSM,r$barcode,"\n");flush.console()
  result<-tryCatch({
    if(file.exists(target)&&identical(readRDS(target)$minfi_version,as.character(packageVersion("minfi")))){
      readRDS(target)$summary
    }else{
      rg<-read.metharray(prefix,force=TRUE)
      rawSNP<-getSnpBeta(rg)
      ms<-MNPpreprocessIllumina(rg)
      bet<-calculate_matrix(ms)
      original.names<-colnames(bet)
      # EPICv2 replicates require an explicit base-CpG projection, recorded below.
      base.names<-sub("_.*$","",original.names)
      projection<-data.frame(array_probe=original.names,base_cpg=base.names,stringsAsFactors=FALSE)
      beta.v<-as.numeric(bet[1,]);names(beta.v)<-base.names
      if(anyDuplicated(base.names))beta.v<-vapply(split(beta.v,base.names),function(x)mean(x,na.rm=TRUE),numeric(1))
      x<-beta.v[cpgs]
      missing<-sum(!is.finite(x))
      if(missing>0){
        probabilities<-setNames(rep(NA_real_,3),levels(ref$group));label<-NA_character_;confidence<-NA_real_;tied<-NA
      }else{
        xx<-as.data.frame(t(x));colnames(xx)<-cpgs
        probabilities<-as.numeric(predict(knn,xx,type="prob")[1,]);names(probabilities)<-levels(ref$group)
        tied<-sum(probabilities==max(probabilities))>1
        label<-if(tied)NA_character_ else names(probabilities)[which.max(probabilities)]
        confidence<-max(probabilities)
      }
      item<-data.frame(cohort=r$cohort,GSM=r$GSM,title=r$title,barcode=r$barcode,array_type=annotation(rg)[[1]],
        full_array_probes=length(original.names),unique_base_probes=length(beta.v),finite_classifier_cpgs=sum(is.finite(x)),
        missing_classifier_cpgs=missing,p_class1=probabilities[1],p_class2=probabilities[2],p_class3=probabilities[3],
        predicted_class=label,confidence=confidence,tied_top=tied,snp_probes=nrow(rawSNP),
        preprocessing="author_control10000_background_illumina_beta_offset100",stringsAsFactors=FALSE,row.names=NULL)
      saveRDS(list(summary=item,classifier_beta=x,snp_beta=rawSNP,probe_projection=projection,
          preprocessMethod=preprocessMethod(ms),minfi_version=as.character(packageVersion("minfi"))),target)
      write.csv(data.frame(cpg=cpgs,beta=as.numeric(x)),file.path(out,paste0(r$GSM,"_403beta.csv")),row.names=FALSE)
      write.csv(data.frame(snp=rownames(rawSNP),beta=as.numeric(rawSNP[,1])),file.path(out,paste0(r$GSM,"_snp.csv")),row.names=FALSE)
      item
    }
  },error=function(e)data.frame(cohort=r$cohort,GSM=r$GSM,title=r$title,barcode=r$barcode,error=conditionMessage(e)))
  rows[[length(rows)+1]]<-result
  # bind with missing fields to preserve failures without dropping successful samples
  fields<-unique(unlist(lapply(rows,names)))
  joined<-do.call(rbind,lapply(rows,function(x){for(z in setdiff(fields,names(x)))x[[z]]<-NA;x[,fields,drop=FALSE]}))
  write.csv(joined,file.path(root,"results",paste0("author_raw_summary_",mode,".csv")),row.names=FALSE)
  print(result);flush.console();gc()
}
capture.output(sessionInfo(),file=file.path(root,"results",paste0("author_raw_session_",mode,".txt")))
