args<-commandArgs(trailingOnly=TRUE)
Sys.setlocale("LC_CTYPE",".UTF-8")
root<-args[1];cohort<-args[2]
batch.size<-if(length(args)>=3)as.integer(args[3]) else 8L
partition<-if(length(args)>=4)as.integer(args[4]) else 1L
partitions<-if(length(args)>=5)as.integer(args[5]) else 1L
.libPaths(c(file.path(root,"R_library"),.libPaths()))
suppressPackageStartupMessages(library(minfi));suppressPackageStartupMessages(library(caret))
stopifnot(as.character(packageVersion("minfi"))=="1.43.1")
functions<-c("MNPpreprocessIllumina","MNPnormalize.illumina.control")
author<-file.path(root,"..","strengthening_20261019/results/official_model")
for(expr in parse(file.path(author,"run_meningioma_risk_analysis.R"))){
 if(is.call(expr)&&identical(expr[[1]],as.name("<-"))&&is.symbol(expr[[2]])&&as.character(expr[[2]])%in%functions)eval(expr)
}
ref<-readRDS(file.path(author,"knn_reference_fc2.rds"));cpgs<-colnames(ref)[colnames(ref)!="group"]
knn<-caret::knn3(x=ref[,cpgs],y=ref$group,k=10)
d<-read.csv(file.path(root,"results",paste0("requested_raw_",cohort,".csv")),check.names=FALSE)
d<-d[((seq_len(nrow(d))-1L)%%partitions)+1L==partition,,drop=FALSE]
out<-file.path(root,"results/raw_preprocessed");dir.create(out,showWarnings=FALSE,recursive=TRUE)
logfile<-file.path(root,"results",paste0("batch_log_",cohort,"_part",partition,".csv"))
logs<-list()
for(start in seq(1L,nrow(d),by=batch.size)){
 chunk<-d[start:min(start+batch.size-1L,nrow(d)),,drop=FALSE]
 already<-vapply(chunk$GSM,function(g)file.exists(file.path(out,paste0(g,"_summary.csv"))),logical(1))
 chunk<-chunk[!already,,drop=FALSE];if(!nrow(chunk))next
 complete<-vapply(chunk$barcode,function(b)all(file.exists(paste0(file.path(root,"raw_idat",cohort,b),c("_Grn.idat","_Red.idat")))),logical(1))
 if(!all(complete)){
   for(g in chunk$GSM[!complete])logs[[length(logs)+1L]]<-data.frame(GSM=g,status="raw_pair_missing",error="download not yet complete")
 }
 chunk<-chunk[complete,,drop=FALSE];if(!nrow(chunk))next
 # Completed single-sample pilots are reprocessed in the first batch to test exact equivalence.
 cat("BATCH",cohort,partition,start,"n",nrow(chunk),"time",as.character(Sys.time()),"\n");flush.console()
 result<-tryCatch({
   rg<-read.metharray(file.path(root,"raw_idat",cohort,chunk$barcode),force=TRUE)
   observed<-basename(sampleNames(rg));stopifnot(identical(as.character(observed),as.character(chunk$barcode)))
   snp<-getSnpBeta(rg)
   ms<-MNPpreprocessIllumina(rg)
   mm<-getMeth(ms);uu<-getUnmeth(ms)
   allnames<-rownames(mm);base<-sub("_.*$","",allnames)
   keep<-which(base%in%cpgs);project<-data.frame(array_probe=allnames[keep],base_cpg=base[keep])
   bb<-mm[keep,,drop=FALSE]/(mm[keep,,drop=FALSE]+uu[keep,,drop=FALSE]+100)
   if(anyDuplicated(project$base_cpg)){
     totals<-rowsum(bb,group=project$base_cpg,reorder=FALSE,na.rm=TRUE)
     counts<-rowsum(is.finite(bb)*1,group=project$base_cpg,reorder=FALSE)
     bb<-totals/counts
   }else rownames(bb)<-project$base_cpg
   roworder<-match(cpgs,rownames(bb));xx<-matrix(NA_real_,nrow=nrow(chunk),ncol=length(cpgs),dimnames=list(chunk$GSM,cpgs))
   good<-!is.na(roworder);xx[,good]<-t(bb[roworder[good],,drop=FALSE])
   for(i in seq_len(nrow(chunk))){
     r<-chunk[i,];x<-xx[i,];missing<-sum(!is.finite(x))
     if(missing){pr<-setNames(rep(NA_real_,3),levels(ref$group));tied<-NA;label<-NA_character_;score<-NA_real_}else{
       dat<-as.data.frame(t(x));colnames(dat)<-cpgs
       pr<-as.numeric(predict(knn,dat,type="prob")[1,]);names(pr)<-levels(ref$group)
       tied<-sum(pr==max(pr))>1;label<-if(tied)NA_character_ else names(pr)[which.max(pr)];score<-max(pr)
     }
     target<-file.path(out,paste0(r$GSM,"_author.rds"))
     equality<-NA_real_
     if(file.exists(target)){
       old<-readRDS(target)
       if(identical(old$minfi_version,"1.43.1"))equality<-max(abs(old$classifier_beta-x),na.rm=TRUE)
       if(is.finite(equality)&&equality>1e-12)stop("Batched vs single-sample preprocessing differs for ",r$GSM,": ",equality)
     }
     item<-data.frame(cohort=cohort,GSM=r$GSM,title=r$title,barcode=r$barcode,array_type=annotation(rg)[[1]],
       full_array_probes=nrow(mm),finite_classifier_cpgs=sum(is.finite(x)),missing_classifier_cpgs=missing,
       p_class1=pr[1],p_class2=pr[2],p_class3=pr[3],predicted_class=label,confidence=score,tied_top=tied,
       snp_probes=nrow(snp),batch_single_max_abs_difference=equality,preprocessing="author_control10000_background_illumina_beta_offset100",
       stringsAsFactors=FALSE,row.names=NULL)
     saveRDS(list(summary=item,classifier_beta=x,snp_beta=snp[,i,drop=FALSE],probe_projection=project,
       preprocessMethod=preprocessMethod(ms),minfi_version="1.43.1"),target)
     write.csv(data.frame(cpg=cpgs,beta=as.numeric(x)),file.path(out,paste0(r$GSM,"_403beta.csv")),row.names=FALSE)
     write.csv(data.frame(snp=rownames(snp),beta=as.numeric(snp[,i])),file.path(out,paste0(r$GSM,"_snp.csv")),row.names=FALSE)
     write.csv(item,file.path(out,paste0(r$GSM,"_summary.csv")),row.names=FALSE)
     logs[[length(logs)+1L]]<-data.frame(GSM=r$GSM,status="success",error="")
   }
   TRUE
 },error=function(e){cat("BATCH_ERROR",conditionMessage(e),"\n");for(g in chunk$GSM)logs[[length(logs)+1L]]<<-data.frame(GSM=g,status="failed",error=conditionMessage(e));FALSE})
 if(length(logs))write.csv(do.call(rbind,logs),logfile,row.names=FALSE)
 rm(list=intersect(ls(),c("rg","ms","mm","uu","bb","xx","snp")));gc()
}
capture.output(sessionInfo(),file=file.path(root,"results",paste0("author_batch_session_",cohort,"_part",partition,".txt")))
