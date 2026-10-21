# Author's published 403-probe KNN reference and caret::knn3(k=10).
# Processed-beta sensitivity only: does NOT reproduce the author's IDAT preprocessing.
args <- commandArgs(trailingOnly=TRUE)
root <- args[1]
suppressPackageStartupMessages(library(caret))
reference <- readRDS(file.path(root,"official_model","knn_reference_fc2.rds"))
cpgs <- colnames(reference)[colnames(reference)!="group"]
model <- caret::knn3(x=reference[,cpgs],y=reference$group,k=10)
means <- colMeans(reference[,cpgs])
coverage <- list()
for (cohort in c("baylor","tubingen","release186","release37","ucsf")) {
  file <- file.path(root,paste0("official_panel_",cohort,".csv.gz"))
  if (!file.exists(file)) next
  raw <- read.csv(gzfile(file),row.names=1,check.names=FALSE)
  raw <- t(as.matrix(raw[cpgs,,drop=FALSE]));storage.mode(raw)<-"numeric"
  finite <- rowSums(is.finite(raw)); eligible <- finite>=ceiling(length(cpgs)*.5)
  dat <- raw[eligible,,drop=FALSE]
  if (!nrow(dat)) {
    cat(cohort,"No samples have >=50% of official probes. No adapted prediction generated.\n")
    next
  }
  for (j in seq_along(cpgs)) dat[!is.finite(dat[,j]),j]<-means[j]
  pr <- predict(model,as.data.frame(dat),type="prob")
  best <- max.col(pr,ties.method="first")
  tie <- rowSums(pr==apply(pr,1,max))>1
  pred <- data.frame(id=rownames(dat),p_low=pr[,1],p_medium=pr[,2],p_high=pr[,3],
      predicted_class=colnames(pr)[best],confidence=apply(pr,1,max),tied_top=tie,
      confidence_ge_075=apply(pr,1,max)>=.75,finite_features=finite[eligible],
      imputed_features=length(cpgs)-finite[eligible],
      in_official_training=rownames(dat)%in%rownames(reference),
      execution="author_KNN_on_processed_beta_training_mean_imputation")
  pred$predicted_class[pred$tied_top] <- NA # no arbitrary tie resolution for outcome models
  write.csv(pred,file.path(root,paste0("official_knn_predictions_",cohort,".csv")),row.names=FALSE)
  cat(cohort,"eligible",nrow(dat),"class counts\n");print(table(pred$predicted_class,useNA="always"))
  # Verify caret probabilities against independent Euclidean KNN calculation, except ties at kth distance.
  check <- min(10,nrow(dat)); maxdiff<-0
  for (i in seq_len(check)) {
    dist <- rowSums((sweep(as.matrix(reference[,cpgs]),2,dat[i,],"-"))^2)
    ord <- order(dist); neighbors <- which(dist<=sort(dist)[10]*(1+1e-12)+1e-14)
    hand <- as.numeric(table(factor(reference$group[neighbors],levels=levels(reference$group))))/length(neighbors)
    maxdiff <- max(maxdiff,max(abs(hand-pr[i,])))
  }
  cat("Independent KNN probability max difference",maxdiff,"\n")
  stopifnot(maxdiff<1e-8)
}
capture.output(sessionInfo(),file=file.path(root,"official_knn_R_sessionInfo.txt"))
