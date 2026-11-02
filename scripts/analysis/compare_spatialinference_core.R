#!/usr/bin/env Rscript
# Optional external-package cross-check, isolated from geoCensor's dependencies.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)) .libPaths(c(args[1],.libPaths()))
stopifnot(requireNamespace('SpatialInference',quietly=TRUE))
f<-'output/conley_reference_checks_2026-10-31';dir.create(f,recursive=TRUE,showWarnings=FALSE)
M<-rbind(c(43,-89),c(43.01,-89),c(44,-89),c(44,-89.01))
S<-matrix(c(1,2,-1,-2,2,-1,3,-4),ncol=2)
rows<-list()
for(kernel in c('bartlett','epanechnikov','uniform','gaussian','parzen','biweight')){
 K<-SpatialInference::DistMat(M,20,kernel,'Haversine')
 reference<-crossprod(S,K%*%S)
 dense<-SpatialInference::XeeXhC(M,20,S,rep(1,4),4L,2L,kernel,'Haversine')
 stopifnot(max(abs(dense-reference))<1e-10,max(abs(K-t(K)))<1e-12,all(diag(K)==1))
 large<-tryCatch(SpatialInference::XeeXhC_Lg(M,20,S,rep(1,4),4L,2L,kernel,'Haversine'),error=function(e)e)
 rows[[kernel]]<-data.frame(kernel=kernel,dense_max_error=max(abs(dense-reference)),
  large_status=if(inherits(large,'error'))conditionMessage(large)else 'completed',
  large_max_error=if(inherits(large,'error'))NA_real_ else max(abs(large-reference)))
}
result<-do.call(rbind,rows)
write.csv(result,file.path(f,'spatialinference_core_comparison.csv'),row.names=FALSE)
writeLines(c(paste('SpatialInference version:',packageVersion('SpatialInference')),
 paste('Remote SHA:',packageDescription('SpatialInference')$RemoteSha),
 'Inputs: 4 latitude/longitude rows, two arbitrary score columns, e=1, n1=4, k=2, cutoff 20 km.',
 'Dense core matches independently assembled S-prime K S using the package kernel matrix.',
 'This checks computational algebra, not statistical coverage or validity for censored data.',
 'The large-data path error is recorded for this exact version and fixture; no general claim is made about other versions.',
 'The package was installed in a temporary reference library, not added as a geoCensor dependency.'
),file.path(f,'spatialinference_comparison_notes.txt'))
print(result)
