#!/usr/bin/env Rscript
suppressPackageStartupMessages(pkgload::load_all('.',quiet=TRUE))
source('scripts/analysis/conley_survreg_reference.R')
source('scripts/analysis/influence_covariogram_reference.R')
set.seed(141)
coords<-as.matrix(expand.grid(x=0:7,y=0:7))+.2
id<-rep(1:64,each=4);value<-exp(2+rep(rnorm(64,sd=.55),each=4)+rnorm(256,sd=.2))
x<-data.frame(detected_val=value,left_cen=FALSE,left_limit=NA_real_,right_cen=FALSE,right_limit=NA_real_)
model<-fit_censored_score_reference(x)
d<-influence_covariogram_reference(model,coords[id,,drop=FALSE])
a<-aggregate_score_locations(model$scores,coords[id,,drop=FALSE])
for(h in c(2,4,8)){
 dx<-abs(outer(coords[,1],coords[,1],'-'))/h;dy<-abs(outer(coords[,2],coords[,2],'-'))/h
 K<-pmax(1-dx,0)*pmax(1-dy,0)
 direct<-drop(crossprod(d$influence,K%*%d$influence))/64^2
 reference<-spatial_hac_reference(model,coords[id,,drop=FALSE],h)$vcov[1,1]
 stopifnot(abs(direct-reference)<1e-10)
}
for(offset in c(100,1000,1e6)){
 shift<-influence_covariogram_reference(model,coords[id,,drop=FALSE]+offset)
 stopifnot(max(abs(d$influence-shift$influence))<1e-12,
  identical(d$bins$n_pairs,shift$bins$n_pairs),
  max(abs(d$bins$covariance-shift$bins$covariance),na.rm=TRUE)<1e-12,
  max(abs(d$bins$distance-shift$bins$distance),na.rm=TRUE)<1e-7)
}
stopifnot(inherits(try(influence_covariogram_reference(model,coords[id,,drop=FALSE],width=0),silent=TRUE),'try-error'))
cat('Joint influence/HAC identity, translation invariance, and lag validation passed.\n')
