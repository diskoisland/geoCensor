#!/usr/bin/env Rscript
suppressPackageStartupMessages(pkgload::load_all('.',quiet=TRUE))
source('scripts/analysis/conley_survreg_reference.R')
f<-'output/conley_hac_validation_2026-10-31';dir.create(f,recursive=TRUE,showWarnings=FALSE)
settings<-data.frame(scenario=c('independent_64','correlated_64','correlated_256','preferential_64_of_144'),
 n_domain=c(64,64,256,144),n_observed=c(64,64,256,64),correlated=c(FALSE,TRUE,TRUE,TRUE),preferential=c(FALSE,FALSE,FALSE,TRUE))
hs<-c(.5,1,2,4,8,16)
kernels<-c('bartlett_product','bartlett_radial','epanechnikov')
readr::write_csv(settings,file.path(f,'settings.csv'))
writeLines(c(
 'Research-only score sandwich; no production package intervals are changed.',
 '200 development and 200 independently seeded confirmation datasets per setting.',
 'Target: superpopulation geometric mean 50. Observation noise sd 0.2; latent effect sd 0.55.',
 'Exponential spatial covariance exp(-distance/0.8); domains 8x8 or 16x16 with coordinates jittered +/-0.15.',
 'Unequal 1,4,12 observations at each site, independent of concentration; declustering cell width 0.5 with two offsets/axis.',
 'Left threshold exp(log(50)-0.3), upper threshold exp(log(50)+1); actual interval-censored lognormal refits.',
 'Preferential negative control selects 64 of 144 sites with probabilities proportional to exp(2*latent/0.55). It supplies no inclusion-weight correction.',
 'Case weights are applied to both location and log-scale likelihood scores. Bread is the inverse weighted observed information.',
 'Co-located score vectors are summed before constructing the kernel-weighted meat S-prime K S.',
 'Fixed bandwidths 0.5,1,2,4,8,16; product Bartlett, radial Bartlett, and radial Epanechnikov kernels. Coordinates are arbitrary metric units.',
 'Normal critical value 1.96; no data-dependent best-bandwidth selection, finite-sample adjustment, eigenvalue clipping, or tuned t degrees of freedom.',
 'Product Bartlett supplies a PSD reference kernel. Radial kernels may yield non-PSD covariance estimates; invalid results are counted and intervals withheld.',
 'Dense pairwise reference is guarded at 2000 unique sites and must not be run on national data.',
 'SpatialInference vignette informs kernel/bandwidth comparisons; its felm wrapper was not substituted for censored estimation.',
 'A new output folder is required if code or design changes. Confirmation scenarios and candidates are fixed before results are inspected.'
),file.path(f,'design_notes.txt'))
writeLines(capture.output(sessionInfo()),file.path(f,'session_info.txt'))
mu<-log(50)
for(phase in c('development','confirmation')){
 base<-if(phase=='development')8100000L else 9100000L
 for(i in 1:200){
  path<-file.path(f,sprintf('%s_%03d.rds',phase,i));if(file.exists(path))next
  rows<-list()
  for(j in seq_len(nrow(settings))){
   s<-settings[j,];set.seed(base+i*10+j)
   side<-as.integer(sqrt(s$n_domain));grid<-expand.grid(x=0:(side-1),y=0:(side-1))
   coords<-as.matrix(grid)+.5+matrix(runif(s$n_domain*2,-.15,.15),ncol=2)
   C<-if(s$correlated).55^2*exp(-as.matrix(dist(coords))/.8)else diag(.55^2,s$n_domain)
   latent<-as.vector(t(chol(C+diag(1e-10,s$n_domain)))%*%rnorm(s$n_domain))
   counts<-sample(c(1,4,12),s$n_domain,replace=TRUE)
   selected<-if(s$preferential)sort(sample.int(s$n_domain,s$n_observed,prob=exp(2*latent/.55)))else seq_len(s$n_domain)
   id<-rep(selected,counts[selected]);logvalue<-mu+latent[id]+rnorm(length(id),sd=.2)
   value<-exp(logvalue);lod<-exp(mu-.3);upper<-exp(mu+1)
   x<-data.frame(detected_val=ifelse(value<lod|value>upper,NA_real_,value),left_cen=value<lod,
    left_limit=ifelse(value<lod,lod,NA_real_),right_cen=value>upper,right_limit=ifelse(value>upper,upper,NA_real_),
    x=coords[id,1],y=coords[id,2])
   w<-cell_declustering_weights(x,'x','y',.5,n_offsets=2)
   model<-tryCatch(fit_censored_score_reference(x,w),error=function(e)e)
   for(kernel in kernels)for(h in hs){
    if(inherits(model,'error')){
     row<-data.frame(theta=NA_real_,se=NA_real_,lo=NA_real_,hi=NA_real_,status=paste('fit_error',conditionMessage(model)))
    }else{
     out<-spatial_hac_reference(model,coords[id,,drop=FALSE],h,kernel)
     good<-out$status=='experimental_spatial_hac'&&is.finite(out$se_log_geomean)
     row<-data.frame(theta=out$theta,se=out$se_log_geomean,
      lo=if(good)out$theta-qnorm(.975)*out$se_log_geomean else NA_real_,
      hi=if(good)out$theta+qnorm(.975)*out$se_log_geomean else NA_real_,status=out$status)
    }
    row$phase<-phase;row$simulation<-i;row$scenario<-s$scenario;row$kernel<-kernel;row$bandwidth<-h
    row$covered<-is.finite(row$lo)&&row$lo<=mu&&row$hi>=mu
    rows[[length(rows)+1]]<-row
   }
  }
  saveRDS(dplyr::bind_rows(rows),path)
  if(i%%20==0){cat(phase,i,'/ 200 Conley datasets complete across all settings\n');flush.console()}
 }
}
paths<-list.files(f,pattern='^(development|confirmation)_[0-9]+.rds$',full.names=TRUE)
x<-dplyr::bind_rows(lapply(paths,readRDS));readr::write_csv(x,file.path(f,'results.csv'))
s<-dplyr::group_modify(dplyr::group_by(x,.data$phase,.data$scenario,.data$kernel,.data$bandwidth),function(z,key){
 valid<-is.finite(z$lo)&is.finite(z$hi);n<-sum(valid);hit<-sum(z$covered[valid]);mc<-if(n)binom.test(hit,n)$conf.int else c(NA,NA)
 data.frame(n_datasets=nrow(z),n_intervals=n,n_invalid=nrow(z)-n,coverage=if(n)hit/n else NA_real_,
 mc_low=mc[1],mc_high=mc[2],coverage_counting_invalid_as_miss=mean(z$covered),
 mean_se=if(n)mean(z$se[valid])else NA_real_,empirical_sd=sd(z$theta,na.rm=TRUE),
 bias_log=mean(z$theta-mu,na.rm=TRUE),mean_log_width=if(n)mean(z$hi[valid]-z$lo[valid])else NA_real_)
})
readr::write_csv(dplyr::ungroup(s),file.path(f,'summary.csv'))
writeLines('Conley score-HAC validation complete; reference remains experimental.',file.path(f,'COMPLETED.txt'))
print(as.data.frame(s))
