#!/usr/bin/env Rscript
# Fixed-design diagnostic extension. No new routine uncertainty method is selected.
suppressPackageStartupMessages(pkgload::load_all('.',quiet=TRUE))
folder<-'output/irregular_spatial_validation_2026-10-31'
dir.create(folder,recursive=TRUE,showWarnings=FALSE)
settings<-data.frame(scenario=c('irregular_independent_uniform','irregular_correlated_uniform',
 'irregular_correlated_uniform_shifted','irregular_correlated_preferential','irregular_correlated_preferential_shifted'),
 correlated=c(FALSE,TRUE,TRUE,TRUE,TRUE),preferential=c(FALSE,FALSE,FALSE,TRUE,TRUE),
 block_size=c(1,3,3,3,3),origin=c(0,0,.5,0,.5))
readr::write_csv(settings,file.path(folder,'settings.csv'))
writeLines(capture.output(sessionInfo()),file.path(folder,'session_info.txt'))
writeLines(c(
 '100 development and 100 independently seeded confirmation datasets per setting; 199 package resamples each.',
 'Candidate domain: 12x12 sites, each jittered independently by +/-0.15 coordinate units around its cell center.',
 'Observe 64 of 144 sites. Uniform selection is independent of values; preferential selection uses probabilities proportional to exp(2 * latent_effect/0.55), without replacement.',
 'Primary target: superpopulation geometric mean 50. Secondary target: realized equal-site mean exp(log(50)+mean(all 144 latent effects)), including unobserved sites.',
 'Site effects sd 0.55; exponential spatial covariance range parameter 0.8 in correlated settings; independent measurement noise sd 0.2.',
 'Each observed site has 1,4,12 observations, selected independently of grade. Declustering cell width 0.5; all repeated records at each site share influence.',
 'Left limit exp(log(50)-0.3), upper limit exp(log(50)+1). Limits and censoring intervals are retained.',
 'Resampling block widths/origins are fixed before results are inspected. Shifted settings share data and seeds with their unshifted counterparts.',
 'Hidden complete-data observed-site means diagnose selection bias; they are not available to a real censored-data analysis.',
 'This study tests preferential omission of sites. No inverse-inclusion-probability correction or spatial selection model is supplied.',
 'All intervals remain experimental. A new output folder is required for changed implementation or design.'
),file.path(folder,'design_notes.txt'))
mu<-log(50);grid<-expand.grid(x=0:11,y=0:11)
for(phase in c('development','confirmation')) {
 base<-if(phase=='development')6100000L else 7100000L
 for(i in 1:100) {
  path<-file.path(folder,sprintf('%s_%03d.rds',phase,i))
  if(file.exists(path))next
  result<-list()
  for(j in seq_len(nrow(settings))) {
   s<-settings[j,]
   set.seed(base+i)
   coords<-as.matrix(grid)+.5+matrix(runif(288,-.15,.15),ncol=2)
   covariance<-if(s$correlated).55^2*exp(-as.matrix(dist(coords))/.8)else diag(.55^2,144)
   latent<-as.vector(t(chol(covariance+diag(1e-10,144)))%*%rnorm(144))
   noise<-matrix(rnorm(144*12,sd=.2),nrow=144)
   counts<-sample(c(1,4,12),144,replace=TRUE)
   probabilities<-if(s$preferential)exp(2*latent/.55)else rep(1,144)
   selected<-sort(sample.int(144,64,replace=FALSE,prob=probabilities))
   id<-rep(selected,counts[selected])
   logvalue<-unlist(lapply(selected,function(h)mu+latent[h]+noise[h,seq_len(counts[h])]))
   value<-exp(logvalue);lod<-exp(mu-.3);upper<-exp(mu+1)
   x<-data.frame(x=coords[id,1],y=coords[id,2],detected_val=ifelse(value<lod|value>upper,NA_real_,value),
    left_cen=value<lod,left_limit=ifelse(value<lod,lod,NA_real_),right_cen=value>upper,right_limit=ifelse(value>upper,upper,NA_real_))
   b<-spatial_block_bootstrap(x,'x','y',.5,s$block_size,n_boot=199,n_offsets=2,
    block_origin=rep(s$origin,2),seed=base+20000+i)
   z<-b$summary;valid<-z$status=='experimental_interval'
   finite<-mu+mean(latent)
   result[[j]]<-data.frame(phase=phase,simulation=i,scenario=s$scenario,geomean=z$geomean,
    lo=z$ci_low,hi=z$ci_high,n_blocks=z$n_blocks,n_success=z$n_success,status=z$status,
    covered=valid&&z$ci_low<=50&&z$ci_high>=50,finite_target=exp(finite),
    finite_covered=valid&&z$ci_low<=exp(finite)&&z$ci_high>=exp(finite),
    log_error_population=log(z$geomean)-mu,log_error_finite=log(z$geomean)-finite,
    complete_observed_site_logmean=mean(as.numeric(tapply(logvalue,id,mean))),
    full_domain_logmean=finite,n_records=nrow(x),fraction_left=mean(x$left_cen))
  }
  saveRDS(dplyr::bind_rows(result),path)
  if(i%%10==0){cat(phase,i,'/ 100 completed for all irregular settings\n');flush.console()}
 }
}
paths<-list.files(folder,pattern='^(development|confirmation)_[0-9]+.rds$',full.names=TRUE)
raw<-dplyr::bind_rows(lapply(paths,readRDS))
readr::write_csv(raw,file.path(folder,'results.csv'))
summary<-dplyr::group_modify(dplyr::group_by(raw,.data$phase,.data$scenario),function(x,key){
 valid<-x$status=='experimental_interval';n<-sum(valid);hit<-sum(x$covered[valid])
 mc<-if(n)binom.test(hit,n)$conf.int else c(NA,NA)
 data.frame(n_datasets=nrow(x),n_intervals=n,coverage=if(n)hit/n else NA_real_,mc_low=mc[1],mc_high=mc[2],
 finite_domain_coverage=if(n)mean(x$finite_covered[valid])else NA_real_,
 mean_log_error_population=mean(x$log_error_population,na.rm=TRUE),
 mean_log_error_finite=mean(x$log_error_finite,na.rm=TRUE),
 complete_data_selection_bias=mean(x$complete_observed_site_logmean-x$full_domain_logmean),
 mean_geomean=mean(x$geomean,na.rm=TRUE),mean_blocks=mean(x$n_blocks),
 fit_success_fraction=sum(x$n_success)/(199*nrow(x)))
})
readr::write_csv(dplyr::ungroup(summary),file.path(folder,'summary.csv'))
writeLines('Irregular sampling validation completed.',file.path(folder,'COMPLETED.txt'))
print(as.data.frame(summary))
