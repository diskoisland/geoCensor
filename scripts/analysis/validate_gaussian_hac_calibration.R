#!/usr/bin/env Rscript
suppressPackageStartupMessages(pkgload::load_all('.',quiet=TRUE))
source('scripts/analysis/conley_survreg_reference.R')
source('scripts/analysis/influence_covariogram_reference.R')
args<-commandArgs(trailingOnly=TRUE)
if(length(args)>1)stop('Supply at most one fresh output directory')
f<-if(length(args))args[1]else 'output/gaussian_hac_calibration_2026-10-31'
dir.create(f,recursive=TRUE,showWarnings=FALSE)
# Never silently reuse checkpoints after source or runtime changes.
sources<-c('scripts/analysis/validate_gaussian_hac_calibration.R',
 'scripts/analysis/influence_covariogram_reference.R','scripts/analysis/conley_survreg_reference.R',
 sort(list.files('R',pattern='\\.[Rr]$',full.names=TRUE)))
signature<-list(format_version=2L,source_md5=tools::md5sum(sources),R=R.version$version.string,platform=R.version$platform,
 survival=as.character(packageVersion('survival')))
manifest<-file.path(f,'study_manifest.rds')
if(file.exists(manifest)){
 if(!identical(signature,readRDS(manifest)))stop('Checkpoint source/runtime differs; use a fresh output directory')
}else{
 if(length(list.files(f,pattern='\\.rds$',recursive=TRUE)))stop('Legacy checkpoints lack a source manifest; preserve them and use a fresh output directory')
 saveRDS(signature,manifest)
}
hs<-c(2,4,8,16);ranges<-c(0,.25,.8,2);sizes<-c(64L,256L)
# 10,000 calibration draws and 5,000 independent Gaussian validation draws/design.
# Experimental critical values are conditional on assumed Gaussian covariance.
critical_rows<-list();gaussian_rows<-list()
make_K<-function(coords,h){dx<-abs(outer(coords[,1],coords[,1],'-'))/h;dy<-abs(outer(coords[,2],coords[,2],'-'))/h;pmax(1-dx,0)*pmax(1-dy,0)}
make_C<-function(coords,range){n<-nrow(coords);if(range==0)diag(.55^2+.2^2/4,n)else .55^2*exp(-as.matrix(dist(coords))/range)+diag(.2^2/4,n)}
for(n in sizes){
 side<-as.integer(sqrt(n));coords<-as.matrix(expand.grid(x=0:(side-1),y=0:(side-1)))+.2
 cp<-file.path(f,paste0('gaussian_design_',n,'.rds'))
 if(file.exists(cp)){saved<-readRDS(cp);critical_rows[[as.character(n)]]<-saved$critical;gaussian_rows[[as.character(n)]]<-saved$validation;next}
 crit<-list();validation<-list()
 for(range in ranges){
  set.seed(10100000+n*10+round(range*10))
  L<-t(chol(make_C(coords,range)))
  Z<-L%*%matrix(rnorm(n*10000),nrow=n)
  means<-colMeans(Z);U<-sweep(Z,2,means,'-')
  for(h in hs){K<-make_K(coords,h);v<-colSums(U*(K%*%U))/n^2
   tstat<-abs(means/sqrt(pmax(v,.Machine$double.eps)))
   crit[[paste(range,h)]]<-data.frame(n_sites=n,assumed_range=range,bandwidth=h,critical=unname(quantile(tstat,.95)),n_calibration=10000)
  }
 }
 crit<-dplyr::bind_rows(crit)
 for(range in ranges){
  set.seed(11100000+n*10+round(range*10))
  Z<-t(chol(make_C(coords,range)))%*%matrix(rnorm(n*5000),nrow=n)
  means<-colMeans(Z);U<-sweep(Z,2,means,'-')
  for(h in hs){v<-colSums(U*(make_K(coords,h)%*%U))/n^2;se<-sqrt(pmax(v,.Machine$double.eps))
   values<-c(normal=qnorm(.975),gaussian_iid=crit$critical[crit$assumed_range==0&crit$bandwidth==h],
    gaussian_reference_0.8=crit$critical[crit$assumed_range==.8&crit$bandwidth==h],
    gaussian_range_grid_max=max(crit$critical[crit$bandwidth==h]))
   for(method in names(values))validation[[paste(range,h,method)]]<-data.frame(n_sites=n,true_range=range,bandwidth=h,method=method,
    critical=values[[method]],n_validation=5000,coverage=mean(abs(means)<=values[[method]]*se))
  }
 }
 validation<-dplyr::bind_rows(validation)
 saveRDS(list(critical=crit,validation=validation),cp)
 critical_rows[[as.character(n)]]<-crit;gaussian_rows[[as.character(n)]]<-validation
 cat('Gaussian calibration/independent validation complete for',n,'sites\n');flush.console()
}
critical<-dplyr::bind_rows(critical_rows);gaussian<-dplyr::bind_rows(gaussian_rows)
readr::write_csv(critical,file.path(f,'critical_values.csv'));readr::write_csv(gaussian,file.path(f,'gaussian_validation.csv'))
writeLines(c(
 'Research sensitivity experiment, not an implementation of a published fixed-b or Muller-Watson method.',
 'Primary target: superpopulation geometric mean 50. Product Bartlett kernel; fixed bandwidths 2,4,8,16; 64 or 256 sites on a fixed regular grid.',
 'Critical values simulated from Gaussian site means using weighted within-site noise variance 0.2^2/4 and latent variance 0.55^2.',
 'Assumed correlation ranges 0 (independent),0.25,0.8,2. Gaussian calibration uses known specified covariance, unavailable in an actual NGDB analysis.',
 '10,000 calibration draws and 5,000 independently seeded Gaussian validation draws per design/range; candidates frozen before censored validation.',
 'Candidates: normal 1.96, Gaussian iid finite-design calibration, Gaussian reference range 0.8, maximum critical across the fixed range grid.',
 'The range-grid maximum covers only this finite assumed class in the Gaussian reference, subject to Monte Carlo error. It is not a bound for arbitrary spatial models.',
 'Censored model validation: 300 development and 500 independently seeded confirmation datasets per each n/range setting; true ranges 0.25,0.8,2.',
 'Four observations/site; latent sd 0.55, observation noise sd 0.2. Left exp(log(50)-0.3), upper exp(log(50)+1). All sites observed; no preferential selection.',
 'Joint weighted censored likelihood scores and bread are used in HAC. Gaussian calibration is deliberately tested for transfer to nonlinear censored estimates, not assumed valid.',
 'Influence covariogram: n times site score multiplied by the location row of the joint bread; bins width 0.5 to lag 6, at least 20 pairs per bin.',
 'First-zero-crossing range is diagnostic only. No interval chooses a bandwidth using the estimated range in this experiment; unavailable ranges are kept missing.',
 'No normality, fixed-range, or covariance assumption in this study is asserted for NGDB data. No production API is changed.'
),file.path(f,'design_notes.txt'))
writeLines(capture.output(sessionInfo()),file.path(f,'session_info.txt'))
for(phase in c('development','confirmation')){
 m<-if(phase=='development')300L else 500L;base<-if(phase=='development')12100000L else 13100000L
 for(i in seq_len(m)){
  path<-file.path(f,sprintf('%s_%03d.rds',phase,i));if(file.exists(path))next
  rows<-list();diag_rows<-list()
  for(n in sizes)for(range in c(.25,.8,2)){
   set.seed(base+i*100+n+round(range*10))
   side<-as.integer(sqrt(n));coords<-as.matrix(expand.grid(x=0:(side-1),y=0:(side-1)))+.2
   C<-.55^2*exp(-as.matrix(dist(coords))/range)+diag(1e-10,n)
   latent<-as.vector(t(chol(C))%*%rnorm(n));id<-rep(seq_len(n),each=4)
   value<-exp(log(50)+latent[id]+rnorm(length(id),sd=.2));lod<-exp(log(50)-.3);upper<-exp(log(50)+1)
   x<-data.frame(detected_val=ifelse(value<lod|value>upper,NA_real_,value),left_cen=value<lod,
    left_limit=ifelse(value<lod,lod,NA_real_),right_cen=value>upper,right_limit=ifelse(value>upper,upper,NA_real_))
   model<-fit_censored_score_reference(x)
   diagnostic<-influence_covariogram_reference(model,coords[id,,drop=FALSE])
   diag_rows[[paste(n,range)]]<-data.frame(phase=phase,simulation=i,n_sites=n,true_range=range,
    estimated_zero_crossing=diagnostic$range,range_status=diagnostic$status)
   for(h in hs){out<-spatial_hac_reference(model,coords[id,,drop=FALSE],h)
    crit<-critical[critical$n_sites==n&critical$bandwidth==h,]
    values<-c(normal=qnorm(.975),gaussian_iid=crit$critical[crit$assumed_range==0],
     gaussian_reference_0.8=crit$critical[crit$assumed_range==.8],gaussian_range_grid_max=max(crit$critical))
    for(method in names(values))rows[[paste(n,range,h,method)]]<-data.frame(phase=phase,simulation=i,n_sites=n,
     true_range=range,bandwidth=h,method=method,critical=values[[method]],theta=out$theta,se=out$se_log_geomean,
     covered=abs(out$theta-log(50))<=values[[method]]*out$se_log_geomean)
   }
  }
  saveRDS(list(results=dplyr::bind_rows(rows),diagnostics=dplyr::bind_rows(diag_rows)),path)
  if(i%%50==0){cat('Censored',phase,i,'/',m,'complete\n');flush.console()}
 }
}
paths<-list.files(f,pattern='^(development|confirmation)_[0-9]+.rds$',full.names=TRUE)
all<-lapply(paths,readRDS);raw<-dplyr::bind_rows(lapply(all,`[[`,'results'));diagnostics<-dplyr::bind_rows(lapply(all,`[[`,'diagnostics'))
readr::write_csv(raw,file.path(f,'censored_results.csv'));readr::write_csv(diagnostics,file.path(f,'influence_ranges.csv'))
summary<-dplyr::group_modify(dplyr::group_by(raw,.data$phase,.data$n_sites,.data$true_range,.data$bandwidth,.data$method),function(x,key){
 hit<-sum(x$covered);n<-nrow(x);mc<-binom.test(hit,n)$conf.int
 data.frame(n_datasets=n,coverage=hit/n,mc_low=mc[1],mc_high=mc[2],bias_log=mean(x$theta-log(50)),mean_se=mean(x$se),critical=x$critical[1])
})
readr::write_csv(dplyr::ungroup(summary),file.path(f,'censored_summary.csv'))
writeLines('Gaussian-reference calibration and independent censored confirmation complete.',file.path(f,'COMPLETED.txt'))
cat('Calibration study complete.\n')
