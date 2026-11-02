#!/usr/bin/env Rscript
# Prespecified benchmark and independent confirmation seeds. No production CI changes.
suppressPackageStartupMessages(pkgload::load_all('.', quiet=TRUE))
folder <- 'output/spatial_bootstrap_diagnostics_2026-10-31'
dir.create(folder,recursive=TRUE,showWarnings=FALSE)
mu <- log(50); k <- 64L
intervals <- function(theta,z) {
  q <- stats::quantile(z,c(.025,.975),names=FALSE)
  list(percentile_log=q,basic_log=2*theta-rev(q),
       normal_log_centered=theta+c(-1,1)*qnorm(.975)*sd(z))
}
record <- function(phase,i,scenario,method,theta,ci,finite_target,se=NA_real_,success=NA_real_) {
  data.frame(phase=phase,simulation=i,scenario=scenario,method=method,
    theta=theta,lo=ci[1],hi=ci[2],covered=is.finite(ci[1])&&ci[1]<=mu&&ci[2]>=mu,
    finite_target=finite_target,finite_covered=is.finite(ci[1])&&ci[1]<=finite_target&&ci[2]>=finite_target,
    boot_se=se,n_success=success)
}
writeLines(c(
 'Primary estimand: superpopulation geometric mean exp(E(log concentration)) = 50.',
 'Secondary estimand: realized equal-site geometric mean exp(mu + mean(latent site effects)); reported separately, never substituted for the primary target.',
 'All 64 sites observed on an 8x8 grid; each site is one equal-area unit. No unsampled-area inference is tested.',
 'Site effect sd 0.55; observation noise sd 0.2. Four observations/site in regular designs.',
 'Uneven independent density: 1, 4, or 12 observations/site, independent of latent value.',
 'Correlated site covariance: 0.55^2 exp(-distance/0.8). Cell size 1, block sizes 1 or 2.',
 'Censoring: left exp(mu-0.3), right exp(mu+1); uncensored counterparts retain the same underlying measurements.',
 'Analytic benchmark: 2000 discovery + 3000 confirmation datasets; 999 resamples each using the exact site-log-mean representation.',
 'Censored diagnostic: 200 discovery + 200 confirmation datasets per setting; 199 actual package bootstrap refits each.',
 'Paired scenarios share site effects and measurements where possible. Confirmation uses disjoint seeds; no settings or methods are selected using discovery results.',
 'Candidate intervals are percentile, basic reflected log-scale, and centered normal log-scale using bootstrap sd. Studentized intervals are evaluated only for the uncensored benchmark with known site-level variance estimates.',
 'Complete-data oracle intervals see hidden uncensored measurements and known covariance; they are diagnostics, not usable methods for censored real data.',
 'All methods and scenarios are fixed before inspecting results. Checkpoints must not be reused for changed code/design.'
),file.path(folder,'design_notes.txt'))
writeLines(capture.output(sessionInfo()),file.path(folder,'session_info.txt'))
# First establish the exact independent Gaussian block-mean benchmark.
for (phase in c('discovery','confirmation')) {
 path <- file.path(folder,paste0('analytic_',phase,'.rds'))
 if(file.exists(path)) next
 n <- if(phase=='discovery') 2000L else 3000L
 seed_base <- if(phase=='discovery') 2100000L else 3100000L
 rows <- vector('list',n)
 for(i in seq_len(n)) {
  set.seed(seed_base+i)
  site <- rnorm(k,mu,sqrt(.55^2+.2^2/4))
  theta<-mean(site); se<-sd(site)/sqrt(k)
  sampled<-matrix(site[sample.int(k,k*999,replace=TRUE)],nrow=k)
  z<-colMeans(sampled)
  boot_variance<-(colSums(sampled^2)-k*z^2)/(k-1)/k
  pivot<-(z-theta)/sqrt(pmax(boot_variance,.Machine$double.eps))
  ci<-intervals(theta,z)
  ci$studentized_log<-theta-rev(quantile(pivot,c(.025,.975),names=FALSE))*se
  ci$analytic_t<-theta+c(-1,1)*qt(.975,k-1)*se
  ci$known_variance_z<-theta+c(-1,1)*qnorm(.975)*sqrt((.55^2+.2^2/4)/k)
  rows[[i]]<-dplyr::bind_rows(lapply(names(ci),function(method)
   record(phase,i,'uncensored_regular',method,theta,ci[[method]],NA_real_,sd(z),999)))
 }
 saveRDS(dplyr::bind_rows(rows),path)
 cat('Analytic benchmark complete:',phase,'\n');flush.console()
}
grid<-expand.grid(x=0:7,y=0:7)
covariance<-.55^2*exp(-as.matrix(dist(grid))/.8)
root<-chol(covariance+diag(1e-10,k))
settings<-data.frame(scenario=c('censored_regular','censored_uneven','correlated_blocks_1','correlated_blocks_2'),
 uneven=c(FALSE,TRUE,TRUE,TRUE),correlated=c(FALSE,FALSE,TRUE,TRUE),block_size=c(1,1,1,2))
readr::write_csv(settings,file.path(folder,'settings.csv'))
for(phase in c('discovery','confirmation')) {
 seed_base<-if(phase=='discovery') 4100000L else 5100000L
 for(i in 1:200) {
  path<-file.path(folder,sprintf('diagnostic_%s_%03d.rds',phase,i))
  if(file.exists(path)) next
  results<-list()
  for(j in seq_len(nrow(settings))) {
   setting<-settings[j,]
   set.seed(seed_base+i)
   latent<-if(setting$correlated)as.vector(t(root)%*%rnorm(k))else rnorm(k,sd=.55)
   # Use common noise arrays before varying density to keep paired comparisons.
   noise<-matrix(rnorm(k*12,sd=.2),nrow=k)
   counts<-if(setting$uneven)sample(c(1,4,12),k,replace=TRUE)else rep(4,k)
   id<-rep(seq_len(k),counts)
   values<-unlist(lapply(seq_len(k),function(h)mu+latent[h]+noise[h,seq_len(counts[h])]))
   finite_target<-mu+mean(latent)
   site_means<-as.numeric(tapply(values,id,mean))
   known_var<-if(setting$correlated)sum(covariance)/k^2 else .55^2/k
   known_var<-known_var+.2^2*sum(1/counts)/k^2
   oracle<-mean(site_means)+c(-1,1)*qnorm(.975)*sqrt(known_var)
   results[[length(results)+1]]<-record(phase,i,setting$scenario,'complete_data_oracle',mean(site_means),oracle,finite_target,sqrt(known_var))
   value<-exp(values);lod<-exp(mu-.3);upper<-exp(mu+1)
   x<-data.frame(x=grid$x[id]+.2,y=grid$y[id]+.2,
    detected_val=ifelse(value<lod|value>upper,NA_real_,value),
    left_cen=value<lod,left_limit=ifelse(value<lod,lod,NA_real_),
    right_cen=value>upper,right_limit=ifelse(value>upper,upper,NA_real_))
   b<-spatial_block_bootstrap(x,'x','y',1,setting$block_size,n_boot=199,n_offsets=2,seed=seed_base+20000+i)
   valid<-b$summary$status=='experimental_interval'
   z<-log(b$replicates$geomean[is.finite(b$replicates$geomean)&b$replicates$geomean>0])
   theta<-log(b$summary$geomean)
   ci<-if(valid)intervals(theta,z)else list(percentile_log=c(NA,NA),basic_log=c(NA,NA),normal_log_centered=c(NA,NA))
   for(method in names(ci))results[[length(results)+1]]<-record(phase,i,setting$scenario,method,theta,ci[[method]],finite_target,if(length(z)>1)sd(z)else NA_real_,b$summary$n_success)
  }
  saveRDS(dplyr::bind_rows(results),path)
  if(i%%20==0){cat('Censored diagnostics',phase,i,'/ 200 complete\n');flush.console()}
 }
}
files<-c(list.files(folder,pattern='^analytic_.*rds$',full.names=TRUE),list.files(folder,pattern='^diagnostic_.*rds$',full.names=TRUE))
raw<-dplyr::bind_rows(lapply(files,readRDS))
readr::write_csv(raw,file.path(folder,'results.csv'))
summary<-raw%>%dplyr::group_by(.data$phase,.data$scenario,.data$method)%>%
 dplyr::group_modify(function(x,key){
  valid<-is.finite(x$lo)&is.finite(x$hi);n<-sum(valid);hit<-sum(x$covered[valid])
  mc<-if(n)binom.test(hit,n)$conf.int else c(NA,NA)
  data.frame(n_datasets=nrow(x),n_intervals=n,coverage=if(n)hit/n else NA_real_,mc_low=mc[1],mc_high=mc[2],
   bias_log=mean(x$theta-mu),empirical_sd_log=sd(x$theta),mean_boot_se=mean(x$boot_se,na.rm=TRUE),
   mean_log_width=mean(x$hi[valid]-x$lo[valid]),
   finite_domain_coverage=if(all(is.na(x$finite_target)))NA_real_ else mean(x$finite_covered[valid]))
 })%>%dplyr::ungroup()
readr::write_csv(summary,file.path(folder,'summary.csv'))
writeLines('Benchmark and confirmation diagnostics complete.',file.path(folder,'COMPLETED.txt'))
print(as.data.frame(summary))
