#!/usr/bin/env Rscript
# Independent direct-survival/character-grid checks on a prespecified real-data panel.
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
out <- normalizePath(args[1],mustWork=TRUE)
config <- read.csv(file.path(out,'config.csv'),colClasses='character')
value <- function(key)config$value[match(key,config$key)]
.libPaths(c(value('library_dir'),.libPaths()))
suppressPackageStartupMessages({library(geoCensor);library(data.table);library(survival)})
setDTthreads(2L)
source_manifest <- read.csv(file.path(out,'source_manifest.csv'),colClasses='character')
stopifnot(identical(unname(tools::md5sum(source_manifest$file)),source_manifest$md5))
widths <- as.numeric(strsplit(value('cell_sizes_m'),',',fixed=TRUE)[[1]])
n_offsets <- as.integer(value('n_offsets'));thresholds <- 10^seq(-2,9,by=.5)
panel <- c('Au','Al','Ag','Cu','Pb','Zn')
intervals_reference <- function(d) {
 left <- !is.na(d$left_cen)&d$left_cen;right <- !is.na(d$right_cen)&d$right_cen
 a <- left&is.finite(d$left_limit)&d$left_limit>0
 b <- right&is.finite(d$right_limit)&d$right_limit>0
 c <- !left&!right&is.finite(d$detected_val)&d$detected_val>0
 lower <- upper <- rep(NA_real_,nrow(d))
 lower[a] <- 0;upper[a] <- d$left_limit[a]
 lower[b] <- d$right_limit[b];upper[b] <- Inf
 lower[c] <- upper[c] <- d$detected_val[c]
 list(ok=a|b|c,lower=lower,upper=upper)
}
fit_reference <- function(d,weights=NULL) {
 iv <- intervals_reference(d);ok <- iv$ok
 if(!is.null(weights))ok <- ok&weights>0
 n <- sum(ok);n_exact <- sum(iv$lower[ok]==iv$upper[ok])
 if(n<6||n_exact<2)return(list(n=n,mean=NA_real_,ci_low=NA_real_,ci_high=NA_real_,status='insufficient'))
 w <- if(is.null(weights))rep(1,n) else weights[ok]
 w <- w/max(w);w <- w/sum(w)*n
 fit <- tryCatch(withCallingHandlers({
  model <- survival::survreg(survival::Surv(iv$lower[ok],iv$upper[ok],type='interval2')~1,
   dist='lognormal',weights=w)
  covariance <- stats::vcov(model)
  valid <- all(is.finite(stats::coef(model)))&&is.finite(model$scale)&&model$scale>0&&
    all(is.finite(covariance))&&covariance[1,1]>0
  if(valid)model else NULL
 },warning=function(w)stop(conditionMessage(w),call.=FALSE)),error=function(e)NULL)
 if(is.null(fit))return(list(n=n,mean=NA_real_,ci_low=NA_real_,ci_high=NA_real_,status='failed'))
 mu <- unname(stats::coef(fit)[1]);se <- sqrt(stats::vcov(fit)[1,1]);z <- stats::qnorm(.975)
 list(n=n,mean=exp(mu),ci_low=exp(mu-z*se),ci_high=exp(mu+z*se),status='fitted')
}
weights_reference <- function(d,size) {
 iv <- intervals_reference(d);spatial <- iv$ok&is.finite(d$x_m)&is.finite(d$y_m)
 w <- numeric(nrow(d));if(!any(spatial))return(w)
 x <- (d$x_m[spatial]-min(d$x_m[spatial]))/size
 y <- (d$y_m[spatial]-min(d$y_m[spatial]))/size
 total <- numeric(length(x))
 for(i in 0:(n_offsets-1L))for(j in 0:(n_offsets-1L)) {
  cells <- paste(sprintf('%.0f',floor(x+i/n_offsets)),sprintf('%.0f',floor(y+j/n_offsets)),sep='|')
  counts <- table(cells)
  total <- total+1/(length(counts)*as.numeric(counts[cells]))
 }
 w[spatial] <- total/sum(total);w
}
select_domain <- function(d,row) {
 keep <- d$state==row$state
 keep <- keep&if(is.na(row$method_collected))is.na(d$method_collected) else
  !is.na(d$method_collected)&d$method_collected==row$method_collected
 d[keep,,drop=FALSE]
}
models <- analytic <- cdfs <- list()
check_fit <- function(row,d,stage,medium,scope,weights=NULL,cell_size=NA_real_) {
 ref <- fit_reference(d,weights)
 stopifnot(row$n_total==ref$n,is.finite(row$geomean)==is.finite(ref$mean))
 error <- NA_real_
 if(is.finite(ref$mean)) {
  error <- abs(log(row$geomean/ref$mean));stopifnot(error<=1e-6)
  if(is.null(weights))stopifnot(abs(log(row$ci_low/ref$ci_low))<=1e-6,
   abs(log(row$ci_high/ref$ci_high))<=1e-6)
 }
 models[[length(models)+1L]] <<- data.frame(stage=stage,medium=medium,element=row$element,
  collection=if(is.na(row$method_collected))NA_character_ else row$method_collected,
  scope=scope,cell_size=cell_size,n_total=ref$n,reference_status=ref$status,
  log_relative_geomean_error=error)
}
check_cdf <- function(d,w,stage,medium,element,scope,cell_size=NA_real_) {
 iv <- intervals_reference(d);ok <- iv$ok
 data <- d
 if(is.null(w)) {weights <- rep(1,nrow(d));name <- NULL} else {
  weights <- w;data$.reference_weight <- w;name <- '.reference_weight'
 }
 ok <- ok&weights>0
 observed <- geoCensor::censored_cdf_bounds(data,thresholds,weight_col=name)$cdf
 if(!any(ok)) {
  stopifnot(all(is.na(observed$cdf_lower)),all(is.na(observed$cdf_upper)))
  error <- NA_real_
 } else {
  weights <- weights[ok]/max(weights[ok]);weights <- weights/sum(weights)
  lo <- iv$lower[ok];hi <- iv$upper[ok];exact <- lo==hi;left <- lo==0;right <- is.infinite(hi)
  low <- vapply(thresholds,function(t)sum(weights[exact&lo<=t])+sum(weights[left&hi<=t]),0.0)
  high <- vapply(thresholds,function(t)sum(weights[exact&lo<=t])+sum(weights[left])+sum(weights[right&lo<t]),0.0)
  error <- max(abs(observed$cdf_lower-low),abs(observed$cdf_upper-high));stopifnot(error<=1e-10)
 }
 cdfs[[length(cdfs)+1L]] <<- data.frame(stage=stage,medium=medium,element=element,
  scope=scope,cell_size=cell_size,n_thresholds=length(thresholds),max_absolute_bound_error=error)
}
for(stage in c('pilot','expansion'))for(medium in c('soil','sed')) {
 if(stage=='pilot'&&medium=='sed')next
 folder <- file.path(out,'results',stage,medium)
 stopifnot(file.exists(file.path(folder,'COMPLETED.txt')))
 paths <- sort(list.files(file.path(folder,'checkpoints'),pattern='[0-9]{3}\\.rds$',full.names=TRUE))
 for(path in paths) {
  checkpoint <- readRDS(path);d <- checkpoint$observations;result <- checkpoint$result
  if(!nrow(d))next
  for(i in seq_len(nrow(result$estimates))) {
   row <- result$estimates[i,];domain <- select_domain(d,row);iv <- intervals_reference(domain)
   exact <- iv$ok&iv$lower==iv$upper;exact[is.na(exact)] <- FALSE
   if(sum(iv$ok)>=6&&all(exact[iv$ok])&&is.finite(row$geomean)) {
    target <- exp(mean(log(domain$detected_val[iv$ok])))
    error <- abs(log(row$geomean/target));stopifnot(error<=1e-6)
    analytic[[length(analytic)+1L]] <- data.frame(stage=stage,medium=medium,element=checkpoint$element,
     scope='unweighted_complete_group',n_total=sum(iv$ok),log_relative_geomean_error=error)
   }
   if(checkpoint$element%in%panel) {
    check_fit(row,domain,stage,medium,'unweighted_group')
    check_cdf(domain,NULL,stage,medium,checkpoint$element,'unweighted_group')
   }
  }
  if(checkpoint$element%in%panel) {
   strata <- result$method_estimates
   for(i in seq_len(nrow(strata))) {
    row <- strata[i,];domain <- select_domain(d,row)
    code <- toupper(trimws(domain$amethod));code[!is.na(code)&code==''] <- NA_character_
    take <- if(is.na(row$method_code))is.na(code) else !is.na(code)&code==row$method_code
    check_fit(row,domain[take,,drop=FALSE],stage,medium,'method_stratum')
   }
   weighted <- result$declustering
   for(i in seq_len(nrow(weighted))) {
    row <- weighted[i,];domain <- select_domain(d,row);w <- weights_reference(domain,row$cell_size)
    iv <- intervals_reference(domain);use <- iv$ok&is.finite(domain$x_m)&is.finite(domain$y_m)
    direct <- geoCensor::cell_declustering_weights(domain[use,,drop=FALSE],'x_m','y_m',row$cell_size,n_offsets)
    if(any(use))stopifnot(max(abs(w[use]-direct))<=1e-12)
    check_fit(row,domain,stage,medium,'declustered_group',w,row$cell_size)
    check_cdf(domain,w,stage,medium,checkpoint$element,'declustered_group',row$cell_size)
    stopifnot(is.na(row$ci_low),is.na(row$ci_high))
   }
  }
 }
 cat(stage,medium,'direct model/weight/CDF panel checked\n');flush.console()
}
folder <- file.path(out,'verification');dir.create(folder,showWarnings=FALSE)
if(length(models))fwrite(rbindlist(models),file.path(folder,'direct_survival_panel.csv'),na='NA')
if(length(analytic))fwrite(rbindlist(analytic),file.path(folder,'analytic_uncensored_checks.csv'),na='NA')
if(length(cdfs))fwrite(rbindlist(cdfs),file.path(folder,'independent_cdf_panel.csv'),na='NA')
stopifnot(identical(unname(tools::md5sum(source_manifest$file)),source_manifest$md5))
writeLines(c('Direct survival/independent character-grid weight/CDF calculations reconcile to the prespecified panel.',
 'Complete uncensored group means match the analytic lognormal MLE where the public fit is available.',
 'Numerical reconciliation does not validate reporting precision, homogeneity, independence or interval coverage.'),
 file.path(folder,'R_VERIFIED.txt'))
