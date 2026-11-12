#!/usr/bin/env Rscript
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
if(!file.exists(script)) script<-gsub('~+~',' ',script,fixed=TRUE)
here <- dirname(normalizePath(script,mustWork=TRUE));source(file.path(here,'functions.R'))
hf_require(c('geoCensor','sf','survival','digest','jsonlite'))
x <- hf_read(here);d <- hf_observations(x$data);m <- x$manifest
stopifnot(!'arrow'%in%loadedNamespaces())
stopifnot(nrow(d)==m$expected$full_reports,sum(d$in_CONUS_map)==m$expected$CONUS_reports)
for(k in names(m$expected$full_methods)) stopifnot(sum(d$method==k)==m$expected$full_methods[[k]])
for(k in names(m$expected$full_censoring))
  stopifnot(sum(x$data$chemistry$parse_status==k)==m$expected$full_censoring[[k]])
# Right-censor numeric markers are checked against the preserved source values.
c <- x$data$chemistry;ii <- which(c$parse_status=='right');q<-c$qualified_value[ii]
decoded <- ifelse(abs(q-.51111)<1e-12,.5,ifelse(abs(q-.25111)<1e-12,.25,floor(q)))
stopifnot(all(abs(q-floor(q)-.11111)<1e-9 | abs(q-.51111)<1e-12 | abs(q-.25111)<1e-12),
          identical(decoded,c$reporting_lower_bound[ii]))
projected <- hf_project(d,m)
hf_validate_membership(projected,x$reference,m)
for(medium in names(x$reference)) {
  a <- if(medium=='pooled') projected else projected[projected$medium==medium,,drop=FALSE]
  b <- x$reference[[medium]];key <- paste(a$ix,a$iy,sep=':');tab <- table(key)
  j <- match(paste(b$ix,b$iy,sep=':'),names(tab));e<-m$expected[[medium]]
  stopifnot(nrow(a)==e$reports,nrow(b)==e$occupied,!anyNA(j),length(tab)==nrow(b),
            all(as.numeric(tab[j])==b$n_reports),sum(is.finite(b$GM_mgkg))==e$coloured,
            sum(!is.finite(b$GM_mgkg))==e$hollow)
}
# Hand-check a likelihood derivative independently with central differences.
z <- log(c(.03,.05,.08,.12,.2,.4));k <- c(1L,0L,0L,0L,2L,2L);theta<-c(-2,.1);eps<-1e-5
ev<-hf_evaluate(theta,z,k)
g<-vapply(1:2,function(j){step<-c(0,0);step[j]<-eps
 (hf_evaluate(theta+step,z,k)$nll-hf_evaluate(theta-step,z,k)$nll)/(2*eps)},numeric(1))
H<-sapply(1:2,function(j){step<-c(0,0);step[j]<-eps
 (hf_evaluate(theta+step,z,k)$score-hf_evaluate(theta-step,z,k)$score)/(2*eps)})
stopifnot(max(abs(ev$score-g))<1e-7,max(abs(ev$hessian-H))<1e-7)
synthetic <- data.frame(ix=0,iy=0,longitude=seq(-110,-109.95,length.out=6),
  latitude=40,method='MS_HF',kind=0L,value_or_actual_bound_mgkg=c(.02,.04,.08,.16,.32,.64))
got <- hf_fit_cell(synthetic,m)
stopifnot(abs(got$GM_mgkg-exp(mean(log(synthetic$value_or_actual_bound_mgkg))))<1e-12,
          abs(got$api_GM_mgkg/got$GM_mgkg-1)<1e-6,!any(grepl('ci_low|ci_high',names(got))))
synthetic$kind <- 1L
stopifnot(is.na(hf_fit_cell(synthetic,m)$GM_mgkg))
synthetic$kind <- c(1L,1L,1L,2L,2L,2L)
stopifnot(is.na(hf_fit_cell(synthetic,m)$GM_mgkg))
synthetic$kind <- 0L;synthetic$value_or_actual_bound_mgkg <- .1
stopifnot(is.na(hf_fit_cell(synthetic,m)$GM_mgkg))
synthetic$value_or_actual_bound_mgkg <- c(.02,.04,.08,.16,.32,.64);synthetic$longitude <- -110
stopifnot(is.na(hf_fit_cell(synthetic,m)$GM_mgkg))
synthetic$longitude <- seq(-110,-109.95,length.out=6);synthetic$kind <- c(1L,0L,0L,0L,0L,2L)
mixed <- hf_fit_cell(synthetic,m)
stopifnot(is.finite(mixed$GM_mgkg),mixed$max_abs_score<=m$score_tolerance,
          mixed$min_curvature_eigenvalue>0,abs(mixed$api_GM_mgkg/mixed$GM_mgkg-1)<1e-5)
# Floor, rather than truncation, is required west/south of the zero origin.
stopifnot(identical(floor(c(-20001,-20000,-1,0,19999,20000)/20000),c(-2,-1,-1,0,0,1)))
bad <- x$data;bad$chemistry$qualifier[which(bad$chemistry$parse_status=='numerical')[1]] <- '<'
stopifnot(inherits(tryCatch(hf_observations(bad),error=identity),'error'))
stopifnot(identical(hf_colours(c(.0001,.05)),rep(hf_colours(.05),2)),
          identical(hf_colours(c(.3,100)),rep(hf_colours(.3),2)))
cat('PASS: bundled hashes, raw qualifiers/bounds, source joins, all report/cell counts, projection,\n',
    'likelihood derivatives, exact/censored fits, unavailable cases and colour tails.\n',sep='')
