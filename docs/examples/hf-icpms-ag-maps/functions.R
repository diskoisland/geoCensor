# Example helpers, MIT license; selected USGS data are CC0-1.0.
# This standalone repository example is not an installed geoCensor command.

hf_require <- function(packages) {
  missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) stop('Install required packages: ',paste(missing,collapse=', '),call.=FALSE)
}

hf_read <- function(directory) {
  hf_require(c('jsonlite','digest'))
  manifest <- jsonlite::read_json(file.path(directory,'manifest.json'),simplifyVector=TRUE)
  for(i in seq_len(nrow(manifest$files))) {
    f <- file.path(directory,manifest$files$file[i])
    if(!file.exists(f) || digest::digest(file=f,algo='sha256')!=manifest$files$sha256[i])
      stop('Missing or changed bundled input: ',basename(f),call.=FALSE)
  }
  list(data=readRDS(file.path(directory,'data.RDS')),
       reference=readRDS(file.path(directory,'reference-cells.RDS')),manifest=manifest)
}

hf_observations <- function(data) {
  d <- data$chemistry; m <- data$metadata
  if(anyDuplicated(d$national_chem_row_ordinal_1based) || anyDuplicated(m$national_geol_row_ordinal_1based))
    stop('Duplicate source row identities',call.=FALSE)
  j <- match(d$metadata_geol_row_ordinal_1based,m$national_geol_row_ordinal_1based)
  if(anyNA(j) || !identical(d$lab_id,m$lab_id[j]) || !identical(d$group_id,m$group_id[j]) ||
     !identical(d$medium_class,m$primary_class[j]) || !identical(d$state_metadata,m$state[j]))
    stop('Source chemistry and metadata disagree',call.=FALSE)
  if(!all(d$species=='Ag' & d$units=='ppm' & d$analytic_method %in% c('MS_HF','MS_HF_UT','MS_HFST')) ||
     !all(d$medium_class %in% c('soil','sediment'))) stop('Unexpected cohort member',call.=FALSE)
  k <- match(d$parse_status,c('numerical','left','right'))-1L
  v <- ifelse(k==0L,d$result_value,ifelse(k==1L,d$reporting_upper_bound,d$reporting_lower_bound))
  if(anyNA(k) || any(!is.finite(v) | v<=0)) stop('Invalid observation or actual censor bound',call.=FALSE)
  blank <- is.na(d$qualifier) | d$qualifier==''
  if(any(k==0L & (!blank | d$qualified_value!=d$result_value)) ||
     any(k==1L & (!d$qualifier %in% c('N','L','<') | d$qualified_value>=0 | abs(d$qualified_value)!=v)) ||
     any(k==2L & !d$qualifier %in% c('G','>'))) stop('Qualifier interpretation disagreement',call.=FALSE)
  lon <- m$longitude[j]; lat <- m$latitude[j]
  located <- is.finite(lon)&is.finite(lat)&abs(lon)<=180&abs(lat)<=90&!(lon==0&lat==0)
  conus <- d$state_metadata %in% c(setdiff(state.abb,c('AK','HI')),'DC') & located
  data.frame(source_chem_row=d$national_chem_row_ordinal_1based,
    source_metadata_row=d$metadata_geol_row_ordinal_1based,lab_id=d$lab_id,group_id=d$group_id,
    medium=d$medium_class,state=d$state_metadata,method=d$analytic_method,source_qualifier=d$qualifier,
    kind=k,value_or_actual_bound_mgkg=v,longitude=lon,latitude=lat,in_CONUS_map=conus)
}

hf_project <- function(d,manifest) {
  hf_require('sf')
  d <- d[d$in_CONUS_map,,drop=FALSE]
  xy <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(d,coords=c('longitude','latitude'),
           crs='OGC:CRS84'),crs=manifest$crs))
  if(any(!is.finite(xy))) stop('Nonfinite projected coordinates')
  d$x_m <- xy[,1];d$y_m <- xy[,2]
  d$ix <- floor(d$x_m/manifest$cell_width_m);d$iy <- floor(d$y_m/manifest$cell_width_m)
  d
}

# Negative log likelihood, score and curvature per report in (mu, log sigma).
hf_evaluate <- function(theta,z,k) {
  if(any(!is.finite(theta)) || abs(theta[1])>200 || theta[2]< -12 || theta[2]>7)
    return(list(nll=Inf,score=c(Inf,Inf),hessian=matrix(NA_real_,2,2)))
  s <- exp(theta[2]); t <- (z-theta[1])/s; point <- k==0L; tp <- t[point]
  sign <- ifelse(k[!point]==1L,1,-1);q <- sign*t[!point]
  lam <- exp(stats::dnorm(q,log=TRUE)-stats::pnorm(q,log.p=TRUE));v <- q+lam
  n <- length(z)
  nll <- -(sum(-theta[2]-z[point]-.5*tp^2-.5*log(2*pi))+sum(stats::pnorm(q,log.p=TRUE)))/n
  score <- c(sum(-tp/s)+sum(sign*lam/s),sum(1-tp^2)+sum(q*lam))/n
  h00 <- length(tp)/s^2+sum(lam*v/s^2)
  h01 <- sum(2*tp/s)+sum(sign*lam*(q*v-1)/s)
  h11 <- sum(2*tp^2)+sum(lam*(q*q*v-q))
  list(nll=nll,score=score,hessian=matrix(c(h00,h01,h01,h11),2,2)/n)
}

hf_fit_cell <- function(d,manifest) {
  v <- d$value_or_actual_bound_mgkg;k <- d$kind;z <- log(v);point <- k==0L
  coords <- nrow(unique(d[,c('longitude','latitude')]))
  distinct <- length(unique(v[point]));n <- nrow(d)
  row <- data.frame(ix=d$ix[1],iy=d$iy[1],x_center_m=(d$ix[1]+.5)*manifest$cell_width_m,
    y_center_m=(d$iy[1]+.5)*manifest$cell_width_m,n_reports=n,n_numerical=sum(point),
    n_left=sum(k==1L),n_right=sum(k==2L),n_coordinate_pairs=coords,
    n_distinct_numerical_values=distinct,GM_mgkg=NA_real_,mu_log=NA_real_,sigma_log=NA_real_,
    normalized_nll=NA_real_,max_abs_score=NA_real_,min_curvature_eigenvalue=NA_real_,
    multistart_objective_gap=NA_real_,api_GM_mgkg=NA_real_,api_method=NA_character_,
    status=if(!any(point)) 'unavailable_no_numerical' else 'unavailable_display_support')
  if(n<manifest$supports$reports || distinct<manifest$supports$numerical_values ||
     coords<manifest$supports$coordinate_pairs) return(row)
  # Public geoCensor API on explicit intervals in mg/kg. Unit weights omit CIs.
  api_data <- data.frame(detected_val=ifelse(point,v,NA_real_),left_cen=k==1L,
    left_limit=ifelse(k==1L,v,NA_real_),right_cen=k==2L,right_limit=ifelse(k==2L,v,NA_real_),
    amethod=d$method,equal_report_weight=1)
  api <- geoCensor::fit_cen_lognorm(api_data,weight_col='equal_report_weight')
  row$api_GM_mgkg <- api$geomean;row$api_method <- api$method
  # Example-level refinement recovers narrow-cluster starts that can fail in
  # the released API. It does not change geoCensor or substitute for censored data.
  candidates <- list()
  if(all(point)) candidates[[1]] <- c(mean(z),log(sqrt(mean((z-mean(z))^2)))) else {
    mu0 <- mean(z[point]);s0 <- sqrt(mean((z[point]-mu0)^2))
    fd <- data.frame(lo=ifelse(k==1L,NA_real_,v),hi=ifelse(k==2L,NA_real_,v))
    for(s in unique(c(s0,.1,1,3))) {
      fit <- tryCatch(withCallingHandlers(survival::survreg(survival::Surv(lo,hi,type='interval2')~1,
        data=fd,dist='lognormal',robust=FALSE,init=c(mu0,log(s)),
        control=survival::survreg.control(maxiter=1000,rel.tolerance=1e-12)),
        warning=function(w)stop(conditionMessage(w),call.=FALSE)),error=function(e)NULL)
      if(!is.null(fit) && is.finite(fit$scale) && fit$scale>0) {
        th <- c(unname(stats::coef(fit)[1]),log(fit$scale));ev <- hf_evaluate(th,z,k)
        if(is.finite(ev$nll) && max(abs(ev$score))>1e-8) {
          polished <- tryCatch(stats::optim(th,function(t)hf_evaluate(t,z,k)$nll,
            gr=function(t)hf_evaluate(t,z,k)$score,method='BFGS',
            control=list(maxit=600,reltol=1e-14)),error=function(e)NULL)
          if(!is.null(polished) && is.finite(polished$value) && polished$value<=ev$nll+1e-12) th <- polished$par
        }
        candidates[[length(candidates)+1L]] <- th
      }
    }
  }
  row$status <- 'unavailable_numerical_guard'
  if(!length(candidates)) return(row)
  evaluations <- lapply(candidates,hf_evaluate,z=z,k=k)
  valid <- vapply(evaluations,function(e) is.finite(e$nll)&&all(is.finite(e$score))&&
    all(is.finite(e$hessian))&&max(abs(e$score))<=manifest$score_tolerance&&
    min(eigen(e$hessian,symmetric=TRUE,only.values=TRUE)$values)>0,logical(1))
  if(!any(valid)) return(row)
  obj <- vapply(evaluations,function(e)e$nll,numeric(1));j <- which(valid)[which.min(obj[valid])]
  ev <- evaluations[[j]];th <- candidates[[j]];gap <- max(0,ev$nll-min(obj))
  if(gap>manifest$objective_gap_tolerance) return(row)
  row$GM_mgkg <- exp(th[1]);row$mu_log <- th[1];row$sigma_log <- exp(th[2])
  row$normalized_nll <- ev$nll;row$max_abs_score <- max(abs(ev$score))
  row$min_curvature_eigenvalue <- min(eigen(ev$hessian,symmetric=TRUE,only.values=TRUE)$values)
  row$multistart_objective_gap <- gap;row$status <- 'conditional_local_lognormal_GM'
  row
}

hf_refit <- function(d,manifest) {
  hf_require(c('geoCensor','survival'))
  results <- list()
  for(medium in c('pooled','soil','sediment')) {
    x <- if(medium=='pooled') d else d[d$medium==medium,,drop=FALSE]
    groups <- split(seq_len(nrow(x)),paste(x$ix,x$iy,sep=':'))
    rows <- vector('list',length(groups))
    for(i in seq_along(groups)) {
      rows[[i]] <- hf_fit_cell(x[groups[[i]],,drop=FALSE],manifest)
      if(i%%1000L==0L) message(medium,': ',i,' / ',length(groups),' occupied cells')
    }
    z <- do.call(rbind,rows);z <- z[order(z$ix,z$iy),,drop=FALSE];row.names(z)<-NULL
    results[[medium]] <- z
    message(medium,': ',sum(is.finite(z$GM_mgkg)),' displayed means')
  }
  results
}

hf_compare <- function(cells,reference,manifest) {
  metrics <- list()
  for(medium in names(reference)) {
    a <- cells[[medium]];b <- reference[[medium]]
    key <- function(x)paste(x$ix,x$iy,sep=':')
    j <- match(key(b),key(a));if(anyNA(j)||nrow(a)!=nrow(b)) stop('Changed grid cells: ',medium)
    a <- a[j,,drop=FALSE]
    for(f in c('ix','iy','n_reports','n_numerical','n_left','n_right','n_coordinate_pairs','n_distinct_numerical_values'))
      if(!identical(as.numeric(a[[f]]),as.numeric(b[[f]]))) stop('Cell count disagreement: ',medium,' / ',f)
    good <- is.finite(b$GM_mgkg)
    if(!identical(is.finite(a$GM_mgkg),good)) stop('Displayed-cell availability changed: ',medium)
    error <- max(abs(a$GM_mgkg[good]/b$GM_mgkg[good]-1))
    if(error>manifest$reference_relative_tolerance) stop('Cell GM mismatch: ',medium)
    metrics[[medium]] <- list(occupied=nrow(a),coloured=sum(good),hollow=sum(!good),
      max_GM_relative_error=error)
  }
  metrics
}

hf_validate_membership <- function(d,reference,manifest) {
  for(medium in names(reference)) {
    x <- if(medium=='pooled') d else d[d$medium==medium,,drop=FALSE]
    b <- reference[[medium]];key <- paste(x$ix,x$iy,sep=':')
    groups <- split(seq_len(nrow(x)),key);expected <- manifest$expected[[medium]]
    ids <- paste(b$ix,b$iy,sep=':')
    if(length(groups)!=nrow(b) || !setequal(ids,names(groups)) ||
       nrow(x)!=expected$reports || nrow(b)!=expected$occupied ||
       sum(is.finite(b$GM_mgkg))!=expected$coloured) stop('Cohort/grid count disagreement: ',medium)
    for(j in seq_len(nrow(b))) {
      a <- x[groups[[ids[j]]],,drop=FALSE];p <- a$kind==0L
      values <- c(nrow(a),sum(p),sum(a$kind==1L),sum(a$kind==2L),
        nrow(unique(a[,c('longitude','latitude')])),length(unique(a$value_or_actual_bound_mgkg[p])))
      fields <- c('n_reports','n_numerical','n_left','n_right','n_coordinate_pairs','n_distinct_numerical_values')
      if(!identical(as.numeric(values),as.numeric(unlist(b[j,fields],use.names=FALSE)))) stop('Source/cell count disagreement: ',medium)
      if(b$x_center_m[j]!=(b$ix[j]+.5)*manifest$cell_width_m ||
         b$y_center_m[j]!=(b$iy[j]+.5)*manifest$cell_width_m) stop('Changed cell centre')
    }
  }
  invisible(TRUE)
}

hf_boundaries <- function(manifest) {
  hf_require(c('sf','maps'))
  m <- maps::map('state',fill=TRUE,plot=FALSE,resolution=0)
  ok <- is.finite(m$x)&is.finite(m$y)
  xy <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(data.frame(lon=m$x[ok],lat=m$y[ok]),
    coords=c('lon','lat'),crs='OGC:CRS84'),crs=manifest$crs))
  xx <- yy <- rep(NA_real_,length(m$x));xx[ok]<-xy[,1];yy[ok]<-xy[,2]
  list(x=xx,y=yy)
}

hf_colours <- function(values) {
  # Same ColorBrewer RdYlBu_r anchors as the earlier maps, in sRGB.
  anchors <- rev(c('#A50026','#D73027','#F46D43','#FDAE61','#FEE090','#FFFFBF',
                   '#E0F3F8','#ABD9E9','#74ADD1','#4575B4','#313695'))
  palette <- grDevices::colorRampPalette(anchors,space='rgb')(256)
  t <- pmin(1,pmax(0,log(values/.05)/log(.30/.05)))
  palette[pmin(256L,1L+floor(t*256))]
}

hf_plot <- function(cells,manifest,borders,medium,western,file) {
  grDevices::png(file,width=if(western)2400 else 3000,height=if(western)2200 else 2000,res=200,
                 type=if(capabilities('cairo'))'cairo' else getOption('bitmapType'))
  on.exit(grDevices::dev.off(),add=TRUE)
  graphics::par(mar=c(0,0,0,0),oma=c(3,1,5,1),family='sans',xaxs='i',yaxs='i')
  graphics::layout(matrix(c(1,2),nrow=1),widths=c(9,1.25))
  if(western) {xr<-c(-2400000,-200000);yr<-c(350000,3170000)} else {
    xr<-range(borders$x,na.rm=TRUE)+c(-45000,45000);yr<-range(borders$y,na.rm=TRUE)+c(-45000,45000)
  }
  graphics::plot.new();graphics::plot.window(xlim=xr,ylim=yr,asp=1)
  graphics::clip(xr[1],xr[2],yr[1],yr[2])
  good <- is.finite(cells$GM_mgkg)&cells$GM_mgkg>0
  if(any(!good)) graphics::symbols(cells$x_center_m[!good],cells$y_center_m[!good],
    circles=rep(manifest$circle_radius_m,sum(!good)),inches=FALSE,add=TRUE,fg='#A8B2BC',bg=NA,lwd=.45)
  if(any(good)) graphics::symbols(cells$x_center_m[good],cells$y_center_m[good],
    circles=rep(manifest$circle_radius_m,sum(good)),inches=FALSE,add=TRUE,fg=NA,bg=hf_colours(cells$GM_mgkg[good]))
  graphics::lines(borders$x,borders$y,col='#3E4852',lwd=.65)
  graphics::plot.new();graphics::plot.window(xlim=c(0,1),ylim=c(0,1))
  edge <- seq(.17,.80,length.out=257);v <- exp(seq(log(.05),log(.30),length.out=256))
  graphics::rect(.12,edge[-257],.31,edge[-1],col=hf_colours(v),border=NA)
  graphics::polygon(c(.12,.215,.31),c(.80,.82,.80),col=hf_colours(.30),border=NA)
  graphics::polygon(c(.12,.215,.31),c(.17,.15,.17),col=hf_colours(.05),border=NA)
  ticks <- c(.05,.075,.10,.15,.20,.30);yy <- .17+.63*log(ticks/.05)/log(.30/.05)
  graphics::segments(.31,yy,.35,yy,col='#333333');graphics::text(.39,yy,
    labels=c('0.05','0.075','0.10','0.15','0.20','0.30'),adj=0,cex=.95)
  graphics::text(.10,.91,'Ag geometric mean',adj=0,cex=.85)
  graphics::text(.10,.875,'mg/kg  |  log scale',adj=0,cex=.85)
  title <- switch(medium,pooled='Soil and sediment pooled',soil='Soil',sediment='Sediment')
  graphics::mtext(paste('Silver in occupied 20 km cells',title,sep=' — '),side=3,outer=TRUE,line=3.2,
                  adj=.03,font=2,cex=1.25)
  graphics::mtext(paste(if(western)'Western United States' else 'Contiguous United States',
    '| Full HF ICP-MS and residue-sintering cohort'),side=3,outer=TRUE,line=1.8,adj=.03,cex=.95)
  graphics::mtext('Hollow: sampled without a displayed mean. Blank: unsampled. Circle centres mark cells, not sites.',
    side=1,outer=TRUE,line=.5,adj=.03,cex=.76)
  graphics::mtext('Local censored-lognormal summaries; no interpolation or spatial confidence intervals. Source: USGS; borders: maps.',
    side=1,outer=TRUE,line=1.7,adj=.03,cex=.72)
  invisible(NULL)
}
