# Research reference implementation, not a production geoCensor API.
# Independently assembled from likelihood derivatives and S' K S; no external
# implementation code is copied. SpatialInference supplies bandwidth/kernel ideas.
# Dense pairwise calculations are deliberately limited to small validation cases.

fit_censored_score_reference <- function(df, weights = rep(1, nrow(df))) {
  z <- geoCensor:::.censored_intervals(df, 'detected_val', 'left_cen', 'left_limit',
    'right_cen', 'right_limit', 'fit_censored_score_reference')
  if (!is.numeric(weights) || length(weights) != nrow(df) || any(!is.finite(weights)) || any(weights < 0)) stop('Invalid weights')
  ok <- z$ok & weights > 0
  if (sum(ok) < 6) stop('Insufficient observations')
  w <- weights[ok] / max(weights[ok]); w <- w / sum(w) * sum(ok)
  lo <- z$lower[ok]; up <- z$upper[ok]
  fit <- withCallingHandlers(survival::survreg(survival::Surv(lo, up, type='interval2')~1,
    dist='lognormal',weights=w,model=TRUE),warning=function(e)stop(conditionMessage(e),call.=FALSE))
  r <- residuals(fit,type='matrix',weighted=FALSE)
  # Both location and log-scale derivatives must use the case weights.
  scores <- cbind(mu=r[,'dg']*w,log_scale=r[,'ds']*w)
  information <- -matrix(c(sum(r[,'ddg']*w),sum(r[,'dsg']*w),
    sum(r[,'dsg']*w),sum(r[,'dds']*w)),2,2)
  if (any(!is.finite(scores)) || any(!is.finite(information)) ||
      min(eigen(information,symmetric=TRUE,only.values=TRUE)$values)<=0) stop('Invalid score information')
  bread <- solve(information)
  if(max(abs(colSums(scores)))/sum(w)>1e-5)stop('Score equations not sufficiently solved')
  list(fit=fit,scores=scores,information=information,bread=bread,
    theta=c(unname(coef(fit)),log(fit$scale)),ok=ok,weights=w)
}

aggregate_score_locations <- function(scores, coordinates) {
  coordinates<-as.matrix(coordinates)
  if(!is.numeric(coordinates)||ncol(coordinates)!=2||nrow(coordinates)!=nrow(scores)||any(!is.finite(coordinates)))stop('Invalid coordinates')
  # 17 significant digits preserve double precision for exact coordinate grouping.
  coordinates[coordinates==0]<-0
  key<-paste(sprintf('%.17g',coordinates[,1]),sprintf('%.17g',coordinates[,2]),sep=':')
  first<-!duplicated(key)
  list(scores=rowsum(scores,key,reorder=FALSE),coordinates=coordinates[first,,drop=FALSE])
}

spatial_hac_reference <- function(model,coordinates,bandwidth,kernel='bartlett_product',max_sites=2000) {
  if(!is.numeric(bandwidth)||length(bandwidth)!=1||!is.finite(bandwidth)||bandwidth<=0)stop('Invalid bandwidth')
  if(!kernel%in%c('bartlett_product','bartlett_radial','epanechnikov','uniform','gaussian'))stop('Unknown kernel')
  a<-aggregate_score_locations(model$scores,coordinates)
  n<-nrow(a$coordinates)
  if(n>max_sites)stop('Research reference site limit exceeded; no dense national run is permitted')
  dx<-abs(outer(a$coordinates[,1],a$coordinates[,1],'-'))/bandwidth
  dy<-abs(outer(a$coordinates[,2],a$coordinates[,2],'-'))/bandwidth
  distance<-sqrt(dx^2+dy^2)
  K<-switch(kernel,
    bartlett_product=pmax(1-dx,0)*pmax(1-dy,0),
    bartlett_radial=pmax(1-distance,0),
    epanechnikov=pmax(1-distance^2,0),
    uniform=1*(distance<=1),
    gaussian=exp(-distance^2/2))
  meat<-crossprod(a$scores,K%*%a$scores)
  V<-model$bread%*%meat%*%model$bread
  V<-(V+t(V))/2
  eigenvalues<-eigen(V,symmetric=TRUE,only.values=TRUE)$values
  positive<-all(is.finite(V)) && min(eigenvalues)>= -1e-10*max(abs(eigenvalues),1e-12)
  se<-if(positive && V[1,1]>0)sqrt(V[1,1])else NA_real_
  list(vcov=V,se_log_geomean=se,n_sites=n,kernel=kernel,bandwidth=bandwidth,
    min_eigenvalue=min(eigenvalues),score_sum=colSums(model$scores),
    status=if(!positive)'non_positive_semidefinite'else if(n<10)'few_locations'else 'experimental_spatial_hac',
    theta=model$theta[1],site_scores=a$scores,coordinates=a$coordinates)
}
