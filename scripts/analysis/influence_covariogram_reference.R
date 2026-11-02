# Research diagnostic only; no automatic production bandwidth selection.
influence_covariogram_reference <- function(model,coordinates,width=.5,max_lag=6,min_pairs=20,max_sites=2000){
 if(!is.numeric(width)||length(width)!=1||!is.finite(width)||width<=0)stop('Invalid lag width')
 if(!is.numeric(max_lag)||length(max_lag)!=1||!is.finite(max_lag)||max_lag<=0)stop('Invalid maximum lag')
 if(!is.numeric(min_pairs)||length(min_pairs)!=1||!is.finite(min_pairs)||min_pairs<1||min_pairs!=floor(min_pairs))stop('Invalid minimum pairs')
 a<-aggregate_score_locations(model$scores,coordinates)
 n<-nrow(a$scores)
 if(n>max_sites)stop('Research covariogram site limit exceeded')
 influence<-as.vector(n*a$scores%*%model$bread[1,])
 d<-as.matrix(dist(a$coordinates));upper<-upper.tri(d)
 distance<-d[upper];products<-outer(influence,influence)[upper]
 coordinate_tolerance<-64*.Machine$double.eps*max(1,abs(a$coordinates))
 if(width<=2*coordinate_tolerance)stop('Lag width is below coordinate precision')
 keep<-is.finite(distance)&distance>0&distance<=max_lag+coordinate_tolerance
 # Right-closed bins need tolerance at exact boundaries to avoid changing bins
 # when coordinates are translated and their subtraction rounds differently.
 ratio<-distance[keep]/width
 bin<-pmax(1L,ceiling(ratio-coordinate_tolerance/width-64*.Machine$double.eps*pmax(1,abs(ratio))))
 products<-products[keep];distance<-distance[keep]
 bins<-data.frame(bin=seq_len(ceiling(max_lag/width)),distance=NA_real_,covariance=NA_real_,n_pairs=0L)
 for(i in seq_len(nrow(bins))){use<-bin==i
  bins$n_pairs[i]<-sum(use)
  if(any(use)){bins$distance[i]<-mean(distance[use]);bins$covariance[i]<-mean(products[use])}
 }
 valid<-bins$n_pairs>=min_pairs&is.finite(bins$covariance)
 z<-bins[valid,,drop=FALSE]
 range<-NA_real_;status<-'insufficient_bins'
 if(nrow(z)>=2){
  if(z$covariance[1]<=0)status<-'no_positive_initial_bin'
  else{
   cross<-which(z$covariance<=0)
   if(length(cross)){j<-cross[1];range<-mean(z$distance[c(j-1,j)]);status<-'first_zero_crossing'}
   else status<-'range_not_identified'
  }
 }
 list(bins=bins,range=range,status=status,n_sites=n,influence=influence)
}
