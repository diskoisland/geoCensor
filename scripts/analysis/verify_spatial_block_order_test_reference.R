#!/usr/bin/env Rscript
source("scripts/analysis/spatial_block_order_test_reference.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply a fresh verification directory")
out <- args[1L]
if (file.exists(file.path(out, "verification.csv"))) stop("Verification already exists")
checks <- list()
record <- function(name, actual, expected, tolerance) {
  error <- max(abs(actual - expected)); stopifnot(is.finite(error), error <= tolerance)
  checks[[length(checks)+1L]] <<- data.frame(check = name, error = error,
    tolerance = tolerance, passed = TRUE)
}
all_permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, 1L, 1L))
  do.call(rbind, lapply(seq_along(x), function(i) cbind(x[i], all_permutations(x[-i]))))
}
xy <- cbind(c(0,.1,.3,.6,.7,1), c(.1,.8,.4,.2,.9,.5))
y <- matrix(c(-.4,.1,0,.5,.8,.2,-.2,.1,.7,.6,-.3,.4), 6L, 2L, byrow = TRUE)
scores <- block_order_scores_reference(y, y)$scores
ranks <- matrix(rank(as.vector(t(y)), ties.method = "average"), 6L, 2L, byrow = TRUE)
record("exact_score_equals_summed_midrank_score", scores, 2*rowSums(ranks)-2*(length(y)+1), 1e-12)
L <- matrix(-.15,6L,2L); U <- matrix(.55,6L,2L)
lo <- ifelse(y < L,-Inf,ifelse(y > U,U,y)); up <- ifelse(y > U,Inf,ifelse(y < L,L,y))
s <- block_order_scores_reference(lo, up)$scores
p <- c(3L,1L,6L,2L,5L,4L)
record("block_permutation_equivariance", block_order_scores_reference(lo[p,],up[p,])$scores, s[p], 0)
record("record_order_invariance", block_order_scores_reference(lo[,2:1],up[,2:1])$scores, s, 0)
set.seed(900010001L); saved <- .Random.seed
plan <- block_order_plan_reference(xy, permutations = 199L, seed = 900010002L)
stopifnot(identical(saved,.Random.seed))
permutations <- all_permutations(1:6)
exact_plan <- plan; exact_plan$index <- permutations
statistics <- block_order_statistics_reference(s, exact_plan)
record("exact_null_Moran_mean_each_scale", colMeans(statistics$per_scale), rep(-1/5,3L), 1e-12)
test <- spatial_block_order_test_reference(lo,up,xy,L,U,TRUE,exact_plan)
tol <- test$numerical_tie_tolerance
exact_p <- vapply(statistics$maximum,function(value) mean(statistics$maximum >= value-tol),0.0)
record("exhaustive_group_p_value", test$p_value, exact_p[1L], 1e-12)
for (alpha in c(.01,.05,.1,.2)) {
  rate <- mean(exact_p <= alpha)
  stopifnot(rate <= alpha+1e-12)
  record(paste0("exhaustive_level_excess_",alpha), max(0,rate-alpha), 0, 1e-12)
}
mc <- spatial_block_order_test_reference(lo,up,xy,L,U,TRUE,plan)
record("plus_one_tail_count",mc$p_value,(1+mc$exceedances)/(plan$random_permutations+1),0)
stopifnot(mc$p_value>=1/200,identical(saved,.Random.seed))
for (j in c(2L,17L,93L)) {
  order <- plan$index[j,]
  direct <- block_order_scores_reference(lo[order,],up[order,])$scores
  direct_plan <- plan; direct_plan$index <- rbind(1:6,1:6)
  expected <- block_order_statistics_reference(direct,direct_plan)$maximum[1L]
  record(paste0("direct_recomputation_matches_score_permutation_",j),
    mc$permutation_maxima[j],expected,1e-12)
}
missing <- lo; missing[1,1] <- NA_real_
stopifnot(spatial_block_order_test_reference(missing,up,xy,L,U,TRUE,plan)$status=="missing_observed_bounds")
heterogeneous <- L; heterogeneous[1,1] <- -.1
guard <- spatial_block_order_test_reference(lo,up,xy,heterogeneous,U,TRUE,plan)
stopifnot(guard$status=="heterogeneous_reporting_design",is.na(guard$p_value))
stopifnot(spatial_block_order_test_reference(lo,up,xy,L,U,FALSE,plan)$status=="exchangeability_not_asserted",
  spatial_block_order_test_reference(lo,up,xy,exchangeable_under_null=TRUE,plan=plan)$status=="reporting_design_not_supplied")
constant <- spatial_block_order_test_reference(matrix(-Inf,6,2),L,xy,L,U,TRUE,plan)
stopifnot(constant$status=="no_order_score_variation",is.na(constant$p_value))
empty <- spatial_block_order_test_reference(matrix(numeric(),0,2),matrix(numeric(),0,2),matrix(numeric(),0,2))
stopifnot(empty$status=="no_data",is.na(empty$p_value))
only_lower <- matrix(-Inf,6L,2L); only_upper <- matrix(-1,6L,2L)
only_lower[4:6,] <- 1; only_upper[4:6,] <- Inf
informative <- spatial_block_order_test_reference(only_lower,only_upper,xy,
  matrix(-1,6,2),matrix(1,6,2),TRUE,plan)
stopifnot(informative$status=="available",is.finite(informative$p_value))
bad_plan <- plan; bad_plan$index[2,1] <- bad_plan$index[2,2]
stopifnot(inherits(tryCatch(block_order_statistics_reference(s,bad_plan),error=function(e)e),"error"))
write.csv(do.call(rbind,checks),file.path(out,"verification.csv"),row.names=FALSE)
write.csv(data.frame(statistic=statistics$maximum,p_value=exact_p),
  file.path(out,"exhaustive_group_audit.csv"),row.names=FALSE)
writeLines(paste(length(checks),"rank/equivariance/exhaustive-level/plus-one comparisons and missing/design/RNG guards passed."),
  file.path(out,"VERIFIED.txt"))
cat(length(checks),"order-score/permutation comparisons and eligibility guards passed.\n")
