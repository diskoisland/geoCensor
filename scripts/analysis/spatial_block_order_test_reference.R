# Research only: positive spatial ordering of exchangeable site blocks.
# No concentration imputation, covariance fit, range estimate or CI procedure.

block_order_scores_reference <- function(lower, upper) {
  n <- nrow(lower); k <- ncol(lower)
  lo <- as.vector(t(lower)); hi <- as.vector(t(upper))
  dominance <- outer(lo, hi, ">") * 1 - outer(hi, lo, "<") * 1
  site <- rep(seq_len(n), each = k)
  dominance[outer(site, site, "==")] <- 0
  scores <- as.vector(rowsum(rowSums(dominance), site, reorder = FALSE))
  stopifnot(sum(scores) == 0)
  list(scores = scores, definite_cross_site_record_pairs = sum(dominance != 0) / 2)
}

block_order_plan_reference <- function(coordinates, scales = c(.05, .15, .4),
    permutations = 999L, seed) {
  xy <- as.matrix(coordinates); n <- nrow(xy)
  if (!is.numeric(xy) || ncol(xy) != 2L || n < 4L || n > 128L || any(!is.finite(xy)) ||
      anyDuplicated(data.frame(xy))) stop("Provide 4 to 128 distinct finite site coordinates")
  if (!is.numeric(scales) || !length(scales) || length(scales) > 8L || any(!is.finite(scales)) ||
      any(scales <= 0) || anyDuplicated(scales)) stop("Invalid fixed geometry scales")
  if (length(permutations) != 1L || !is.finite(permutations) || permutations < 1L ||
      permutations != as.integer(permutations) || !is.finite(seed)) stop("Invalid permutation count/seed")
  distance <- as.matrix(dist(xy))
  weights <- lapply(scales, function(scale) {
    w <- exp(-distance/scale); diag(w) <- 0
    if (!is.finite(sum(w)) || sum(w) <= 0) stop("Distance scales have no numerical kernel support")
    w * n / sum(w)
  })
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had) saved <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had) assign(".Random.seed", saved, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv))
  set.seed(seed)
  index <- rbind(seq_len(n), t(replicate(permutations, sample.int(n))))
  list(coordinates = xy, scales = scales, weights = weights, index = index,
    random_permutations = permutations, seed = seed)
}

block_order_statistics_reference <- function(scores, plan) {
  index <- plan$index; n <- length(scores)
  if (ncol(index) != n || !all(index[1L, ] == seq_len(n)) ||
      any(!is.finite(index)) || any(index < 1 | index > n | index != floor(index)) ||
      any(apply(index, 1L, anyDuplicated) > 0)) stop("Invalid complete-site permutation plan")
  centered <- scores - mean(scores)
  denominator <- sum(centered^2)
  if (!is.finite(denominator) || denominator <= 0) stop("No site-order score variation")
  values <- matrix(centered[as.vector(index)], nrow = nrow(index), ncol = n)
  statistics <- vapply(plan$weights, function(w) rowSums((values %*% w) * values) / denominator,
    numeric(nrow(values)))
  if (is.null(dim(statistics))) statistics <- matrix(statistics, ncol = 1L)
  colnames(statistics) <- paste0("scale_", plan$scales)
  list(per_scale = statistics, maximum = apply(statistics, 1L, max))
}

spatial_block_order_test_reference <- function(lower, upper, coordinates,
    lower_reporting = NULL, upper_reporting = NULL, exchangeable_under_null = FALSE,
    plan) {
  lower <- as.matrix(lower); upper <- as.matrix(upper); xy <- as.matrix(coordinates)
  n <- nrow(lower); k <- ncol(lower)
  unavailable <- function(status, scores = NULL, definite_pairs = NA_real_) list(status = status,
    p_value = NA_real_, statistic = NA_real_, n_sites = n, records_per_site = k,
    scores = scores, definite_cross_site_record_pairs = definite_pairs)
  if (!length(lower) && !length(upper)) return(unavailable("no_data"))
  if (!is.numeric(lower) || !is.numeric(upper) || !identical(dim(lower), dim(upper)) ||
      n < 4L || k < 1L || k > 4L || !identical(dim(xy), c(n, 2L))) stop("Invalid balanced site-block inputs")
  if (anyNA(lower) || anyNA(upper)) return(unavailable("missing_observed_bounds"))
  if (any(lower > upper | lower == Inf | upper == -Inf | (lower == -Inf & upper == Inf)))
    return(unavailable("invalid_observed_bounds"))
  if (!isTRUE(exchangeable_under_null)) return(unavailable("exchangeability_not_asserted"))
  if (is.null(lower_reporting) || is.null(upper_reporting)) return(unavailable("reporting_design_not_supplied"))
  L <- as.matrix(lower_reporting); U <- as.matrix(upper_reporting)
  if (!is.numeric(L) || !is.numeric(U) || !identical(dim(L), dim(lower)) || !identical(dim(U), dim(lower)))
    stop("Reporting design must describe every record, including exact observations")
  if (anyNA(L) || anyNA(U)) return(unavailable("missing_reporting_design"))
  if (any(L >= U | L == Inf | U == -Inf)) return(unavailable("invalid_reporting_design"))
  common <- all(vapply(seq_len(k), function(j) all(L[, j] == L[1L, j]) && all(U[, j] == U[1L, j]), TRUE))
  if (!common) return(unavailable("heterogeneous_reporting_design"))
  exact <- lower == upper
  left <- lower == -Inf & is.finite(upper) & upper == L
  right <- upper == Inf & is.finite(lower) & lower == U
  if (any(!(exact | left | right)) || any(exact & (lower < L | upper > U)))
    return(unavailable("bounds_reporting_design_mismatch"))
  if (!identical(xy, plan$coordinates)) stop("Plan coordinates do not match observed site blocks")
  score <- block_order_scores_reference(lower, upper)
  if (sum((score$scores - mean(score$scores))^2) == 0)
    return(unavailable("no_order_score_variation", score$scores, score$definite_cross_site_record_pairs))
  statistic <- block_order_statistics_reference(score$scores, plan)
  values <- statistic$maximum; observed <- values[1L]
  tolerance <- 64 * .Machine$double.eps * max(1, abs(values))
  exceedances <- sum(values[-1L] >= observed - tolerance)
  p <- (1 + exceedances) / length(values)
  stopifnot(is.finite(p), p > 0, p <= 1)
  list(status = "available", p_value = p, statistic = observed,
    n_sites = n, records_per_site = k, scores = score$scores,
    definite_cross_site_record_pairs = score$definite_cross_site_record_pairs,
    observed_per_scale = statistic$per_scale[1L, ], random_permutations = length(values) - 1L,
    exceedances = exceedances, numerical_tie_tolerance = tolerance,
    permutation_maxima = values, permutation_per_scale = statistic$per_scale)
}
