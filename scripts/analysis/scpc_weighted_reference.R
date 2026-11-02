# Conditional-on-locations weighted Gaussian research extension. This is not a
# published SCPC guarantee for estimated spatial weights or a production API.
source("scripts/analysis/scpc_envelope_reference.R")

scpc_weighted_basis_reference <- function(coordinates, weights, c0, q) {
  s <- as.matrix(coordinates); n <- nrow(s)
  if (!is.numeric(s) || ncol(s) != 2L || n < 4L || n > 500L ||
      any(!is.finite(s)) || nrow(unique(s)) != n) stop("Invalid coordinates")
  if (!is.numeric(weights) || length(weights) != n || any(!is.finite(weights)) ||
      any(weights <= 0)) stop("Strictly positive site weights required")
  if (length(q) != 1L || !is.finite(q) || q != floor(q) || q < 1 || q >= n) stop("Invalid q")
  if (length(c0) != 1L || !is.finite(c0) || c0 <= 0) stop("Invalid c0")
  a <- weights / max(weights); a <- a / sum(a)
  distance <- as.matrix(dist(s))
  B <- diag(n * a, n) %*% (diag(n) - outer(rep(1, n), a))
  covariance <- B %*% exp(-c0 * distance) %*% t(B)
  e <- eigen((covariance + t(covariance)) / 2, symmetric = TRUE)
  if (e$values[q] <= 0 || e$values[q] - e$values[q + 1L] < 1e-10 * e$values[1L]) {
    stop("Nonpositive or tied weighted eigenspace")
  }
  V <- e$vectors[, seq_len(q), drop = FALSE]
  # theta_hat=a'Y, psi=n*diag(a)*(Y-1*theta_hat). Hence (sqrt(n)*theta_hat,
  # V'psi) = W'Y. The weighted residual operator is not the ordinary M matrix.
  W <- cbind(sqrt(n) * a, t(B) %*% V)
  stopifnot(max(abs(colSums(B))) < 1e-10, max(abs(B %*% rep(1, n))) < 1e-10)
  list(n = n, q = q, c0 = c0, coordinates = s, distance = distance,
       weights = a, B = B, V = V, W = W, iid_Omega = crossprod(W))
}

scpc_gaussian_norm_mean_reference <- function(covariance) {
  e <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(e)) || any(e <= 0)) stop("Invalid Gaussian norm covariance")
  f <- function(v) vapply(v, function(z) {
    if (z == 0) return(sum(e))
    if (z == 1) return(1)
    -expm1(-0.5 * sum(log1p(2 * e * (z / (1 - z))^2))) / z^2
  }, numeric(1))
  # Laplace identity for E sqrt(Z'AZ), transformed to the unit interval.
  answer <- integrate(f, 0, 1, rel.tol = 1e-9, abs.tol = 1e-11, subdivisions = 1000L)
  if (!identical(answer$message, "OK")) stop(answer$message)
  answer$value / sqrt(pi)
}

scpc_weighted_envelope_reference <- function(basis, grid_points = 31L,
                                             alpha = 0.05) {
  if (grid_points < 5L || grid_points != floor(grid_points)) stop("Invalid grid_points")
  n <- basis$n; delta <- 1e-10
  tail_tv <- delta * sqrt(n * (n - 1) / (8 * (1 - (n - 1) * delta)))
  cmax <- max(basis$c0, -log(delta) / min(basis$distance[upper.tri(basis$distance)]))
  top <- log(cmax / basis$c0)
  cache <- new.env(parent = emptyenv())
  point <- function(u) {
    key <- sprintf("%.17g", u)
    if (!exists(key, cache, inherits = FALSE)) {
      S <- exp(-basis$c0 * exp(u) * basis$distance)
      Omega <- crossprod(basis$W, S %*% basis$W)
      assign(key, scpc_gaussian_critical_reference(Omega, alpha), cache)
    }
    get(key, cache, inherits = FALSE)
  }
  u <- seq(0, top, length.out = grid_points)
  cv <- vapply(u, point, numeric(1))
  peaks <- which(cv >= c(-Inf, head(cv, -1L)) & cv >= c(tail(cv, -1L), -Inf))
  all_cv <- cv
  for (i in peaks) {
    bounds <- u[c(max(1L, i - 1L), min(length(u), i + 1L))]
    if (diff(bounds) > 0) all_cv <- c(all_cv,
      optimize(point, bounds, maximum = TRUE, tol = 1e-5)$objective)
  }
  # Even under iid observations, unequal weights do not yield a Student-t
  # statistic. Use the actual transformed covariance at the independent limit.
  iid_cv <- scpc_gaussian_critical_reference(basis$iid_Omega, alpha - tail_tv)
  critical <- max(all_cv, iid_cv) * (1 + 1e-6)
  expected_iid_length <- 2 * critical / sqrt(basis$q * n) *
    scpc_gaussian_norm_mean_reference(basis$iid_Omega[-1L, -1L, drop = FALSE])
  list(critical = critical, iid_critical = iid_cv, tail_tv = tail_tv,
       cmax = cmax, expected_iid_length = expected_iid_length,
       grid = data.frame(decay = basis$c0 * exp(u), point_critical = cv),
       status = "weighted_numerical_envelope_not_global_certificate")
}

scpc_weighted_select_reference <- function(coordinates, weights, c0,
                                           grid_points = 31L, alpha = 0.05) {
  # Search all nonzero centered components; no arbitrary q truncation.
  fits <- lapply(seq_len(nrow(coordinates) - 1L), function(q) {
    basis <- scpc_weighted_basis_reference(coordinates, weights, c0, q)
    cal <- scpc_weighted_envelope_reference(basis, grid_points, alpha)
    list(basis = basis, calibration = cal)
  })
  tab <- do.call(rbind, lapply(fits, function(x) data.frame(q = x$basis$q,
    critical = x$calibration$critical, expected_iid_length = x$calibration$expected_iid_length)))
  list(selected = fits[[which.min(tab$expected_iid_length)]], candidates = tab)
}

scpc_weighted_audit_reference <- function(fit, grid_points = 241L) {
  b <- fit$basis; cal <- fit$calibration
  top <- log(cal$cmax / b$c0)
  u <- sort(unique(c(0, top, (seq_len(grid_points) - 0.5) / grid_points * top)))
  point <- function(u) {
    S <- exp(-b$c0 * exp(u) * b$distance)
    scpc_gaussian_rejection_reference(crossprod(b$W, S %*% b$W), cal$critical)
  }
  rp <- vapply(u, point, numeric(1))
  peaks <- which(rp >= c(-Inf, head(rp, -1L)) & rp >= c(tail(rp, -1L), -Inf))
  values <- rp
  for (i in peaks) {
    bounds <- u[c(max(1L, i - 1L), min(length(u), i + 1L))]
    if (diff(bounds) > 0) values <- c(values,
      optimize(point, bounds, maximum = TRUE, tol = 1e-6)$objective)
  }
  list(max_rejection = max(values),
       iid_tail_bound = scpc_gaussian_rejection_reference(b$iid_Omega, cal$critical) + cal$tail_tv,
       grid = data.frame(decay = b$c0 * exp(u), rejection = rp))
}
