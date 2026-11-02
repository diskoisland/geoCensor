# Numerical SCPC reference for small Gaussian designs. Original implementation
# of Muller & Watson (2022) equations, not copied author code. The compact
# interval uses a finite grid and local searches: no global continuum guarantee
# is asserted. The far-decay tail has an explicit Gaussian TV bound.
source("scripts/analysis/scpc_gaussian_reference.R")

scpc_decay_from_average_reference <- function(coordinates, average_correlation = 0.03) {
  d <- as.numeric(dist(as.matrix(coordinates)))
  if (!length(d) || any(!is.finite(d)) || min(d) <= 0 ||
      length(average_correlation) != 1L || !is.finite(average_correlation) ||
      average_correlation <= 0 || average_correlation >= 1) stop("Invalid design or average correlation")
  upper <- 1 / min(d)
  while (mean(exp(-upper * d)) > average_correlation) upper <- upper * 2
  uniroot(function(c) mean(exp(-c * d)) - average_correlation,
          c(0, upper), tol = 1e-10)$root
}

scpc_envelope_critical_reference <- function(basis, grid_points = 31L,
                                            alpha = 0.05, tail_delta = 1e-10) {
  if (grid_points < 5L || grid_points != floor(grid_points)) stop("Invalid grid_points")
  n <- basis$n
  distance <- basis$distance
  dmin <- min(distance[upper.tri(distance)])
  epsilon <- (n - 1) * tail_delta
  if (!is.finite(epsilon) || epsilon >= 1 || tail_delta <= 0) stop("Invalid tail_delta")
  # For c >= cmax, all off-diagonals are <= tail_delta. Pinsker plus the
  # Gaussian KL bound implies TV <= delta*sqrt(n*(n-1)/(8*(1-epsilon))).
  tail_tv <- tail_delta * sqrt(n * (n - 1) / (8 * (1 - epsilon)))
  if (tail_tv >= alpha) stop("Tail approximation too loose")
  cmax <- max(basis$c0, -log(tail_delta) / dmin)
  top <- log(cmax / basis$c0)
  cache <- new.env(parent = emptyenv())
  critical_at <- function(u) {
    key <- sprintf("%.17g", u)
    if (!exists(key, cache, inherits = FALSE)) {
      Sigma <- exp(-basis$c0 * exp(u) * distance)
      Omega <- crossprod(basis$W, Sigma %*% basis$W)
      cv <- scpc_gaussian_critical_reference(Omega, alpha)
      assign(key, cv, cache)
    }
    get(key, cache, inherits = FALSE)
  }
  u <- seq(0, top, length.out = grid_points)
  values <- vapply(u, critical_at, numeric(1))
  peaks <- which(values >= c(-Inf, head(values, -1L)) &
                 values >= c(tail(values, -1L), -Inf))
  candidates_u <- u
  candidates_cv <- values
  for (i in peaks) {
    a <- u[max(1L, i - 1L)]; b <- u[min(length(u), i + 1L)]
    if (b > a) {
      opt <- optimize(critical_at, c(a, b), maximum = TRUE, tol = 1e-5)
      candidates_u <- c(candidates_u, opt$maximum)
      candidates_cv <- c(candidates_cv, opt$objective)
    }
  }
  # A small numerical margin is fixed before validation. The tail uses an iid
  # critical value at alpha-tail_tv; the compact search remains approximate.
  iid_cv <- qt(1 - (alpha - tail_tv) / 2, basis$q)
  best <- which.max(candidates_cv)
  cv <- max(candidates_cv, iid_cv) * (1 + 1e-6)
  list(critical = cv, q = basis$q, c0 = basis$c0, cmax = cmax,
       tail_delta = tail_delta, tail_tv = tail_tv,
       worst_decay = basis$c0 * exp(candidates_u[best]),
       grid = data.frame(decay = basis$c0 * exp(u), point_critical = values),
       iid_critical = iid_cv, grid_points = grid_points,
       status = "numerical_envelope_not_global_certificate")
}

scpc_select_reference <- function(coordinates, c0, q_candidates = 1:12,
                                  grid_points = 31L, alpha = 0.05) {
  if (!length(q_candidates) || anyDuplicated(q_candidates)) stop("Invalid candidate set")
  fits <- lapply(q_candidates, function(q) {
    basis <- scpc_basis_reference(coordinates, c0, q)
    calibration <- scpc_envelope_critical_reference(basis, grid_points, alpha)
    expected_length <- 2 * calibration$critical * sqrt(2 / (q * basis$n)) *
      exp(lgamma((q + 1) / 2) - lgamma(q / 2))
    list(basis = basis, calibration = calibration, expected_length = expected_length)
  })
  tab <- do.call(rbind, lapply(fits, function(x) data.frame(q = x$basis$q,
    critical = x$calibration$critical, expected_iid_interval_length = x$expected_length,
    worst_decay = x$calibration$worst_decay,
    tail_tv = x$calibration$tail_tv, cmax = x$calibration$cmax)))
  best <- which.min(tab$expected_iid_interval_length)
  list(selected = fits[[best]], candidates = tab,
       upper_candidate_boundary = tab$q[best] == max(q_candidates))
}

scpc_envelope_audit_reference <- function(fit, grid_points = 241L) {
  b <- fit$basis; cal <- fit$calibration
  # Midpoint mesh differs from the calibration mesh; endpoints are added.
  top <- log(cal$cmax / b$c0)
  u <- sort(unique(c(0, top, (seq_len(grid_points) - 0.5) / grid_points * top)))
  rejection <- function(u) {
    S <- exp(-b$c0 * exp(u) * b$distance)
    scpc_gaussian_rejection_reference(crossprod(b$W, S %*% b$W), cal$critical)
  }
  rp <- vapply(u, rejection, numeric(1))
  peaks <- which(rp >= c(-Inf, head(rp, -1L)) & rp >= c(tail(rp, -1L), -Inf))
  local_rp <- rp
  for (i in peaks) {
    a <- u[max(1L, i - 1L)]; z <- u[min(length(u), i + 1L)]
    if (z > a) local_rp <- c(local_rp,
      optimize(rejection, c(a, z), maximum = TRUE, tol = 1e-6)$objective)
  }
  list(max_rejection = max(local_rp),
       iid_rejection = 2 * pt(-cal$critical, b$q),
       tail_upper_bound = 2 * pt(-cal$critical, b$q) + cal$tail_tv,
       grid = data.frame(decay = b$c0 * exp(u), rejection = rp))
}
