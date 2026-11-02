# Research-only checks of fixed-q SCPC Gaussian formulas, not a package API.
# Equations (3)-(4) and Lemma 1 of Muller & Watson (2022).
# This is original R code from the mathematical formulas, not a translation of
# the supplied Fortran program. No covariance envelope or component selection
# is performed here, so this is not a complete published SCPC implementation.

scpc_basis_reference <- function(coordinates, c0, q) {
  s <- as.matrix(coordinates)
  if (!is.numeric(s) || ncol(s) != 2L || nrow(s) < 4L ||
      any(!is.finite(s)) || nrow(unique(s)) != nrow(s)) stop("Invalid coordinates")
  n <- nrow(s)
  if (n > 500L) stop("Small-design reference limited to 500 locations")
  if (length(c0) != 1L || !is.finite(c0) || c0 <= 0) stop("Invalid c0")
  if (length(q) != 1L || !is.finite(q) || q != floor(q) || q < 1 || q > n - 2L) {
    stop("q must leave at least one unselected centered component")
  }
  distance <- as.matrix(dist(s))
  Sigma <- exp(-c0 * distance)
  M <- diag(n) - matrix(1 / n, n, n)
  centered <- M %*% Sigma %*% M
  decomposition <- eigen((centered + t(centered)) / 2, symmetric = TRUE)
  if (decomposition$values[q] <= 0 ||
      decomposition$values[q] - decomposition$values[q + 1L] <
        1e-10 * decomposition$values[1L]) stop("Nonpositive or tied retained eigenspace")
  V <- decomposition$vectors[, seq_len(q), drop = FALSE]
  list(W = cbind(rep(1 / sqrt(n), n), V), q = q, n = n,
       c0 = c0, distance = distance, coordinates = s)
}

scpc_gaussian_rejection_reference <- function(Omega, critical) {
  Omega <- as.matrix(Omega)
  q <- nrow(Omega) - 1L
  if (!is.numeric(Omega) || ncol(Omega) != nrow(Omega) || q < 1L ||
      any(!is.finite(Omega)) || max(abs(Omega - t(Omega))) > 1e-10) {
    stop("Invalid projected covariance")
  }
  if (length(critical) != 1L || !is.finite(critical) || critical <= 0) stop("Invalid critical value")
  L <- t(chol(Omega))
  # Multiplying each row by the diagonal avoids an unnecessary dense product.
  diagonal <- c(1, rep(-critical^2 / q, q))
  e <- eigen(crossprod(L, L * diagonal), symmetric = TRUE, only.values = TRUE)$values
  if (sum(e > 0) != 1L || any(e[-1L] >= 0)) stop("Invalid quadratic-form inertia")
  a <- -e[-1L] / e[1L]
  # Substitute x=cos(theta)^2 in the one-positive-root integral. Logarithms
  # avoid products overflowing near the endpoint; no eigenvalues are clipped.
  integrand <- function(theta) {
    exp(-0.5 * rowSums(log1p(outer(1 / cos(theta)^2, a))))
  }
  answer <- integrate(integrand, 0, pi / 2, rel.tol = 1e-10,
                      abs.tol = 1e-12, subdivisions = 1000L)
  if (!identical(answer$message, "OK")) stop(answer$message)
  2 / pi * answer$value
}

scpc_gaussian_critical_reference <- function(Omega, alpha = 0.05) {
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("Invalid alpha")
  upper <- 2
  while (scpc_gaussian_rejection_reference(Omega, upper) > alpha) {
    upper <- upper * 2
    if (upper > 1e6) stop("Critical-value bracket exceeded")
  }
  uniroot(function(cv) scpc_gaussian_rejection_reference(Omega, cv) - alpha,
          c(1e-6, upper), tol = 1e-9)$root
}
