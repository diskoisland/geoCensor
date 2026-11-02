source("scripts/analysis/scpc_weighted_reference.R")
set.seed(11030100)
s <- matrix(runif(32), ncol = 2L)
c0 <- scpc_decay_from_average_reference(s)
for (q in c(1L, 2L, 5L, 10L)) {
  got <- scpc_gaussian_norm_mean_reference(diag(q))
  expected <- sqrt(2) * exp(lgamma((q + 1) / 2) - lgamma(q / 2))
  stopifnot(abs(got - expected) < 1e-9)
}
u <- scpc_basis_reference(s, c0, 3L)
w <- scpc_weighted_basis_reference(s, rep(1, nrow(s)), c0, 3L)
stopifnot(max(abs(tcrossprod(u$W[, -1L]) - tcrossprod(w$V))) < 1e-9,
          max(abs(w$iid_Omega - diag(4L))) < 1e-9)
a <- scpc_envelope_critical_reference(u)
b <- scpc_weighted_envelope_reference(w)
stopifnot(abs(a$critical - b$critical) < 1e-7)
weights <- seq_len(nrow(s))^1.5
w <- scpc_weighted_basis_reference(s, weights, c0, 3L)
rescaled <- scpc_weighted_basis_reference(s, weights * 123, c0, 3L)
stopifnot(max(abs(w$iid_Omega - rescaled$iid_Omega)) < 1e-9)
y <- rnorm(nrow(s)); theta <- sum(w$weights * y)
psi <- nrow(s) * w$weights * (y - theta)
transformed <- as.vector(crossprod(w$W, y))
stopifnot(abs(transformed[1L] - sqrt(nrow(s)) * theta) < 1e-10,
          max(abs(transformed[-1L] - crossprod(w$V, psi))) < 1e-10,
          abs(sum(psi)) < 1e-10)
L <- chol(w$iid_Omega[-1L, -1L, drop = FALSE])
z <- matrix(rnorm(300000L), ncol = 3L) %*% L
norms <- sqrt(rowSums(z^2))
stopifnot(abs(mean(norms) - scpc_gaussian_norm_mean_reference(crossprod(L))) <
            5 * sd(norms) / sqrt(nrow(z)))
cat("Weighted residual algebra, equal-weight reduction, weight-scale invariance,\n",
    "Gaussian norm integration, and independent norm simulation checks passed.\n")
