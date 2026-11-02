#!/usr/bin/env Rscript
# Original, deterministic counterchecks of generic method interpretations.
# These do not reproduce either author's software or simulation experiments.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Supply an optional output directory")
out <- if (length(args)) args[1L] else "output/additional_methods_review_2026-11-01"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# A one-dimensional special case of randomized quadratic projection:
# beta ~ N(0,1), theta = argmin (theta-beta)^2 subject to theta >= 0.
# theta=max(beta,0) is feasible but is not N(0,1) conditioned on theta>=0.
projected_mean <- integrate(function(z) z * dnorm(z), 0, Inf)$value
projected_second <- integrate(function(z) z^2 * dnorm(z), 0, Inf)$value
truncated_mean <- projected_mean / pnorm(0, lower.tail = FALSE)
truncated_second <- projected_second / pnorm(0, lower.tail = FALSE)
stopifnot(abs(projected_mean - dnorm(0)) < 1e-10,
          abs(truncated_mean - sqrt(2 / pi)) < 1e-10)
projection <- data.frame(
  distribution = c("quadratic_projection_max_Z_0", "normal_conditioned_Z_ge_0"),
  mass_at_zero = c(pnorm(0), 0),
  mean = c(projected_mean, truncated_mean),
  variance = c(projected_second - projected_mean^2,
               truncated_second - truncated_mean^2),
  q025 = c(0, qnorm(0.5 + 0.5 * 0.025)),
  q975 = c(qnorm(0.975), qnorm(0.5 + 0.5 * 0.975)))
stopifnot(projection$mass_at_zero[1L] == 0.5,
          projection$mean[1L] < projection$mean[2L])
write.csv(projection, file.path(out, "projection_vs_truncation.csv"), row.names = FALSE)

# Spherical covariance on an expanding 2D unit grid, with true support 1.2.
# Only the four cardinal neighbors have nonzero covariance. Even an oracle
# cutoff equal to the true support leaves attenuation from a fixed taper.
radius <- 1.2
spherical <- function(distance) {
  u <- distance / radius
  ifelse(u < 1, 1 - 1.5 * u + 0.5 * u^3, 0)
}
rho <- spherical(1)
neighbor_weight <- c(bartlett = 1 - 1 / radius,
                     epanechnikov = 1 - (1 / radius)^2,
                     uniform = 1)
finite <- do.call(rbind, lapply(c(3L, 7L, 12L), function(m) {
  s <- as.matrix(expand.grid(x = seq_len(m), y = seq_len(m)))
  d <- as.matrix(dist(s)); n <- nrow(s)
  S <- spherical(d)
  stopifnot(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values) > 0)
  true <- sum(S) / n
  stopifnot(abs(true - (1 + 4 * (1 - 1 / m) * rho)) < 1e-12)
  do.call(rbind, lapply(names(neighbor_weight), function(kernel) {
    K <- switch(kernel, bartlett = pmax(0, 1 - d / radius),
                epanechnikov = pmax(0, 1 - (d / radius)^2),
                uniform = 1 * (d < radius))
    oracle <- sum(K * S) / n
    formula <- 1 + 4 * (1 - 1 / m) * rho * neighbor_weight[kernel]
    stopifnot(abs(oracle - formula) < 1e-12)
    M <- diag(n) - matrix(1 / n, n, n)
    residual <- sum(K * (M %*% S %*% M)) / n
    data.frame(grid_side = m, n = n, kernel = kernel,
               variance_sqrt_n_mean = true,
               oracle_taper_expectation = oracle,
               demeaned_taper_expectation = residual,
               oracle_ratio_to_true = oracle / true)
  }))
}))
write.csv(finite, file.path(out, "fixed_range_taper_finite.csv"), row.names = FALSE)
true_limit <- 1 + 4 * rho
taper_limit <- 1 + 4 * rho * neighbor_weight
limits <- data.frame(kernel = names(neighbor_weight),
  radius = radius, nearest_neighbor_covariance = rho,
  true_variance_limit = true_limit, taper_expectation_limit = unname(taper_limit),
  ratio_to_true = unname(taper_limit / true_limit))
stopifnot(all(limits$ratio_to_true[limits$kernel != "uniform"] < 1),
          abs(limits$ratio_to_true[limits$kernel == "uniform"] - 1) < 1e-14)
write.csv(limits, file.path(out, "fixed_range_taper_limits.csv"), row.names = FALSE)
writeLines(c(
  "Deterministic literature counterchecks — November 1, 2026",
  "Scalar constrained Gaussian projection has a boundary atom; truncation does not.",
  "Direct integration verified analytic means and second moments.",
  "Spherical covariance matrices checked positive definite on 3x3, 7x7 and 12x12 grids.",
  "Direct matrix sums matched exact finite-grid covariance and tapered expectations.",
  "Fixed Bartlett/Epanechnikov tapers attenuate nonzero interior covariance even with an oracle support cutoff.",
  "These are original generic method checks, not author-code audits, full theorem audits or simulations of geoCensor.",
  "No production API or estimator was changed."), file.path(out, "method_counterchecks.txt"))
cat("Deterministic projection and fixed-range taper counterchecks passed.\n")
