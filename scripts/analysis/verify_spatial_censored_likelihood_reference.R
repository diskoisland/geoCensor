#!/usr/bin/env Rscript
source("scripts/analysis/spatial_censored_likelihood_reference.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Supply an optional output directory")
out <- if (length(args)) args[1L] else "output/spatial_likelihood_reference_2026-11-01"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checks <- list()
record <- function(name, actual, expected, tolerance) {
  error <- max(abs(actual - expected))
  if (!is.finite(error) || error > tolerance) stop(name, ": discrepancy ", error)
  checks[[length(checks) + 1L]] <<- data.frame(check = name, error = error,
    tolerance = tolerance, passed = TRUE)
}
must_fail <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL},
  error = function(e) e), "error"))

# Complete observations: independently assembled Gaussian density and GLS mean.
S <- matrix(c(1.2, .3, .1, .3, .8, .2, .1, .2, .6), 3L)
y <- c(-.5, .1, .9)
m <- spatial_censored_setup_reference(y, y, S)
for (mu in c(-1, 0, 1)) record(paste("all_exact", mu),
  spatial_censored_loglik_reference(m, mu),
  mvtnorm::dmvnorm(y, rep(mu, 3L), S, log = TRUE), 1e-12)
gls <- sum(solve(S, y)) / sum(solve(S, rep(1, 3L)))
record("all_exact_location_GLS", spatial_censored_fit_mu_reference(m)$theta,
       gls, 1e-6)

# Diagonal covariance: independent exact densities and tail probabilities.
lo <- c(.2, -Inf, .8, -.7); up <- c(.2, -.1, Inf, .9)
sd <- c(.6, .8, 1.1, .9)
m <- spatial_censored_setup_reference(lo, up, diag(sd^2))
for (mu in c(-.4, 0, .5)) {
  expected <- dnorm(.2, mu, sd[1L], log = TRUE) +
    pnorm(-.1, mu, sd[2L], log.p = TRUE) +
    pnorm(.8, mu, sd[3L], lower.tail = FALSE, log.p = TRUE) +
    log(pnorm(.9, mu, sd[4L]) - pnorm(-.7, mu, sd[4L]))
  u <- spatial_censored_uniforms_reference(m, 4096L, 321L)
  record(paste("diagonal_reduction", mu),
    spatial_censored_loglik_reference(m, mu, u), expected, 1e-11)
}

# One censored dimension: exact one-dimensional conditional integration.
S2 <- matrix(c(1, .45, .45, .7), 2L)
for (bounds in list(c(-Inf, .6), c(-.5, Inf), c(-.5, .6))) {
  m <- spatial_censored_setup_reference(c(.3, bounds[1L]), c(.3, bounds[2L]), S2)
  mu <- .1
  integrand <- function(z) vapply(z, function(v)
    mvtnorm::dmvnorm(c(.3, v), rep(mu, 2L), S2), 0.0)
  expected <- log(integrate(integrand, bounds[1L], bounds[2L],
                            rel.tol = 1e-11)$value)
  record(paste("conditional_integral", paste(bounds, collapse = ":")),
    spatial_censored_loglik_reference(m, mu), expected, 1e-10)
}

# Check a correlated four-variate mixed exact/left/right likelihood, using an
# independent deterministic orthant integration plus permutation and reflection.
S4 <- .6^abs(outer(1:4, 1:4, "-")) + diag(.15, 4L)
lo <- c(.3, -Inf, .7, -Inf); up <- c(.3, .2, Inf, .4)
m <- spatial_censored_setup_reference(lo, up, S4)
u <- spatial_censored_uniforms_reference(m, 524288L, 456L)
reference <- spatial_censored_loglik_reference(m, .1, engine = "miwa", miwa_steps = 256L)
record("simulated_vs_Miwa", spatial_censored_loglik_reference(m, .1, u),
       reference, 0.003)
p <- c(4L, 2L, 1L, 3L)
permuted <- spatial_censored_setup_reference(lo[p], up[p], S4[p, p])
record("partition_permutation", spatial_censored_loglik_reference(permuted, .1,
  engine = "miwa", miwa_steps = 256L), reference, 1e-8)
reflected <- spatial_censored_setup_reference(-up, -lo, S4)
record("reflection", spatial_censored_loglik_reference(reflected, -.1,
  engine = "miwa", miwa_steps = 256L), reference, 1e-8)

# A remote finite interval must not lose its probability through CDF cancellation.
tail_model <- spatial_censored_setup_reference(9, 10, matrix(1, 1L, 1L))
tail <- pnorm(9, lower.tail = FALSE, log.p = TRUE) +
  log1p(-exp(pnorm(10, lower.tail = FALSE, log.p = TRUE) -
             pnorm(9, lower.tail = FALSE, log.p = TRUE)))
record("remote_tail_interval", spatial_censored_loglik_reference(tail_model, 0),
       tail, 1e-12)

# Multivariate upper tails previously lost their small probability through CDF
# cancellation. Check both diagonal and correlated cases against analytic/1D
# integrals; all-zero sampled probability must not become a finite floor.
tail_model <- spatial_censored_setup_reference(c(-Inf, 9), c(10, Inf), diag(2L))
tail_u <- spatial_censored_uniforms_reference(tail_model, 4096L, 771L)
record("multivariate_remote_right_tail", spatial_censored_loglik_reference(tail_model, 0, tail_u),
  pnorm(10, log.p = TRUE) + pnorm(9, lower.tail = FALSE, log.p = TRUE), 1e-10)
rho <- .1; S_tail <- matrix(c(1, rho, rho, 1), 2L)
tail_model <- spatial_censored_setup_reference(c(-Inf, 9), c(0, Inf), S_tail)
tail_u <- spatial_censored_uniforms_reference(tail_model, 1048576L, 772L)
target <- log(integrate(function(z) dnorm(z) * pnorm(-rho * z / sqrt(1 - rho^2)),
  9, Inf, abs.tol = 0, rel.tol = 1e-10)$value)
record("correlated_remote_right_tail", spatial_censored_loglik_reference(tail_model, 0, tail_u),
       target, .003)
tail_model <- spatial_censored_setup_reference(c(-.3, 9), c(.4, 10), S_tail)
tail_u <- spatial_censored_uniforms_reference(tail_model, 1048576L, 773L)
target <- log(integrate(function(z) dnorm(z) *
  (pnorm((.4 - rho * z) / sqrt(1 - rho^2)) - pnorm((-.3 - rho * z) / sqrt(1 - rho^2))),
  9, 10, abs.tol = 0, rel.tol = 1e-10)$value)
record("correlated_remote_finite_box", spatial_censored_loglik_reference(tail_model, 0, tail_u),
       target, .003)
impossible_numerically <- spatial_censored_setup_reference(c(1000, 1000), c(Inf, Inf), S_tail)
must_fail(spatial_censored_loglik_reference(impossible_numerically, 0,
  spatial_censored_uniforms_reference(impossible_numerically, 4096L, 774L)))

# Reproducibility without perturbing the simulation generator.
set.seed(987L); old <- .Random.seed
u1 <- spatial_censored_uniforms_reference(m, 128L, 111L)
stopifnot(identical(old, .Random.seed))
u2 <- spatial_censored_uniforms_reference(m, 128L, 111L)
stopifnot(identical(u1, u2), identical(old, .Random.seed))
record("fixed_uniform_reproducibility",
  spatial_censored_loglik_reference(m, .1, u1),
  spatial_censored_loglik_reference(m, .1, u2), 0)
saved <- .Random.seed
rm(".Random.seed", envir = .GlobalEnv)
invisible(spatial_censored_loglik_reference(m, .1, u2))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
assign(".Random.seed", saved, envir = .GlobalEnv)

all_left <- spatial_censored_setup_reference(rep(-Inf, 2L), c(.6, .6), S2)
all_right <- spatial_censored_setup_reference(c(.6, .6), rep(Inf, 2L), S2)
stopifnot(spatial_censored_fit_mu_reference(all_left)$status == "no_finite_location_mle",
          spatial_censored_fit_mu_reference(all_right)$status == "no_finite_location_mle")
must_fail(spatial_censored_setup_reference(NA_real_, 1, matrix(1, 1L, 1L)))
must_fail(spatial_censored_setup_reference(-Inf, Inf, matrix(1, 1L, 1L)))
must_fail(spatial_censored_setup_reference(c(0, 0), c(0, 0), matrix(1, 2L, 2L)))
must_fail(spatial_censored_loglik_reference(m, 0, matrix(0, 2L, 10L)))

write.csv(do.call(rbind, checks), file.path(out, "verification.csv"), row.names = FALSE)
writeLines(c("Known-covariance spatial censored likelihood reference — November 1, 2026",
  "Gaussian/GLS, diagonal reduction, direct conditional integrals, Miwa comparison,",
  "permutation/reflection, remote tails, fixed uniforms/RNG preservation and invalid-input guards passed.",
  "No likelihood-based confidence interval or covariance-parameter estimator validated."),
  file.path(out, "VERIFIED.txt"))
cat(nrow(do.call(rbind, checks)), "numerical equalities passed, plus RNG and failure guards.\n")
