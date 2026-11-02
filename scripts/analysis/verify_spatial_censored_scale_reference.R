#!/usr/bin/env Rscript
source("scripts/analysis/spatial_censored_scale_reference.R")
out <- "output/spatial_scale_reference_2026-11-01"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checks <- list()
record <- function(name, actual, expected, tolerance) {
  error <- max(abs(actual - expected))
  stopifnot(is.finite(error), error <= tolerance)
  checks[[length(checks) + 1L]] <<- data.frame(check = name, error = error,
    tolerance = tolerance, passed = TRUE)
}
y <- c(-.5, .1, .9, .3)
S <- .4^abs(outer(1:4, 1:4, "-")) + diag(.2, 4L)
model <- spatial_censored_setup_reference(y, y, S)
for (multiplier in c(.4, 1, 2.5)) {
  mu <- .2
  record(paste("Gaussian_scaled_density", multiplier),
    spatial_censored_scaled_loglik_reference(model, mu / multiplier, 1 / multiplier),
    mvtnorm::dmvnorm(y, rep(mu, 4L), S * multiplier^2, log = TRUE), 1e-11)
}
gls <- sum(solve(S, y)) / sum(solve(S, rep(1, 4L)))
v <- drop(crossprod(y - gls, solve(S, y - gls))) / length(y)
fit <- spatial_censored_fit_scale_reference(model)
record("all_exact_GLS_location", fit$theta, gls, 1e-5)
record("all_exact_variance_MLE", fit$multiplier^2, v, 1e-5)

# Compare parameter rescaling with rebuilding the original covariance from scratch.
model <- spatial_censored_setup_reference(c(.2, -Inf, .7), c(.2, .3, Inf), S[1:3, 1:3])
for (multiplier in c(.6, 1, 1.8)) {
  mu <- .15
  rebuilt <- spatial_censored_setup_reference(model$lower, model$upper, model$S * multiplier^2)
  record(paste("conditional_scaled_covariance", multiplier),
    spatial_censored_scaled_loglik_reference(model, mu / multiplier, 1 / multiplier,
      engine = "miwa"),
    spatial_censored_loglik_reference(rebuilt, mu, engine = "miwa"), 1e-10)
}

# Independent censored MLE agrees with the existing survival likelihood.
z <- c(-.8, -.4, -.2, .1, .3, .5, .7, 1.1, 1.2, 1.5)
lo <- up <- z; lo[1:2] <- -Inf; up[1:2] <- -.3
lo[10L] <- 1.3; up[10L] <- Inf
model <- spatial_censored_setup_reference(lo, up, diag(length(z)))
u <- spatial_censored_uniforms_reference(model, 4096L, 481L)
fit <- spatial_censored_fit_scale_reference(model, u)
reference <- survival::survreg(survival::Surv(lo, up, type = "interval2") ~ 1, dist = "gaussian")
record("independent_censored_location", fit$theta, unname(coef(reference)), 1e-4)
record("independent_censored_scale", fit$multiplier, reference$scale, 1e-4)
model <- spatial_censored_setup_reference(c(-Inf, 0), c(-1, Inf), diag(2L))
stopifnot(spatial_censored_fit_scale_reference(model)$status == "no_finite_common_window_scale_mle")
unsupported <- spatial_censored_setup_reference(c(0, 3), c(1, 4), diag(2L))
stopifnot(inherits(tryCatch({spatial_censored_fit_scale_reference(unsupported); NULL},
  error = function(e) e), "error"))
write.csv(do.call(rbind, checks), file.path(out, "verification.csv"), row.names = FALSE)
writeLines("Fixed-shape/estimated-scale parameterization checks passed; no interval calibration.",
  file.path(out, "VERIFIED.txt"))
cat(length(checks), "fixed-shape/estimated-scale checks passed.\n")
