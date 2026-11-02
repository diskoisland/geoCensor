#!/usr/bin/env Rscript
source("scripts/analysis/spatial_censored_covariance_reference.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply a fresh verification directory")
out <- args[1L]
if (file.exists(file.path(out, "verification.csv"))) stop("Verification output already exists")
checks <- list()
record <- function(name, actual, expected, tolerance) {
  error <- max(abs(actual - expected))
  stopifnot(is.finite(error), error <= tolerance)
  checks[[length(checks) + 1L]] <<- data.frame(check = name, error = error,
    tolerance = tolerance, passed = TRUE)
}
xy <- rbind(c(0, 0), c(0, 0), c(.5, 0), c(.5, 0), c(.2, .4), c(.2, .4))
S <- spatial_covariance_matrix_reference(xy, .7, .2, .25)
record("marginal_variance", diag(S), rep(.7^2, 6L), 1e-12)
record("colocated_covariance", S[1, 2], .7^2 * .75, 1e-12)
record("separated_covariance", S[1, 3], .7^2 * .75 * exp(-.5/.2), 1e-12)
y <- c(-.4, .1, .3, 1, -.1, .5)
for (j in 1:3) {
  sigma <- c(.3, .8, 1.5)[j]; range <- c(.03, .2, 1)[j]; nu <- c(.1, .4, .9)[j]
  record(paste0("full_exact_density_", j),
    spatial_covariance_loglik_reference(y, y, xy, .2, sigma, range, nu),
    mvtnorm::dmvnorm(y, rep(.2, 6L),
      spatial_covariance_matrix_reference(xy, sigma, range, nu), log = TRUE), 1e-10)
}
lo <- c(.2, -.1, -Inf, .9); up <- c(.2, -.1, .1, Inf)
expected <- sum(dnorm(c(.2, -.1), .3, .8, log = TRUE)) +
  pnorm(.1, .3, .8, log.p = TRUE) + pnorm(.9, .3, .8, lower.tail = FALSE, log.p = TRUE)
record("independent_censored_density", spatial_covariance_loglik_reference(lo, up,
  xy[1:4, ], .3, .8, .2, 1, engine = "miwa"), expected, 1e-10)
record("no_signal_range_unidentified", spatial_covariance_matrix_reference(xy, .8, .01, 1),
  spatial_covariance_matrix_reference(xy, .8, 2, 1), 1e-12)
model <- spatial_censored_setup_reference(lo, up, diag(4L))
u <- spatial_censored_uniforms_reference(model, 131072L, 700010003L)
set.seed(700010004L); rng <- .Random.seed
simulation <- spatial_covariance_loglik_reference(lo, up, xy[1:4, ], .3, .8, .2, .3, u)
stopifnot(identical(rng, .Random.seed))
record("mixed_censoring_simulation_vs_Miwa", simulation,
  spatial_covariance_loglik_reference(lo, up, xy[1:4, ], .3, .8, .2, .3, engine = "miwa"), .002)

# Preserve the stronger boundary/unconverged candidate rather than selecting a
# lower-likelihood interior/converged candidate to make availability look good.
bounds_lo <- c(-8, -4, -4, -4); bounds_hi <- c(8, 4, 4, 4)
rows <- data.frame(start = 1:2, mu = c(0, 0), log_sigma = c(0, 0),
  log_range = c(-4, 0), logit_nugget = c(0, 0), loglik = c(10, 9), convergence = c(0L, 0L))
boundary <- spatial_covariance_select_reference(rows, bounds_lo, bounds_hi)
stopifnot(boundary$status == "search_boundary", is.na(boundary$theta),
  !is.null(boundary$candidate), boundary$candidate$loglik == 10)
rows$log_range[1L] <- 0; rows$convergence[1L] <- 1L
unresolved <- spatial_covariance_select_reference(rows, bounds_lo, bounds_hi)
stopifnot(unresolved$status == "optimizer_unresolved", is.na(unresolved$theta))
empty <- spatial_censored_fit_covariance_reference(numeric(), numeric(), matrix(numeric(), 0L, 2L))
stopifnot(empty$status == "no_data", is.na(empty$theta), is.na(empty$sigma))
none <- spatial_censored_fit_covariance_reference(c(-Inf, 1), c(0, Inf), xy[c(1, 3), ])
stopifnot(none$status == "no_exact_scale_identification", is.na(none$theta))
stopifnot(inherits(tryCatch(spatial_covariance_matrix_reference(xy, 1, .2, 0),
  error = function(e) e), "error"))

# Independent reduction of the complete-data optimization: profile location and
# scale analytically, leaving only range and nugget to optimize by L-BFGS-B.
# This verification field is development data and never a fresh pilot response.
set.seed(700010001L)
sites <- matrix(runif(24L), ncol = 2L); coordinates <- sites[rep(1:12, each = 2L), ]
covariance <- spatial_covariance_matrix_reference(coordinates, .8, .15, .3)
y <- .7 + as.vector(t(chol(covariance)) %*% rnorm(nrow(coordinates)))
fit <- spatial_censored_fit_covariance_reference(y, y, coordinates,
  maxit = 4000L, reltol = 1e-9)
profile <- function(par) {
  C <- spatial_covariance_matrix_reference(coordinates, 1, exp(par[1L]), plogis(par[2L]))
  inv <- chol2inv(chol(C)); mu <- sum(inv %*% y) / sum(inv)
  sigma <- sqrt(drop(crossprod(y - mu, inv %*% (y - mu))) / length(y))
  if (sigma < .03 || sigma > 10) stop("Verification scale outside numerical search")
  -mvtnorm::dmvnorm(y, rep(mu, length(y)), C * sigma^2, log = TRUE)
}
reference <- lapply(1:3, function(i) optim(c(log(c(.03, .2, .8)[i]), qlogis(c(.05, .3, .8)[i])),
  profile, method = "L-BFGS-B", lower = c(log(.01), qlogis(.02)),
  upper = c(log(2), qlogis(.98)), control = list(maxit = 300L, factr = 1e7)))
best <- max(vapply(reference, function(x) -x$value, 0.0))
stopifnot(!is.null(fit$candidate))
record("complete_data_four_parameter_vs_analytic_profile_max", fit$candidate$loglik, best, .005)
saveRDS(list(fit = fit, reference = reference, y = y, coordinates = coordinates,
  development_seed = 700010001L), file.path(out, "verification_development_fit.rds"))
write.csv(do.call(rbind, checks), file.path(out, "verification.csv"), row.names = FALSE)
writeLines(paste(length(checks), "numerical comparisons plus no-data, no-exact, boundary, unresolved,
singular-covariance and RNG guards passed. Development field is not validation data."),
  file.path(out, "VERIFIED.txt"))
cat(length(checks), "covariance comparisons and missing/unresolved guards passed.\n")
