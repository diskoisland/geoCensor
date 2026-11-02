#!/usr/bin/env Rscript
source("scripts/analysis/spatial_covariance_endpoints_reference.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply fresh verification directory")
out <- args[1L]
if (file.exists(file.path(out, "verification.csv"))) stop("Verification already exists")
checks <- list()
record <- function(name, actual, expected, tolerance) {
  error <- max(abs(actual - expected)); stopifnot(is.finite(error), error <= tolerance)
  checks[[length(checks) + 1L]] <<- data.frame(check = name, error = error,
    tolerance = tolerance, passed = TRUE)
}
xy <- matrix(0, 2L, 2L)
for (nu in c(.3, .0002)) {
  left <- spatial_endpoint_setup_reference(c(-Inf, -Inf), c(0, 0), xy)
  mixed <- spatial_endpoint_setup_reference(c(0, -Inf), c(Inf, 0), xy)
  rho <- 1 - nu
  record(paste0("left_pair_closed_form_", nu),
    spatial_endpoint_site_loglik_reference(left, 0, 1, nu),
    log(.25 + asin(rho) / (2 * pi)), 1e-6)
  record(paste0("opposite_pair_closed_form_", nu),
    spatial_endpoint_site_loglik_reference(mixed, 0, 1, nu),
    log(.25 - asin(rho) / (2 * pi)), 1e-6)
}
coordinates <- rbind(c(0, 0), c(0, 0), c(1, 0), c(1, 0))
lo <- c(.2, -Inf, -.1, .9); up <- c(.2, .1, -.1, Inf)
s <- spatial_endpoint_setup_reference(lo, up, coordinates)
S <- .8^2 * ((1 - .3) * outer(c(1, 1, 2, 2), c(1, 1, 2, 2), "==") + diag(.3, 4L))
model <- spatial_censored_setup_reference(lo, up, S)
record("conditional_site_factorization_vs_global_Miwa",
  spatial_endpoint_site_loglik_reference(s, .3, .8, .3),
  spatial_censored_loglik_reference(model, .3, engine = "miwa"), 1e-8)
record("independent_endpoint_vs_site_nugget_one",
  spatial_endpoint_site_loglik_reference(s, .3, .8, 1),
  spatial_endpoint_independent_loglik_reference(s, .3, .8), 1e-12)
record("stable_univariate_right_tail_interval",
  spatial_endpoint_univariate_logmass_reference(8, 9, 0, 1),
  log(pnorm(8, lower.tail = FALSE) - pnorm(9, lower.tail = FALSE)), 1e-12)
y <- c(-.4, -.2, .7, .8)
exact <- spatial_endpoint_setup_reference(y, y, coordinates)
record("exact_site_factorization_vs_multivariate_density",
  spatial_endpoint_site_loglik_reference(exact, .3, .8, .3),
  mvtnorm::dmvnorm(y, rep(.3, 4L), S, log = TRUE), 1e-10)

# Balanced Gaussian site random effects: independent analytic variance MLEs.
means <- c(-1, -.6, .1, .5, .8, 1.4); deviation <- c(.1, .2, .15, .1, .25, .1)
id <- rep(1:6, each = 2L)
y <- means[id] + rep(c(-1, 1), 6L) * deviation[id]
coordinates <- cbind(id/6, rep(0, length(id)))
setup <- spatial_endpoint_setup_reference(y, y, coordinates)
within_var <- sum((y - means[id])^2) / 6
between_var <- 2 * sum((means - mean(means))^2) / 6
site_var <- (between_var - within_var) / 2
expected_sigma <- sqrt(site_var + within_var); expected_nu <- within_var / expected_sigma^2
fit <- spatial_endpoint_optimize_reference(function(par)
  spatial_endpoint_site_loglik_reference(setup, par[1L], exp(par[2L]), plogis(par[3L])),
  rbind(c(mean(y), log(.5), qlogis(.2)), c(mean(y), log(1), qlogis(.8))),
  c(setup$center - 8, log(.03), qlogis(.0001)), c(setup$center + 8, log(10), qlogis(.9999)),
  maxit = 3000L, reltol = 1e-10)
stopifnot(fit$status == "finite")
record("site_MLE_location", fit$candidate$par[1L], mean(y), 1e-5)
record("site_MLE_marginal_SD", exp(fit$candidate$par[2L]), expected_sigma, 1e-5)
record("site_MLE_nugget_fraction", plogis(fit$candidate$par[3L]), expected_nu, 1e-5)
saveRDS(list(fit = fit, y = y, coordinates = coordinates), file.path(out, "verification_site_fit.rds"))

branch <- function(mu, sigma, range, nugget) list(status = "finite",
  candidate = list(mu = mu, sigma = sigma, range = range, nugget = nugget))
branches <- list(spatial = branch(.1, .8, .2, .3), site_only = branch(.2, .8, NA_real_, .3),
  independent = branch(.3, .8, NA_real_, 1))
audit <- matrix(rep(c(-10, -11, -9), each = 8L), nrow = 8L)
selected <- spatial_endpoint_compare_reference(branches, audit, rep(0, 3L))
stopifnot(selected$status == "finite_numerical_winner", selected$structure == "independent",
  is.na(selected$range), selected$nugget_fraction == 1)
audit[, 2L] <- -8
selected <- spatial_endpoint_compare_reference(branches, audit, rep(0, 3L))
stopifnot(selected$structure == "site_only", is.na(selected$range))
audit[, 1L] <- -7
selected <- spatial_endpoint_compare_reference(branches, audit, rep(0, 3L))
stopifnot(selected$structure == "spatial", selected$range == .2)
branches$spatial$status <- "artificial_search_boundary"
selected <- spatial_endpoint_compare_reference(branches, audit, rep(0, 3L))
stopifnot(selected$status == "winning_branch_unresolved", is.na(selected$theta), is.na(selected$range))
branches$spatial$status <- "finite"
audit[,] <- -9
tied <- spatial_endpoint_compare_reference(branches, audit, rep(0, 3L))
stopifnot(tied$status == "comparison_unresolved", is.na(tied$range), is.na(tied$structure))
none <- spatial_censored_fit_endpoints_reference(numeric(), numeric(), matrix(numeric(), 0L, 2L), 1L, 2L, 3L)
stopifnot(none$status == "no_data", is.na(none$theta), is.na(none$range))
no_exact <- spatial_endpoint_setup_reference(c(-Inf, 1), c(0, Inf), rbind(c(0, 0), c(1, 0)))
stopifnot(no_exact$status == "no_exact_scale_identification")
stopifnot(inherits(tryCatch(spatial_endpoint_site_loglik_reference(s, 0, 1, 0), error = function(e)e), "error"),
  inherits(tryCatch(spatial_endpoint_setup_reference(rep(1, 3L), rep(1, 3L), matrix(0, 3L, 2L)),
    error = function(e)e), "error"))
write.csv(do.call(rbind, checks), file.path(out, "verification.csv"), row.names = FALSE)
writeLines(paste(length(checks), "analytic/factorization comparisons and endpoint/tie/missing/boundary guards passed."),
  file.path(out, "VERIFIED.txt"))
cat(length(checks), "endpoint comparisons and missing/ambiguity guards passed.\n")
