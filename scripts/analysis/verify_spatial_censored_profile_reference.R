#!/usr/bin/env Rscript
source("scripts/analysis/spatial_censored_profile_reference.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the fresh verification output directory")
out <- args[1L]
if (file.exists(file.path(out, "verification.csv"))) stop("Verification output already exists")
checks <- list()
record <- function(name, actual, expected, tolerance) {
  error <- max(abs(actual - expected))
  stopifnot(is.finite(error), error <= tolerance)
  checks[[length(checks) + 1L]] <<- data.frame(check = name, error = error,
    tolerance = tolerance, passed = TRUE)
}
y <- c(-.6, .1, .4, 1.2, .7)
S <- .5^abs(outer(1:5, 1:5, "-")) + diag(.15, 5L)
model <- spatial_censored_setup_reference(y, y, S)
for (mu in c(-3, -.2, .5, 4)) {
  expected <- sqrt(drop(crossprod(y - mu, solve(S, y - mu))) / length(y))
  result <- spatial_censored_profile_scale_reference(model, mu)
  stopifnot(result$status == "finite")
  record(paste0("fixed_location_scale_MLE_", mu), result$multiplier, expected, 1e-5)
  record(paste0("fixed_location_profile_density_", mu), result$loglik,
    mvtnorm::dmvnorm(y, rep(mu, length(y)), S * expected^2, log = TRUE), 1e-9)
}
boundary <- spatial_censored_profile_scale_reference(model, 1000)
stopifnot(boundary$status == "profile_search_boundary")
constant <- spatial_censored_setup_reference(rep(1, 3L), rep(1, 3L), diag(3L))
stopifnot(spatial_censored_profile_scale_reference(constant, 1)$status == "profile_search_boundary")
stopifnot(inherits(tryCatch(spatial_censored_profile_scale_reference(model, NA_real_),
  error = function(e) e), "error"))

# Two censored coordinates allow an independent deterministic Miwa rectangle
# check of the simulation profile, with common uniforms preserved across probes.
mixed <- spatial_censored_setup_reference(c(.2, -.1, -Inf, .9),
  c(.2, -.1, .1, Inf), S[1:4, 1:4])
u <- spatial_censored_uniforms_reference(mixed, 131072L, 600010001L)
set.seed(600010002L); rng <- .Random.seed
a <- spatial_censored_profile_scale_reference(mixed, .3, u)
stopifnot(identical(rng, .Random.seed))
b <- spatial_censored_profile_scale_reference(mixed, .3, engine = "miwa")
stopifnot(a$status == "finite", b$status == "finite")
record("censored_simulation_vs_Miwa_scale", a$multiplier, b$multiplier, .002)
record("censored_simulation_vs_Miwa_loglik", a$loglik, b$loglik, .002)
write.csv(do.call(rbind, checks), file.path(out, "verification.csv"), row.names = FALSE)
writeLines("10 analytic/deterministic profile comparisons and search-boundary, argument and RNG guards passed.",
  file.path(out, "VERIFIED.txt"))
cat(length(checks), "profile comparisons plus boundary/argument/RNG guards passed.\n")
