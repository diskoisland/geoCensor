#!/usr/bin/env Rscript
# Post-hoc mechanism diagnostic, not independent statistical validation.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Supply the completed pilot directory")
pilot <- if (length(args)) args[1L] else "output/spatial_likelihood_reference_2026-11-01"
if (!file.exists(file.path(pilot, "PILOT_COMPLETED.txt"))) stop("Pilot incomplete")
out <- "output/censored_scale_identification_2026-11-01"
if (dir.exists(out) && length(list.files(out))) stop("Use a fresh diagnostic directory")
dir.create(out, recursive = TRUE)
cases <- lapply(list.files(file.path(pilot, "cases"), pattern = "\\.rds$", full.names = TRUE), readRDS)
scale <- vapply(cases, function(x) if (identical(x$working$status, "finite"))
  x$working$log_scale else -Inf, 0.0)
selected <- which.max(scale)
case <- cases[[selected]]
key <- case$row$key
field_key <- sub("_(moderate|heavy)$", "", key)
field <- readRDS(file.path(pilot, "fields", paste0(field_key, ".rds")))
w <- field$weights / sum(field$weights) * length(field$weights)
lo <- case$lower; up <- case$upper; exact <- lo == up
left <- is.infinite(lo); right <- is.infinite(up)
if (!any(left) || !any(right) || !any(exact)) stop("Diagnostic requires all three observation types")
likelihood <- function(mu, eta) {
  sigma <- exp(eta); a <- (lo - mu) / sigma; b <- (up - mu) / sigma
  values <- dnorm(lo, mu, sigma, log = TRUE)
  censored <- !exact
  positive <- a[censored] > 0
  x <- ifelse(positive, pnorm(a[censored], lower.tail = FALSE, log.p = TRUE),
              pnorm(b[censored], log.p = TRUE))
  y <- ifelse(positive, pnorm(b[censored], lower.tail = FALSE, log.p = TRUE),
              pnorm(a[censored], log.p = TRUE))
  values[censored] <- x + log(-expm1(y - x))
  sum(w * values)
}
mu_hat <- case$working$theta; eta_hat <- case$working$log_scale
maximum <- likelihood(mu_hat, eta_hat)
profile_one <- function(mu) {
  fit <- optimize(function(eta) likelihood(mu, eta), c(-6, 6), maximum = TRUE, tol = 1e-9)
  if (min(fit$maximum + 6, 6 - fit$maximum) < .01) stop("Scale profile hit numerical bound")
  c(loglik = fit$objective, sigma = exp(fit$maximum))
}
truth <- readRDS(file.path(pilot, "pilot_manifest.rds"))$config$truth
grid <- sort(unique(c(seq(mu_hat - 4, max(up[is.finite(up)]) + 1, length.out = 240L),
                      mu_hat, truth, case$row$weighted_known_theta)))
profile <- t(vapply(grid, profile_one, c(loglik = 0, sigma = 0)))
stopifnot(abs(profile_one(mu_hat)[1L] - maximum) < 1e-5,
          max(profile[, "loglik"]) <= maximum + 1e-5)
profile <- data.frame(mu = grid, loglik = profile[, "loglik"],
  relative_loglik = profile[, "loglik"] - maximum, sigma = profile[, "sigma"])
write.csv(profile, file.path(out, "location_profile.csv"), row.names = FALSE)
lower_limit <- unique(up[left]); upper_limit <- unique(lo[right])
stopifnot(length(lower_limit) == 1L, length(upper_limit) == 1L)
weighted_fraction <- c(sum(w[left]), sum(w[exact]), sum(w[right])) / sum(w)
raw_fraction <- c(mean(left), mean(exact), mean(right))
predicted <- c(pnorm(lower_limit, mu_hat, exp(eta_hat)),
  pnorm(upper_limit, mu_hat, exp(eta_hat)) - pnorm(lower_limit, mu_hat, exp(eta_hat)),
  pnorm(upper_limit, mu_hat, exp(eta_hat), lower.tail = FALSE))
counts <- data.frame(type = c("left", "exact_window", "right"),
  count = c(sum(left), sum(exact), sum(right)),
  record_fraction = raw_fraction, weight_fraction = weighted_fraction,
  fitted_normal_window_fraction = predicted)
write.csv(counts, file.path(out, "observed_and_fitted_fractions.csv"), row.names = FALSE)
q_low <- qnorm(weighted_fraction[1L]); q_high <- qnorm(1 - weighted_fraction[3L])
bin_sigma <- (upper_limit - lower_limit) / (q_high - q_low)
bin_mu <- lower_limit - bin_sigma * q_low
known_sigma <- sqrt(.55^2 + .2^2)
target_profile <- profile_one(truth)
summary <- data.frame(key = key, selection = "largest fitted scale among finite pilot working-independence fits",
  n_records = length(w), n_exact = sum(exact), weight_fraction_exact = weighted_fraction[2L],
  fitted_mu = mu_hat, fitted_sigma = exp(eta_hat), known_marginal_sigma = known_sigma,
  common_target = truth, profile_loglik_at_target_minus_max = target_profile[1L] - maximum,
  profiled_sigma_at_target = target_profile[2L], count_only_mu = bin_mu, count_only_sigma = bin_sigma)
write.csv(summary, file.path(out, "selected_case.csv"), row.names = FALSE)

png(file.path(out, "scale_identification.png"), width = 1900, height = 1000, res = 150)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
plot(profile$mu, profile$relative_loglik, type = "l", lwd = 2, col = "#2166ac",
  xlab = "Common log location (scale profiled)", ylab = "Log likelihood minus maximum",
  main = "Working-independence likelihood profile")
abline(v = mu_hat, col = "#b35806", lty = 2)
abline(v = truth, col = "grey35", lty = 3)
legend("bottomleft", legend = c("Likelihood profile", "Fitted location", "Generating target"),
  col = c("#2166ac", "#b35806", "grey35"), lty = c(1, 2, 3), bty = "n", cex = .85)
barplot(t(as.matrix(counts[, c("record_fraction", "weight_fraction", "fitted_normal_window_fraction")])),
  beside = TRUE, names.arg = c("Below lower limit", "Exact window", "Above upper limit"),
  col = c("grey65", "#2166ac", "#b35806"), ylim = c(0, 1),
  ylab = "Fraction", main = "Few exact values after declustering")
legend("topright", legend = c("Raw records", "Cell weights", "Fitted normal window mass"),
  fill = c("grey65", "#2166ac", "#b35806"), bty = "n", cex = .8)
invisible(dev.off())
writeLines(c(
  "# Heavy censoring and scale identification — November 1, 2026", "",
  "This is an exploratory mechanism check selected after the first feasibility pilot, not new confirmation. The case was selected by the largest fitted lognormal scale among finite working-independence fits, not by a favorable spatial result.", "",
  paste("Selected case:", key), "",
  sprintf("There are %d records and %d exact values. Exact values carry only %.3f%% of the cell weights. The fitted log location is %.4f and log-scale SD is %.4f, compared with the generating marginal SD %.4f.",
    length(w), sum(exact), 100 * weighted_fraction[2L], mu_hat, exp(eta_hat), known_sigma), "",
  "The common lower and upper censoring limits define a narrow exact-observation window. Once almost all weighted mass is in the two tails, the fit can reproduce those tail fractions with a very large scale and a low location. Counts alone approximately imply sigma=(upper-lower)/(qnorm(1-right_mass)-qnorm(left_mass)); a tiny exact-window mass makes that denominator small. Exact-value densities still belong in the likelihood: the count-only calculation is explanatory, not a replacement estimator.", "",
  sprintf("The count-only approximation gives location %.4f and SD %.4f. Independently profiling the actual weighted likelihood reproduces the numerical maximum to tolerance. At the generating location, the profiled log likelihood is %.4f below the maximum, with SD %.4f.",
    bin_mu, bin_sigma, target_profile[1L] - maximum, target_profile[2L]), "",
  "This shows why a finite, converged censored estimate can remain scientifically weak. It is a sampling/identification issue, separate from the spatial probability-integrator bug. Knowing the covariance or fixing scale changes the problem; an oracle pilot cannot establish that estimating those nuisance parameters will solve it.", "",
  "No likelihood-ratio cutoff, confidence interval, identifiability threshold or production data filter is proposed. Profile values and observed/fitted fractions are retained. The current low-exact-information review flag remains appropriate; no observations were deleted and no extreme estimate was clipped.", "",
  "Next: assess the location/scale interaction under a declared estimated-scale spatial benchmark before considering routine interval changes. Preserve this selected-case diagnostic as exploratory evidence."),
  file.path(out, "README.md"))
saveRDS(list(pilot = pilot, key = key, original_case = case, field = field), file.path(out, "selected_source.rds"))
file.copy("scripts/analysis/diagnose_censored_scale_profile.R", out)
cat("Selected-case scale profile and censoring-mass diagnostic complete.\n")
