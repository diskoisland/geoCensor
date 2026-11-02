#!/usr/bin/env Rscript
out <- "output/spatial_scale_development_2026-11-01"
if (!file.exists(file.path(out, "COMPLETED.txt"))) stop("Selected-case development incomplete")
x <- read.csv(file.path(out, "selected_case_results.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(x) == 4L, !anyDuplicated(x$key))
cases <- lapply(file.path(out, "cases", paste0(x$key, ".rds")), readRDS)
candidates <- do.call(rbind, lapply(seq_along(cases), function(i) {
  do.call(rbind, lapply(c("coarse", "fine"), function(stage) {
    a <- cases[[i]][[stage]]$candidates
    if (!nrow(a)) return(NULL)
    cbind(key = x$key[i], stage = stage, a)
  }))
}))
write.csv(candidates, file.path(out, "optimization_starts.csv"), row.names = FALSE)
audits <- do.call(rbind, lapply(seq_along(cases), function(i) {
  a <- cases[[i]]$audit
  if (is.null(a)) return(NULL)
  delta <- sweep(a[, c(1L, 3L, 4L, 5L), drop = FALSE], 1L, a[, 2L])
  data.frame(key = x$key[i], perturbation = c("mu_minus_0.02", "mu_plus_0.02",
    "log_multiplier_minus_0.02", "log_multiplier_plus_0.02"),
    mean_delta_loglik = colMeans(delta), mc_se = apply(delta, 2L, sd) / sqrt(nrow(delta)),
    central_loglik_mc_se = sd(a[, 2L]) / sqrt(nrow(a)))
}))
if (is.null(audits)) audits <- data.frame(key = character(), perturbation = character(),
  mean_delta_loglik = numeric(), mc_se = numeric(), central_loglik_mc_se = numeric())
write.csv(audits, file.path(out, "local_integration_audit.csv"), row.names = FALSE)
fmt <- function(v) if (is.finite(v)) sprintf("%.4f", v) else "NA"
table <- vapply(seq_len(nrow(x)), function(i) {
  a <- x[i, ]
  sprintf("| %s | %d | %s | %s | %s | %s | %s |", a$key, a$n_exact, a$status,
    fmt(a$oracle_scale_theta), fmt(a$estimated_scale_theta),
    fmt(a$estimated_marginal_sigma), fmt(a$weighted_independent_theta))
}, "")
moderate <- x[x$key == "clustered_03_decay20_moderate", ]
heavy <- x[x$key == "clustered_03_decay20_heavy", ]
one <- x[x$key == "clustered_03_decay10_heavy", ]
writeLines(c("# Fixed spatial shape, estimated scale: selected-case development — November 1, 2026", "",
  "This removes supplied scale from four deliberately selected completed-pilot cases. Spatial decay and the nugget/signal ratio remain supplied at their generating values. This is neither fully estimated covariance nor independent confirmation; the cases were chosen after seeing numerical or identification problems.", "",
  "## What this development check found", "",
  sprintf("All four selected fits returned finite numerical solutions. In the paired clustered field, moderate censoring gave location %s and marginal SD %s, whereas heavy censoring gave location %s and SD %s. The generating location/SD were 3.912023/0.585235. The one-exact case instead estimated marginal SD %s. Removing supplied scale therefore leaves pronounced censoring sensitivity in these selected cases; numerical success is not a solution to the statistical failure.",
    fmt(moderate$estimated_scale_theta), fmt(moderate$estimated_marginal_sigma),
    fmt(heavy$estimated_scale_theta), fmt(heavy$estimated_marginal_sigma),
    fmt(one$estimated_marginal_sigma)), "",
  "The private reference maximizes the joint censored Gaussian likelihood in delta=mu/multiplier and kappa=1/multiplier, where covariance is the supplied shape divided by kappa squared. Complete Gaussian density/GLS/variance MLE, rescaled conditional likelihood and independent censored Gaussian fits were checked before these cases. No code is promoted to the public API and no confidence interval is calculated.", "",
  "## Retained case results", "",
  "The common generating log-location is log(50)=3.912023. Generating marginal log SD is 0.585235. Oracle-scale and working-independence columns are saved prior results on the same responses. Selection was not randomized, so these four rows must not be summarized as an unbiased performance/coverage estimate.", "",
  "| Case | Exact | Status | Known-scale spatial mu | Estimated-scale spatial mu | Estimated marginal SD | Weighted-independent mu |",
  "| --- | ---: | --- | ---: | ---: | ---: | ---: |", table, "",
  "One case has only one exact value, which the private research harness permits. The public package still requires two; this is not a change to its fitting rules.", "",
  "## Numerical evidence and limits", "",
  "optimization_starts.csv retains every coarse/refined start, convergence status, invalid objective evaluation and likelihood. Boundary or unsuccessful solutions are not accepted. local_integration_audit.csv evaluates paired likelihood changes from location +/-0.02 and log-scale multiplier +/-0.02, using eight new batches. Central loglik Monte Carlo error and paired errors are numerical diagnostics, not sampling uncertainty or confidence limits.", "",
  "A finite best start does not certify global numerical accuracy or statistical reliability. Compare coarse/refined parameters and the agreement of starts before designing a larger experiment. The observation covariance still benefits from true shape/nugget knowledge. Unknown spatial range, unknown nugget ratio, trends, mixed sample support, preferential selection and a calibrated confidence procedure remain untested.", "",
  "## Next gate", "",
  "Use this selected-case work to identify implementation or model failures only. Any method choice, regularization or covariance fit informed by these responses requires a fresh, prespecified statistical development/confirmation study. Retain the earlier heavy-censoring failure evidence and all real-USGS results unchanged. The newly added censoring-weight fractions expose observation support without claiming to solve inference."),
  file.path(out, "README.md"))
cat("Four selected cases and all optimization/integration statuses reconciled.\n")
