#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Supply an optional output directory")
out <- if (length(args)) args[1L] else "output/spatial_likelihood_reference_2026-11-01"
if (!file.exists(file.path(out, "PILOT_COMPLETED.txt"))) stop("Pilot is incomplete")
protocol <- readRDS(file.path(out, "pilot_manifest.rds"))$config
x <- read.csv(file.path(out, "pilot_results.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(x) == 32L, !anyDuplicated(x$key),
  length(list.files(file.path(out, "layouts"), pattern = "\\.rds$")) == 8L,
  length(list.files(file.path(out, "fields"), pattern = "\\.rds$")) == 16L,
  all(x$n_exact + x$n_left + x$n_right == 96L))
methods <- c("spatial", "equal_known", "weighted_known", "weighted_unknown")
summary <- do.call(rbind, lapply(methods, function(method) {
  do.call(rbind, lapply(c("moderate", "heavy"), function(severity) {
    a <- x[x$censoring == severity, ]; theta <- a[[paste0(method, "_theta")]]
    ok <- a[[paste0(method, "_status")]] == "finite" & is.finite(theta)
    errors <- theta[ok] - protocol$truth
    data.frame(method = method, censoring = severity, attempted = nrow(a),
      available = sum(ok), unavailable = sum(!ok),
      mean_error_available = if (any(ok)) mean(errors) else NA_real_,
      rmse_available = if (any(ok)) sqrt(mean(errors^2)) else NA_real_,
      mean_abs_error_available = if (any(ok)) mean(abs(errors)) else NA_real_)
  }))
}))
write.csv(summary, file.path(out, "pilot_summary.csv"), row.names = FALSE)
paired <- do.call(rbind, lapply(methods[-1L], function(method) {
  do.call(rbind, lapply(c("moderate", "heavy"), function(severity) {
    a <- x[x$censoring == severity, ]
    ok <- a$spatial_status == "finite" & a[[paste0(method, "_status")]] == "finite" &
      is.finite(a$spatial_theta) & is.finite(a[[paste0(method, "_theta")]])
    spatial <- a$spatial_theta[ok] - protocol$truth
    other <- a[[paste0(method, "_theta")]][ok] - protocol$truth
    data.frame(comparator = method, censoring = severity, paired_cases = sum(ok),
      spatial_mean_error = if (any(ok)) mean(spatial) else NA_real_,
      comparator_mean_error = if (any(ok)) mean(other) else NA_real_,
      mean_absolute_error_difference = if (any(ok)) mean(abs(spatial) - abs(other)) else NA_real_)
  }))
}))
write.csv(paired, file.path(out, "pilot_paired_errors.csv"), row.names = FALSE)
finite <- x$spatial_status == "finite" & is.finite(x$spatial_theta)
failures <- x[!finite, c("key", "n_exact", "n_left", "n_right", "spatial_status", "message")]
write.csv(failures, file.path(out, "pilot_failures.csv"), row.names = FALSE)

png(file.path(out, "pilot_point_errors.png"), width = 1800, height = 1050, res = 150)
layout(matrix(c(1, 2, 3, 3), 2L, 2L, byrow = TRUE), heights = c(1, .13))
par(mar = c(5, 5, 4, 1))
values <- c(x$weighted_unknown_theta, x$spatial_theta) - protocol$truth
lim <- range(c(-.05, .05, values[is.finite(values)]))
lim <- lim + c(-.03, .03)
for (severity in c("moderate", "heavy")) {
  a <- x[x$censoring == severity, ]
  valid <- a$spatial_status == "finite" & a$weighted_unknown_status == "finite"
  colors <- ifelse(a$regime == "uniform", "#2166ac", "#b35806")
  plot(a$weighted_unknown_theta - protocol$truth,
       a$spatial_theta - protocol$truth, type = "n", xlim = lim, ylim = lim,
       xlab = "Weighted independence estimate: log-location error",
       ylab = "Oracle spatial estimate: log-location error",
       main = paste(severity, "censoring"))
  abline(h = 0, v = 0, col = "grey75", lty = 2)
  abline(0, 1, col = "grey50", lty = 3)
  points((a$weighted_unknown_theta - protocol$truth)[valid],
         (a$spatial_theta - protocol$truth)[valid], col = colors[valid],
         pch = ifelse(a$decay[valid] == 10, 16, 17))
  mtext(sprintf("%d/16 paired finite fits; development only", sum(valid)), side = 3,
        line = .15, cex = .8)
}
par(mar = rep(0, 4)); plot.new()
legend("center", legend = c("Uniform layout", "Clustered layout", "Decay 10", "Decay 20"),
  col = c("#2166ac", "#b35806", "grey30", "grey30"), pch = c(16, 16, 16, 17),
  horiz = TRUE, bty = "n", cex = 1)
invisible(dev.off())

fmt <- function(value) if (is.finite(value)) sprintf("%.4f", value) else "NA"
table <- vapply(seq_len(nrow(summary)), function(i) {
  a <- summary[i, ]
  sprintf("| %s | %s | %d/%d | %s | %s |", a$method, a$censoring,
    a$available, a$attempted, fmt(a$mean_error_available), fmt(a$rmse_available))
}, "")
message <- if (all(finite)) "All 32 oracle spatial fits were finite." else
  sprintf("%d of 32 oracle spatial fits were finite; %d failures remain in the record.", sum(finite), sum(!finite))
writeLines(c(
  "# Known-covariance spatial likelihood feasibility pilot — November 1, 2026", "",
  message, "",
  if (file.exists(file.path(out, "NUMERICAL_FINDING.md")))
    "**Numerical limitation identified:** the original tail calculation can return nonfinite or tolerance-floor artifacts. Read [NUMERICAL_FINDING.md](NUMERICAL_FINDING.md). This run is preserved development evidence, not an accepted likelihood benchmark." else character(),
  "",
  "This is a small development experiment, not independent confirmation or a calibrated confidence-interval method. The spatial model is given the true covariance. No confidence limits are calculated.", "",
  "## Target and design", "",
  "Eight fresh exogenous layouts (four uniform, four clustered), 24 sites and four records per site. Sixteen new fields at decays 10/20 are paired across moderate/heavy censoring, giving 32 attempted analyses. The target is log(50), the common superpopulation log location, not the realized field average. The complete protocol, seeds, runtime and frozen source snapshots are retained.", "",
  "The first comparator also knows the true marginal variance but treats records as independent. The next adds cell weights. The final comparator is the earlier weighted working-independence fit with unknown scale. Differences from that final comparator include nuisance-parameter knowledge as well as spatial modeling; they must not all be credited to dependence modeling.", "",
  "## Descriptive point errors", "",
  "| Method | Censoring | Finite/attempted | Mean log-location error | RMSE |",
  "| --- | --- | --- | --- | --- |", table, "",
  "Available-only summaries are subject to selection; all failures are retained in pilot_failures.csv and raw case files. Paired comparisons use the same available cases. There are only eight layouts: repeated fields/severities are not independent surveys. These summaries have no confirmation Monte Carlo confidence bars and do not establish superiority.", "",
  "## Numerical integration", "",
  "The observed likelihood factors the exact Gaussian density and a conditional Gaussian rectangle probability. Fixed ordinary Monte Carlo uniforms are reused across candidate locations, avoiding an independently changing objective at each optimizer call. Coarse/refined fits use 16,384/65,536 uniforms from separate seeds. Eight independent 32,768-draw batches then evaluate the refined location and offsets -0.05/0/+0.05. These audits concern numerical integration, not sampling uncertainty in the estimate.", "",
  sprintf("Across finite fits the largest absolute coarse/refined location shift was %s; the largest audit Monte Carlo SE of log likelihood was %s.",
    fmt(if (any(finite)) max(abs(x$refit_shift[finite])) else NA_real_),
    fmt(if (any(finite)) max(x$audit_loglik_mc_se[finite]) else NA_real_)), "",
  "Inspect local paired likelihood differences and their integration uncertainty in pilot_results.csv. A finite optimizer result is not a certificate of adequate integration, covariance specification or statistical calibration. Nonfinite probability evaluations are failures, not silently clipped probabilities.", "",
  "## Interpretation and next gate", "",
  "The oracle covariance isolates location estimation and numerical feasibility. Real covariance parameters remain unknown. Before any larger statistical benchmark, resolve recorded numerical failures, verify a converged probability calculation, and compare known versus estimated covariance/scale. Reuse of these pilot responses for numerical development must be labeled development; new confirmation requires fresh data/seeds. Preserve the completed weighted-SCPC failure study unchanged.", "",
  "The reference has passed exact Gaussian/GLS reduction, independent censoring reduction, direct conditional integrations, deterministic low-dimensional Miwa comparison, permutation/reflection checks, remote tails, fixed-uniform reproducibility, RNG preservation and input/failure guards. See verification.csv. These checks support the calculation; they do not validate 95% coverage.", "",
  "No public API, production confidence limits, package version, library dependencies or USGS estimates changed."),
  file.path(out, "README.md"))
cat("Pilot counts reconciled; descriptive and numerical-audit reports written.\n")
