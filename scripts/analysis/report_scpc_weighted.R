#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop("Supply optional weighted Gaussian and censored directories")
g <- if (length(args)) args[1L] else "output/scpc_weighted_gaussian_2026-11-01"
f <- if (length(args) == 2L) args[2L] else "output/scpc_weighted_censored_2026-11-01"
stopifnot(file.exists(file.path(g, "COMPLETED.txt")))
cal <- read.csv(file.path(g, "calibration_summary.csv"))
x <- read.csv(file.path(g, "gaussian_validation.csv"))
stopifnot(nrow(cal) == 4L, nrow(x) == 96L, sum(x$draws) / 3 == 1920000L)
methods <- c("weight_aware", "geometry_only", "site_cluster_t")
labels <- c("Weight-aware calibration", "Geometry-only calibration", "Site-cluster t")
cols <- c("#1768AC", "#C56718", "#62676D")
plot_coverage <- function(data, file, title, gaussian = TRUE) {
  png(file, width = 2000, height = 1600, res = 160)
  layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3L, byrow = TRUE), heights = c(1, 1, 0.12))
  par(mar = c(4.5, 4.5, 3, 1), las = 1)
  lower <- max(0, min(data$mc_low, 0.9) - 0.025)
  for (design in cal$design) {
    z <- data[data$design == design, ]
    steps <- if (gaussian) 4L else 3L
    plot(NA, xlim = c(0.7, steps + 0.3), ylim = c(lower, 1), xaxt = "n",
         xlab = "Decay / c0 (smaller = stronger dependence)", ylab = title,
         main = gsub("_", " ", design))
    axis(1, at = seq_len(steps), labels = if (gaussian) c("0.5\noutside", "1", "2", "iid") else c("1", "2", "iid"))
    abline(h = 0.95, lty = 2L, col = "#777777")
    for (m in seq_along(methods)) {
      a <- z[z$method == methods[m], ]; a <- a[order(a$actual_decay_over_c0), ]
      xx <- seq_len(nrow(a)) + (m - 2L) * 0.06
      lines(xx, a$coverage, type = "b", col = cols[m], pch = 14L + m, lwd = 1.5)
      segments(xx, a$mc_low, xx, a$mc_high, col = cols[m])
    }
  }
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", legend = labels, col = cols, pch = 15:17, horiz = TRUE, bty = "n")
  dev.off()
}
z <- x[x$phase == "confirmation", ]
plot_coverage(z, file.path(g, "confirmation_coverage.png"), "Gaussian coverage")
tab <- vapply(seq_len(nrow(cal)), function(i) sprintf("| %s | %.2f | %d | %.3f | %.3f | %.6f |",
  cal$design[i], cal$kish_weight_n[i], cal$selected_q[i], cal$critical[i],
  cal$iid_point_critical[i], cal$max_audited_rejection[i]), character(1))
inside <- z[z$method == "weight_aware" & z$inside_assumed_envelope, ]
outside <- z[z$method == "weight_aware" & !z$inside_assumed_envelope, ]
writeLines(c("# Conditional weighted Gaussian calibration — November 1, 2026", "",
  "This original research extension calibrates a weighted mean conditional on its observed locations and strictly positive weights. It is not a published SCPC guarantee for estimated weights, an inference method for preferential sampling, or a production geoCensor API.", "",
  "## Why weighting changes the calculation", "",
  "For normalized site weights a and n sites, the Gaussian estimate is theta = a'Y. Its fitted estimating-equation influences are psi = B*Y, where B = n*diag(a)*(I - 1*a'). The test numerator is sqrt(n)*a'Y, and denominator components are V'B*Y. V contains leading eigenvectors of B*Sigma(c0)*B'. The transformed covariance for calibration is W'*Sigma(c)*W, with W = [sqrt(n)*a, B'*V]. Both the numerator and weighted centering matter.", "",
  "At independence, W'*W is generally not identity, so ordinary Student-t critical values and the equal-weight expected-length formula cannot simply be reused. The independent endpoint is calibrated from its actual Gaussian covariance. Expected iid length uses numerical integration of the Gaussian quadratic norm. Algebra, weight-scale invariance, equal-weight reduction, chi-norm identities, and independent norm-simulation checks passed before the study.", "",
  "## Prespecified designs and calibration", "",
  "Fresh jittered-grid and exogenously clustered patterns contain 36 or 64 sites. Coordinate-only cell weights use width 0.25 in unit-domain coordinates and 100 origins. Aggregating weights computed for four colocated records per site was verified to reproduce the site weights. Kish weight n describes imbalance, not spatially independent replication. Weights and locations remain fixed within each simulated design.", "",
  "The baseline c0 implies 3% unweighted average pairwise exponential correlation. Every feasible component count, 1 through n-1, is considered. Expected iid interval length selects the count before response simulation. Calibration uses 31 then 61 log-decay nodes plus local searches; the selected statistic is audited on a separate 241-midpoint grid plus local searches. The Gaussian independence tail uses a total-variation bound. The compact envelope has a numerical search, not a global supremum certificate.", "",
  "| Design | Kish weight n | Selected q | Envelope critical | Iid point critical | Audited rejection maximum |",
  "| --- | ---: | ---: | ---: | ---: | ---: |", tab, "",
  "Critical values alone are not comparable across these differently scaled projection statistics. Expected iid interval length includes the scale of each denominator and is the quantity used for component selection.", "",
  "Component choices were stable under grid refinement. The transformed Gaussian probabilities and expected interval length were then independently checked using full-vector draws. All bases, critical values, weights, and comparisons were frozen before development and confirmation.", "",
  "## Independent confirmation", "",
  "Each design and covariance case has 20,000 development and 100,000 separately seeded confirmation draws: 1.92 million full-vector draws and 5.76 million interval evaluations. Actual decay is 0.5, 1, or 2 times c0, or independence. The 0.5 setting deliberately violates the assumed dependence envelope. No simulation response selected weights, tuning parameters, or covariance class.", "",
  sprintf("For the weight-aware method, confirmation coverage inside the tested envelope was %.2f–%.2f%%. Under stronger-than-allowed dependence it was %.2f–%.2f%%. These are conditional Gaussian results for these designs and assumptions.",
    100 * min(inside$coverage), 100 * max(inside$coverage), 100 * min(outside$coverage), 100 * max(outside$coverage)), "",
  "The geometry-only comparator still estimates the same weighted mean and uses weighted fitted influences, but chooses/calibrates its spatial components as if weights were equal. The site-cluster t comparator allows within-site repeats and weight imbalance in its sandwich estimate but treats distinct sites as independent. All methods are compared on the same draws.", "",
  "![Gaussian confirmation coverage](confirmation_coverage.png)", "",
  "## Limits and reproduction", "",
  "There is no trend, chemistry-dependent location selection, missing terrain, or sample-support mixture in this experiment. Conditional success does not validate uncertainty from changing site patterns or preferential exploration. The finite-domain weighted mean and the stationary population location parameter are different targets; this study estimates the latter.", "",
  "Reproduce with scripts/analysis/validate_scpc_weighted_gaussian.R and report_scpc_weighted.R. The runner refuses a nonempty output directory. Source snapshots, runtime information, component candidates, calibration/audit tables, source/input hashes, and frozen weights/bases are preserved. No production package code or full-data USGS estimate changed."), file.path(g, "README.md"))
if (!file.exists(file.path(f, "COMPLETED.txt"))) {
  cat("Weighted Gaussian report complete; censored study not yet complete.\n"); quit(status = 0L)
}
raw <- read.csv(file.path(f, "censored_results.csv"))
x <- read.csv(file.path(f, "censored_summary.csv"))
stopifnot(nrow(raw) == 8400L, nrow(x) == 72L)
z <- x[x$phase == "confirmation", ]
plot_coverage(z, file.path(f, "confirmation_coverage.png"), "Coverage over all attempts", FALSE)
a <- z[z$method == "weight_aware", ]
tab <- vapply(seq_len(nrow(a)), function(i) sprintf("| %s | %s | %.1f%% | %.1f–%.1f%% | %d/%d | %.4f |",
  a$design[i], a$actual_decay_over_c0[i], 100 * a$coverage[i], 100 * a$mc_low[i], 100 * a$mc_high[i],
  a$available[i], a$attempted[i], a$bias_log[i]), character(1))
writeLines(c("# Weighted censored-score transfer — November 1, 2026", "",
  "This independent transfer test applies the frozen weighted Gaussian calibration to the full-joint estimated influences of a weighted censored lognormal likelihood. It conditions on fixed exogenous location patterns and coordinate-derived weights. It does not validate a production USGS confidence interval.", "",
  "## Prespecified design", "",
  "The four fresh weighted Gaussian designs are reused, with weights and component counts frozen before responses. Each site has four observations. Log concentration is log(50) plus a spatial Gaussian latent component (SD 0.55) and independent observation noise (SD 0.2). Left and upper censoring thresholds are exp(log(50)-0.3) and exp(log(50)+1); the likelihood retains their limits. Latent decays are c0, 2*c0, or independence. The target is the common stationary superpopulation log location log(50), not a finite-domain spatial average.", "",
  "Each design/decay setting has 200 development and 500 independently seeded confirmation datasets: 8,400 attempted fits, 25,200 interval evaluations, and 72 summaries. Observation weights are site weights divided by four; their scale is normalized in the likelihood. The location influence includes joint location/log-scale nuisance effects. No calibration, weight, or component count was retuned using development or confirmation.", "",
  "## Independent confirmation", "",
  sprintf("Available fits: %d/%d. Failures are retained and counted as unavailable intervals rather than silently dropped.", sum(raw$status == "available"), nrow(raw)), "",
  sprintf("Weight-aware confirmation coverage across all tested censored settings was %.1f–%.1f%%. The upper end is conservative relative to nominal 95%%, not evidence of exact calibration.",
    100 * min(a$coverage), 100 * max(a$coverage)), "",
  "| Design | Latent decay/c0 | Weight-aware coverage | 95% Monte Carlo interval | Available | Log-location bias |",
  "| --- | ---: | ---: | ---: | ---: | ---: |", tab, "",
  "![Weighted censored confirmation coverage](confirmation_coverage.png)", "",
  "The geometry-only comparator uses the weighted estimate/influences but an equal-weight spatial calibration; the site-cluster t comparator treats distinct sites as independent. Identical simulated data are used for all three comparisons. Monte Carlo intervals describe uncertainty in measured coverage, not intervals for a geochemical sample. Conditional coverage among available intervals is also retained in censored_summary.csv.", "",
  "## Interpretation boundary", "",
  "The Gaussian transformation is exact for its linear weighted mean. Transferring it to nonlinear fitted censored influences is an empirical experiment, not an inherited finite-Gaussian guarantee. Weights are geometry-derived but fixed conditional on each of four location patterns; no claim about unconditional random-location uncertainty follows. Preferential sampling, spatial trends, support/method mixtures, and unobserved terrain remain separate problems. A favorable result cannot establish the amount of bias in NGDB.", "",
  "Reproduce with scripts/analysis/validate_scpc_weighted_censored.R using the completed Gaussian input and a fresh output directory. Source snapshots, source/calibration hashes, runtime information, censoring/exact-weight fractions, all fit statuses, and summaries are preserved. Reports are generated by scripts/analysis/report_scpc_weighted.R. Package R code, dependencies, version, routine intervals, and full-data USGS estimates are unchanged."), file.path(f, "README.md"))
cat("Weighted Gaussian and censored reports reconciled.\n")
