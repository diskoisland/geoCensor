#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop("Supply optional Gaussian and censored result directories")
g <- if (length(args)) args[1L] else "output/scpc_envelope_validation_2026-11-01"
c <- if (length(args) == 2L) args[2L] else "output/scpc_censored_transfer_site_control_2026-11-01"
stopifnot(file.exists(file.path(g, "COMPLETED.txt")))
cal <- read.csv(file.path(g, "calibration_summary.csv"))
x <- read.csv(file.path(g, "gaussian_validation.csv"))
stopifnot(nrow(cal) == 4L, nrow(x) == 120L)
confirm <- x[x$phase == "confirmation", ]
method_names <- c("scpc_numerical_envelope", "projection_student", "iid_mean_t")
labels <- c("Calibrated projection", "Projection + ordinary t", "Ordinary mean t")
colors <- c("#1768AC", "#C56718", "#62676D")
png(file.path(g, "confirmation_coverage.png"), width = 2000, height = 1600, res = 160)
layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE), heights = c(1, 1, 0.12))
par(mar = c(5, 4.5, 3, 1), las = 1)
for (design in cal$design) {
  z <- confirm[confirm$design == design, ]
  plot(NA, xlim = c(0.7, 5.3), ylim = c(min(z$mc_low, 0.9) - 0.025, 1),
       xaxt = "n", xlab = "Decay / c0 (smaller = stronger dependence)",
       ylab = "Coverage", main = gsub("_", " ", design))
  axis(1, at = 1:5, labels = c("0.5\noutside", "1", "2", "10", "iid"))
  abline(h = 0.95, lty = 2, col = "#777777")
  for (m in seq_along(method_names)) {
    a <- z[z$method == method_names[m], ]; a <- a[order(a$actual_decay_over_c0), ]
    xx <- seq_len(nrow(a)) + (m - 2) * 0.06
    lines(xx, a$empirical_coverage, type = "b", pch = 14 + m, col = colors[m], lwd = 1.5)
    segments(xx, a$mc_low, xx, a$mc_high, col = colors[m])
  }
}
par(mar = c(0, 0, 0, 0))
plot.new(); legend("center", legend = labels, col = colors, pch = 15:17, horiz = TRUE, bty = "n")
dev.off()
inside <- confirm[confirm$method == method_names[1L] & confirm$inside_assumed_envelope, ]
outside <- confirm[confirm$method == method_names[1L] & !confirm$inside_assumed_envelope, ]
tab <- vapply(seq_len(nrow(cal)), function(i) sprintf("| %s | %d | %.3f | %d | %.3f | %.6f |",
  cal$design[i], cal$n_sites[i], cal$c0[i], cal$selected_q[i], cal$critical[i], cal$audit_max_rejection[i]), character(1))
writeLines(c("# Numerical SCPC envelope benchmark — November 1, 2026", "",
  "This original R research reference extends the earlier fixed-component formula checks to component selection and numerical covariance-envelope calibration. No production interval changed.", "",
  "## Prespecified design", "",
  "Four fixed designs use 36 or 64 jittered-grid or exogenously clustered locations. The benchmark decay c0 is chosen from locations alone to imply 3% average pairwise exponential correlation. Candidate component counts are 1–20; expected iid interval length selects q. No response data selects the covariance class or component count.", "",
  "Each candidate is calibrated over exponential decays c >= c0. The compact interval uses 31 then 61 log-spaced nodes plus local peak searches. Selected statistics are separately audited on 241 midpoint nodes plus local searches. The near-independence tail uses an explicit Gaussian total-variation bound. The compact search is numerical, not a rigorous global supremum certificate. No censored, preferential-selection, arbitrary-covariance, or estimated-weight guarantee follows.", "",
  "| Design | Locations | c0 | Selected q | Critical | Largest audited rejection |",
  "| --- | ---: | ---: | ---: | ---: | ---: |", tab, "",
  "Component choices were stable under grid refinement; no choice reached the upper component boundary. Scale/rotation invariance, calibration refinement, tail bounds, and independently computed Gaussian rejection probabilities passed. Critical values and bases were frozen before response simulation.", "",
  "## Independent simulation", "",
  "Each of four designs and five decay settings has 20,000 development and 100,000 independently seeded confirmation draws: 2.4 million full-vector Gaussian draws and 7.2 million interval evaluations across three methods. The tested decays are 0.5, 1, 2, 10 times c0 and independence. The 0.5 setting deliberately violates the dependence envelope. Development results were not used to retune the method. Analytic rejection and simulated rejection agree within the prespecified Monte Carlo tolerances; iid expected interval length also agrees.", "",
  sprintf("Confirmation coverage for the numerical-envelope method within the tested envelope ranged from %.2f%% to %.2f%%. Under the stronger-than-allowed setting it ranged from %.2f%% to %.2f%%. These are simulated Gaussian results conditional on the designs and assumptions.",
    100 * min(inside$empirical_coverage), 100 * max(inside$empirical_coverage),
    100 * min(outside$empirical_coverage), 100 * max(outside$empirical_coverage)), "",
  "![Independent confirmation coverage](confirmation_coverage.png)", "",
  "The narrower iid and ordinary-t projection comparators show what is lost by ignoring calibration. Broadening the dependence class could increase interval width; the class cannot be justified by selecting one that gives convenient significance.", "",
  "## Limits and reproduction", "",
  "The Intel/IMSL author program was not run. This is an independent numerical construction from the published formulas, not a claim of byte-for-byte replication. The source snapshot, source_manifest.csv, frozen_calibrations.rds, calibration_summary.csv, component_candidates.csv, dense_envelope_audit.csv, and gaussian_validation.csv preserve the calculation. The finite component candidate set is explicit.", "",
  "Run scripts/analysis/validate_scpc_envelope.R with an optional fresh output directory; it refuses nonempty output directories. Report with scripts/analysis/report_scpc_validation.R. Estimated censored influences require a separate transfer study, and weighted/preferential settings remain separate questions."), file.path(g, "README.md"))

if (!file.exists(file.path(c, "COMPLETED.txt"))) {
  cat("Gaussian report complete; censored study is not yet complete.\n")
  quit(status = 0L)
}
raw <- read.csv(file.path(c, "censored_results.csv"))
y <- read.csv(file.path(c, "censored_summary.csv"))
stopifnot(nrow(raw) == 8400L, nrow(y) == 72L)
z <- y[y$phase == "confirmation", ]
png(file.path(c, "confirmation_coverage.png"), width = 2000, height = 1600, res = 160)
layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE), heights = c(1, 1, 0.12))
par(mar = c(4.5, 4.5, 3, 1), las = 1)
for (design in cal$design) {
  a <- z[z$design == design, ]
  plot(NA, xlim = c(0.8, 3.2), ylim = c(min(a$mc_low, 0.85) - 0.025, 1), xaxt = "n",
       xlab = "Latent decay / c0", ylab = "Coverage over all attempts", main = gsub("_", " ", design))
  axis(1, at = 1:3, labels = c("1", "2", "iid")); abline(h = 0.95, lty = 2, col = "#777777")
  for (i in 1:3) {
    b <- a[a$method == c("scpc_equal_weight_transfer", "site_cluster_t", "working_independence")[i], ]
    b <- b[order(b$actual_decay_over_c0), ]; xx <- 1:3 + (i - 2) * 0.06
    lines(xx, b$coverage_all_attempts, type = "b", col = colors[i], pch = 14 + i, lwd = 1.5)
    segments(xx, b$mc_low, xx, b$mc_high, col = colors[i])
  }
}
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", legend = c("Frozen SCPC transfer", "Site-cluster t", "Working independence"),
       col = colors, pch = 15:17, horiz = TRUE, bty = "n")
dev.off()
transfer <- z[z$method == "scpc_equal_weight_transfer", ]
tab <- vapply(seq_len(nrow(transfer)), function(i) {
  a <- transfer[i, ]; sprintf("| %s | %s | %.1f%% | %.1f–%.1f%% | %d/%d | %.4f |",
    a$design, a$actual_decay_over_c0, 100 * a$coverage_all_attempts,
    100 * a$mc_low, 100 * a$mc_high, a$available, a$attempted, a$bias_log)
}, character(1))
writeLines(c("# Equal-weight censored-score transfer — November 1, 2026", "",
  "This first transfer study uses fixed observed locations and equal observation weights. No sites are preferentially omitted, and no estimated declustering weights are tested. Gaussian calibration and component counts remain frozen from the separate envelope benchmark.", "",
  "## Design", "",
  "The same four 36/64-site geometries each have four observations per site. Log concentration is log(50) plus a spatial Gaussian latent term (SD 0.55) and independent observation noise (SD 0.2). Latent decays are c0, 2*c0, or independence. Observations below exp(log(50)-0.3) or above exp(log(50)+1) are left/right censored; their limits remain in the likelihood. The target is the common superpopulation log location log(50), not the realized finite-domain spatial average.", "",
  "Each setting has 200 development and 500 independently seeded confirmation datasets: 8,400 attempted censored fits and 25,200 interval evaluations. The location influence uses the full joint location/log-scale bread and aggregated site scores. The numerical SCPC standard error projects these estimated site influences onto the frozen spatial components. Neither SCPC calibration nor method was retuned after either phase.", "",
  "A working-independence likelihood interval ignores both within-site and between-site dependence. A site-cluster sandwich comparator, with an n/(n-1) site correction and Student-t critical value on n-1 degrees of freedom, permits within-site dependence but assumes distinct sites are independent. This third comparator was added after inspecting the initial two-method study; its comparison is secondary. A new output directory preserves that initial study, and the original SCPC fit/interval columns are checked unchanged. Even with independent sites, four measurements at a site share a latent component, which explains undercoverage of the working-independence comparator.", "",
  sprintf("Available fits: %d/%d. The Gaussian covariance-envelope result is not assumed to transfer to nonlinear estimated influences; coverage here is measured empirically.", sum(raw$status == "available"), nrow(raw)), "",
  "## Independent confirmation", "",
  "| Design | Latent decay/c0 | Coverage | 95% Monte Carlo interval | Available | Log-location bias |",
  "| --- | ---: | ---: | ---: | ---: | ---: |", tab, "",
  "![Censored confirmation coverage](confirmation_coverage.png)", "",
  "Coverage over all attempts counts an unavailable interval as a noncoverage; conditional coverage is separately retained in censored_summary.csv. A nominal 95% method need not give exactly 95% in 500 datasets. These Monte Carlo intervals quantify simulation error, not uncertainty in a geochemical concentration estimate.", "",
  "## Decision boundary", "",
  "This is a restricted equal-weight transfer experiment. It does not justify default USGS confidence intervals, score stationarity under estimated declustering weights, or repair of preferential-sampling bias. Review the observed successes and failures without redefining the covariance class using confirmation data. Weighted and changing-location designs require fresh independent experiments.", "",
  "Reproduce with scripts/analysis/validate_scpc_censored_transfer.R (optional Gaussian input and fresh output directories) and scripts/analysis/report_scpc_validation.R. Source snapshots, input calibration hashes, runtime information, raw fit status, censoring counts, and summaries are preserved. No package API, dependency, version, or full-data estimate changed."), file.path(c, "README.md"))
cat("Gaussian and censored reports reconciled and written.\n")
