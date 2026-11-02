#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Supply an optional completed study directory")
out <- if (length(args)) args[1L] else "output/scpc_random_designs_2026-11-01"
if (!file.exists(file.path(out, "COMPLETED.txt"))) stop("Study is not complete; do not report partial confirmation as final")
protocol <- readRDS(file.path(out, "study_manifest.rds"))$config
g <- read.csv(file.path(out, "gaussian_results.csv"))
c <- read.csv(file.path(out, "censored_results.csv"))
l <- read.csv(file.path(out, "layout_summary.csv"))
stopifnot(nrow(l) == 200L, nrow(g) == 800L, nrow(c) == 30000L,
          sum(g$draws) == 800000L, length(unique(c$key)) == 200L)

# A layout, not an individual field, is the independent Monte Carlo unit.
layout_ci <- function(value) {
  n <- length(value); se <- sd(value) / sqrt(n); half <- qt(0.975, n - 1L) * se
  c(mean = mean(value), se = se, low = max(0, mean(value) - half), high = min(1, mean(value) + half))
}
gaussian_layout <- do.call(rbind, lapply(c("weight_aware", "site_cluster"), function(method) {
  available <- if (method == "weight_aware") g$calibration_available else rep(TRUE, nrow(g))
  analytic <- g[[paste0("analytic_", method, "_coverage")]]; analytic[!available] <- 0
  data.frame(g[c("key", "phase", "regime", "layout", "decay", "inside_envelope", "draws")],
    method = method, coverage = g[[paste0(method, "_coverage")]], analytic_coverage = analytic,
    available_fraction = as.numeric(available), width = g[[paste0(method, "_width")]])
}))
gaussian_groups <- interaction(gaussian_layout$phase, gaussian_layout$regime,
                              gaussian_layout$decay, gaussian_layout$method, drop = TRUE)
gaussian_summary <- do.call(rbind, lapply(split(gaussian_layout, gaussian_groups), function(x) {
  ci <- layout_ci(x$coverage)
  data.frame(phase = x$phase[1L], regime = x$regime[1L], decay = x$decay[1L],
    inside_envelope = x$inside_envelope[1L], method = x$method[1L], layouts = nrow(x), draws = sum(x$draws),
    coverage = ci["mean"], analytic_coverage = mean(x$analytic_coverage),
    layout_mc_se = ci["se"], mc_low = ci["low"], mc_high = ci["high"],
    available_layouts = sum(x$available_fraction), mean_width = mean(x$width, na.rm = TRUE))
}))
write.csv(gaussian_summary, file.path(out, "gaussian_summary.csv"), row.names = FALSE)

censored_groups <- interaction(c$key, c$decay, c$censoring, drop = TRUE)
censored_layout <- do.call(rbind, lapply(split(c, censored_groups), function(x) {
  do.call(rbind, lapply(c("weight_aware", "site_cluster"), function(method) {
    available <- x[[paste0(method, "_available")]]
    covered <- x[[paste0(method, "_covered")]]
    data.frame(key = x$key[1L], phase = x$phase[1L], regime = x$regime[1L], layout = x$layout[1L],
      decay = x$decay[1L], censoring = x$censoring[1L], method = method,
      attempted = nrow(x), available = sum(available), covered = sum(covered),
      coverage = mean(covered), available_fraction = mean(available),
      mean_exact = mean(x$n_exact), mean_weighted_exact_fraction = mean(x$weighted_exact_fraction),
      mean_point_bias_log = mean(x$theta - protocol$truth, na.rm = TRUE),
      mean_width = mean(x[[paste0(method, "_width")]], na.rm = TRUE))
  }))
}))
write.csv(censored_layout, file.path(out, "censored_layout_rates.csv"), row.names = FALSE)
groups <- interaction(censored_layout$phase, censored_layout$regime, censored_layout$decay,
                      censored_layout$censoring, censored_layout$method, drop = TRUE)
censored_summary <- do.call(rbind, lapply(split(censored_layout, groups), function(x) {
  ci <- layout_ci(x$coverage)
  data.frame(phase = x$phase[1L], regime = x$regime[1L], decay = x$decay[1L],
    censoring = x$censoring[1L], method = x$method[1L], layouts = nrow(x), attempted = sum(x$attempted),
    available = sum(x$available), unavailable = sum(x$attempted - x$available), coverage = ci["mean"],
    coverage_available = sum(x$covered) / sum(x$available), layout_mc_se = ci["se"],
    mc_low = ci["low"], mc_high = ci["high"], mean_exact = mean(x$mean_exact),
    mean_weighted_exact_fraction = mean(x$mean_weighted_exact_fraction),
    mean_layout_point_bias_log = mean(x$mean_point_bias_log, na.rm = TRUE),
    mean_layout_interval_width = mean(x$mean_width, na.rm = TRUE))
}))
write.csv(censored_summary, file.path(out, "censored_summary.csv"), row.names = FALSE)

# Paired differences retain the common layout and field used by both methods.
w <- censored_layout[censored_layout$method == "weight_aware", ]
s <- censored_layout[censored_layout$method == "site_cluster", ]
pair_keys <- c("key", "phase", "regime", "layout", "decay", "censoring")
paired <- merge(w[c(pair_keys, "coverage")], s[c(pair_keys, "coverage")], by = pair_keys,
                suffixes = c("_spatial", "_site"))
paired$difference <- paired$coverage_spatial - paired$coverage_site
paired_groups <- interaction(paired$phase, paired$regime, paired$decay, paired$censoring, drop = TRUE)
paired_summary <- do.call(rbind, lapply(split(paired, paired_groups), function(x) {
  se <- sd(x$difference) / sqrt(nrow(x)); half <- qt(0.975, nrow(x) - 1L) * se
  data.frame(phase = x$phase[1L], regime = x$regime[1L], decay = x$decay[1L],
    censoring = x$censoring[1L], layouts = nrow(x), coverage_difference = mean(x$difference),
    layout_mc_se = se, mc_low = mean(x$difference) - half, mc_high = mean(x$difference) + half)
}))
write.csv(paired_summary, file.path(out, "paired_coverage_difference.csv"), row.names = FALSE)
moderate <- censored_layout[censored_layout$censoring == "moderate", ]
heavy <- censored_layout[censored_layout$censoring == "heavy", ]
severity_keys <- c("key", "phase", "regime", "layout", "decay", "method")
severity_pairs <- merge(heavy[c(severity_keys, "coverage")], moderate[c(severity_keys, "coverage")],
                        by = severity_keys, suffixes = c("_heavy", "_moderate"))
severity_pairs$difference <- severity_pairs$coverage_heavy - severity_pairs$coverage_moderate
severity_groups <- interaction(severity_pairs$phase, severity_pairs$regime,
                               severity_pairs$decay, severity_pairs$method, drop = TRUE)
severity_summary <- do.call(rbind, lapply(split(severity_pairs, severity_groups), function(x) {
  se <- sd(x$difference) / sqrt(nrow(x)); half <- qt(0.975, nrow(x) - 1L) * se
  data.frame(phase = x$phase[1L], regime = x$regime[1L], decay = x$decay[1L], method = x$method[1L],
    layouts = nrow(x), heavy_minus_moderate_coverage = mean(x$difference), layout_mc_se = se,
    mc_low = mean(x$difference) - half, mc_high = mean(x$difference) + half)
}))
write.csv(severity_summary, file.path(out, "paired_censoring_difference.csv"), row.names = FALSE)
failure <- c[c$fit_status != "available", ]
if (nrow(failure)) {
  failure_table <- aggregate(list(count = failure$key), failure[c("phase", "regime", "censoring", "message")], length)
} else failure_table <- data.frame(phase = character(), regime = character(), censoring = character(), message = character(), count = integer())
write.csv(failure_table, file.path(out, "fit_failures.csv"), row.names = FALSE)

methods <- c("weight_aware", "site_cluster")
colors <- c("#1768AC", "#C56718")
legend_labels <- c("Weight-aware spatial interval", "Site-cluster t")
confirmation_g <- gaussian_summary[gaussian_summary$phase == "confirmation", ]
png(file.path(out, "gaussian_confirmation.png"), width = 1800, height = 950, res = 160)
layout(matrix(c(1, 2, 3, 3), nrow = 2L, byrow = TRUE), heights = c(1, 0.14))
par(mar = c(4.5, 4.5, 3, 1), las = 1)
lower <- max(0, min(confirmation_g$mc_low, 0.9) - 0.025)
for (regime in protocol$regimes) {
  z <- confirmation_g[confirmation_g$regime == regime, ]
  plot(NA, xlim = c(0.8, 4.2), ylim = c(lower, 1), xaxt = "n", main = paste(regime, "random surveys"),
       xlab = "Physical decay (smaller = stronger dependence)", ylab = "Gaussian coverage")
  axis(1, at = 1:4, labels = c("5\noutside", "10", "20", "iid")); abline(h = 0.95, lty = 2L, col = "#777777")
  for (m in seq_along(methods)) {
    a <- z[z$method == methods[m], ]; a <- a[order(a$decay), ]; xx <- 1:4 + (m - 1.5) * 0.06
    lines(xx, a$coverage, type = "b", col = colors[m], pch = 14L + m, lwd = 1.5)
    segments(xx, a$mc_low, xx, a$mc_high, col = colors[m])
  }
}
par(mar = c(0, 0, 0, 0)); plot.new(); legend("center", legend = legend_labels, col = colors, pch = 15:16, horiz = TRUE, bty = "n")
dev.off()
confirmation_c <- censored_summary[censored_summary$phase == "confirmation", ]
png(file.path(out, "censored_confirmation.png"), width = 2000, height = 1600, res = 160)
layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3L, byrow = TRUE), heights = c(1, 1, 0.12))
par(mar = c(4.5, 4.5, 3, 1), las = 1)
lower <- max(0, min(confirmation_c$mc_low, 0.9) - 0.025)
for (censoring in names(protocol$left_offsets)) for (regime in protocol$regimes) {
  z <- confirmation_c[confirmation_c$regime == regime & confirmation_c$censoring == censoring, ]
  plot(NA, xlim = c(0.8, 3.2), ylim = c(lower, 1), xaxt = "n", main = paste(regime, ":", censoring, "censoring"),
       xlab = "Physical latent decay", ylab = "Coverage over all attempts")
  axis(1, at = 1:3, labels = c("10", "20", "iid")); abline(h = 0.95, lty = 2L, col = "#777777")
  for (m in seq_along(methods)) {
    a <- z[z$method == methods[m], ]; a <- a[order(a$decay), ]; xx <- 1:3 + (m - 1.5) * 0.06
    lines(xx, a$coverage, type = "b", col = colors[m], pch = 14L + m, lwd = 1.5)
    segments(xx, a$mc_low, xx, a$mc_high, col = colors[m])
  }
}
par(mar = c(0, 0, 0, 0)); plot.new(); legend("center", legend = legend_labels, col = colors, pch = 15:16, horiz = TRUE, bty = "n")
dev.off()

a <- confirmation_c[confirmation_c$method == "weight_aware", ]
a <- a[order(match(a$censoring, c("moderate", "heavy")), a$regime, a$decay), ]
moderate_results <- a[a$censoring == "moderate", ]
heavy_results <- a[a$censoring == "heavy", ]
gaussian_inside <- confirmation_g[confirmation_g$method == "weight_aware" & confirmation_g$inside_envelope, ]
table_lines <- vapply(seq_len(nrow(a)), function(i) sprintf("| %s | %s | %s | %.1f%% | %.1f–%.1f%% | %d/%d |",
  a$regime[i], a$censoring[i], a$decay[i], 100 * a$coverage[i], 100 * a$mc_low[i], 100 * a$mc_high[i],
  a$available[i], a$attempted[i]), character(1))
writeLines(c("# Changing survey locations and declustering weights — November 1, 2026", "",
  "This research study tests the weighted spatial interval across new exogenous site layouts, rather than only a few fixed geometries. It adds a heavy-censoring stress test. It does not test response-dependent preferential sampling or establish a production interval for NGDB.", "",
  "## Main finding", "",
  sprintf("Across independent confirmation layouts, interval coverage was %.2f–%.2f%% under moderate censoring and %.2f–%.2f%% under heavy censoring. Under heavy censoring, coverage among available intervals was still only %.2f–%.2f%%; unavailable fits do not explain the whole shortfall. The method remains experimental, with heavy censoring an identified failure mode.",
    100 * min(moderate_results$coverage), 100 * max(moderate_results$coverage),
    100 * min(heavy_results$coverage), 100 * max(heavy_results$coverage),
    100 * min(heavy_results$coverage_available), 100 * max(heavy_results$coverage_available)), "",
  sprintf("Heavy-censoring mean layout-level log-location bias ranged from %.3f to %.3f among available fits, versus %.3f to %.3f under moderate censoring. These are log-scale biases, not arithmetic-concentration bias estimates. Spatial standard-error calibration does not itself correct this point-estimation bias.",
    min(heavy_results$mean_layout_point_bias_log), max(heavy_results$mean_layout_point_bias_log),
    min(moderate_results$mean_layout_point_bias_log), max(moderate_results$mean_layout_point_bias_log)), "",
  "## Design frozen before chemistry", "",
  "There are 20 development and 80 independent confirmation layouts for each of two sampling regimes: 200 layouts total, each with 24 sites. Uniform layouts sample the unit square. Clustered layouts have 18 random sites in [0,0.25] x [0,0.25] and six random sites across the whole window; every part of the window has positive sampling support. All locations are independent of the latent chemistry process.", "",
  "Every layout recalculates cell-declustering weights (width 0.25, 100 origins), spatial components, and numerical calibration. All feasible component counts 1–23 are considered. Each geometry/calibration is saved before response generation. Physical benchmark decay c0=10 is fixed across layouts. Actual decays are also fixed physical values, not chosen after seeing the geometry or responses. The compact envelope uses grid refinement/local searches and an independent midpoint/local-search audit; this is numerical evidence, not a global compact-interval certificate.", "",
  "## Independent Gaussian check", "",
  sprintf("There are 800,000 full-vector Gaussian draws across 200 layouts and four covariance cases, producing 1.6 million interval evaluations. Calibration was available for %d/%d layouts. The decay-5 case deliberately violates the assumed decay>=10 envelope.", sum(l$calibration_status == "available"), nrow(l)), "",
  "![Gaussian confirmation across random layouts](gaussian_confirmation.png)", "",
  sprintf("Gaussian confirmation coverage within the tested envelope ranged from %.3f%% to %.3f%%. This supports the weighted linear-Gaussian calculation across changing exogenous layouts, while the heavy-censoring transfer requires further work.",
    100 * min(gaussian_inside$coverage), 100 * max(gaussian_inside$coverage)), "",
  "gaussian_summary.csv reports both exact conditional Gaussian probabilities averaged over layouts and empirical coverage. A numerical conditional-calibration failure contributes an unavailable interval and zero coverage, not a replacement layout.", "",
  "## Censored transfer and heavier censoring", "",
  "For every layout and each latent decay (10, 20, iid), 25 new fields are generated. Each field has four measurements per site with latent log SD 0.55 and independent record-noise SD 0.2. The same field is paired across moderate/heavy censoring: left thresholds exp(log(50)-0.3) and exp(log(50)+0.6), with the common upper threshold exp(log(50)+1). Both censoring limits remain in the likelihood. There are 15,000 distinct underlying fields and 30,000 attempted censored analyses, producing 60,000 interval evaluations across two methods.", "",
  "The target is the common superpopulation geometric-mean parameter 50 (log location log(50)), not the realized spatial mean of an unsampled domain. Observation weights are site weights divided by four. No response selects the weights, covariance class, components, or calibration.", "",
  "Datasets with no exact observations at these two common bounds are marked unavailable: this particular two-parameter likelihood has no finite maximum. Other fit/calibration failures are retained. Coverage over all attempts counts every unavailable interval as noncoverage; available-interval coverage is separately reported and is subject to selection. Heavy censoring can therefore reveal both estimation failures and undercoverage.", "",
  "The research score harness can fit a finite likelihood with one exact observation; the public fit_cen_lognorm function requires at least two. This experiment is not a test of the public function's eligibility rule. No observation-count rule was tuned using confirmation results.", "",
  sprintf("Available likelihood fits: %d/%d. Calibration and fit failure details are retained in layout_summary.csv and fit_failures.csv.", sum(c$fit_status == "available"), nrow(c)), "",
  "| Sampling regime | Censoring | Latent decay | Interval coverage | 95% layout-based Monte Carlo interval | Available intervals |",
  "| --- | --- | ---: | ---: | ---: | ---: |", table_lines, "",
  "![Censored confirmation across random layouts](censored_confirmation.png)", "",
  "## Counting uncertainty correctly", "",
  "The independent Monte Carlo unit is the survey layout. First compute coverage within each layout; then average those rates over the independent layouts in each regime/decay/censoring setting. Monte Carlo standard errors are the between-layout SD divided by sqrt(number of layouts), with t-based intervals on number-of-layouts minus one degrees of freedom. Repeated fields at one layout are not treated as independent surveys. The confidence bars describe uncertainty in measured coverage, not uncertainty in a geochemical concentration estimate.", "",
  "censored_layout_rates.csv preserves layout rates, availability, exact-observation support, widths, and point bias. paired_coverage_difference.csv compares methods on the same layouts/fields with a layout-based standard error. Unavailable intervals remain counted.", "",
  "paired_censoring_difference.csv compares heavy with moderate censoring on the same underlying fields and layouts. This isolates the detection-limit stress within this experiment. Direct numerical comparisons with the earlier four-layout study also involve a different covariance-envelope specification and site count, so they should not be attributed solely to changing layouts.", "",
  "## Reproduction and limits", "",
  "Run scripts/analysis/validate_scpc_random_designs.R, then report_scpc_random_designs.R. Identical-source/runtime/configuration checkpoints can resume; changed provenance is rejected and requires a fresh directory. Source snapshots, fingerprints, seeds/protocol, coordinates, weights, calibration candidates, audits, and fit statuses are preserved. The original partial run and a rejected search-speed experiment are retained; the original search and statistical specification were kept. A matrix-operation optimization was independently checked against the original probability/critical-value source snapshot before this study.", "",
  "This is a small-design stationary lognormal experiment. It does not model variable supports, laboratory effects, geological trends, changing site counts, zero-probability inaccessible terrain, or chemistry-dependent location selection. An exogenous cluster process is not a preferential process. No package R API, dependency, version, real-USGS estimate, or routine confidence interval was changed."), file.path(out, "README.md"))
stopifnot(nrow(gaussian_summary) == 32L, nrow(censored_summary) == 48L,
          nrow(censored_layout) == 2400L, nrow(paired_summary) == 24L, nrow(severity_summary) == 24L)
cat("Random-design reports reconciled with layout-based Monte Carlo uncertainty.\n")
