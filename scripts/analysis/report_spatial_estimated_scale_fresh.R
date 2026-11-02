#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Supply an optional study directory")
out <- if (length(args)) args[1L] else "output/spatial_scale_fresh_2026-11-01"
if (!file.exists(file.path(out, "COMPLETED.txt"))) stop("Fresh pilot incomplete")
manifest <- readRDS(file.path(out, "study_manifest.rds"))
stopifnot(identical(tools::md5sum(names(manifest$sources)), manifest$sources),
  identical(tools::md5sum(names(manifest$protocol)), manifest$protocol))
x <- read.csv(file.path(out, "fresh_results.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(x) == 64L, length(unique(x$layout_key)) == 16L,
  !anyDuplicated(x$key), all(x$n == 48L))
config <- manifest$config
cases <- lapply(file.path(out, "cases", paste0(x$key, ".rds")), readRDS)
stopifnot(identical(do.call(rbind, lapply(cases, `[[`, "row"))$key, x$key))

# Ratio-of-cluster-totals delta-method MCSE allows unequal availability.
# All independently generated layouts remain in the denominator calculation;
# two responses at one layout are never counted as independent survey layouts.
cluster_mcse <- function(value, cluster, all_clusters) {
  keep <- is.finite(value)
  numerator <- denominator <- setNames(rep(0, length(all_clusters)), all_clusters)
  for (k in all_clusters) {
    at <- keep & cluster == k
    numerator[k] <- sum(value[at]); denominator[k] <- sum(at)
  }
  m <- length(all_clusters)
  if (m < 2L || sum(denominator > 0) < 2L) return(NA_real_)
  estimate <- sum(numerator) / sum(denominator)
  sqrt(m / (m - 1L) * sum((numerator - estimate * denominator)^2)) / sum(denominator)
}
stopifnot(abs(cluster_mcse(c(1, 1, 3, 3), c("a", "a", "b", "b"),
  c("a", "b")) - 1) < 1e-12,
  abs(cluster_mcse(c(1, 1, 4), c("a", "a", "b"), c("a", "b")) - 4/3) < 1e-12,
  is.na(cluster_mcse(c(1, NA_real_), c("a", "b"), c("a", "b"))))
methods <- c(spatial_estimated_scale = "spatial", spatial_known_scale = "oracle",
  independent_equal = "equal", independent_cell = "weighted", fully_observed = "complete")
long <- do.call(rbind, lapply(names(methods), function(method) {
  prefix <- methods[[method]]
  status <- if (prefix == "complete") rep("finite", nrow(x)) else x[[paste0(prefix, "_status")]]
  theta <- x[[paste0(prefix, "_theta")]]
  sigma <- if (prefix == "oracle") rep(sqrt(config$latent_sd^2 + config$noise_sd^2), nrow(x)) else
    x[[paste0(prefix, "_sigma")]]
  available <- status == "finite" & is.finite(theta)
  data.frame(x[, c("key", "layout_key", "regime", "decay", "censoring", "n_exact",
    "public_count_eligible", "weighted_exact_fraction", "numerical_review")],
    method = method, status = status, available = available,
    theta = ifelse(available, theta, NA_real_), sigma = ifelse(available, sigma, NA_real_),
    error = ifelse(available, theta - config$truth, NA_real_))
}))
write.csv(long, file.path(out, "method_results.csv"), row.names = FALSE)
summaries <- list(); pairs <- list(); k <- 0L; q <- 0L
for (regime in config$regimes) for (severity in names(config$left_offset)) {
  subset <- x[x$regime == regime & x$censoring == severity, ]
  clusters <- sort(unique(subset$layout_key))
  stopifnot(nrow(subset) == 16L, length(clusters) == 8L)
  for (method in names(methods)) {
    z <- long[long$regime == regime & long$censoring == severity & long$method == method, ]
    error <- z$error; available <- is.finite(error)
    k <- k + 1L
    summaries[[k]] <- data.frame(regime = regime, censoring = severity, method = method,
      attempted = nrow(z), available = sum(available), unavailable = sum(!available),
      layout_clusters = length(clusters), available_layouts = length(unique(z$layout_key[available])),
      bias = if (any(available)) mean(error[available]) else NA_real_,
      bias_layout_mcse = cluster_mcse(error, z$layout_key, clusters),
      rmse = if (any(available)) sqrt(mean(error[available]^2)) else NA_real_,
      mean_absolute_error = if (any(available)) mean(abs(error[available])) else NA_real_,
      median_absolute_error = if (any(available)) median(abs(error[available])) else NA_real_,
      mean_sigma = if (any(available)) mean(z$sigma[available]) else NA_real_,
      count_ineligible_available = sum(available & !z$public_count_eligible),
      spatial_numerical_review = if (method == "spatial_estimated_scale")
        sum(available & z$numerical_review) else NA_integer_)
  }
  for (comparator in c("oracle", "equal", "weighted")) {
    ok <- subset$spatial_status == "finite" & subset[[paste0(comparator, "_status")]] == "finite"
    delta_abs <- ifelse(ok, abs(subset$spatial_theta - config$truth) -
      abs(subset[[paste0(comparator, "_theta")]] - config$truth), NA_real_)
    delta_squared <- ifelse(ok, (subset$spatial_theta - config$truth)^2 -
      (subset[[paste0(comparator, "_theta")]] - config$truth)^2, NA_real_)
    q <- q + 1L
    pairs[[q]] <- data.frame(regime = regime, censoring = severity,
      comparison = paste0("estimated_spatial_minus_", comparator),
      attempted = nrow(subset), common_finite = sum(ok),
      mean_abs_error_difference = if (any(ok)) mean(delta_abs[ok]) else NA_real_,
      abs_error_difference_layout_mcse = cluster_mcse(delta_abs, subset$layout_key, clusters),
      mean_squared_error_difference = if (any(ok)) mean(delta_squared[ok]) else NA_real_,
      squared_error_difference_layout_mcse = cluster_mcse(delta_squared, subset$layout_key, clusters))
  }
}
summary <- do.call(rbind, summaries); paired <- do.call(rbind, pairs)
write.csv(summary, file.path(out, "performance_summary.csv"), row.names = FALSE)
write.csv(paired, file.path(out, "paired_comparisons.csv"), row.names = FALSE)
# Secondary descriptive restriction uses the existing public count rule.
# It does not replace the frozen all-case summaries or certify eligibility.
eligible <- long[long$public_count_eligible, ]
eligible_groups <- split(eligible, interaction(eligible$regime, eligible$censoring,
  eligible$method, drop = TRUE))
eligible_summary <- do.call(rbind, lapply(eligible_groups, function(z) {
  ok <- is.finite(z$error)
  clusters <- sort(unique(x$layout_key[x$regime == z$regime[1L]]))
  data.frame(regime = z$regime[1L], censoring = z$censoring[1L], method = z$method[1L],
    attempted_count_eligible = nrow(z), available = sum(ok),
    bias = if (any(ok)) mean(z$error[ok]) else NA_real_,
    bias_layout_mcse = cluster_mcse(z$error, z$layout_key, clusters),
    rmse = if (any(ok)) sqrt(mean(z$error[ok]^2)) else NA_real_)
}))
write.csv(eligible_summary, file.path(out, "secondary_count_eligible_summary.csv"), row.names = FALSE)
candidates <- do.call(rbind, lapply(seq_along(cases), function(i) {
  do.call(rbind, lapply(c("coarse", "fine"), function(stage) {
    if (stage == "fine" && cases[[i]]$coarse$status != "finite") return(NULL)
    z <- cases[[i]][[stage]]$candidates
    if (!nrow(z)) return(NULL)
    cbind(key = x$key[i], stage = stage, z)
  }))
}))
write.csv(candidates, file.path(out, "optimization_starts.csv"), row.names = FALSE)
audits <- do.call(rbind, lapply(seq_along(cases), function(i) {
  audit <- cases[[i]]$audit
  if (is.null(audit) || audit$status != "finite") return(NULL)
  a <- audit$values; delta <- sweep(a[, c(1L, 3L, 4L, 5L), drop = FALSE], 1L, a[, 2L])
  data.frame(key = x$key[i], perturbation = c("mu_minus_0.02", "mu_plus_0.02",
    "log_multiplier_minus_0.02", "log_multiplier_plus_0.02"),
    mean_delta_loglik = colMeans(delta), mc_se = apply(delta, 2L, sd) / sqrt(nrow(delta)),
    central_loglik_mc_se = sd(a[, 2L]) / sqrt(nrow(a)))
}))
write.csv(audits, file.path(out, "local_integration_audit.csv"), row.names = FALSE)

long$method <- factor(long$method, levels = names(methods), labels = c(
  "Spatial\nestimated scale", "Spatial\nknown scale", "Independent\nequal weights",
  "Independent\ncell weights", "Fully\nobserved"))
long$regime <- factor(long$regime, levels = config$regimes)
long$censoring <- factor(long$censoring, levels = names(config$left_offset))
p <- ggplot2::ggplot(long[long$available, ], ggplot2::aes(x = .data$method, y = .data$error)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey45", linetype = "dashed") +
  ggplot2::geom_boxplot(outlier.shape = NA, width = .5, fill = "grey92") +
  ggplot2::geom_point(ggplot2::aes(color = factor(.data$decay)),
    position = ggplot2::position_jitter(width = .12, height = 0, seed = 560010001L), size = 1.8, alpha = .8) +
  ggplot2::facet_wrap(ggplot2::vars(regime, censoring), ncol = 2, scales = "free_y") +
  ggplot2::scale_color_manual(values = c("10" = "#2166ac", "20" = "#b35806"), name = "Spatial decay") +
  ggplot2::labs(x = NULL, y = "Estimated log location minus generating value",
    title = "Fresh survey layouts: estimating spatial scale under censoring",
    subtitle = "16 layouts; 32 fields; 48 records per field. Spatial covariance shape is still known.",
    caption = "Available fits only; every numerical-review case retained. Panels have different vertical scales.\nFully observed benchmarks repeat for paired censoring limits. No confidence interval or coverage calculation.") +
  ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out, "location_errors.png"), p, width = 12, height = 8, dpi = 180)
saveRDS(p, file.path(out, "plot.rds"))
scale_data <- long[long$available, ]
stopifnot(all(is.finite(scale_data$sigma)), all(scale_data$sigma > 0))
scale_data$error <- log(scale_data$sigma / sqrt(config$latent_sd^2 + config$noise_sd^2))
p_scale <- p + scale_data + ggplot2::labs(
  y = "Log ratio of fitted to generating marginal SD",
  title = "Fresh survey layouts: sensitivity of fitted scale",
  caption = "Zero means the fitted SD equals the generating SD. Every available fit is retained.\nKnown-scale values are supplied; fully observed benchmarks repeat across paired censoring limits.")
ggplot2::ggsave(file.path(out, "scale_errors.png"), p_scale, width = 12, height = 8, dpi = 180)
fmt <- function(z) if (is.finite(z)) sprintf("%.3f", z) else "NA"
table <- vapply(seq_len(nrow(summary)), function(i) {
  a <- summary[i, ]
  sprintf("| %s | %s | %s | %d/%d | %s | %s | %s |", a$regime, a$censoring, a$method,
    a$available, a$attempted, fmt(a$bias), fmt(a$bias_layout_mcse), fmt(a$rmse))
}, "")
writeLines(c("# Fresh estimated-scale point-estimation pilot — November 1, 2026", "",
  "The protocol was frozen before generating these responses. Sixteen fresh exogenous layouts (eight uniform/eight clustered), two spatial decays and paired censoring limits give 32 fields and 64 analyses, with 24 sites and two records/site. This is a small fresh pilot, not a powered coverage study. The prior selected-case pilot had four records/site; these are different designs.", "",
  sprintf("Estimated-scale spatial fits were finite in %d/64 attempts; %d of those meet a prespecified numerical-review flag. %d cases have fewer than two exact values and remain explicit private research cases. All failures and flagged cases are retained; unavailable fits have no finite error to plot and are counted separately below.",
    sum(x$spatial_status == "finite"), sum(x$spatial_status == "finite" & x$numerical_review), sum(!x$public_count_eligible)), "",
  "The method still receives true spatial decay and nugget/signal ratio, fitting only an overall covariance scale. It is not a fully estimated spatial model or a solution to preferential sampling. Geometry is selected independently of chemistry. The target is the generating common log location, log(50). No public API, estimate eligibility, confidence interval or package version changes follow from this pilot.", "",
  "## Availability and available-fit errors", "",
  "Bias and RMSE use available fits, whose subsets may differ by method. All attempted cases appear in availability counts. MCSE is sampling Monte Carlo uncertainty for the bias, clustered across eight independent layouts per regime; both decay responses share their layout. It is a small-sample delta-method summary, not a geochemical population confidence interval. Likelihood integration error is recorded separately.", "",
  "| Layout | Censoring | Method | Available | Bias | Bias layout MCSE | RMSE |",
  "| --- | --- | --- | ---: | ---: | ---: | ---: |", table, "",
  "The complete-data benchmark is analytic GLS plus the scale MLE on the same uncensored field; it is repeated for each censoring partner. Pairwise absolute/squared-error differences use common finite cases in paired_comparisons.csv. Negative differences favor estimated-scale spatial fitting. Point estimates with eight layouts/regime cannot establish superiority or reliable interval coverage.", "",
  "secondary_count_eligible_summary.csv is a secondary descriptive restriction to cases with at least two exact records, matching the existing public count rule. It was added to separate the two private one-exact cases from public relevance, after the frozen all-case run. It does not replace the primary table, create a new threshold, or establish validity for count-eligible samples.", "",
  "![Log-location errors](location_errors.png)", "",
  "![Marginal-scale errors](scale_errors.png)", "",
  "## Numerical audit and reproducibility", "",
  "optimization_starts.csv preserves every coarse/fine start, convergence and invalid probe. local_integration_audit.csv retains new-batch likelihood perturbations and numerical MCSEs. Numerical-review flags in fresh_results.csv do not exclude cases. Detailed coarse/refined fits and eight-batch audits are retained per case. The frozen protocol defines all integration draws, seeds and review thresholds; fingerprints guard source/protocol/runtime reuse.", "",
  "Read numerical discrepancies before interpreting an apparent statistical improvement. Numerical convergence is not accuracy or statistical validation. Any targeted refinement must preserve these original results in a separate output directory. Unknown range/nugget, response-dependent sampling, sample-support mixtures and calibrated intervals remain open. All previous simulations and original USGS data remain untouched."),
  file.path(out, "README.md"))
invisible(file.copy("scripts/analysis/report_spatial_estimated_scale_fresh.R",
  file.path(out, "source_snapshot"), overwrite = FALSE))
cat("Fresh pilot summaries, all statuses and layout-grouped MCSE reconciled.\n")
