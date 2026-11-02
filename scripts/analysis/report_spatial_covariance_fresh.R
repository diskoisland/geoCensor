#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Supply an optional pilot directory")
out <- if (length(args)) args[1L] else "output/spatial_covariance_fresh_2026-11-01"
if (!file.exists(file.path(out, "COMPLETED.txt"))) stop("Estimated-covariance pilot incomplete")
manifest <- readRDS(file.path(out, "study_manifest.rds"))
stopifnot(identical(tools::md5sum(names(manifest$sources)), manifest$sources),
  identical(tools::md5sum(names(manifest$inputs)), manifest$inputs))
x <- read.csv(file.path(out, "covariance_results.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(x) == 16L, length(unique(x$layout_key)) == 8L, !anyDuplicated(x$key))
unresolved <- x$fit_status != "finite_interior"
stopifnot(all(is.na(x$theta[unresolved])), all(is.na(x$sigma[unresolved])),
  all(is.na(x$range[unresolved])), all(is.na(x$nugget_fraction[unresolved])))
cases <- lapply(file.path(out, "cases", paste0(x$key, ".rds")), readRDS)
status_counts <- as.data.frame(with(x, table(regime, censoring, fit_status)), stringsAsFactors = FALSE)
status_counts <- status_counts[status_counts$Freq > 0, ]
names(status_counts)[4L] <- "attempts"
write.csv(status_counts, file.path(out, "fit_status_counts.csv"), row.names = FALSE)

evaluations <- do.call(rbind, lapply(seq_along(cases), function(i) {
  stages <- if (x$fine_attempted[i]) c("coarse", "fine") else "coarse"
  do.call(rbind, lapply(stages, function(stage) {
    z <- cases[[i]][[stage]]$evaluations
    if (!nrow(z)) return(NULL)
    cbind(key = x$key[i], stage = stage, z)
  }))
}))
candidates <- do.call(rbind, lapply(seq_along(cases), function(i) {
  stages <- if (x$fine_attempted[i]) c("coarse", "fine") else "coarse"
  do.call(rbind, lapply(stages, function(stage) {
    z <- cases[[i]][[stage]]$candidates
    if (!nrow(z)) return(NULL)
    cbind(key = x$key[i], stage = stage, z)
  }))
}))
write.csv(evaluations, file.path(out, "all_likelihood_evaluations.csv"), row.names = FALSE)
write.csv(candidates, file.path(out, "optimization_candidates.csv"), row.names = FALSE)
audits <- do.call(rbind, lapply(seq_along(cases), function(i) {
  a <- cases[[i]]$audit
  if (is.null(a) || a$status != "finite") return(NULL)
  delta <- sweep(a$values[, -1L, drop = FALSE], 1L, a$values[, 1L])
  data.frame(key = x$key[i], perturbation = c("mu_minus_0.02", "mu_plus_0.02",
    "log_sigma_minus_0.02", "log_sigma_plus_0.02", "log_range_minus_0.05",
    "log_range_plus_0.05", "logit_nugget_minus_0.05", "logit_nugget_plus_0.05"),
    mean_delta_loglik = colMeans(delta), delta_mc_se = apply(delta, 2L, sd) / sqrt(nrow(delta)),
    center_loglik_mc_se = sd(a$values[, 1L]) / sqrt(nrow(delta)))
}))
write.csv(audits, file.path(out, "independent_numerical_audit.csv"), row.names = FALSE)

# Sampling Monte Carlo summaries group responses by independent survey layout.
# With only four layouts/regime they are especially imprecise, even when finite.
cluster_mcse <- function(value, cluster, all_clusters) {
  good <- is.finite(value)
  numerator <- denominator <- setNames(rep(0, length(all_clusters)), all_clusters)
  for (k in all_clusters) {
    at <- good & cluster == k
    numerator[k] <- sum(value[at]); denominator[k] <- sum(at)
  }
  m <- length(all_clusters)
  if (m < 2L || sum(denominator > 0) < 2L) return(NA_real_)
  estimate <- sum(numerator) / sum(denominator)
  sqrt(m / (m - 1L) * sum((numerator - estimate * denominator)^2)) / sum(denominator)
}
stopifnot(abs(cluster_mcse(c(1, 1, 3, 3), c("a", "a", "b", "b"), c("a", "b")) - 1) < 1e-12,
  abs(cluster_mcse(c(1, 1, 4), c("a", "a", "b"), c("a", "b")) - 4/3) < 1e-12,
  is.na(cluster_mcse(c(1, NA_real_), c("a", "b"), c("a", "b"))))
methods <- c("estimated_covariance", "known_shape_control", "independent")
long <- do.call(rbind, lapply(methods, function(method) {
  if (method == "estimated_covariance") {
    status <- x$fit_status; theta <- x$theta; sigma <- x$sigma
  } else if (method == "known_shape_control") {
    status <- x$oracle_status; theta <- x$oracle_theta; sigma <- x$oracle_sigma
  } else {
    status <- x$independent_status; theta <- x$independent_theta; sigma <- x$independent_sigma
  }
  available <- if (method == "estimated_covariance") status == "finite_interior" else status == "finite"
  stopifnot(all(is.finite(theta[available])))
  data.frame(x[, c("key", "layout_key", "regime", "censoring", "n_exact", "public_count_eligible")],
    method = method, status = status, available = available, theta = theta, sigma = sigma,
    point_status = ifelse(method == "estimated_covariance" & !is.na(x$numerical_review) &
      x$numerical_review, "Flagged primary fit", "Other shown fits"),
    error = ifelse(available, theta - manifest$config$truth, NA_real_))
}))
write.csv(long, file.path(out, "accepted_method_results.csv"), row.names = FALSE)
summaries <- list(); pairwise <- list(); k <- q <- 0L
for (regime in manifest$config$regimes) for (severity in names(manifest$config$left_offset)) {
  selected <- x$regime == regime & x$censoring == severity
  group <- x[selected, ]; clusters <- sort(unique(group$layout_key))
  stopifnot(nrow(group) == 4L, length(clusters) == 4L)
  for (method in methods) {
    z <- long[long$regime == regime & long$censoring == severity & long$method == method, ]
    good <- is.finite(z$error); k <- k + 1L
    summaries[[k]] <- data.frame(regime = regime, censoring = severity, method = method,
      attempted = nrow(z), available = sum(good), unavailable = sum(!good),
      layout_clusters = length(clusters),
      bias_available = if (any(good)) mean(z$error[good]) else NA_real_,
      bias_layout_mcse = cluster_mcse(z$error, z$layout_key, clusters),
      rmse_available = if (any(good)) sqrt(mean(z$error[good]^2)) else NA_real_,
      accepted_count_ineligible = sum(good & !z$public_count_eligible))
  }
  for (control in c("known_shape_control", "independent")) {
    primary <- long[long$regime == regime & long$censoring == severity &
      long$method == "estimated_covariance", ]
    other <- long[long$regime == regime & long$censoring == severity & long$method == control, ]
    other <- other[match(primary$key, other$key), ]
    common <- primary$available & other$available
    delta <- ifelse(common, abs(primary$error) - abs(other$error), NA_real_)
    q <- q + 1L
    pairwise[[q]] <- data.frame(regime = regime, censoring = severity,
      comparison = paste0("estimated_covariance_minus_", control), attempted = nrow(primary),
      common_finite = sum(common), mean_absolute_error_difference = if (any(common)) mean(delta[common]) else NA_real_,
      difference_layout_mcse = cluster_mcse(delta, primary$layout_key, clusters))
  }
}
summary <- do.call(rbind, summaries)
write.csv(summary, file.path(out, "availability_and_errors.csv"), row.names = FALSE)
write.csv(do.call(rbind, pairwise), file.path(out, "common_case_comparisons.csv"), row.names = FALSE)

# The parameter figure deliberately uses separately labeled raw candidates,
# including boundaries/unconverged results; it never fills accepted NA fields.
context <- x[, c("key", "regime", "censoring", "fit_status", "numerical_review")]
raw <- do.call(rbind, list(
  data.frame(context, parameter = "Log location error",
    value = x$candidate_theta - manifest$config$truth),
  data.frame(context, parameter = "Log marginal SD ratio",
    value = log(x$candidate_sigma / x$true_sigma)),
  data.frame(context, parameter = "Log spatial range ratio",
    value = log(x$candidate_range / x$true_range)),
  data.frame(context, parameter = "Nugget fraction error",
    value = x$candidate_nugget_fraction - x$true_nugget)))
write.csv(raw, file.path(out, "raw_candidate_parameter_errors.csv"), row.names = FALSE)
raw$regime <- factor(raw$regime, levels = manifest$config$regimes)
raw$censoring <- factor(raw$censoring, levels = names(manifest$config$left_offset))
raw$parameter <- factor(raw$parameter, levels = c("Log location error", "Log marginal SD ratio",
  "Log spatial range ratio", "Nugget fraction error"))
raw$display_status <- ifelse(raw$fit_status == "finite_interior" &
  !is.na(raw$numerical_review) & raw$numerical_review, "interior_review_flag", raw$fit_status)
p <- ggplot2::ggplot(raw[is.finite(raw$value), ],
  ggplot2::aes(.data$regime, .data$value, color = .data$censoring, shape = .data$display_status)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey45", linetype = "dashed") +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = .14, height = 0, seed = 760010001L),
    size = 2.4, alpha = .9) +
  ggplot2::facet_wrap(ggplot2::vars(parameter), scales = "free_y", ncol = 2) +
  ggplot2::scale_color_manual(values = c(moderate = "#2166ac", heavy = "#b35806"), name = "Censoring") +
  ggplot2::scale_shape_manual(values = c(finite_interior = 16, search_boundary = 17,
    interior_review_flag = 4, optimizer_unresolved = 3), name = "Raw candidate status") +
  ggplot2::labs(x = "Survey geometry", y = NULL,
    title = "Fresh pilot: covariance estimated from observed bounds and coordinates",
    subtitle = "Every finite raw candidate is shown. Boundary/unresolved candidates are not accepted estimates.",
    caption = "Eight independent layouts; paired censoring; four layouts/regime. Numerical review flags remain in the tables.\nDifferent vertical scales; zero is the generating value. Confidence limits remain experimental.") +
  ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
ggplot2::ggsave(file.path(out, "covariance_candidate_errors.png"), p, width = 11, height = 8, dpi = 180)
long$method <- factor(long$method, levels = methods, labels = c("Estimated covariance", "Known-shape control", "Independent"))
long$regime <- factor(long$regime, levels = manifest$config$regimes)
long$censoring <- factor(long$censoring, levels = names(manifest$config$left_offset))
qplot <- ggplot2::ggplot(long, ggplot2::aes(.data$method, .data$error, group = .data$layout_key)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey45", linetype = "dashed") +
  ggplot2::geom_line(color = "grey65", alpha = .7, na.rm = TRUE) +
  ggplot2::geom_point(ggplot2::aes(shape = .data$point_status), color = "#2166ac", size = 2, na.rm = TRUE) +
  ggplot2::scale_shape_manual(values = c("Flagged primary fit" = 4, "Other shown fits" = 16), name = NULL) +
  ggplot2::facet_wrap(ggplot2::vars(regime, censoring), ncol = 2, scales = "free_y") +
  ggplot2::labs(x = NULL, y = "Estimated common log location minus generating value",
    title = "Location errors among computationally accepted fits",
    subtitle = "Lines pair methods on the same field. Missing estimates remain missing.",
    caption = "Availability differs by method; read availability_and_errors.csv before comparing errors.\nFinite interior fitting is not statistical validation. Four independent layouts/regime; no interval calculation.") +
  ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out, "accepted_location_errors.png"), qplot, width = 11, height = 8, dpi = 180)
saveRDS(list(candidates = p, accepted_errors = qplot), file.path(out, "plots.rds"))
fmt <- function(z) if (is.finite(z)) sprintf("%.3f", z) else "NA"
table <- vapply(seq_len(nrow(summary)), function(i) {
  z <- summary[i, ]
  sprintf("| %s | %s | %s | %d/%d | %s | %s | %s |", z$regime, z$censoring, z$method,
    z$available, z$attempted, fmt(z$bias_available), fmt(z$bias_layout_mcse), fmt(z$rmse_available))
}, "")
counts <- vapply(seq_len(nrow(status_counts)), function(i) {
  z <- status_counts[i, ]; sprintf("| %s | %s | %s | %d |", z$regime, z$censoring, z$fit_status, z$attempts)
}, "")
writeLines(c("# Fresh estimated spatial covariance pilot — November 1, 2026", "",
  "The protocol and source were frozen before generating validation responses. Eight fresh exogenous layouts/fields (four uniform/four clustered), 24 sites and two records/site give sixteen paired moderate/heavy analyses. There are eight independent surveys, not sixteen. This is a small numerical/point-estimation pilot, not validation of calibrated uncertainty or general covariance recovery.", "",
  "The primary fit estimates common log location, marginal SD, exponential range and nugget fraction jointly from observed bounds and coordinates. It receives no generating covariance, parameters, latent responses or oracle estimates. The covariance family is still specified (stationary isotropic exponential plus independent record noise); this is not unrestricted covariance estimation or a correction for preferential sampling.", "",
  "A separate [higher-precision numerical follow-up](../spatial_covariance_precision_2026-11-01/README.md) moved the sole flagged-interior case to a search-boundary outcome. The original primary categories and tables below are preserved, with flagged fits explicitly marked. Do not interpret the original thirteen interior outcomes as thirteen fully reviewed reliable estimates.", "",
  "## Missing and unresolved results", "",
  sprintf("There are %d/16 finite interior outcomes, %d search-boundary outcomes, %d optimizer-unresolved outcomes and %d other outcomes. %d cases have fewer than two exact records and remain explicit private research. No unresolved candidate was converted into an accepted estimate; accepted location/SD/range/nugget fields remain NA.",
    sum(x$fit_status == "finite_interior"), sum(x$fit_status == "search_boundary"),
    sum(x$fit_status == "optimizer_unresolved"),
    sum(!x$fit_status %in% c("finite_interior", "search_boundary", "optimizer_unresolved")),
    sum(!x$public_count_eligible)), "",
  "| Geometry | Censoring | Outcome | Attempts |", "| --- | --- | --- | ---: |", counts, "",
  "Numerical limits are not physical constraints or regularization claims. The highest finite candidate controls status before convergence/boundary classification; a weaker interior/converged candidate is never substituted for a stronger unresolved one. Raw candidate fields and all probes are retained separately. An unresolved search boundary need not imply no statistical MLE; this numerical search has not resolved it. No-data/no-exact scale failures have their own statuses. Missing observed concentrations and unsampled terrain are not filled.", "",
  "![Raw covariance candidate errors](covariance_candidate_errors.png)", "",
  "## Availability and conditional error summaries", "",
  "Bias/RMSE below use computationally accepted fits only, with explicit denominators. Finite-interior status is not scientific validity. Numerical-review flags remain in the primary rows and are not discarded. Different methods' available samples may differ; common_case_comparisons.csv pairs only common finite cases. Sampling Monte Carlo errors group by layout and are especially imprecise with four layouts/regime. They are separate from integration MCSE and are not geochemical confidence intervals.", "",
  "| Geometry | Censoring | Method | Available | Bias | Bias layout MCSE | RMSE |",
  "| --- | --- | --- | ---: | ---: | ---: | ---: |", table, "",
  "The known-shape control receives the generating covariance shape and estimates scale/location on the same censored responses; it does not initialize the primary fit. The independent control ignores spatial dependence. Neither control is automatically a valid confidence method. Raw boundary/unconverged parameter errors are diagnostics, not accepted accuracy estimates.", "",
  "![Accepted location errors](accepted_location_errors.png)", "",
  "## Numerical audit and provenance", "",
  sprintf("%d finite raw candidates receive independent numerical audits; %d audits are available. %d attempted cases carry a numerical-review flag. All finite/unavailable/outside-search trial evaluations are retained in all_likelihood_evaluations.csv, and every final optimization candidate is in optimization_candidates.csv.",
    sum(is.finite(x$candidate_theta)), sum(x$audit_status == "finite"), sum(x$numerical_review, na.rm = TRUE)), "",
  "Coarse fitting uses three data-only starts and 2,048 fixed uniforms; refinement uses 8,192 new uniforms and two starts. Eight fresh 16,384-draw batches evaluate frozen candidates and perturbations, including directions beyond search limits where mathematically valid. independent_numerical_audit.csv reports paired likelihood differences and numerical MCSEs. These checks do not prove a global maximum or sampling reliability.", "",
  "The protocol's numerical budgets refer to optim maxit controls (800 coarse, 1,200 refined); they are not a guarantee of exactly that many function evaluations. Every actual evaluation is saved.", "",
  "Ten covariance/density/optimizer comparisons and no-data/no-exact, boundary/unresolved selection, singular-covariance and RNG guards passed before fresh responses. The separate verification field is development only. Source/protocol/runtime fingerprints and snapshots guard checkpoint reuse. All original papers, USGS data and previous research outputs remain untouched.", "",
  "Confidence limits remain experimental. No CI, p-value, production function, new data exclusion, missing-cell concentration, package version or release change follows from this pilot. Any rule informed by these outcomes needs another prespecified validation with fresh responses."),
  file.path(out, "README.md"))
invisible(file.copy("scripts/analysis/report_spatial_covariance_fresh.R", file.path(out, "source_snapshot"), overwrite = FALSE))
cat("Fresh covariance statuses, missing fields, method denominators and audits reported.\n")
