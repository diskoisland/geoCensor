#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1L] else "output/spatial_endpoints_fresh_2026-11-01"
stopifnot(file.exists(file.path(out, "COMPLETED.txt")))
manifest <- readRDS(file.path(out, "study_manifest.rds"))
stopifnot(identical(tools::md5sum(names(manifest$sources)), manifest$sources),
  identical(tools::md5sum(names(manifest$inputs)), manifest$inputs))
design <- read.csv(file.path(out, "design_manifest.csv"), stringsAsFactors = FALSE)
stopifnot(identical(unname(tools::md5sum(design$file)), design$md5))
x <- read.csv(file.path(out, "endpoint_results.csv"), stringsAsFactors = FALSE)
cases <- lapply(file.path(out, "cases", paste0(x$key, ".rds")), readRDS)
stopifnot(nrow(x) == 24L, length(unique(x$layout_key)) == 4L,
  length(unique(x$field_key)) == 12L, !anyDuplicated(x$key))
unresolved <- x$status != "finite_numerical_winner"
stopifnot(all(is.na(x$theta[unresolved])), all(is.na(x$sigma[unresolved])),
  all(is.na(x$range[unresolved])), all(is.na(x$nugget_fraction[unresolved])))
branches <- do.call(rbind, lapply(seq_along(cases), function(i) {
  fit <- cases[[i]]$fit
  do.call(rbind, lapply(c("spatial", "site_only", "independent"), function(name) {
    branch <- fit$branches[[name]]; candidate <- branch$candidate
    value <- function(key) if (is.null(candidate)) NA_real_ else candidate[[key]]
    audit <- if (is.null(fit$audit)) rep(NA_real_, 8L) else fit$audit[, name]
    data.frame(key = x$key[i], field_key = x$field_key[i], regime = x$regime[i],
      truth_structure = x$truth_structure[i], censoring = x$censoring[i],
      branch = name, status = if (is.null(branch)) "not_attempted" else branch$status,
      candidate_mu = value("mu"), candidate_sigma = value("sigma"), candidate_range = value("range"),
      candidate_nugget = value("nugget"), training_loglik = value("loglik"),
      audit_mean_loglik = if (all(is.finite(audit))) mean(audit) else NA_real_,
      numerical_mc_se = if (all(is.finite(audit))) sd(audit)/sqrt(length(audit)) else NA_real_,
      deterministic_crosscheck = if (is.null(fit$deterministic_error)) NA_real_ else fit$deterministic_error[name])
  }))
}))
stopifnot(nrow(branches) == 72L, all(is.na(branches$candidate_range[branches$branch != "spatial"])))
write.csv(branches, file.path(out, "branch_candidates.csv"), row.names = FALSE)
comparisons <- do.call(rbind, lapply(seq_along(cases), function(i) {
  z <- cases[[i]]$fit$comparison
  if (is.null(z)) return(NULL)
  cbind(key = x$key[i], z)
}))
write.csv(comparisons, file.path(out, "numerical_comparisons.csv"), row.names = FALSE)
traces <- list(); candidates <- list(); k <- 0L
for (i in seq_along(cases)) {
  fit <- cases[[i]]$fit
  for (stage in c("spatial_coarse", "spatial_refined", "site_only")) {
    branch <- switch(stage, spatial_coarse = fit$coarse,
      spatial_refined = fit$branches$spatial, site_only = fit$branches$site_only)
    if (is.null(branch) || is.null(branch$evaluations)) next
    k <- k + 1L; traces[[k]] <- cbind(key = x$key[i], stage = stage, branch$evaluations)
    candidates[[k]] <- do.call(rbind, lapply(branch$candidates, function(z)
      data.frame(key = x$key[i], stage = stage, start = z$start,
        parameter = paste(sprintf("%.17g", z$par), collapse = ";"), loglik = z$loglik,
        convergence = z$convergence, message = z$message)))
  }
}
write.csv(do.call(rbind, traces), file.path(out, "all_optimization_evaluations.csv"), row.names = FALSE)
write.csv(do.call(rbind, candidates), file.path(out, "optimizer_candidates.csv"), row.names = FALSE)
counts <- as.data.frame(with(x, table(regime, truth_structure, censoring, status)), stringsAsFactors = FALSE)
counts <- counts[counts$Freq > 0, ]; names(counts)[5L] <- "attempts"
write.csv(counts, file.path(out, "outcome_counts.csv"), row.names = FALSE)
range_policy <- data.frame(x[, c("key", "truth_structure", "status", "selected_structure", "range")],
  explanation = ifelse(x$status != "finite_numerical_winner", "overall comparison/search unresolved",
    ifelse(x$selected_structure == "spatial", "positive range in numerical spatial winner", "range not a parameter of endpoint")))
write.csv(range_policy, file.path(out, "range_reporting_audit.csv"), row.names = FALSE)
summary <- do.call(rbind, lapply(manifest$config$structures, function(structure) {
  do.call(rbind, lapply(names(manifest$config$offsets), function(severity) {
    z <- x[x$truth_structure == structure & x$censoring == severity, ]
    finite <- z$status == "finite_numerical_winner"
    data.frame(truth_structure = structure, censoring = severity, attempted = nrow(z),
      finite_numerical_winner = sum(finite), unresolved = sum(!finite),
      endpoint_winner = sum(finite & z$selected_structure != "spatial", na.rm = TRUE),
      positive_range_output = sum(is.finite(z$range)),
      available_bias = if (any(finite)) mean(z$theta[finite] - manifest$config$truth) else NA_real_,
      available_rmse = if (any(finite)) sqrt(mean((z$theta[finite] - manifest$config$truth)^2)) else NA_real_)
  }))
}))
write.csv(summary, file.path(out, "availability_summary.csv"), row.names = FALSE)

counts$regime <- factor(counts$regime, levels = manifest$config$regimes)
counts$censoring <- factor(counts$censoring, levels = names(manifest$config$offsets))
counts$truth_structure <- factor(counts$truth_structure, levels = manifest$config$structures)
p <- ggplot2::ggplot(counts, ggplot2::aes(.data$truth_structure, .data$attempts, fill = .data$status)) +
  ggplot2::geom_col(width = .7) +
  ggplot2::facet_wrap(ggplot2::vars(regime, censoring), ncol = 2) +
  ggplot2::scale_fill_manual(values = c(finite_numerical_winner = "#2166ac", comparison_unresolved = "#b35806",
    winning_branch_unresolved = "#9970ab", branch_search_unresolved = "#bdbdbd",
    no_exact_scale_identification = "#636363", numerical_failure = "#d73027", audit_unavailable = "#e08214"), name = "Outcome") +
  ggplot2::scale_y_continuous(breaks = 0:2, limits = c(0, 2)) +
  ggplot2::labs(x = "Generating covariance structure (not supplied to fit)", y = "Attempts",
    title = "Fresh endpoint pilot: keep unresolved comparisons explicit",
    subtitle = "Two layouts per geometry; paired censoring. Numerical winner is not a statistical family test.",
    caption = "Four layouts, twelve response fields and twenty-four analyses. No data are filled; missing outputs stay NA.\nConfidence limits remain experimental.") +
  ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
ggplot2::ggsave(file.path(out, "endpoint_outcomes.png"), p, width = 11, height = 8, dpi = 180)

plot_data <- rbind(data.frame(branches[, c("key", "field_key", "truth_structure", "censoring")],
  method = branches$branch, error = branches$candidate_mu - manifest$config$truth),
  data.frame(x[, c("key", "field_key", "truth_structure", "censoring")],
    method = "overall", error = x$theta - manifest$config$truth))
plot_data$method <- factor(plot_data$method, levels = c("spatial", "site_only", "independent", "overall"),
  labels = c("Spatial\ncandidate", "Site-only\ncandidate", "Noise\ncandidate", "Overall\nresult"))
plot_data$truth_structure <- factor(plot_data$truth_structure, levels = manifest$config$structures)
plot_data$censoring <- factor(plot_data$censoring, levels = names(manifest$config$offsets))
plot_data$point_type <- ifelse(x$public_count_eligible[match(plot_data$key, x$key)],
  "At least two exact records", "Private: fewer than two exact")
q <- ggplot2::ggplot(plot_data, ggplot2::aes(.data$method, .data$error, group = .data$field_key)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey45", linetype = "dashed") +
  ggplot2::geom_line(color = "grey65", alpha = .7, na.rm = TRUE) +
  ggplot2::geom_point(ggplot2::aes(shape = .data$point_type), color = "#2166ac", size = 1.8, na.rm = TRUE) +
  ggplot2::scale_shape_manual(values = c("At least two exact records" = 16,
    "Private: fewer than two exact" = 4), name = NULL) +
  ggplot2::facet_grid(truth_structure ~ censoring, scales = "free_y") +
  ggplot2::labs(x = NULL, y = "Common log location minus generating value",
    title = "Branch-specific candidates and the overall numerical result",
    subtitle = "Raw branch candidates include unresolved fits. Overall values are missing when comparison/search is unresolved.",
    caption = "Different vertical scales by row. Lines pair the same field across candidates; no population intervals.\nAvailable-only errors are conditional and do not establish superiority or reliable family detection.") +
  ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out, "endpoint_location_candidates.png"), q, width = 11, height = 10, dpi = 180)
saveRDS(list(outcomes = p, candidates = q), file.path(out, "plots.rds"))
table <- vapply(seq_len(nrow(summary)), function(i) {
  z <- summary[i, ]
  sprintf("| %s | %s | %d | %d | %d | %d |", z$truth_structure, z$censoring,
    z$attempted, z$finite_numerical_winner, z$unresolved, z$positive_range_output)
}, "")
writeLines(c("# Fresh explicit covariance endpoint pilot — November 1, 2026", "",
  "The revision and verification were frozen before fresh validation responses. Four exogenous layouts (two uniform/two clustered) carry twelve independent conditional response fields: spatial, site-only and independent-noise on each layout, paired across censoring severity. Twenty-four analyses are not twenty-four independent surveys. With only two layouts/regime this is an endpoint-behavior/numerical pilot, not calibrated confidence-limit or general performance validation.", "",
  "The fitter sees only observed bounds, coordinates and fixed integration seeds. No generating labels, true covariance or latent responses are supplied. Independent noise has no range parameter. The site-only limit permits within-site repeat correlation but has no positive-distance range parameter. Both branch ranges are NA; there is no numerical stand-in for unsupported length scales.", "",
  "Read the [scientific interpretation and limits](INTERPRETATION.md) alongside these numerical summaries.", "",
  "## Outcomes and missing ranges", "",
  sprintf("There are %d finite numerical winners, %d unresolved overall outcomes and %d finite positive-range outputs. All %d endpoint branch ranges (including unavailable candidates) are missing. Overall unresolved parameter tuples are NA. %d cases have fewer than two exact records and remain private research, not a change to public eligibility.",
    sum(x$status == "finite_numerical_winner"), sum(x$status != "finite_numerical_winner"),
    sum(is.finite(x$range)), sum(branches$branch != "spatial"), sum(!x$public_count_eligible)), "",
  "| Generating structure | Censoring | Attempts | Numerical winner | Unresolved | Positive range output |",
  "| --- | --- | ---: | ---: | ---: | ---: |", table, "",
  sprintf("A secondary description using the existing public count rule has %d eligible attempts: %d numerical winners and %d unresolved comparisons. The primary table retains all cases. Private fewer-than-two-exact cases are marked in the candidate plot, not silently removed or admitted to the public fitting rule.",
    sum(x$public_count_eligible), sum(x$public_count_eligible & x$status == "finite_numerical_winner"),
    sum(x$public_count_eligible & x$status != "finite_numerical_winner")), "",
  "![Endpoint outcomes](endpoint_outcomes.png)", "",
  "A finite numerical winner only separates the assessed candidate likelihoods at this calculation's precision. It does not establish the true covariance family, prove a global maximum, certify parameter identification or supply a confidence interval. A positive fitted range under a noise-only generating model is possible through sampling/optimization variation; no hypothesis-test interpretation is assigned. Nested endpoint equivalence and unresolved comparisons are expected.", "",
  sprintf("In this pilot, %d endpoint-generated analyses still yield positive numerical range outputs. These are descriptive counts on shared layouts, not a calibrated false-positive rate. A separate statistical evidence check is needed before treating fitted range as detected cross-site signal.",
    sum(x$truth_structure != "spatial" & is.finite(x$range))), "",
  "## What is retained", "",
  "branch_candidates.csv contains all three branch rows for every attempt, with raw parameters, training/audit loglikelihoods, numerical MCSEs and deterministic cross-checks. numerical_comparisons.csv retains gaps and numerical margins. all_optimization_evaluations.csv and optimizer_candidates.csv retain every trial, rejected search point and final start. range_reporting_audit.csv explains every range output/NA. No weaker resolved branch replaces a stronger unresolved candidate, and no missing concentration or unsampled terrain is filled.", "",
  "![Location candidates](endpoint_location_candidates.png)", "",
  "The candidate plot includes branch-specific means even when overall comparison is unresolved; those are explicitly candidates, not accepted overall estimates. Available-only errors in availability_summary.csv retain their denominators and cannot establish superiority. Deterministic endpoints can have zero integration MCSE while retaining substantial sampling uncertainty. This numerical comparison is not a covariance-family test.", "",
  "## Numerical development and limits", "",
  "Eleven analytic/factorization comparisons and endpoint/tie/missing/boundary guards passed before validation. The initial near-perfect-correlation Miwa128 failure and source are preserved in NUMERICAL_DEVELOPMENT.md and development_snapshot. The corrected site-only likelihood uses exact/conditional Gaussian formulas and bivariate TVPACK for one-sided doubly-censored pairs, checked independently with Miwa4096. General finite doubly-censored boxes and more than two records per exact coordinate are unsupported in this research reference.", "",
  "The finite spatial branch uses expanded numerical range/nugget searches, with every artificial boundary still labeled unresolved. Zero-noise singular limits and arbitrary long-range upper boundaries are not silently regularized. Branch means are compared in eight new 32,768-draw batches; endpoint scores are deterministic, and their cross-engine disagreement contributes to the numerical margin. The margin is a resolution heuristic, not an inferential cutoff. No Monte Carlo probability floor is used.", "",
  "Source/protocol/runtime hashes and separate design hashes guard initialization, parallel workers and resume. Distinct field checkpoints prevent output races. Original records, papers and prior studies are preserved. No production API, count rule, package version, prediction in unsampled terrain or confidence limit changed. Any outcome-informed revision requires another frozen study with fresh responses. Confidence limits remain experimental."),
  file.path(out, "README.md"))
invisible(file.copy("scripts/analysis/report_spatial_endpoints_fresh.R", file.path(out, "source_snapshot"), overwrite = TRUE))
cat("Endpoint branch outputs, comparison ambiguity, missing ranges and plots reported.\n")
