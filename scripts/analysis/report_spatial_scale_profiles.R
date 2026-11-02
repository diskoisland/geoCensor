#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Supply an optional diagnostic directory")
out <- if (length(args)) args[1L] else "output/spatial_scale_profiles_2026-11-01"
if (!file.exists(file.path(out, "COMPLETED.txt"))) stop("Selected profiles incomplete")
manifest <- readRDS(file.path(out, "diagnostic_manifest.rds"))
stopifnot(identical(tools::md5sum(names(manifest$sources)), manifest$sources),
  identical(tools::md5sum(names(manifest$inputs)), manifest$inputs))
s <- read.csv(file.path(out, "case_summary.csv"))
g <- read.csv(file.path(out, "location_profiles.csv"))
a <- read.csv(file.path(out, "independent_anchor_audit.csv"))
refined <- read.csv(file.path(out, "anchor_refinement.csv"))
stopifnot(nrow(s) == 4L, nrow(a) == 48L, nrow(refined) == 24L,
  !anyDuplicated(s$key), all(s$n_exact >= 2L), all(s$n == 48L))
cases <- lapply(file.path(out, "cases", paste0(s$key, ".rds")), readRDS)
pilot <- dirname(names(manifest$inputs)[basename(names(manifest$inputs)) == "fresh_results.csv"])
stopifnot(length(pilot) == 1L)
support <- do.call(rbind, lapply(s$key, function(key) {
  original <- readRDS(file.path(pilot, "cases", paste0(key, ".rds")))
  field <- readRDS(file.path(pilot, "fields", paste0(sub("_(heavy|moderate)$", "", key), ".rds")))
  exact <- original$lower == original$upper
  data.frame(key = key, n_exact_records = sum(exact),
    n_exact_sites = length(unique(field$id[exact])), n_sites = length(unique(field$id)),
    weighted_exact_fraction = sum(field$weights[exact]))
}))
write.csv(support, file.path(out, "exact_site_support.csv"), row.names = FALSE)
evaluations <- do.call(rbind, lapply(seq_along(cases), function(i) {
  do.call(rbind, lapply(seq_along(cases[[i]]$profiles), function(j) {
    cbind(key = s$key[i], mu = cases[[i]]$grid$mu[j], cases[[i]]$profiles[[j]]$evaluations)
  }))
}))
write.csv(evaluations, file.path(out, "all_profile_evaluations.csv"), row.names = FALSE)
starts <- do.call(rbind, lapply(seq_along(cases), function(i)
  cbind(key = s$key[i], cases[[i]]$reference$candidates)))
write.csv(starts, file.path(out, "reference_optimization_starts.csv"), row.names = FALSE)

# Each curve uses the same free-model reference. Fixed-scale likelihoods are
# deliberately not re-centered at their own maxima.
plot_data <- rbind(
  data.frame(g[, c("key", "regime", "censoring", "mu")],
    model = "Scale profiled", delta = ifelse(g$status == "finite", g$loglik_minus_reference, NA_real_)),
  data.frame(g[, c("key", "regime", "censoring", "mu")],
    model = "True scale fixed", delta = ifelse(g$fixed_scale_status == "finite", g$fixed_scale_minus_reference, NA_real_)))
plot_data$regime <- factor(plot_data$regime, levels = c("uniform", "clustered"))
plot_data$censoring <- factor(plot_data$censoring, levels = c("moderate", "heavy"))
s$regime <- factor(s$regime, levels = c("uniform", "clustered"))
s$censoring <- factor(s$censoring, levels = c("moderate", "heavy"))
colors <- c("Scale profiled" = "#2166ac", "True scale fixed" = "#b35806")
p <- ggplot2::ggplot(plot_data, ggplot2::aes(.data$mu, .data$delta, color = .data$model)) +
  ggplot2::geom_line(linewidth = .8, na.rm = TRUE) +
  ggplot2::geom_vline(xintercept = manifest$config$truth, linetype = "dotted", color = "grey30") +
  ggplot2::geom_vline(data = s, ggplot2::aes(xintercept = .data$reference_mu),
    linetype = "dashed", color = "#2166ac", inherit.aes = FALSE) +
  ggplot2::facet_wrap(ggplot2::vars(regime, censoring), ncol = 2) +
  ggplot2::scale_color_manual(values = colors, name = NULL) +
  ggplot2::coord_cartesian(ylim = c(-8, .1)) +
  ggplot2::labs(x = "Proposed common log location", y = "Log likelihood minus free-model reference",
    title = "Selected survey fields: location/scale likelihood tradeoff",
    subtitle = "Solid curves share one reference maximum. Dotted: generating location; dashed: fitted location.",
    caption = "Post hoc selected cases; covariance shape is known. No confidence cutoff or interval.\nVertical view is limited to eight loglik units; full unclipped evaluations and search statuses are saved.") +
  ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out, "location_profiles.png"), p, width = 11, height = 7.5, dpi = 180)

g$regime <- factor(g$regime, levels = c("uniform", "clustered"))
g$censoring <- factor(g$censoring, levels = c("moderate", "heavy"))
g$plot_sigma <- ifelse(g$status == "finite", g$marginal_sigma, NA_real_)
q <- ggplot2::ggplot(g, ggplot2::aes(.data$mu, .data$plot_sigma)) +
  ggplot2::geom_line(color = "#2166ac", linewidth = .8, na.rm = TRUE) +
  ggplot2::geom_hline(yintercept = sqrt(.55^2 + .2^2), linetype = "dotted", color = "grey35") +
  ggplot2::geom_vline(xintercept = manifest$config$truth, linetype = "dotted", color = "grey35") +
  ggplot2::geom_point(data = s, ggplot2::aes(x = .data$reference_mu,
    y = .data$reference_marginal_sigma), inherit.aes = FALSE, color = "#b35806", size = 2) +
  ggplot2::facet_wrap(ggplot2::vars(regime, censoring), ncol = 2) +
  ggplot2::scale_y_log10() +
  ggplot2::labs(x = "Proposed common log location", y = "Profiled marginal log-response SD (log axis)",
    title = "How the fitted scale compensates for a changed location",
    subtitle = "Orange: joint fit. Dotted lines: generating location and marginal SD.",
    caption = "All finite main-grid profiles shown, including very poor likelihoods.\nThese are selected-case diagnostics, not population intervals or a filtering rule.") +
  ggplot2::theme_minimal(base_size = 11)
ggplot2::ggsave(file.path(out, "profiled_scale.png"), q, width = 11, height = 7.5, dpi = 180)
saveRDS(list(location = p, scale = q), file.path(out, "plots.rds"))

fmt <- function(x) sprintf("%.4f", x)
table <- vapply(seq_len(nrow(s)), function(i) {
  x <- s[i, ]; at_truth <- a[a$key == x$key & a$anchor == "generating_location" & a$mode == "profiled_scale", ]
  stopifnot(nrow(at_truth) == 1L)
  sprintf("| %s | %d | %d | %s | %s | %s | %s |", x$key, x$n_exact,
    support$n_exact_sites[match(x$key, support$key)],
    fmt(x$reference_mu), fmt(x$reference_marginal_sigma),
    fmt(at_truth$mean_loglik_delta), fmt(at_truth$delta_mc_se))
}, "")
writeLines(c("# Selected spatial location/scale profiles — November 1, 2026", "",
  "This is post hoc mechanism development on two completed fresh-pilot fields, with both moderate/heavy censoring versions retained. In each geometry regime we chose the largest absolute heavy-censoring location error among cases with at least two exact records. Selection uses the known generating value: it cannot estimate the frequency of weak support or serve as independent confirmation.", "",
  "A finite fit can have broad, asymmetric finite-sample likelihood support. That is distinct from proving structural non-identifiability of the Gaussian model. At each proposed common log location, scale is re-optimized. The spatial covariance shape and nugget/signal ratio still use generating values. Spatial profiles use the joint covariance likelihood; cell-weight support shares are context only, not likelihood weights.", "",
  "## Selected fits and independent numerical audit", "",
  "The generating location is log(50)=3.912023 and marginal log-response SD is 0.585235. The final two columns evaluate frozen, refined profile parameters in eight new integration batches. Differences are relative to the profiled-scale reference location, and MCSE concerns numerical integration only. They are not likelihood-ratio p-values, confidence intervals or sampling standard errors.", "",
  "| Case | Exact records | Exact sites | Refit location | Marginal SD | Loglik delta at generating location | Paired numerical MCSE |",
  "| --- | ---: | ---: | ---: | ---: | ---: | ---: |", table, "",
  "Site IDs here are known simulation locations. The two selected heavy cases have two and four distinct exact sites, so their sparse information is not solely duplicate exact records. These descriptive counts introduce no new threshold. Real-data coordinate matches would require separate interpretation of rounding and sample support.", "",
  "Both heavy cases also fall below the existing default ten-exact-observation review convention. They pass the minimum-two-exact fitting count but would still receive a few-exact-information warning. These profiles investigate that warning; neither passing nor failing a count convention is a calibrated inference test.", "",
  "![Location likelihood profiles](location_profiles.png)", "",
  "The fixed-scale and profiled curves share the free-model reference maximum; the fixed-scale curve is not independently re-centered. Values outside the displayed vertical range remain in location_profiles.csv. Lines show approximate numerical profiles over the evaluated grid, not confidence bands or calibrated intervals.", "",
  "![Profiled covariance scale](profiled_scale.png)", "",
  "The scale plot includes poor-likelihood locations too. Its curves show how scale changes to accommodate each location; they do not assert that every plotted combination is equally plausible. Read it together with the likelihood-loss plot.", "",
  "## Checks and retained evidence", "",
  sprintf("The main grid has %d points across four cases: %d finite profiles, %d search-boundary outcomes, %d other outcomes. %d unavailable profile probes remain in all_profile_evaluations.csv. Probe failure is distinct from an unavailable final profile; no likelihood was clipped to a finite floor.",
    nrow(g), sum(s$n_finite_profiles), sum(s$n_boundary_profiles), sum(s$n_other_profiles), sum(s$invalid_profile_probes)), "",
  sprintf("The fixed-scale comparison has %d unavailable evaluations, recorded explicitly in location_profiles.csv rather than replaced with finite losses.",
    sum(g$fixed_scale_status != "finite")), "",
  "Ten analytic/deterministic comparisons and boundary/argument/RNG guards passed before calculating these cases. The main curves use 32,768 common uniforms per case. Six anchors are re-profiled with 65,536 new uniforms, and their parameters are frozen before eight independent 32,768-draw audits. All starts, scalar evaluations, source snapshots, source/input fingerprints and search statuses are retained. A reported numerical search boundary is not a regularized or statistically bounded estimate.", "",
  "case_summary.csv retains source/refit changes. anchor_refinement.csv separates parameter changes from likelihood gains evaluated under the same new uniforms. independent_anchor_audit.csv retains all anchor differences and numerical errors. No asymptotic likelihood cutoff, CI, identification threshold, data deletion or public fitting change follows from this diagnostic.", "",
  "Any new rule informed by these selected responses requires fresh validation, including estimated covariance parameters and sampling uncertainty. Original inputs, earlier failed numerical evidence and USGS analysis outputs remain unchanged."),
  file.path(out, "README.md"))
invisible(file.copy("scripts/analysis/report_spatial_scale_profiles.R",
  file.path(out, "source_snapshot"), overwrite = FALSE))
cat("Selected profiles, full evaluations and independent numerical audits reported.\n")
