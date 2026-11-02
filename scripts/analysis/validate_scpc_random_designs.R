#!/usr/bin/env Rscript
# Random exogenous survey layouts. No response selects the spatial calibration.
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
source("scripts/analysis/scpc_weighted_reference.R")
source("scripts/analysis/conley_survreg_reference.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Supply an optional output directory")
out <- if (length(args)) args[1L] else "output/scpc_random_designs_2026-11-01"
config <- list(n = 24L, cell_size = 0.25, origins_per_axis = 10L, c0 = 10,
  development_layouts_per_regime = 20L, confirmation_layouts_per_regime = 80L,
  gaussian_draws_per_layout_case = 1000L, censored_fields_per_layout_case = 25L,
  regimes = c("uniform", "clustered"), gaussian_decay = c(5, 10, 20, Inf),
  censored_decay = c(10, 20, Inf), left_offsets = c(moderate = -0.3, heavy = 0.6),
  right_offset = 1, truth = log(50), latent_sd = 0.55, noise_sd = 0.2)
sources <- c("scripts/analysis/validate_scpc_random_designs.R",
  "scripts/analysis/scpc_gaussian_reference.R", "scripts/analysis/scpc_envelope_reference.R",
  "scripts/analysis/scpc_weighted_reference.R", "scripts/analysis/conley_survreg_reference.R",
  sort(list.files("R", pattern = "\\.[Rr]$", full.names = TRUE)))
signature <- list(format = 1L, config = config, sources = tools::md5sum(sources),
  R = R.version.string, platform = R.version$platform, survival = as.character(packageVersion("survival")))
manifest <- file.path(out, "study_manifest.rds")
if (file.exists(manifest)) {
  if (!identical(readRDS(manifest), signature)) stop("Source/runtime/configuration differs; use a fresh output directory")
} else {
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Unverified results exist; use a fresh output directory")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "source_snapshot"))
  invisible(file.copy(sources, file.path(out, "source_snapshot")))
  saveRDS(signature, manifest)
  write.csv(data.frame(file = sources, md5 = unname(signature$sources)),
            file.path(out, "source_manifest.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
}
if (file.exists(file.path(out, "COMPLETED.txt"))) {
  cat("Completed study has matching source/runtime/configuration; no results changed.\n")
  quit(status = 0L)
}
dir.create(file.path(out, "layouts"), showWarnings = FALSE)
dir.create(file.path(out, "checkpoints"), showWarnings = FALSE)
writeLines(c("Random exogenous layout study, frozen before response simulation.",
  "20 development and 80 confirmation layouts per regime; 24 sites; uniform or 18 corner+6 whole-window locations.",
  "Sampling locations are independent of the chemistry process. All locations in the unit square have positive sampling support.",
  "Cell width 0.25, 100 origins, four records/site. Weights, components, and calibration recalculated for every layout.",
  "Physical benchmark decay c0=10 and actual decays are fixed across layouts; they are not selected to preserve an average correlation.",
  "31/61-node calibration refinement and 241-midpoint/local-peak audit; all q=1,...,23 are considered.",
  "Gaussian: 1000 fields per layout at actual decays 5(outside),10,20,iid.",
  "Censored: 25 independent fields per layout at decays10,20,iid, each paired across moderate/heavy censoring.",
  "Left offsets -0.3/+0.6 and upper offset+1 from log(50); latent SD0.55; record noise SD0.2.",
  "No exact observations at these common censoring bounds is recorded as unavailable (no finite two-parameter MLE).",
  "Other numerical/fit/calibration failures are retained. Unavailable intervals count as noncoverage.",
  "Coverage uncertainty is estimated across independent layouts; fields at one layout are not counted as independent surveys.",
  "Source/runtime/configuration manifest guards checkpoint reuse. No confirmation tuning or production interval claim."),
  file.path(out, "design_protocol.txt"))

for (phase in c("development", "confirmation")) {
  layouts <- if (phase == "development") config$development_layouts_per_regime else config$confirmation_layouts_per_regime
  for (r in seq_along(config$regimes)) for (i in seq_len(layouts)) {
    regime <- config$regimes[r]
    key <- sprintf("%s_%s_%03d", phase, regime, i)
    path <- file.path(out, "checkpoints", paste0(key, ".rds"))
    if (file.exists(path)) {
      stopifnot(identical(readRDS(path)$key, key)); next
    }
    layout_path <- file.path(out, "layouts", paste0(key, ".rds"))
    if (file.exists(layout_path)) {
      layout <- readRDS(layout_path); stopifnot(identical(layout$key, key))
    } else {
      seed <- if (phase == "development") 210010000L else 220010000L
      set.seed(seed + r * 1000000L + i)
      s <- if (regime == "uniform") matrix(runif(48L), ncol = 2L) else
        rbind(matrix(runif(36L, 0, 0.25), ncol = 2L), matrix(runif(12L), ncol = 2L))
      a <- cell_declustering_weights(data.frame(x = s[, 1L], y = s[, 2L]), "x", "y", config$cell_size, config$origins_per_axis)
      calibration <- tryCatch({
        coarse <- scpc_weighted_select_reference(s, a, config$c0, 31L)
        fine <- scpc_weighted_select_reference(s, a, config$c0, 61L)
        audit <- scpc_weighted_audit_reference(fine$selected, 241L)
        relative <- max(abs(coarse$candidates$critical / fine$candidates$critical - 1))
        if (relative >= 5e-5 || coarse$selected$basis$q != fine$selected$basis$q ||
            audit$max_rejection > 0.05 + 1e-7 || audit$iid_tail_bound > 0.05) stop("Calibration refinement/audit failed")
        list(status = "available", message = "", fit = fine$selected, candidates = fine$candidates,
             relative = relative, audit = audit)
      }, error = function(e) list(status = "unavailable", message = conditionMessage(e), fit = NULL))
      layout <- list(key = key, phase = phase, regime = regime, index = i, s = s, weights = a,
                     kish_n = 1 / sum(a^2), calibration = calibration)
      # Geometry, weights, and calibration are saved before generating responses.
      saveRDS(layout, layout_path)
    }
    s <- layout$s; a <- layout$weights; n <- config$n
    B <- diag(n * a, n) %*% (diag(n) - outer(rep(1, n), a))
    available <- layout$calibration$status == "available"
    fit <- layout$calibration$fit
    Vfull <- qr.Q(qr(cbind(rep(1 / sqrt(n), n), diag(n)[, -1L])))[, -1L, drop = FALSE]
    siteW <- cbind(sqrt(n) * a, t(B) %*% Vfull)
    distance <- as.matrix(dist(s)); gaussian <- list()
    for (j in seq_along(config$gaussian_decay)) {
      decay <- config$gaussian_decay[j]
      S <- if (is.infinite(decay)) diag(n) else exp(-decay * distance)
      seed <- if (phase == "development") 230010000L else 240010000L
      set.seed(seed + r * 1000000L + i * 1000L + j)
      y <- matrix(rnorm(config$gaussian_draws_per_layout_case * n), ncol = n) %*% chol(S)
      theta <- as.vector(y %*% a); psi <- (y - theta) %*% diag(n * a, n)
      site_se <- sqrt(rowSums(psi^2) / (n * (n - 1L))); site_cv <- qt(0.975, n - 1L)
      site_probability <- 1 - scpc_gaussian_rejection_reference(crossprod(siteW, S %*% siteW), site_cv)
      empirical_site <- mean(abs(theta) <= site_cv * site_se)
      p <- NA_real_; empirical <- 0; width <- NA_real_
      if (available) {
        b <- fit$basis; cv <- fit$calibration$critical
        se <- sqrt(rowSums((psi %*% b$V)^2) / (b$q * n))
        p <- 1 - scpc_gaussian_rejection_reference(crossprod(b$W, S %*% b$W), cv)
        empirical <- mean(abs(theta) <= cv * se); width <- mean(2 * cv * se)
        if (abs(empirical - p) > 6 * sqrt(p * (1 - p) / nrow(y)) + 0.001) stop("Gaussian simulation disagrees with analytic probability")
      }
      gaussian[[j]] <- data.frame(key = key, phase = phase, regime = regime, layout = i,
        decay = decay, inside_envelope = decay >= config$c0, draws = nrow(y),
        calibration_available = available, weight_aware_coverage = empirical, analytic_weight_aware_coverage = p,
        site_cluster_coverage = empirical_site, analytic_site_cluster_coverage = site_probability,
        weight_aware_width = width, site_cluster_width = mean(2 * site_cv * site_se))
    }
    records <- list(); k <- 0L; id <- rep(seq_len(n), each = 4L); w <- a[id] / 4
    for (j in seq_along(config$censored_decay)) {
      decay <- config$censored_decay[j]
      S <- if (is.infinite(decay)) diag(n) else exp(-decay * distance)
      C <- chol(S)
      for (h in seq_len(config$censored_fields_per_layout_case)) {
        seed <- if (phase == "development") 250010000L else 260010000L
        set.seed(seed + r * 1000000L + i * 1000L + j * 100L + h)
        latent <- as.vector(t(C) %*% rnorm(n)) * config$latent_sd
        logvalue <- config$truth + latent[id] + rnorm(length(id), sd = config$noise_sd)
        for (censoring in names(config$left_offsets)) {
          left_offset <- config$left_offsets[[censoring]]
          left <- logvalue < config$truth + left_offset
          right <- logvalue > config$truth + config$right_offset
          n_exact <- sum(!left & !right)
          result <- tryCatch({
            if (n_exact == 0L) stop("no_exact_observations_with_common_bounds")
            x <- data.frame(detected_val = ifelse(left | right, NA_real_, exp(logvalue)),
              left_cen = left, left_limit = ifelse(left, exp(config$truth + left_offset), NA_real_),
              right_cen = right, right_limit = ifelse(right, exp(config$truth + config$right_offset), NA_real_))
            model <- fit_censored_score_reference(x, w)
            influence <- as.vector(n * rowsum(model$scores, id, reorder = FALSE) %*% model$bread[1L, ])
            if (abs(mean(influence)) > 1e-5) stop("Influence equations not solved")
            site_se <- sqrt(sum(influence^2) / (n * (n - 1L)))
            if (!is.finite(site_se) || site_se <= 0) stop("Invalid site standard error")
            se <- if (available) sqrt(sum(crossprod(fit$basis$V, influence)^2) / (fit$basis$q * n)) else NA_real_
            if (available && (!is.finite(se) || se <= 0)) stop("Invalid spatial standard error")
            list(theta = model$theta[1L], se = se, site_se = site_se, status = "available", message = "")
          }, error = function(e) list(theta = NA_real_, se = NA_real_, site_se = NA_real_, status = "unavailable", message = conditionMessage(e)))
          ok <- result$status == "available"; spatial_ok <- ok && available
          spatial_cv <- if (available) fit$calibration$critical else NA_real_
          k <- k + 1L
          records[[k]] <- data.frame(key = key, phase = phase, regime = regime, layout = i, decay = decay,
            field = h, censoring = censoring, n_exact = n_exact, n_left = sum(left), n_right = sum(right),
            weighted_exact_fraction = sum(w[!left & !right]), theta = result$theta,
            fit_status = result$status, message = result$message,
            weight_aware_available = spatial_ok, site_cluster_available = ok,
            weight_aware_covered = if (spatial_ok) abs(result$theta - config$truth) <= spatial_cv * result$se else FALSE,
            site_cluster_covered = if (ok) abs(result$theta - config$truth) <= site_cv * result$site_se else FALSE,
            weight_aware_width = if (spatial_ok) 2 * spatial_cv * result$se else NA_real_,
            site_cluster_width = if (ok) 2 * site_cv * result$site_se else NA_real_)
        }
      }
    }
    saveRDS(list(key = key, gaussian = do.call(rbind, gaussian), censored = do.call(rbind, records)), path)
    if (i %% 10L == 0L) {
      cat(phase, regime, i, "/", layouts, "layouts complete\n"); flush.console()
    }
  }
}
paths <- list.files(file.path(out, "checkpoints"), pattern = "\\.rds$", full.names = TRUE)
all <- lapply(paths, readRDS)
gaussian <- do.call(rbind, lapply(all, `[[`, "gaussian"))
censored <- do.call(rbind, lapply(all, `[[`, "censored"))
layouts <- lapply(list.files(file.path(out, "layouts"), pattern = "\\.rds$", full.names = TRUE), readRDS)
layout_summary <- do.call(rbind, lapply(layouts, function(x) {
  cal <- x$calibration
  data.frame(key = x$key, phase = x$phase, regime = x$regime, layout = x$index,
    kish_weight_n = x$kish_n, min_weight = min(x$weights), max_weight = max(x$weights),
    calibration_status = cal$status, message = cal$message,
    q = if (cal$status == "available") cal$fit$basis$q else NA_integer_,
    critical = if (cal$status == "available") cal$fit$calibration$critical else NA_real_,
    audit_max_rejection = if (cal$status == "available") cal$audit$max_rejection else NA_real_)
}))
stopifnot(length(layouts) == 200L, length(all) == 200L, nrow(gaussian) == 800L, nrow(censored) == 30000L)
stopifnot(identical(tools::md5sum(sources), signature$sources))
write.csv(layout_summary, file.path(out, "layout_summary.csv"), row.names = FALSE)
write.csv(gaussian, file.path(out, "gaussian_results.csv"), row.names = FALSE)
write.csv(censored, file.path(out, "censored_results.csv"), row.names = FALSE)
writeLines(c("Random exogenous survey-layout study complete.",
  "200 layouts, 800,000 Gaussian vectors, 30,000 attempted censored fits.",
  "Moderate/heavy censoring paired on common fields; coverage uncertainty must be grouped by layout.",
  "All calibration/fit failures retained. No preferential-sampling or production interval claim."),
  file.path(out, "COMPLETED.txt"))
