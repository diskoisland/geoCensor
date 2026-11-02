#!/usr/bin/env Rscript
source("scripts/analysis/spatial_censored_profile_reference.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Supply an optional output directory")
out <- if (length(args)) args[1L] else "output/spatial_scale_profiles_2026-11-01"
pilot <- "output/spatial_scale_fresh_2026-11-01"
keys <- c("uniform_06_decay10_heavy", "uniform_06_decay10_moderate",
  "clustered_07_decay20_heavy", "clustered_07_decay20_moderate")
config <- list(keys = keys, main_draws = 32768L, refined_draws = 65536L,
  audit_draws = 32768L, audit_batches = 8L,
  location_grid = seq(-4, 6, length.out = 61L),
  offsets = c(.02, .05, .1, .2, .5, 1, 2),
  log_multiplier_grid = seq(-4, 4, length.out = 17L), truth = log(50))
sources <- c("scripts/analysis/diagnose_spatial_scale_profiles.R",
  "scripts/analysis/spatial_censored_profile_reference.R",
  "scripts/analysis/verify_spatial_censored_profile_reference.R",
  "scripts/analysis/spatial_censored_scale_reference.R",
  "scripts/analysis/spatial_censored_likelihood_reference.R")
inputs <- c(file.path(out, c("PROTOCOL.md", "verification.csv", "VERIFIED.txt")),
  file.path(pilot, c("study_manifest.rds", "fresh_results.csv")),
  file.path(pilot, "cases", paste0(keys, ".rds")),
  file.path(pilot, "fields", paste0(unique(sub("_(moderate|heavy)$", "", keys)), ".rds")))
stopifnot(all(file.exists(inputs)))
signature <- list(config = config, sources = tools::md5sum(sources),
  inputs = tools::md5sum(inputs), R = R.version.string, platform = R.version$platform,
  RNG = RNGkind(), mvtnorm = as.character(packageVersion("mvtnorm")))
manifest <- file.path(out, "diagnostic_manifest.rds")
if (file.exists(manifest)) {
  if (!identical(readRDS(manifest), signature)) stop("Diagnostic source/input/runtime changed")
} else {
  if (length(setdiff(list.files(out), c("PROTOCOL.md", "verification.csv", "VERIFIED.txt"))))
    stop("Unverified diagnostic outputs exist")
  dir.create(file.path(out, "source_snapshot"))
  invisible(file.copy(sources, file.path(out, "source_snapshot")))
  saveRDS(signature, manifest)
  write.csv(data.frame(file = c(sources, inputs),
    md5 = unname(c(signature$sources, signature$inputs))),
    file.path(out, "source_input_manifest.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
}
if (file.exists(file.path(out, "COMPLETED.txt"))) {
  cat("Completed selected profiles match provenance; no results changed.\n")
  quit(status = 0L)
}
dir.create(file.path(out, "cases"), showWarnings = FALSE)
pilot_results <- read.csv(file.path(pilot, "fresh_results.csv"))
selected <- unlist(lapply(c("uniform", "clustered"), function(regime) {
  x <- pilot_results[pilot_results$regime == regime &
    pilot_results$censoring == "heavy" & pilot_results$public_count_eligible, ]
  key <- x$key[which.max(abs(x$spatial_theta - config$truth))]
  c(key, sub("_heavy$", "_moderate", key))
}))
stopifnot(identical(keys, selected))

for (i in seq_along(keys)) {
  path <- file.path(out, "cases", paste0(keys[i], ".rds"))
  if (file.exists(path)) next
  started <- proc.time()[[3L]]
  original <- readRDS(file.path(pilot, "cases", paste0(keys[i], ".rds")))
  field <- readRDS(file.path(pilot, "fields", paste0(sub("_(moderate|heavy)$", "", keys[i]), ".rds")))
  model <- spatial_censored_setup_reference(original$lower, original$upper, field$covariance)
  stopifnot(model$no >= 2L)
  main_u <- spatial_censored_uniforms_reference(model, config$main_draws, 610010000L + i)
  reference <- spatial_censored_fit_scale_reference(model, main_u, starts = rbind(
    c(original$fine$theta / original$fine$multiplier, 1 / original$fine$multiplier),
    c(original$oracle$theta, 1)))
  if (reference$status != "finite") stop("New numerical reference fit unavailable: ", keys[i])
  mu_hat <- reference$theta
  locations <- sort(unique(c(config$location_grid, mu_hat, config$truth,
    original$fine$theta, original$oracle$theta,
    original$upper[is.infinite(original$lower)],
    original$lower[is.infinite(original$upper)],
    mu_hat + config$offsets, mu_hat - config$offsets)))
  profiles <- vector("list", length(locations)); rows <- vector("list", length(locations))
  fixed_evaluation <- function(mu, u) tryCatch(list(status = "finite",
    loglik = spatial_censored_loglik_reference(model, mu, u), message = ""),
    error = function(e) list(status = "unavailable", loglik = NA_real_, message = conditionMessage(e)))
  for (j in seq_along(locations)) {
    profile <- spatial_censored_profile_scale_reference(model, locations[j], main_u,
      log_multiplier_grid = config$log_multiplier_grid)
    fixed <- fixed_evaluation(locations[j], main_u)
    profiles[[j]] <- profile
    rows[[j]] <- data.frame(key = keys[i], regime = original$row$regime,
      censoring = original$row$censoring, mu = locations[j], status = profile$status,
      multiplier = profile$multiplier,
      marginal_sigma = sqrt(field$covariance[1L, 1L]) * profile$multiplier,
      loglik = profile$loglik, loglik_minus_reference = profile$loglik - reference$loglik,
      fixed_scale_status = fixed$status, fixed_scale_loglik = fixed$loglik,
      fixed_scale_minus_reference = fixed$loglik - reference$loglik,
      fixed_scale_message = fixed$message,
      invalid_evaluations = sum(!is.finite(profile$evaluations$loglik)))
    if (j %% 20L == 0L) {
      cat(keys[i], j, "/", length(locations), "profile points\n"); flush.console()
    }
  }
  grid <- do.call(rbind, rows)
  center <- match(mu_hat, grid$mu)
  # A constrained fixed-scale point must not exceed the profiled value; all use
  # common integration uniforms. Do not silently clamp numerical discrepancies.
  finite <- grid$status == "finite" & grid$fixed_scale_status == "finite"
  stopifnot(abs(grid$loglik_minus_reference[center]) < 1e-4,
    max(grid$loglik_minus_reference[grid$status == "finite"]) < 1e-3,
    all(grid$loglik[finite] >= grid$fixed_scale_loglik[finite] - 1e-8))
  anchors <- c(mu_hat, config$truth, mu_hat - .5, mu_hat + .5, mu_hat - 1, mu_hat + 1)
  anchor_labels <- c("reference_location", "generating_location", "mu_minus_0.5",
    "mu_plus_0.5", "mu_minus_1", "mu_plus_1")
  refined_u <- spatial_censored_uniforms_reference(model, config$refined_draws, 620010000L + i)
  refined <- lapply(anchors, function(mu) spatial_censored_profile_scale_reference(model,
    mu, refined_u, log_multiplier_grid = config$log_multiplier_grid))
  if (!all(vapply(refined, function(x) x$status == "finite", TRUE)))
    stop("Anchor profile unavailable or search-bound: ", keys[i])
  audit <- matrix(NA_real_, config$audit_batches, length(anchors) * 2L)
  for (b in seq_len(config$audit_batches)) {
    u <- spatial_censored_uniforms_reference(model, config$audit_draws, 630010000L + i * 100L + b)
    for (a in seq_along(anchors)) {
      lambda <- refined[[a]]$multiplier
      audit[b, a] <- spatial_censored_scaled_loglik_reference(model, anchors[a] / lambda, 1 / lambda, u)
      audit[b, length(anchors) + a] <- spatial_censored_loglik_reference(model, anchors[a], u)
    }
  }
  anchor_rows <- do.call(rbind, lapply(seq_along(anchors), function(a) {
    main_index <- match(anchors[a], grid$mu)
    stopifnot(!is.na(main_index))
    main_lambda <- grid$multiplier[main_index]
    main_under_fine <- spatial_censored_scaled_loglik_reference(model,
      anchors[a] / main_lambda, 1 / main_lambda, refined_u)
    data.frame(key = keys[i], anchor = anchor_labels[a], mu = anchors[a],
      main_multiplier = main_lambda, refined_multiplier = refined[[a]]$multiplier,
      log_multiplier_shift = log(refined[[a]]$multiplier / main_lambda),
      refined_loglik_gain = refined[[a]]$loglik - main_under_fine)
  }))
  summary <- data.frame(key = keys[i], regime = original$row$regime,
    censoring = original$row$censoring, n = model$n, n_exact = model$no,
    n_left = sum(is.infinite(original$lower)), n_right = sum(is.infinite(original$upper)),
    weighted_exact_fraction = original$row$weighted_exact_fraction,
    source_mu = original$fine$theta, reference_mu = mu_hat,
    reference_marginal_sigma = sqrt(field$covariance[1L, 1L]) * reference$multiplier,
    reference_refit_shift = mu_hat - original$fine$theta,
    n_profile_points = nrow(grid), n_finite_profiles = sum(grid$status == "finite"),
    n_boundary_profiles = sum(grid$status == "profile_search_boundary"),
    n_other_profiles = sum(!grid$status %in% c("finite", "profile_search_boundary")),
    invalid_profile_probes = sum(grid$invalid_evaluations),
    elapsed_seconds = proc.time()[[3L]] - started)
  saveRDS(list(summary = summary, grid = grid, profiles = profiles, reference = reference,
    anchor_labels = anchor_labels, anchors = anchors, refined = refined,
    anchor_rows = anchor_rows, audit = audit, original_key = keys[i],
    seeds = c(main = 610010000L + i, refined = 620010000L + i,
      setNames(630010000L + i * 100L + seq_len(config$audit_batches),
        paste0("audit", seq_len(config$audit_batches))))), path)
  cat(i, "/ 4", keys[i], "saved; seconds", round(summary$elapsed_seconds, 1), "\n")
  flush.console()
}
cases <- lapply(file.path(out, "cases", paste0(keys, ".rds")), readRDS)
write.csv(do.call(rbind, lapply(cases, `[[`, "summary")), file.path(out, "case_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(cases, `[[`, "grid")), file.path(out, "location_profiles.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(cases, `[[`, "anchor_rows")), file.path(out, "anchor_refinement.csv"), row.names = FALSE)
audits <- do.call(rbind, lapply(seq_along(cases), function(i) {
  a <- cases[[i]]; paired <- sweep(a$audit, 1L, a$audit[, 1L])
  data.frame(key = keys[i], anchor = rep(a$anchor_labels, 2L),
    mode = rep(c("profiled_scale", "fixed_true_scale"), each = length(a$anchor_labels)),
    mu = rep(a$anchors, 2L), mean_loglik_delta = colMeans(paired),
    delta_mc_se = apply(paired, 2L, sd) / sqrt(nrow(paired)),
    absolute_loglik_mc_se = apply(a$audit, 2L, sd) / sqrt(nrow(paired)))
}))
write.csv(audits, file.path(out, "independent_anchor_audit.csv"), row.names = FALSE)
stopifnot(identical(tools::md5sum(sources), signature$sources),
  identical(tools::md5sum(inputs), signature$inputs))
writeLines(c("Selected-case spatial scale profiles complete — November 1, 2026.",
  "Post hoc mechanism diagnostic; no confidence cutoff, production threshold or validation claim.",
  "All grid evaluations, unavailable probes and independent anchor audits retained."),
  file.path(out, "COMPLETED.txt"))
cat("Four selected profiles and anchor audits complete.\n")
