#!/usr/bin/env Rscript
# Equal record counts per site; cell weights fixed conditional on fresh exogenous
# locations. Weighted Gaussian calibration is frozen before censored responses.
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
source("scripts/analysis/conley_survreg_reference.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop("Supply optional Gaussian input and fresh output directories")
input <- if (length(args)) args[1L] else "output/scpc_weighted_gaussian_2026-11-01"
out <- if (length(args) == 2L) args[2L] else "output/scpc_weighted_censored_2026-11-01"
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Use a fresh output directory")
stopifnot(file.exists(file.path(input, "COMPLETED.txt")))
frozen <- readRDS(file.path(input, "frozen_calibrations.rds"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
sources <- c("scripts/analysis/validate_scpc_weighted_censored.R", "scripts/analysis/conley_survreg_reference.R",
             sort(list.files("R", pattern = "\\.[Rr]$", full.names = TRUE)),
             file.path(input, "frozen_calibrations.rds"))
dir.create(file.path(out, "source_snapshot"))
invisible(file.copy(sources, file.path(out, "source_snapshot")))
write.csv(data.frame(file = sources, md5 = unname(tools::md5sum(sources))),
          file.path(out, "source_manifest.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
truth <- log(50); results <- list(); count <- 0L
for (phase in c("development", "confirmation")) {
  simulations <- if (phase == "development") 200L else 500L
  base_seed <- if (phase == "development") 170010000L else 180010000L
  for (d in seq_along(frozen$designs)) {
    design <- names(frozen$designs)[d]; pair <- frozen$fits[[design]]
    b <- pair$weighted$basis; old <- pair$geometry_only
    n <- b$n; id <- rep(seq_len(n), each = 4L)
    weights <- b$weights[id] / 4
    for (j in seq_along(c(1, 2, Inf))) {
      multiple <- c(1, 2, Inf)[j]
      Sigma <- if (is.infinite(multiple)) diag(n) else exp(-multiple * b$c0 * b$distance)
      C <- chol(Sigma)
      rows <- vector("list", simulations)
      for (i in seq_len(simulations)) {
        set.seed(base_seed + d * 100000L + j * 10000L + i)
        latent <- as.vector(t(C) %*% rnorm(n)) * 0.55
        logvalue <- truth + latent[id] + rnorm(length(id), sd = 0.2)
        left <- logvalue < truth - 0.3; right <- logvalue > truth + 1
        x <- data.frame(detected_val = ifelse(left | right, NA_real_, exp(logvalue)),
          left_cen = left, left_limit = ifelse(left, exp(truth - 0.3), NA_real_),
          right_cen = right, right_limit = ifelse(right, exp(truth + 1), NA_real_))
        answer <- tryCatch({
          fit <- fit_censored_score_reference(x, weights)
          scores <- rowsum(fit$scores, id, reorder = FALSE)
          psi <- as.vector(n * scores %*% fit$bread[1L, ])
          if (abs(mean(psi)) > 1e-5) stop("Influence equations not solved")
          se <- sqrt(sum(crossprod(b$V, psi)^2) / (b$q * n))
          alt <- sqrt(sum(crossprod(old$basis$W[, -1L, drop = FALSE], psi)^2) / (old$basis$q * n))
          site <- sqrt(sum(psi^2) / (n * (n - 1L)))
          if (any(!is.finite(c(se, alt, site))) || min(se, alt, site) <= 0) stop("Invalid standard error")
          list(theta = fit$theta[1L], se = se, alternative_se = alt, site_se = site,
               status = "available", message = "")
        }, error = function(e) list(theta = NA_real_, se = NA_real_, alternative_se = NA_real_, site_se = NA_real_,
                                    status = "unavailable", message = conditionMessage(e)))
        rows[[i]] <- data.frame(phase = phase, design = design, n_sites = n,
          actual_decay_over_c0 = multiple, simulation = i, n_observations = length(id),
          n_left = sum(left), n_right = sum(right), n_exact = sum(!left & !right),
          weighted_exact_fraction = sum(weights[!left & !right]),
          theta = answer$theta, weighted_se = answer$se, geometry_se = answer$alternative_se,
          site_se = answer$site_se, weighted_critical = pair$weighted$calibration$critical,
          geometry_critical = old$calibration$critical, q = b$q,
          status = answer$status, message = answer$message)
      }
      count <- count + 1L; results[[count]] <- do.call(rbind, rows)
    }
    cat("Weighted censored", phase, design, "complete\n"); flush.console()
  }
}
raw <- do.call(rbind, results)
write.csv(raw, file.path(out, "censored_results.csv"), row.names = FALSE)
keys <- interaction(raw$phase, raw$design, raw$actual_decay_over_c0, drop = TRUE)
summary <- lapply(split(raw, keys), function(x) {
  do.call(rbind, lapply(c("weight_aware", "geometry_only", "site_cluster_t"), function(method) {
    se <- switch(method, weight_aware = x$weighted_se, geometry_only = x$geometry_se, x$site_se)
    critical <- switch(method, weight_aware = x$weighted_critical[1L],
                       geometry_only = x$geometry_critical[1L], qt(0.975, x$n_sites[1L] - 1L))
    ok <- x$status == "available" & is.finite(se) & is.finite(x$theta)
    hits <- sum(ok & abs(x$theta - truth) <= critical * se, na.rm = TRUE)
    ci <- binom.test(hits, nrow(x))$conf.int
    data.frame(phase = x$phase[1L], design = x$design[1L], n_sites = x$n_sites[1L],
      actual_decay_over_c0 = x$actual_decay_over_c0[1L], method = method,
      attempted = nrow(x), available = sum(ok), failed = sum(!ok), coverage = hits / nrow(x),
      coverage_available = if (any(ok)) hits / sum(ok) else NA_real_, mc_low = ci[1L], mc_high = ci[2L],
      bias_log = if (any(ok)) mean(x$theta[ok] - truth) else NA_real_,
      mean_log_width = if (any(ok)) mean(2 * critical * se[ok]) else NA_real_,
      mean_weighted_exact_fraction = mean(x$weighted_exact_fraction))
  }))
})
write.csv(do.call(rbind, summary), file.path(out, "censored_summary.csv"), row.names = FALSE)
stopifnot(nrow(raw) == 8400L)
writeLines(c("Weighted censored transfer study complete; 8,400 attempted fits, three interval methods.",
  "Conditional on four fixed exogenous geometries and coordinate-derived cell weights.",
  "Fresh independent development and confirmation seeds; no calibration retuning.",
  "No variable-site selection, chemistry-dependent sampling, or real-data coverage guarantee."),
  file.path(out, "COMPLETED.txt"))
