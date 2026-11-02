#!/usr/bin/env Rscript
# First transfer test only: equal weights, fixed observed locations, no omitted
# sites or preferential sampling. Gaussian calibration is frozen independently.
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
source("scripts/analysis/conley_survreg_reference.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop("Supply optional Gaussian input and fresh output directories")
input <- if (length(args)) args[1L] else "output/scpc_envelope_validation_2026-11-01"
out <- if (length(args) == 2L) args[2L] else "output/scpc_censored_transfer_site_control_2026-11-01"
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Use a fresh output directory")
stopifnot(file.exists(file.path(input, "COMPLETED.txt")))
designs <- readRDS(file.path(input, "designs.rds"))
fits <- readRDS(file.path(input, "frozen_calibrations.rds"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
sources <- c("scripts/analysis/validate_scpc_censored_transfer.R",
             "scripts/analysis/conley_survreg_reference.R",
             sort(list.files("R", pattern = "\\.[Rr]$", full.names = TRUE)),
             file.path(input, c("designs.rds", "frozen_calibrations.rds")))
write.csv(data.frame(file = sources, md5 = unname(tools::md5sum(sources))),
          file.path(out, "source_manifest.csv"), row.names = FALSE)
snapshot <- file.path(out, "source_snapshot"); dir.create(snapshot)
for (p in sources) file.copy(p, file.path(snapshot, basename(p)))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
truth <- log(50)
all <- list(); total <- 0L
for (phase in c("development", "confirmation")) {
  simulations <- if (phase == "development") 200L else 500L
  base_seed <- if (phase == "development") 11024000L else 12024000L
  for (d in seq_along(designs)) {
    design <- names(designs)[d]; b <- fits[[design]]$basis
    n <- b$n; id <- rep(seq_len(n), each = 4L)
    for (j in seq_along(c(1, 2, Inf))) {
      multiple <- c(1, 2, Inf)[j]
      Sigma <- if (is.infinite(multiple)) diag(n) else exp(-multiple * b$c0 * b$distance)
      C <- chol(Sigma)
      rows <- vector("list", simulations)
      for (i in seq_len(simulations)) {
        set.seed(base_seed + d * 100000L + j * 10000L + i)
        latent <- as.vector(t(C) %*% rnorm(n)) * 0.55
        logvalue <- truth + latent[id] + rnorm(length(id), sd = 0.2)
        value <- exp(logvalue)
        left <- logvalue < truth - 0.3
        right <- logvalue > truth + 1
        x <- data.frame(detected_val = ifelse(left | right, NA_real_, value),
                        left_cen = left, left_limit = ifelse(left, exp(truth - 0.3), NA_real_),
                        right_cen = right, right_limit = ifelse(right, exp(truth + 1), NA_real_))
        result <- tryCatch({
          model <- fit_censored_score_reference(x)
          site_score <- rowsum(model$scores, id, reorder = FALSE)
          influence <- as.vector(n * site_score %*% model$bread[1L, ])
          if (abs(mean(influence)) > 1e-5) stop("Influence mean not sufficiently zero")
          projection <- crossprod(b$W[, -1L, drop = FALSE], influence)
          se <- sqrt(sum(projection^2) / (b$q * n))
          site_se <- sqrt(sum(influence^2) / (n * (n - 1L)))
          working_se <- sqrt(model$fit$var[1L, 1L])
          if (!is.finite(se) || se <= 0 || !is.finite(site_se) || site_se <= 0 ||
              !is.finite(working_se) || working_se <= 0) stop("Invalid standard error")
          list(theta = model$theta[1L], scpc_se = se, site_se = site_se, working_se = working_se,
               status = "available", message = "")
        }, error = function(e) list(theta = NA_real_, scpc_se = NA_real_, site_se = NA_real_, working_se = NA_real_,
                                    status = "unavailable", message = conditionMessage(e)))
        rows[[i]] <- data.frame(phase = phase, design = design, n_sites = n,
          actual_decay_over_c0 = multiple, simulation = i,
          n_observations = length(id), n_left = sum(left), n_right = sum(right),
          n_exact = sum(!left & !right), theta = result$theta,
          scpc_se = result$scpc_se, site_se = result$site_se, working_se = result$working_se,
          scpc_critical = fits[[design]]$calibration$critical,
          q = b$q, status = result$status, message = result$message)
      }
      total <- total + 1L; all[[total]] <- do.call(rbind, rows)
    }
    cat("Censored", phase, design, "complete\n"); flush.console()
  }
}
raw <- do.call(rbind, all)
write.csv(raw, file.path(out, "censored_results.csv"), row.names = FALSE)
key <- interaction(raw$phase, raw$design, raw$actual_decay_over_c0, drop = TRUE)
summary <- lapply(split(raw, key), function(x) {
  do.call(rbind, lapply(c("scpc_equal_weight_transfer", "site_cluster_t", "working_independence"), function(method) {
    se <- switch(method, working_independence = x$working_se, site_cluster_t = x$site_se, x$scpc_se)
    critical <- switch(method, working_independence = qnorm(0.975),
                       site_cluster_t = qt(0.975, x$n_sites[1L] - 1L), x$scpc_critical[1L])
    ok <- x$status == "available" & is.finite(se) & is.finite(x$theta)
    covered <- ok & abs(x$theta - truth) <= critical * se
    hits <- sum(covered, na.rm = TRUE)
    ci <- binom.test(hits, nrow(x))$conf.int
    data.frame(phase = x$phase[1L], design = x$design[1L], n_sites = x$n_sites[1L],
      actual_decay_over_c0 = x$actual_decay_over_c0[1L], method = method,
      attempted = nrow(x), available = sum(ok), failed = sum(!ok),
      coverage_all_attempts = hits / nrow(x), coverage_available = if (any(ok)) hits / sum(ok) else NA_real_,
      mc_low = ci[1L], mc_high = ci[2L], bias_log = if (any(ok)) mean(x$theta[ok] - truth) else NA_real_,
      mean_log_interval_width = if (any(ok)) mean(2 * critical * se[ok]) else NA_real_,
      mean_left_fraction = mean(x$n_left / x$n_observations),
      mean_right_fraction = mean(x$n_right / x$n_observations))
  }))
})
write.csv(do.call(rbind, summary), file.path(out, "censored_summary.csv"), row.names = FALSE)
stopifnot(nrow(raw) == 8400L)
writeLines(c("Equal-weight censored transfer study complete; 8,400 attempted fits; three interval methods.",
  "Critical values, bases, and component counts frozen from the independent Gaussian study.",
  "No retuning after development or confirmation. No estimated declustering weights or omitted sites.",
  "Site-cluster comparator added after initial two-method results; original frozen SCPC calculations are unchanged.",
  "No claim of a calibrated production interval for NGDB."), file.path(out, "COMPLETED.txt"))
