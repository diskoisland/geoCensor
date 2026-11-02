#!/usr/bin/env Rscript
# Post hoc numerical audit; does not replace primary validation outcomes.
source("scripts/analysis/spatial_censored_covariance_reference.R")
pilot <- "output/spatial_covariance_fresh_2026-11-01"
out <- "output/spatial_covariance_precision_2026-11-01"
stopifnot(file.exists(file.path(pilot, "COMPLETED.txt")))
original <- read.csv(file.path(pilot, "covariance_results.csv"))
keys <- original$key[original$fit_status == "finite_interior" &
  !is.na(original$numerical_review) & original$numerical_review]
config <- list(keys = keys, draws = 32768L, maxit = 2400L, reltol = 1e-8,
  audit_draws = 32768L, audit_batches = 8L)
sources <- c("scripts/analysis/refine_spatial_covariance_precision.R",
  "scripts/analysis/spatial_censored_covariance_reference.R",
  "scripts/analysis/spatial_censored_likelihood_reference.R")
inputs <- c(file.path(out, "PROTOCOL.md"), file.path(pilot, c("study_manifest.rds", "covariance_results.csv")),
  file.path(pilot, "cases", paste0(keys, ".rds")),
  file.path(pilot, "fields", paste0(sub("_(moderate|heavy)$", "", keys), ".rds")))
signature <- list(config = config, sources = tools::md5sum(sources), inputs = tools::md5sum(inputs),
  R = R.version.string, RNG = RNGkind(), mvtnorm = as.character(packageVersion("mvtnorm")))
manifest <- file.path(out, "refinement_manifest.rds")
if (file.exists(manifest)) {
  if (!identical(readRDS(manifest), signature)) stop("Refinement provenance changed")
} else {
  if (length(setdiff(list.files(out), "PROTOCOL.md"))) stop("Unverified refinement directory")
  dir.create(file.path(out, "source_snapshot"))
  invisible(file.copy(sources, file.path(out, "source_snapshot")))
  saveRDS(signature, manifest)
  write.csv(data.frame(file = c(sources, inputs), md5 = unname(c(signature$sources, signature$inputs))),
    file.path(out, "source_input_manifest.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
}
if (file.exists(file.path(out, "COMPLETED.txt"))) {
  cat("Completed numerical refinement matches provenance; no results changed.\n"); quit(status = 0L)
}
dir.create(file.path(out, "cases"), showWarnings = FALSE)
rows <- list(); audits <- list()
for (i in seq_along(keys)) {
  path <- file.path(out, "cases", paste0(keys[i], ".rds"))
  if (file.exists(path)) {
    saved <- readRDS(path); rows[[i]] <- saved$row; audits[[i]] <- saved$audit_rows; next
  }
  started <- proc.time()[[3L]]
  a <- readRDS(file.path(pilot, "cases", paste0(keys[i], ".rds")))
  field <- readRDS(file.path(pilot, "fields", paste0(sub("_(moderate|heavy)$", "", keys[i]), ".rds")))
  unit <- spatial_censored_setup_reference(a$lower, a$upper, diag(length(a$lower)))
  u <- spatial_censored_uniforms_reference(unit, config$draws, 770010000L + i)
  data_start <- c(median(a$lower[a$lower == a$upper]), log(.3), log(.03), qlogis(.05))
  fit <- spatial_censored_fit_covariance_reference(a$lower, a$upper, field$coordinates, u,
    starts = rbind(a$fine$candidate$par, data_start), maxit = config$maxit, reltol = config$reltol)
  raw <- fit$candidate
  audit <- NULL; audit_rows <- NULL; old_ll <- NA_real_
  if (!is.null(raw)) {
    old <- a$fine$candidate
    old_ll <- spatial_covariance_loglik_reference(a$lower, a$upper, field$coordinates,
      old$theta, old$sigma, old$range, old$nugget_fraction, u)
    offsets <- rbind(c(0, 0, 0, 0), c(-.02, 0, 0, 0), c(.02, 0, 0, 0),
      c(0, -.02, 0, 0), c(0, .02, 0, 0), c(0, 0, -.05, 0), c(0, 0, .05, 0),
      c(0, 0, 0, -.05), c(0, 0, 0, .05))
    audit <- matrix(NA_real_, config$audit_batches, nrow(offsets))
    for (b in seq_len(config$audit_batches)) {
      uniform <- spatial_censored_uniforms_reference(unit, config$audit_draws, 780010000L + i * 100L + b)
      for (j in seq_len(nrow(offsets))) {
        par <- raw$par + offsets[j, ]
        audit[b, j] <- spatial_covariance_loglik_reference(a$lower, a$upper, field$coordinates,
          par[1L], exp(par[2L]), exp(par[3L]), plogis(par[4L]), uniform)
      }
    }
    delta <- sweep(audit[, -1L, drop = FALSE], 1L, audit[, 1L])
    audit_rows <- data.frame(key = keys[i], perturbation = c("mu_minus_0.02", "mu_plus_0.02",
      "log_sigma_minus_0.02", "log_sigma_plus_0.02", "log_range_minus_0.05", "log_range_plus_0.05",
      "logit_nugget_minus_0.05", "logit_nugget_plus_0.05"), mean_delta_loglik = colMeans(delta),
      delta_mc_se = apply(delta, 2L, sd) / sqrt(nrow(delta)),
      center_loglik_mc_se = sd(audit[, 1L]) / sqrt(nrow(delta)))
  }
  value <- function(name) if (is.null(raw)) NA_real_ else raw[[name]]
  row <- data.frame(key = keys[i], primary_status = a$row$fit_status, refined_status = fit$status,
    primary_candidate_theta = a$row$candidate_theta, refined_theta = fit$theta,
    refined_sigma = fit$sigma, refined_range = fit$range, refined_nugget_fraction = fit$nugget_fraction,
    raw_theta = value("theta"), raw_sigma = value("sigma"), raw_range = value("range"),
    raw_nugget_fraction = value("nugget_fraction"),
    raw_location_shift = value("theta") - a$row$candidate_theta,
    raw_loglik_gain_on_new_uniforms = value("loglik") - old_ll,
    elapsed_seconds = proc.time()[[3L]] - started)
  saveRDS(list(row = row, fit = fit, audit = audit, audit_rows = audit_rows), path)
  rows[[i]] <- row; audits[[i]] <- audit_rows
  cat(keys[i], fit$status, "raw location", value("theta"), "raw nugget", value("nugget_fraction"),
    "seconds", round(row$elapsed_seconds, 1), "\n"); flush.console()
}
if (length(rows)) write.csv(do.call(rbind, rows), file.path(out, "refinement_results.csv"), row.names = FALSE)
if (length(Filter(Negate(is.null), audits))) write.csv(do.call(rbind, audits),
  file.path(out, "independent_refinement_audit.csv"), row.names = FALSE)
stopifnot(identical(tools::md5sum(sources), signature$sources),
  identical(tools::md5sum(inputs), signature$inputs))
writeLines(paste(length(keys), "flagged-interior numerical follow-ups complete; primary pilot unchanged."),
  file.path(out, "COMPLETED.txt"))
cat("Selected numerical follow-up complete.\n")
