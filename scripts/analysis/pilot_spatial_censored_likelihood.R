#!/usr/bin/env Rscript
# Development/feasibility pilot only; no confidence interval calibration.
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
source("scripts/analysis/spatial_censored_likelihood_reference.R")
source("scripts/analysis/conley_survreg_reference.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Supply an optional output directory")
out <- if (length(args)) args[1L] else "output/spatial_likelihood_reference_2026-11-01"
protocol <- file.path(out, "PILOT_PROTOCOL.md")
if (!file.exists(protocol)) stop("The pilot protocol must exist before generating responses")
config <- list(n_sites = 24L, repeats = 4L, layouts_per_regime = 4L,
  regimes = c("uniform", "clustered"), decay = c(10, 20),
  truth = log(50), latent_sd = .55, noise_sd = .2,
  left_offset = c(moderate = -.3, heavy = .6), right_offset = 1,
  coarse_draws = 16384L, fine_draws = 65536L,
  audit_draws = 32768L, audit_batches = 8L, audit_offsets = c(-.05, 0, .05))
sources <- c("scripts/analysis/pilot_spatial_censored_likelihood.R",
  "scripts/analysis/spatial_censored_likelihood_reference.R",
  "scripts/analysis/verify_spatial_censored_likelihood_reference.R",
  "scripts/analysis/conley_survreg_reference.R",
  sort(list.files("R", pattern = "\\.[Rr]$", full.names = TRUE)))
signature <- list(config = config, sources = tools::md5sum(sources),
  protocol = tools::md5sum(protocol), R = R.version.string,
  platform = R.version$platform, mvtnorm = as.character(packageVersion("mvtnorm")),
  survival = as.character(packageVersion("survival")))
manifest <- file.path(out, "pilot_manifest.rds")
if (file.exists(manifest)) {
  if (!identical(readRDS(manifest), signature)) stop("Pilot provenance changed; use a fresh declared output directory")
} else {
  allowed <- c("PILOT_PROTOCOL.md", "verification.csv", "VERIFIED.txt")
  if (length(setdiff(list.files(out, all.files = TRUE, no.. = TRUE), allowed)))
    stop("Unverified pilot files already exist")
  dir.create(file.path(out, "source_snapshot"), showWarnings = FALSE)
  invisible(file.copy(sources, file.path(out, "source_snapshot")))
  saveRDS(signature, manifest)
  write.csv(data.frame(file = sources, md5 = unname(signature$sources)),
    file.path(out, "source_manifest.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
}
if (file.exists(file.path(out, "PILOT_COMPLETED.txt"))) {
  cat("Completed feasibility pilot matches provenance; no results changed.\n")
  quit(status = 0L)
}
for (directory in c("layouts", "fields", "cases"))
  dir.create(file.path(out, directory), showWarnings = FALSE)

independent_known <- function(lower, upper, sd, weights) {
  exact <- lower == upper
  if (!any(exact) && (all(lower == -Inf) || all(upper == Inf)))
    return(list(status = "no_finite_location_mle", theta = NA_real_))
  weights <- weights / sum(weights)
  loglik <- function(mu) {
    values <- numeric(length(lower))
    values[exact] <- dnorm(lower[exact], mu, sd, log = TRUE)
    if (any(!exact)) {
      lo <- (lower[!exact] - mu) / sd; up <- (upper[!exact] - mu) / sd
      positive <- lo > 0
      x <- ifelse(positive, pnorm(lo, lower.tail = FALSE, log.p = TRUE),
                  pnorm(up, log.p = TRUE))
      y <- ifelse(positive, pnorm(up, lower.tail = FALSE, log.p = TRUE),
                  pnorm(lo, log.p = TRUE))
      values[!exact] <- x + log(-expm1(y - x))
    }
    sum(weights * values)
  }
  bounds <- c(lower[is.finite(lower)], upper[is.finite(upper)])
  interval <- median(bounds) + c(-1, 1) * sd * 16
  fit <- optimize(loglik, interval, maximum = TRUE, tol = 1e-7)
  if (min(fit$maximum - interval[1L], interval[2L] - fit$maximum) < .01)
    stop("Independent comparator bracket exhausted")
  list(status = "finite", theta = fit$maximum)
}

for (r in seq_along(config$regimes)) for (i in seq_len(config$layouts_per_regime)) {
  index <- (r - 1L) * config$layouts_per_regime + i
  key <- sprintf("%s_%02d", config$regimes[r], i)
  location_file <- file.path(out, "layouts", paste0(key, ".rds"))
  if (file.exists(location_file)) layout <- readRDS(location_file) else {
    layout_seed <- 310010000L + index; set.seed(layout_seed)
    s <- if (r == 1L) matrix(runif(48L), ncol = 2L) else
      rbind(matrix(runif(36L, 0, .25), ncol = 2L), matrix(runif(12L), ncol = 2L))
    a <- cell_declustering_weights(data.frame(x = s[, 1L], y = s[, 2L]),
                                    "x", "y", .25, 10L)
    layout <- list(key = key, seed = layout_seed, s = s, weights = a)
    saveRDS(layout, location_file)
  }
  for (j in seq_along(config$decay)) {
    field_index <- (index - 1L) * length(config$decay) + j
    field_key <- paste0(key, "_decay", config$decay[j])
    field_file <- file.path(out, "fields", paste0(field_key, ".rds"))
    if (file.exists(field_file)) field <- readRDS(field_file) else {
      field_seed <- 320010000L + field_index; set.seed(field_seed)
      R <- exp(-config$decay[j] * as.matrix(dist(layout$s)))
      id <- rep(seq_len(config$n_sites), each = config$repeats)
      latent <- as.vector(t(chol(R)) %*% rnorm(config$n_sites)) * config$latent_sd
      z <- config$truth + latent[id] + rnorm(length(id), sd = config$noise_sd)
      S <- config$latent_sd^2 * R[id, id] + diag(config$noise_sd^2, length(id))
      field <- list(key = field_key, seed = field_seed, id = id, log_values = z,
        covariance = S, weights = layout$weights[id] / config$repeats)
      saveRDS(field, field_file)
    }
    for (severity in seq_along(config$left_offset)) {
      case_index <- (field_index - 1L) * length(config$left_offset) + severity
      case_key <- paste0(field_key, "_", names(config$left_offset)[severity])
      case_file <- file.path(out, "cases", paste0(case_key, ".rds"))
      if (file.exists(case_file)) next
      left_limit <- config$truth + config$left_offset[severity]
      right_limit <- config$truth + config$right_offset
      left <- field$log_values < left_limit
      right <- field$log_values > right_limit
      lower <- ifelse(left, -Inf, ifelse(right, right_limit, field$log_values))
      upper <- ifelse(right, Inf, ifelse(left, left_limit, field$log_values))
      model <- spatial_censored_setup_reference(lower, upper, field$covariance)
      started <- proc.time()[[3L]]
      spatial <- tryCatch({
        coarse_uniform <- spatial_censored_uniforms_reference(model, config$coarse_draws,
          330010000L + case_index)
        coarse <- spatial_censored_fit_mu_reference(model, coarse_uniform)
        if (coarse$status != "finite") list(status = coarse$status,
          theta = NA_real_, coarse_theta = NA_real_, audit = NULL,
          message = "Known-covariance location has no finite maximum") else {
          fine_uniform <- spatial_censored_uniforms_reference(model, config$fine_draws,
            340010000L + case_index)
          fine <- spatial_censored_fit_mu_reference(model, fine_uniform)
          if (fine$status != "finite") stop("Refined fit unavailable")
          audit <- matrix(NA_real_, config$audit_batches, length(config$audit_offsets))
          for (b in seq_len(config$audit_batches)) {
            uniforms <- spatial_censored_uniforms_reference(model, config$audit_draws,
              350010000L + case_index * 100L + b)
            audit[b, ] <- vapply(fine$theta + config$audit_offsets, function(mu)
              spatial_censored_loglik_reference(model, mu, uniforms), 0.0)
          }
          list(status = "finite", theta = fine$theta, coarse_theta = coarse$theta,
            coarse_loglik = coarse$loglik, fine_loglik = fine$loglik,
            audit = audit, message = "")
        }
      }, error = function(e) list(status = "numerical_failure", theta = NA_real_,
        coarse_theta = NA_real_, audit = NULL, message = conditionMessage(e)))
      sd <- sqrt(config$latent_sd^2 + config$noise_sd^2)
      comparator <- function(weights) tryCatch(independent_known(lower, upper, sd, weights),
        error = function(e) list(status = conditionMessage(e), theta = NA_real_))
      equal <- comparator(rep(1, model$n)); weighted <- comparator(field$weights)
      working <- tryCatch({
        if (!model$no) stop("unknown_scale_no_exact_observations")
        df <- data.frame(detected_val = ifelse(left | right, NA_real_, exp(field$log_values)),
          left_cen = left, left_limit = ifelse(left, exp(left_limit), NA_real_),
          right_cen = right, right_limit = ifelse(right, exp(right_limit), NA_real_))
        fit <- fit_censored_score_reference(df, field$weights)
        list(status = "finite", theta = fit$theta[1L], log_scale = fit$theta[2L])
      }, error = function(e) list(status = conditionMessage(e), theta = NA_real_))
      audit <- spatial$audit
      center_se <- if (!is.null(audit)) sd(audit[, 2L]) / sqrt(nrow(audit)) else NA_real_
      delta <- if (!is.null(audit)) sweep(audit[, c(1L, 3L), drop = FALSE], 1L, audit[, 2L]) else NULL
      row <- data.frame(key = case_key, regime = config$regimes[r], layout = i,
        decay = config$decay[j], censoring = names(config$left_offset)[severity],
        n_exact = model$no, n_left = sum(left), n_right = sum(right),
        spatial_status = spatial$status, spatial_theta = spatial$theta,
        coarse_spatial_theta = spatial$coarse_theta,
        refit_shift = spatial$theta - spatial$coarse_theta,
        audit_loglik_mc_se = center_se,
        audit_lower_delta = if (!is.null(delta)) mean(delta[, 1L]) else NA_real_,
        audit_lower_delta_mc_se = if (!is.null(delta)) sd(delta[, 1L]) / sqrt(nrow(delta)) else NA_real_,
        audit_upper_delta = if (!is.null(delta)) mean(delta[, 2L]) else NA_real_,
        audit_upper_delta_mc_se = if (!is.null(delta)) sd(delta[, 2L]) / sqrt(nrow(delta)) else NA_real_,
        equal_known_status = equal$status, equal_known_theta = equal$theta,
        weighted_known_status = weighted$status, weighted_known_theta = weighted$theta,
        weighted_unknown_status = working$status, weighted_unknown_theta = working$theta,
        elapsed_seconds = proc.time()[[3L]] - started, message = spatial$message)
      saveRDS(list(row = row, lower = lower, upper = upper, spatial = spatial,
        equal = equal, weighted = weighted, working = working,
        integration_seeds = list(coarse = 330010000L + case_index,
          fine = 340010000L + case_index,
          audit = 350010000L + case_index * 100L + seq_len(config$audit_batches))), case_file)
      cat(case_index, "/ 32", case_key, spatial$status, "elapsed",
          round(row$elapsed_seconds, 1), "seconds\n"); flush.console()
    }
  }
}
paths <- list.files(file.path(out, "cases"), pattern = "\\.rds$", full.names = TRUE)
cases <- lapply(paths, readRDS)
stopifnot(length(cases) == 32L,
  identical(tools::md5sum(sources), signature$sources),
  identical(tools::md5sum(protocol), signature$protocol))
results <- do.call(rbind, lapply(cases, `[[`, "row"))
write.csv(results, file.path(out, "pilot_results.csv"), row.names = FALSE)
writeLines(c("Feasibility/development pilot complete — November 1, 2026",
  "Eight new exogenous layouts; 16 unique underlying fields; 32 paired-severity analyses.",
  "Oracle covariance supplied. No confidence interval or independent confirmation claim.",
  "All statuses, fit shifts and independent integration audits retained."),
  file.path(out, "PILOT_COMPLETED.txt"))
cat("All 32 feasibility cases recorded.\n")
