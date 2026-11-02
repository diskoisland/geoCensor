#!/usr/bin/env Rscript
# Frozen fresh-response point-estimation pilot; not a coverage experiment.
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
source("scripts/analysis/spatial_censored_scale_reference.R")
source("scripts/analysis/conley_survreg_reference.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Supply an optional output directory")
out <- if (length(args)) args[1L] else "output/spatial_scale_fresh_2026-11-01"
protocol <- file.path(out, "PROTOCOL.md")
if (!file.exists(protocol)) stop("Freeze the protocol before generating responses")
config <- list(n_sites = 24L, repeats = 2L, layouts_per_regime = 8L,
  regimes = c("uniform", "clustered"), decays = c(10, 20),
  truth = log(50), latent_sd = .55, noise_sd = .2,
  left_offset = c(moderate = -.3, heavy = .6), right_offset = 1,
  cell_size = .25, origins_per_axis = 10L,
  coarse_draws = 8192L, fine_draws = 32768L,
  audit_draws = 8192L, audit_batches = 8L,
  location_review = .02, log_multiplier_review = .05, audit_mc_se_review = .05)
sources <- c("scripts/analysis/validate_spatial_estimated_scale_fresh.R",
  "scripts/analysis/spatial_censored_scale_reference.R",
  "scripts/analysis/spatial_censored_likelihood_reference.R",
  "scripts/analysis/conley_survreg_reference.R",
  sort(list.files("R", pattern = "\\.[Rr]$", full.names = TRUE)))
signature <- list(config = config, sources = tools::md5sum(sources),
  protocol = tools::md5sum(protocol), R = R.version.string,
  platform = R.version$platform, mvtnorm = as.character(packageVersion("mvtnorm")),
  survival = as.character(packageVersion("survival")))
manifest <- file.path(out, "study_manifest.rds")
if (file.exists(manifest)) {
  if (!identical(readRDS(manifest), signature)) stop("Study provenance changed; use a fresh directory")
} else {
  if (length(setdiff(list.files(out, all.files = TRUE, no.. = TRUE), "PROTOCOL.md")))
    stop("Unverified study directory")
  dir.create(file.path(out, "source_snapshot"))
  invisible(file.copy(sources, file.path(out, "source_snapshot")))
  saveRDS(signature, manifest)
  write.csv(data.frame(file = c(sources, protocol),
    md5 = unname(c(signature$sources, signature$protocol))),
    file.path(out, "source_input_manifest.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
}
if (file.exists(file.path(out, "COMPLETED.txt"))) {
  cat("Fresh pilot completed; provenance matches and results unchanged.\n")
  quit(status = 0L)
}
for (directory in c("layouts", "fields", "cases"))
  dir.create(file.path(out, directory), showWarnings = FALSE)

failure <- function(status, message = "") list(status = status, theta = NA_real_,
  multiplier = NA_real_, loglik = NA_real_, candidates = data.frame(), message = message)
safe <- function(expr) tryCatch(expr, error = function(e)
  failure("numerical_failure", conditionMessage(e)))
independent <- function(lower, upper, weights, n_exact) {
  if (!n_exact) return(failure("no_finite_common_window_scale_mle"))
  safe({
    left <- is.infinite(lower); right <- is.infinite(upper)
    df <- data.frame(detected_val = ifelse(left | right, NA_real_, exp(lower)),
      left_cen = left, left_limit = ifelse(left, exp(upper), NA_real_),
      right_cen = right, right_limit = ifelse(right, exp(lower), NA_real_))
    fit <- fit_censored_score_reference(df, weights)
    list(status = "finite", theta = fit$theta[1L],
      multiplier = exp(fit$theta[2L]), message = "")
  })
}
total <- length(config$regimes) * config$layouts_per_regime *
  length(config$decays) * length(config$left_offset)
for (r in seq_along(config$regimes)) for (i in seq_len(config$layouts_per_regime)) {
  index <- (r - 1L) * config$layouts_per_regime + i
  key <- sprintf("%s_%02d", config$regimes[r], i)
  location_file <- file.path(out, "layouts", paste0(key, ".rds"))
  if (file.exists(location_file)) layout <- readRDS(location_file) else {
    set.seed(510010000L + index)
    s <- if (r == 1L) matrix(runif(48L), ncol = 2L) else
      rbind(matrix(runif(36L, 0, .25), ncol = 2L), matrix(runif(12L), ncol = 2L))
    weights <- cell_declustering_weights(data.frame(x = s[, 1L], y = s[, 2L]),
      "x", "y", config$cell_size, config$origins_per_axis)
    layout <- list(key = key, seed = 510010000L + index, s = s, weights = weights)
    saveRDS(layout, location_file)
  }
  for (j in seq_along(config$decays)) {
    field_index <- (index - 1L) * length(config$decays) + j
    field_key <- paste0(key, "_decay", config$decays[j])
    field_file <- file.path(out, "fields", paste0(field_key, ".rds"))
    if (file.exists(field_file)) field <- readRDS(field_file) else {
      set.seed(520010000L + field_index)
      correlation <- exp(-config$decays[j] * as.matrix(dist(layout$s)))
      id <- rep(seq_len(config$n_sites), each = config$repeats)
      latent <- as.vector(t(chol(correlation)) %*% rnorm(config$n_sites)) * config$latent_sd
      y <- config$truth + latent[id] + rnorm(length(id), sd = config$noise_sd)
      covariance <- config$latent_sd^2 * correlation[id, id] +
        diag(config$noise_sd^2, length(id))
      inv <- chol2inv(chol(covariance))
      complete_mu <- sum(inv %*% y) / sum(inv)
      complete_multiplier <- sqrt(drop(crossprod(y - complete_mu,
        inv %*% (y - complete_mu))) / length(y))
      field <- list(key = field_key, seed = 520010000L + field_index, id = id,
        log_values = y, covariance = covariance,
        weights = layout$weights[id] / config$repeats,
        complete_mu = complete_mu, complete_multiplier = complete_multiplier)
      saveRDS(field, field_file)
    }
    for (severity in seq_along(config$left_offset)) {
      case_index <- (field_index - 1L) * length(config$left_offset) + severity
      case_key <- paste0(field_key, "_", names(config$left_offset)[severity])
      case_file <- file.path(out, "cases", paste0(case_key, ".rds"))
      if (file.exists(case_file)) next
      started <- proc.time()[[3L]]
      left_limit <- config$truth + config$left_offset[severity]
      right_limit <- config$truth + config$right_offset
      left <- field$log_values < left_limit; right <- field$log_values > right_limit
      lower <- ifelse(left, -Inf, ifelse(right, right_limit, field$log_values))
      upper <- ifelse(right, Inf, ifelse(left, left_limit, field$log_values))
      model <- spatial_censored_setup_reference(lower, upper, field$covariance)
      coarse_u <- spatial_censored_uniforms_reference(model, config$coarse_draws,
        530010000L + case_index)
      fine_u <- spatial_censored_uniforms_reference(model, config$fine_draws,
        540010000L + case_index)
      oracle <- safe(spatial_censored_fit_mu_reference(model, fine_u))
      if (!model$no) {
        coarse <- fine <- spatial_censored_fit_scale_reference(model)
      } else {
        # Location-only oracle is a starting device, not the fit target or selection criterion.
        center <- if (oracle$status == "finite") oracle$theta else median(field$log_values[!left & !right])
        coarse <- safe(spatial_censored_fit_scale_reference(model, coarse_u,
          starts = rbind(c(center, 1), c(center / 5, .2))))
        if (coarse$status == "finite") {
          fine <- safe(spatial_censored_fit_scale_reference(model, fine_u,
            starts = rbind(c(coarse$theta / coarse$multiplier, 1 / coarse$multiplier), c(center, 1))))
        } else fine <- coarse
      }
      audit <- NULL
      if (fine$status == "finite") {
        offsets <- rbind(c(-.02, 0), c(0, 0), c(.02, 0), c(0, -.02), c(0, .02))
        audit <- safe({
          values <- matrix(NA_real_, config$audit_batches, nrow(offsets))
          for (b in seq_len(config$audit_batches)) {
            u <- spatial_censored_uniforms_reference(model, config$audit_draws,
              550010000L + case_index * 100L + b)
            values[b, ] <- vapply(seq_len(nrow(offsets)), function(k) {
              multiplier <- fine$multiplier * exp(offsets[k, 2L])
              spatial_censored_scaled_loglik_reference(model,
                (fine$theta + offsets[k, 1L]) / multiplier, 1 / multiplier, u)
            }, 0.0)
          }
          list(status = "finite", values = values)
        })
      }
      equal <- independent(lower, upper, rep(1, model$n), model$no)
      weighted <- independent(lower, upper, field$weights, model$no)
      shift <- fine$theta - coarse$theta
      scale_shift <- log(fine$multiplier / coarse$multiplier)
      audited <- !is.null(audit) && audit$status == "finite"
      center_se <- if (audited) sd(audit$values[, 2L]) / sqrt(config$audit_batches) else NA_real_
      directional <- FALSE
      if (audited) {
        deltas <- sweep(audit$values[, c(1L, 3L, 4L, 5L), drop = FALSE], 1L, audit$values[, 2L])
        directional <- any(colMeans(deltas) > 3 * apply(deltas, 2L, sd) / sqrt(config$audit_batches))
      }
      review <- fine$status == "finite" && (!audited ||
        abs(shift) > config$location_review || abs(scale_shift) > config$log_multiplier_review ||
        center_se > config$audit_mc_se_review || directional)
      row <- data.frame(key = case_key, layout_key = key, regime = config$regimes[r],
        layout = i, decay = config$decays[j], censoring = names(config$left_offset)[severity],
        n = model$n, n_exact = model$no, n_left = sum(left), n_right = sum(right),
        public_count_eligible = model$no >= 2L,
        weighted_exact_fraction = sum(field$weights[!left & !right]),
        spatial_status = fine$status, spatial_theta = fine$theta,
        spatial_multiplier = fine$multiplier,
        spatial_sigma = sqrt(field$covariance[1L, 1L]) * fine$multiplier,
        oracle_status = oracle$status, oracle_theta = if (oracle$status == "finite") oracle$theta else NA_real_,
        equal_status = equal$status, equal_theta = equal$theta, equal_sigma = equal$multiplier,
        weighted_status = weighted$status, weighted_theta = weighted$theta, weighted_sigma = weighted$multiplier,
        complete_theta = field$complete_mu,
        complete_sigma = sqrt(field$covariance[1L, 1L]) * field$complete_multiplier,
        refit_theta_shift = shift, refit_log_multiplier_shift = scale_shift,
        audit_status = if (is.null(audit)) "not_attempted" else audit$status,
        audit_loglik_mc_se = center_se, audit_positive_direction = directional,
        numerical_review = review,
        elapsed_seconds = proc.time()[[3L]] - started,
        message = if (is.null(fine$message)) "" else fine$message)
      saveRDS(list(row = row, lower = lower, upper = upper,
        coarse = coarse, fine = fine, oracle = oracle, audit = audit,
        equal = equal, weighted = weighted,
        integration_seeds = c(coarse = 530010000L + case_index,
          fine = 540010000L + case_index,
          setNames(550010000L + case_index * 100L + seq_len(config$audit_batches),
            paste0("audit", seq_len(config$audit_batches))))), case_file)
      cat(case_index, "/", total, case_key, fine$status,
        "exact", model$no, "review", review, "seconds", round(row$elapsed_seconds, 1), "\n")
      flush.console()
    }
  }
}
paths <- list.files(file.path(out, "cases"), pattern = "\\.rds$", full.names = TRUE)
cases <- lapply(paths, readRDS)
stopifnot(length(cases) == total,
  identical(tools::md5sum(sources), signature$sources),
  identical(tools::md5sum(protocol), signature$protocol))
results <- do.call(rbind, lapply(cases, `[[`, "row"))
stopifnot(!anyDuplicated(results$key), all(results$n_exact + results$n_left + results$n_right == results$n))
write.csv(results, file.path(out, "fresh_results.csv"), row.names = FALSE)
writeLines(c("Fresh estimated-scale point-estimation pilot complete — November 1, 2026.",
  "16 independent layouts; 32 fields; 64 paired censoring analyses, 48 records each.",
  "No coverage, population CI, preferential-sampling correction or fully estimated covariance claim.",
  "All unavailable and numerical-review cases retained."), file.path(out, "COMPLETED.txt"))
cat("All", total, "fresh cases recorded.\n")
