#!/usr/bin/env Rscript
source("scripts/analysis/spatial_censored_covariance_reference.R")
source("scripts/analysis/spatial_censored_scale_reference.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Supply an optional output directory")
out <- if (length(args)) args[1L] else "output/spatial_covariance_fresh_2026-11-01"
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
config <- list(n_sites = 24L, repeats = 2L, layouts_per_regime = 4L,
  regimes = c("uniform", "clustered"), decay_by_layout = c(10, 20, 10, 20),
  truth = log(50), latent_sd = .55, noise_sd = .2,
  left_offset = c(moderate = -.3, heavy = .6), right_offset = 1,
  coarse_draws = 2048L, fine_draws = 8192L,
  coarse_maxit = 800L, fine_maxit = 1200L, reltol = 1e-6,
  audit_draws = 16384L, audit_batches = 8L)
sources <- c("scripts/analysis/validate_spatial_covariance_fresh.R",
  "scripts/analysis/spatial_censored_covariance_reference.R",
  "scripts/analysis/verify_spatial_censored_covariance_reference.R",
  "scripts/analysis/spatial_censored_likelihood_reference.R",
  "scripts/analysis/spatial_censored_scale_reference.R")
inputs <- file.path(out, c("PROTOCOL.md", "verification.csv", "VERIFIED.txt",
  "verification_development_fit.rds"))
stopifnot(all(file.exists(inputs)))
signature <- list(config = config, sources = tools::md5sum(sources), inputs = tools::md5sum(inputs),
  R = R.version.string, platform = R.version$platform, RNG = RNGkind(),
  mvtnorm = as.character(packageVersion("mvtnorm")), survival = as.character(packageVersion("survival")))
manifest <- file.path(out, "study_manifest.rds")
if (file.exists(manifest)) {
  if (!identical(readRDS(manifest), signature)) stop("Covariance pilot source/protocol/runtime changed")
} else {
  allowed <- basename(inputs)
  if (length(setdiff(list.files(out, all.files = TRUE, no.. = TRUE), allowed))) stop("Unverified pilot files")
  dir.create(file.path(out, "source_snapshot"))
  invisible(file.copy(sources, file.path(out, "source_snapshot")))
  saveRDS(signature, manifest)
  write.csv(data.frame(file = c(sources, inputs), md5 = unname(c(signature$sources, signature$inputs))),
    file.path(out, "source_input_manifest.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
}
if (file.exists(file.path(out, "COMPLETED.txt"))) {
  cat("Completed covariance pilot matches provenance; no results changed.\n")
  quit(status = 0L)
}
for (directory in c("layouts", "fields", "cases"))
  dir.create(file.path(out, directory), showWarnings = FALSE)
missing_fit <- function(status, message = "") list(status = status, theta = NA_real_,
  sigma = NA_real_, range = NA_real_, nugget_fraction = NA_real_,
  candidate = NULL, candidates = data.frame(), evaluations = data.frame(), message = message)
safe <- function(expr) tryCatch(expr, error = function(e) missing_fit("numerical_failure", conditionMessage(e)))
independent <- function(lower, upper) {
  if (!any(lower == upper)) return(missing_fit("no_exact_scale_identification"))
  safe({
    fit <- withCallingHandlers(survival::survreg(
      survival::Surv(lower, upper, type = "interval2") ~ 1, dist = "gaussian"),
      warning = function(w) stop(conditionMessage(w)))
    residual <- residuals(fit, type = "matrix", weighted = FALSE)
    info <- -matrix(c(sum(residual[, "ddg"]), sum(residual[, "dsg"]),
      sum(residual[, "dsg"]), sum(residual[, "dds"])), 2L, 2L)
    if (!is.finite(coef(fit)) || !is.finite(fit$scale) || fit$scale <= 0 ||
        any(!is.finite(info)) || min(eigen(info, symmetric = TRUE, only.values = TRUE)$values) <= 0 ||
        max(abs(colSums(residual[, c("dg", "ds"), drop = FALSE]))) / length(lower) > 1e-5)
      stop("Independent fit failed finite/score/information checks")
    list(status = "finite", theta = unname(coef(fit)), sigma = fit$scale)
  })
}
oracle_control <- function(lower, upper, covariance, u) {
  safe({
    model <- spatial_censored_setup_reference(lower, upper, covariance)
    fit <- spatial_censored_fit_scale_reference(model, u)
    if (!nrow(fit$candidates)) return(list(status = fit$status, theta = NA_real_, sigma = NA_real_, fit = fit))
    finite <- which(is.finite(fit$candidates$loglik))
    if (!length(finite)) return(list(status = "no_finite_candidate", theta = NA_real_, sigma = NA_real_, fit = fit))
    best <- fit$candidates[finite[which.max(fit$candidates$loglik[finite])], ]
    status <- if (best$convergence != 0L) "optimizer_unresolved" else if (
      best$precision <= 1.01e-5 || best$precision >= .99e3) "search_boundary" else "finite"
    list(status = status, theta = if (status == "finite") best$theta else NA_real_,
      sigma = if (status == "finite") sqrt(covariance[1L, 1L]) * best$multiplier else NA_real_, fit = fit)
  })
}

for (r in seq_along(config$regimes)) for (i in seq_len(config$layouts_per_regime)) {
  index <- (r - 1L) * config$layouts_per_regime + i
  key <- sprintf("%s_%02d", config$regimes[r], i)
  layout_file <- file.path(out, "layouts", paste0(key, ".rds"))
  if (file.exists(layout_file)) layout <- readRDS(layout_file) else {
    set.seed(710010000L + index)
    sites <- if (r == 1L) matrix(runif(48L), ncol = 2L) else
      rbind(matrix(runif(36L, 0, .25), ncol = 2L), matrix(runif(12L), ncol = 2L))
    layout <- list(key = key, seed = 710010000L + index, sites = sites)
    saveRDS(layout, layout_file)
  }
  field_file <- file.path(out, "fields", paste0(key, ".rds"))
  if (file.exists(field_file)) field <- readRDS(field_file) else {
    set.seed(720010000L + index)
    decay <- config$decay_by_layout[i]
    correlation <- exp(-decay * as.matrix(dist(layout$sites)))
    id <- rep(seq_len(config$n_sites), each = config$repeats)
    latent <- as.vector(t(chol(correlation)) %*% rnorm(config$n_sites)) * config$latent_sd
    values <- config$truth + latent[id] + rnorm(length(id), sd = config$noise_sd)
    S <- config$latent_sd^2 * correlation[id, id] + diag(config$noise_sd^2, length(id))
    field <- list(key = key, seed = 720010000L + index, site_id = id,
      coordinates = layout$sites[id, ], log_values = values, covariance = S,
      true_range = 1 / decay,
      true_sigma = sqrt(config$latent_sd^2 + config$noise_sd^2),
      true_nugget = config$noise_sd^2 / (config$latent_sd^2 + config$noise_sd^2))
    saveRDS(field, field_file)
  }
  for (severity in seq_along(config$left_offset)) {
    case_index <- (index - 1L) * length(config$left_offset) + severity
    case_key <- paste0(key, "_", names(config$left_offset)[severity])
    path <- file.path(out, "cases", paste0(case_key, ".rds"))
    if (file.exists(path)) next
    started <- proc.time()[[3L]]
    low <- config$truth + config$left_offset[severity]; high <- config$truth + config$right_offset
    left <- field$log_values < low; right <- field$log_values > high
    lower <- ifelse(left, -Inf, ifelse(right, high, field$log_values))
    upper <- ifelse(right, Inf, ifelse(left, low, field$log_values))
    exact <- lower == upper
    unit_model <- spatial_censored_setup_reference(lower, upper, diag(length(lower)))
    coarse_u <- spatial_censored_uniforms_reference(unit_model, config$coarse_draws, 730010000L + case_index)
    fine_u <- spatial_censored_uniforms_reference(unit_model, config$fine_draws, 740010000L + case_index)
    # Primary fitting receives only observed bounds, coordinates and integration
    # uniforms. No generating covariance/parameter or latent response is passed.
    coarse <- safe(spatial_censored_fit_covariance_reference(lower, upper, field$coordinates,
      coarse_u, maxit = config$coarse_maxit, reltol = config$reltol))
    fine_attempted <- !is.null(coarse$candidate)
    if (fine_attempted) {
      data_start <- c(median(lower[exact]), log(.3), log(.03), qlogis(.05))
      fine <- safe(spatial_censored_fit_covariance_reference(lower, upper, field$coordinates,
        fine_u, starts = rbind(coarse$candidate$par, data_start),
        maxit = config$fine_maxit, reltol = config$reltol))
    } else fine <- coarse
    audit <- NULL; review <- NA
    shift <- rep(NA_real_, 4L); center_se <- NA_real_; directional <- NA
    if (!is.null(fine$candidate)) {
      raw <- fine$candidate
      if (!is.null(coarse$candidate)) shift <- c(raw$theta - coarse$candidate$theta,
        log(raw$sigma / coarse$candidate$sigma), log(raw$range / coarse$candidate$range),
        raw$nugget_fraction - coarse$candidate$nugget_fraction)
      offsets <- rbind(c(0, 0, 0, 0), c(-.02, 0, 0, 0), c(.02, 0, 0, 0),
        c(0, -.02, 0, 0), c(0, .02, 0, 0), c(0, 0, -.05, 0), c(0, 0, .05, 0),
        c(0, 0, 0, -.05), c(0, 0, 0, .05))
      audit <- tryCatch({
        values <- matrix(NA_real_, config$audit_batches, nrow(offsets))
        for (b in seq_len(config$audit_batches)) {
          u <- spatial_censored_uniforms_reference(unit_model, config$audit_draws,
            750010000L + case_index * 100L + b)
          for (j in seq_len(nrow(offsets))) {
            par <- raw$par + offsets[j, ]
            values[b, j] <- spatial_covariance_loglik_reference(lower, upper, field$coordinates,
              par[1L], exp(par[2L]), exp(par[3L]), plogis(par[4L]), u)
          }
        }
        list(status = "finite", values = values)
      }, error = function(e) list(status = "unavailable", message = conditionMessage(e)))
      audited <- audit$status == "finite"
      if (audited) {
        center_se <- sd(audit$values[, 1L]) / sqrt(config$audit_batches)
        delta <- sweep(audit$values[, -1L, drop = FALSE], 1L, audit$values[, 1L])
        directional <- any(colMeans(delta) > 3 * apply(delta, 2L, sd) / sqrt(config$audit_batches))
      }
      review <- !audited || any(abs(shift) > c(.05, .1, .2, .05)) ||
        center_se > .05 || isTRUE(directional)
    }
    # Truth is supplied only to this separately labeled comparator.
    control <- oracle_control(lower, upper, field$covariance, fine_u)
    equal <- independent(lower, upper)
    raw <- fine$candidate
    value <- function(name) if (is.null(raw)) NA_real_ else raw[[name]]
    row <- data.frame(key = case_key, layout_key = key, regime = config$regimes[r],
      layout = i, decay = config$decay_by_layout[i], censoring = names(config$left_offset)[severity],
      n = length(lower), n_exact = sum(exact), n_left = sum(left), n_right = sum(right),
      n_exact_sites = length(unique(field$site_id[exact])), public_count_eligible = sum(exact) >= 2L,
      fit_status = fine$status, fine_attempted = fine_attempted,
      theta = fine$theta, sigma = fine$sigma, range = fine$range, nugget_fraction = fine$nugget_fraction,
      candidate_theta = value("theta"), candidate_sigma = value("sigma"),
      candidate_range = value("range"), candidate_nugget_fraction = value("nugget_fraction"),
      true_range = field$true_range, true_sigma = field$true_sigma, true_nugget = field$true_nugget,
      oracle_status = control$status, oracle_theta = control$theta, oracle_sigma = control$sigma,
      independent_status = equal$status, independent_theta = equal$theta, independent_sigma = equal$sigma,
      theta_shift = shift[1L], log_sigma_shift = shift[2L], log_range_shift = shift[3L],
      nugget_shift = shift[4L], audit_status = if (is.null(audit)) "not_attempted" else audit$status,
      audit_loglik_mc_se = center_se, audit_improving_direction = directional, numerical_review = review,
      elapsed_seconds = proc.time()[[3L]] - started,
      message = if (is.null(fine$message)) "" else fine$message)
    saveRDS(list(row = row, lower = lower, upper = upper, coarse = coarse, fine = fine,
      audit = audit, oracle = control, independent = equal,
      seeds = c(coarse = 730010000L + case_index, fine = 740010000L + case_index,
        setNames(750010000L + case_index * 100L + seq_len(config$audit_batches),
          paste0("audit", seq_len(config$audit_batches))))), path)
    cat(case_index, "/ 16", case_key, fine$status, "exact", sum(exact),
      "numerical review", review, "seconds", round(row$elapsed_seconds, 1), "\n")
    flush.console()
  }
}
cases <- lapply(list.files(file.path(out, "cases"), pattern = "\\.rds$", full.names = TRUE), readRDS)
stopifnot(length(cases) == 16L, identical(tools::md5sum(sources), signature$sources),
  identical(tools::md5sum(inputs), signature$inputs))
results <- do.call(rbind, lapply(cases, `[[`, "row"))
stopifnot(!anyDuplicated(results$key), all(results$n_exact + results$n_left + results$n_right == results$n))
unresolved <- results$fit_status != "finite_interior"
stopifnot(all(is.na(results$theta[unresolved])), all(is.na(results$range[unresolved])),
  all(is.na(results$sigma[unresolved])), all(is.na(results$nugget_fraction[unresolved])))
write.csv(results, file.path(out, "covariance_results.csv"), row.names = FALSE)
writeLines(c("Fresh estimated-covariance pilot complete — November 1, 2026.",
  "Eight layouts/fields, sixteen paired-severity analyses, 48 records each.",
  "Missing and unresolved accepted estimates remain NA; raw candidates are separate.",
  "No population interval or production validity claim; confidence limits remain experimental."),
  file.path(out, "COMPLETED.txt"))
cat("All sixteen fresh covariance attempts retained.\n")
