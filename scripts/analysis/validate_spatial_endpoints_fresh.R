#!/usr/bin/env Rscript
source("scripts/analysis/spatial_covariance_endpoints_reference.R")
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1L] else "output/spatial_endpoints_fresh_2026-11-01"
mode <- if (length(args) >= 2L) args[2L] else "init"
if (length(args) > 2L || !mode %in% c("init", "worker1", "worker2", "finish")) stop("Use init, worker1, worker2 or finish")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
config <- list(sites = 24L, repeats = 2L, layouts_per_regime = 2L, workers = 2L,
  regimes = c("uniform", "clustered"), structures = c("spatial", "site_only", "independent"),
  truth = log(50), sigma = sqrt(.55^2 + .2^2), nugget = .2^2/(.55^2 + .2^2), range = .08,
  offsets = c(moderate = -.3, heavy = .6), upper_offset = 1,
  coarse_draws = 2048L, fine_draws = 8192L, audit_draws = 32768L, audit_batches = 8L)
sources <- c("scripts/analysis/validate_spatial_endpoints_fresh.R",
  "scripts/analysis/spatial_covariance_endpoints_reference.R",
  "scripts/analysis/verify_spatial_covariance_endpoints_reference.R",
  "scripts/analysis/spatial_censored_covariance_reference.R",
  "scripts/analysis/spatial_censored_likelihood_reference.R")
inputs <- c(file.path(out, c("PROTOCOL.md", "verification.csv", "VERIFIED.txt", "verification_site_fit.rds",
  "NUMERICAL_DEVELOPMENT.md")), list.files(file.path(out, "development_snapshot"), full.names = TRUE))
stopifnot(all(file.exists(inputs)))
signature <- list(config = config, sources = tools::md5sum(sources), inputs = tools::md5sum(inputs),
  R = R.version.string, platform = R.version$platform, RNG = RNGkind(),
  mvtnorm = as.character(packageVersion("mvtnorm")), survival = as.character(packageVersion("survival")))
manifest <- file.path(out, "study_manifest.rds")
if (file.exists(manifest)) {
  if (!identical(readRDS(manifest), signature)) stop("Endpoint study source/protocol/runtime changed")
} else {
  if (mode != "init") stop("Initialize the frozen study first")
  allowed <- c(basename(inputs), "development_snapshot")
  if (length(setdiff(list.files(out), allowed))) stop("Unverified endpoint study outputs exist")
  dir.create(file.path(out, "source_snapshot"))
  invisible(file.copy(sources, file.path(out, "source_snapshot")))
  saveRDS(signature, manifest)
  write.csv(data.frame(file = c(sources, inputs), md5 = unname(c(signature$sources, signature$inputs))),
    file.path(out, "source_input_manifest.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
}
verify_design <- function() {
  inventory <- read.csv(file.path(out, "design_manifest.csv"), stringsAsFactors = FALSE)
  stopifnot(identical(unname(tools::md5sum(inventory$file)), inventory$md5))
}
if (file.exists(file.path(out, "INITIALIZED.txt"))) verify_design()
if (file.exists(file.path(out, "COMPLETED.txt"))) {
  cat("Completed endpoint study matches source/runtime/design fingerprints; no results changed.\n")
  quit(status = 0L)
}
for (directory in c("layouts", "fields", "cases")) dir.create(file.path(out, directory), showWarnings = FALSE)
if (mode == "init") {
  if (file.exists(file.path(out, "INITIALIZED.txt"))) {
    cat("Frozen endpoint layouts/fields already initialized.\n"); quit(status = 0L)
  }
  jobs <- list(); field_index <- case_index <- 0L
  for (r in seq_along(config$regimes)) for (i in seq_len(config$layouts_per_regime)) {
    index <- (r - 1L) * config$layouts_per_regime + i
    key <- sprintf("%s_%02d", config$regimes[r], i)
    set.seed(810010000L + index)
    sites <- if (r == 1L) matrix(runif(48L), ncol = 2L) else
      rbind(matrix(runif(36L, 0, .25), ncol = 2L), matrix(runif(12L), ncol = 2L))
    saveRDS(list(key = key, seed = 810010000L + index, sites = sites), file.path(out, "layouts", paste0(key, ".rds")))
    id <- rep(seq_len(config$sites), each = config$repeats)
    coordinates <- sites[id, ]; distance <- as.matrix(dist(coordinates))
    for (structure in config$structures) {
      field_index <- field_index + 1L
      field_key <- paste(key, structure, sep = "_")
      kernel <- switch(structure, spatial = exp(-distance/config$range),
        site_only = outer(id, id, "==") * 1, independent = diag(length(id)))
      nu <- if (structure == "independent") 1 else config$nugget
      S <- config$sigma^2 * ((1 - nu) * kernel + diag(nu, length(id)))
      set.seed(820010000L + field_index)
      log_values <- config$truth + as.vector(t(chol(S)) %*% rnorm(length(id)))
      field <- list(key = field_key, layout_key = key, seed = 820010000L + field_index,
        coordinates = coordinates, site_id = id, covariance = S, log_values = log_values,
        truth_structure = structure, truth_mu = config$truth, truth_sigma = config$sigma,
        truth_range = if (structure == "spatial") config$range else NA_real_, truth_nugget = nu)
      saveRDS(field, file.path(out, "fields", paste0(field_key, ".rds")))
      for (severity in names(config$offsets)) {
        case_index <- case_index + 1L
        jobs[[case_index]] <- data.frame(key = paste(field_key, severity, sep = "_"),
          field_key = field_key, layout_key = key, regime = config$regimes[r],
          truth_structure = structure, censoring = severity, case_index = case_index,
          worker = (field_index - 1L) %% config$workers + 1L)
      }
    }
  }
  write.csv(do.call(rbind, jobs), file.path(out, "jobs.csv"), row.names = FALSE)
  files <- c(file.path(out, "jobs.csv"), list.files(file.path(out, "layouts"), full.names = TRUE),
    list.files(file.path(out, "fields"), full.names = TRUE))
  write.csv(data.frame(file = files, md5 = unname(tools::md5sum(files))), file.path(out, "design_manifest.csv"), row.names = FALSE)
  writeLines("Four layouts, twelve independent conditional response fields and twenty-four paired analyses frozen.",
    file.path(out, "INITIALIZED.txt"))
  cat("Four fresh layouts/twelve fields initialized; worker1 and worker2 can run independently.\n")
  quit(status = 0L)
}
stopifnot(file.exists(file.path(out, "INITIALIZED.txt")))
jobs <- read.csv(file.path(out, "jobs.csv"), stringsAsFactors = FALSE)
if (mode == "finish") {
  paths <- file.path(out, "cases", paste0(jobs$key, ".rds"))
  if (!all(file.exists(paths))) stop("Some endpoint cases are incomplete; do not fill missing checkpoints")
  cases <- lapply(paths, readRDS); rows <- do.call(rbind, lapply(cases, `[[`, "row"))
  stopifnot(nrow(rows) == 24L, identical(rows$key, jobs$key), !anyDuplicated(rows$key))
  unresolved <- rows$status != "finite_numerical_winner"
  stopifnot(all(is.na(rows$theta[unresolved])), all(is.na(rows$range[unresolved])),
    all(is.na(rows$sigma[unresolved])), all(is.na(rows$nugget_fraction[unresolved])),
    all(is.na(rows$range[!is.na(rows$selected_structure) & rows$selected_structure != "spatial"])))
  write.csv(rows, file.path(out, "endpoint_results.csv"), row.names = FALSE)
  stopifnot(identical(tools::md5sum(sources), signature$sources), identical(tools::md5sum(inputs), signature$inputs))
  verify_design()
  writeLines(c("Fresh endpoint pilot complete — November 1, 2026.",
    "Four layouts, twelve fields, twenty-four analyses; unsupported/ambiguous parameters remain NA.",
    "Numerical candidate comparison is not a statistical family test or confidence interval."), file.path(out, "COMPLETED.txt"))
  cat("All twenty-four endpoint attempts retained and finalized.\n"); quit(status = 0L)
}
worker <- as.integer(sub("worker", "", mode))
for (j in which(jobs$worker == worker)) {
  job <- jobs[j, ]; path <- file.path(out, "cases", paste0(job$key, ".rds"))
  if (file.exists(path)) {stopifnot(identical(readRDS(path)$row$key, job$key)); next}
  started <- proc.time()[[3L]]
  field <- readRDS(file.path(out, "fields", paste0(job$field_key, ".rds")))
  low <- config$truth + config$offsets[[job$censoring]]; high <- config$truth + config$upper_offset
  left <- field$log_values < low; right <- field$log_values > high
  lower <- ifelse(left, -Inf, ifelse(right, high, field$log_values))
  upper <- ifelse(right, Inf, ifelse(left, low, field$log_values))
  exact <- lower == upper; i <- job$case_index
  # The primary estimator receives only observed bounds/coordinates and seeds.
  fit <- tryCatch(spatial_censored_fit_endpoints_reference(lower, upper, field$coordinates,
    830010000L + i, 840010000L + i, 850010000L + 100L * i,
    coarse_draws = config$coarse_draws, fine_draws = config$fine_draws,
    audit_draws = config$audit_draws, audit_batches = config$audit_batches),
    error = function(e) list(status = "numerical_failure", structure = NA_character_,
      theta = NA_real_, sigma = NA_real_, range = NA_real_, nugget_fraction = NA_real_,
      branches = list(), message = conditionMessage(e)))
  value <- function(name, default = NA_real_) if (is.null(fit[[name]])) default else fit[[name]]
  candidate_structure <- value("candidate_structure", NA_character_)
  candidate <- if (!is.na(candidate_structure)) fit$branches[[candidate_structure]]$candidate else NULL
  raw <- function(name) if (is.null(candidate)) NA_real_ else candidate[[name]]
  branch_status <- function(name) if (is.null(fit$branches[[name]])) "not_attempted" else fit$branches[[name]]$status
  row <- data.frame(key = job$key, field_key = job$field_key, layout_key = job$layout_key,
    regime = job$regime, truth_structure = job$truth_structure, censoring = job$censoring,
    n = length(lower), n_exact = sum(exact), n_left = sum(left), n_right = sum(right),
    n_exact_sites = length(unique(field$site_id[exact])), public_count_eligible = sum(exact) >= 2L,
    status = fit$status, selected_structure = fit$structure,
    theta = fit$theta, sigma = fit$sigma, range = fit$range, nugget_fraction = fit$nugget_fraction,
    candidate_structure = candidate_structure, candidate_theta = raw("mu"), candidate_sigma = raw("sigma"),
    candidate_range = raw("range"), candidate_nugget = raw("nugget"),
    spatial_status = branch_status("spatial"), site_only_status = branch_status("site_only"),
    independent_status = branch_status("independent"),
    min_audited_gap = if (is.null(fit$comparison)) NA_real_ else min(fit$comparison$mean_gap),
    max_numerical_margin = if (is.null(fit$comparison)) NA_real_ else max(fit$comparison$numerical_margin),
    deterministic_crosscheck_error = if (is.null(fit$deterministic_error)) NA_real_ else max(fit$deterministic_error),
    elapsed_seconds = proc.time()[[3L]] - started,
    message = if (is.null(fit$message)) "" else fit$message)
  saveRDS(list(row = row, lower = lower, upper = upper, fit = fit,
    seeds = c(coarse = 830010000L + i, fine = 840010000L + i, audit_base = 850010000L + 100L * i)), path)
  cat("worker", worker, i, "/24", job$key, fit$status, "exact", sum(exact),
    "seconds", round(row$elapsed_seconds, 1), "\n"); flush.console()
}
verify_design()
stopifnot(identical(tools::md5sum(sources), signature$sources), identical(tools::md5sum(inputs), signature$inputs))
cat("Endpoint worker", worker, "complete.\n")
