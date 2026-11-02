#!/usr/bin/env Rscript
source("scripts/analysis/spatial_censored_scale_reference.R")
pilot <- "output/spatial_likelihood_reflected_2026-11-01"
out <- "output/spatial_scale_development_2026-11-01"
keys <- c("uniform_01_decay10_heavy", "clustered_03_decay20_heavy",
          "clustered_03_decay20_moderate", "clustered_03_decay10_heavy")
sources <- c("scripts/analysis/diagnose_spatial_estimated_scale.R",
  "scripts/analysis/spatial_censored_scale_reference.R",
  "scripts/analysis/spatial_censored_likelihood_reference.R")
inputs <- c(file.path(out, "PROTOCOL.md"), file.path(pilot, "cases", paste0(keys, ".rds")),
  unique(file.path(pilot, "fields", paste0(sub("_(moderate|heavy)$", "", keys), ".rds"))))
signature <- list(sources = tools::md5sum(sources), inputs = tools::md5sum(inputs),
  R = R.version.string, mvtnorm = as.character(packageVersion("mvtnorm")),
  survival = as.character(packageVersion("survival")))
manifest <- file.path(out, "manifest.rds")
if (file.exists(manifest)) {
  if (!identical(readRDS(manifest), signature)) stop("Development provenance changed")
} else {
  if (length(setdiff(list.files(out), "PROTOCOL.md"))) stop("Unverified development directory")
  dir.create(file.path(out, "source_snapshot"))
  invisible(file.copy(sources, file.path(out, "source_snapshot")))
  saveRDS(signature, manifest)
  write.csv(data.frame(file = c(sources, inputs), md5 = unname(c(signature$sources, signature$inputs))),
    file.path(out, "source_input_manifest.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
}
if (file.exists(file.path(out, "COMPLETED.txt"))) {
  cat("Completed selected-case scale development matches provenance.\n"); quit(status = 0L)
}
dir.create(file.path(out, "cases"), showWarnings = FALSE)
all_rows <- list()
for (i in seq_along(keys)) {
  path <- file.path(out, "cases", paste0(keys[i], ".rds"))
  if (file.exists(path)) {all_rows[[i]] <- readRDS(path)$row; next}
  case <- readRDS(file.path(pilot, "cases", paste0(keys[i], ".rds")))
  field <- readRDS(file.path(pilot, "fields", paste0(sub("_(moderate|heavy)$", "", keys[i]), ".rds")))
  model <- spatial_censored_setup_reference(case$lower, case$upper, field$covariance)
  center <- case$row$spatial_theta
  started <- proc.time()[[3L]]
  coarse_u <- spatial_censored_uniforms_reference(model, 16384L, 380010000L + i)
  coarse <- spatial_censored_fit_scale_reference(model, coarse_u,
    starts = rbind(c(center, 1), c(center / 5, .2)))
  if (coarse$status != "finite") {
    fine <- coarse; audit <- NULL
  } else {
    fine_u <- spatial_censored_uniforms_reference(model, 65536L, 390010000L + i)
    fine <- spatial_censored_fit_scale_reference(model, fine_u,
      starts = rbind(c(coarse$theta / coarse$multiplier, 1 / coarse$multiplier), c(center, 1)))
    audit <- NULL
    if (fine$status == "finite") {
      offsets <- rbind(c(-.02, 0), c(0, 0), c(.02, 0), c(0, -.02), c(0, .02))
      audit <- matrix(NA_real_, 8L, nrow(offsets))
      for (b in 1:8) {
        u <- spatial_censored_uniforms_reference(model, 32768L, 400010000L + i * 100L + b)
        audit[b, ] <- vapply(seq_len(nrow(offsets)), function(j) {
          mu <- fine$theta + offsets[j, 1L]
          multiplier <- fine$multiplier * exp(offsets[j, 2L])
          spatial_censored_scaled_loglik_reference(model, mu / multiplier, 1 / multiplier, u)
        }, 0.0)
      }
    }
  }
  equal <- tryCatch({
    fit <- withCallingHandlers(survival::survreg(
      survival::Surv(case$lower, case$upper, type = "interval2") ~ 1,
      dist = "gaussian"), warning = function(w) stop(conditionMessage(w)))
    if (any(!is.finite(coef(fit))) || !is.finite(fit$scale) || fit$scale <= 0)
      stop("Invalid independent estimate")
    list(status = "finite", theta = unname(coef(fit)), sigma = fit$scale)
  }, error = function(e) list(status = "unavailable", theta = NA_real_, sigma = NA_real_, message = conditionMessage(e)))
  row <- data.frame(key = keys[i], n_exact = model$no, status = fine$status,
    oracle_scale_theta = center, estimated_scale_theta = fine$theta,
    scale_multiplier = fine$multiplier,
    estimated_marginal_sigma = sqrt(field$covariance[1L, 1L]) * fine$multiplier,
    coarse_theta = coarse$theta, coarse_multiplier = coarse$multiplier,
    weighted_independent_theta = case$working$theta,
    weighted_independent_sigma = exp(case$working$log_scale),
    equal_independent_status = equal$status,
    equal_independent_theta = equal$theta, equal_independent_sigma = equal$sigma,
    elapsed_seconds = proc.time()[[3L]] - started)
  saveRDS(list(row = row, coarse = coarse, fine = fine, audit = audit, independent = equal), path)
  all_rows[[i]] <- row
  cat(i, "/ 4", keys[i], fine$status, "location", fine$theta,
      "scale multiplier", fine$multiplier, "elapsed", round(row$elapsed_seconds, 1), "seconds\n")
  flush.console()
}
stopifnot(identical(tools::md5sum(sources), signature$sources),
          identical(tools::md5sum(inputs), signature$inputs))
write.csv(do.call(rbind, all_rows), file.path(out, "selected_case_results.csv"), row.names = FALSE)
writeLines("Four selected-case estimated-scale development attempts recorded; not confirmation.",
  file.path(out, "COMPLETED.txt"))
