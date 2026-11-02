#!/usr/bin/env Rscript
# Post hoc endpoint diagnostic, not a replacement estimator or model test.
source("scripts/analysis/spatial_censored_covariance_reference.R")
pilot <- "output/spatial_covariance_fresh_2026-11-01"
out <- "output/spatial_covariance_precision_2026-11-01"
stopifnot(file.exists(file.path(out, "COMPLETED.txt")))
if (file.exists(file.path(out, "independent_endpoint.csv"))) stop("Endpoint diagnostic already exists")
manifest <- readRDS(file.path(out, "refinement_manifest.rds"))
keys <- manifest$config$keys
rows <- list(); files <- c("scripts/analysis/diagnose_covariance_independent_endpoint.R",
  "scripts/analysis/spatial_censored_covariance_reference.R",
  "scripts/analysis/spatial_censored_likelihood_reference.R")
for (i in seq_along(keys)) {
  paths <- c(file.path(pilot, "cases", paste0(keys[i], ".rds")),
    file.path(pilot, "fields", paste0(sub("_(moderate|heavy)$", "", keys[i]), ".rds")),
    file.path(out, "cases", paste0(keys[i], ".rds")))
  files <- c(files, paths)
  original <- readRDS(paths[1L]); field <- readRDS(paths[2L]); refined <- readRDS(paths[3L])
  candidate <- refined$fit$candidate
  independent <- original$independent
  if (is.null(candidate) || independent$status != "finite") stop("Endpoint/candidate unavailable")
  exact <- original$lower == original$upper
  left <- is.infinite(original$lower); right <- is.infinite(original$upper)
  stopifnot(all(exact | left | right))
  # At nugget fraction one the observations are independent and range drops
  # out exactly. No range value is supplied or estimated for this endpoint.
  ll <- sum(dnorm(original$lower[exact], independent$theta, independent$sigma, log = TRUE)) +
    sum(pnorm(original$upper[left], independent$theta, independent$sigma, log.p = TRUE)) +
    sum(pnorm(original$lower[right], independent$theta, independent$sigma, lower.tail = FALSE, log.p = TRUE))
  unit <- spatial_censored_setup_reference(original$lower, original$upper, diag(length(original$lower)))
  deltas <- numeric(8L)
  for (b in seq_len(8L)) {
    u <- spatial_censored_uniforms_reference(unit, 32768L, 790010000L + i * 100L + b)
    value <- spatial_covariance_loglik_reference(original$lower, original$upper, field$coordinates,
      candidate$theta, candidate$sigma, candidate$range, candidate$nugget_fraction, u)
    deltas[b] <- ll - value
  }
  rows[[i]] <- data.frame(key = keys[i], refined_status = refined$row$refined_status,
    endpoint_theta = independent$theta, endpoint_sigma = independent$sigma,
    endpoint_nugget_fraction = 1, endpoint_range = NA_real_,
    endpoint_analytic_loglik = ll, mean_endpoint_minus_candidate_loglik = mean(deltas),
    paired_numerical_mc_se = sd(deltas) / sqrt(length(deltas)))
}
hash <- tools::md5sum(files)
write.csv(data.frame(file = files, md5 = unname(hash)), file.path(out, "endpoint_source_input_manifest.csv"), row.names = FALSE)
write.csv(do.call(rbind, rows), file.path(out, "independent_endpoint.csv"), row.names = FALSE)
invisible(file.copy("scripts/analysis/diagnose_covariance_independent_endpoint.R", file.path(out, "source_snapshot")))
stopifnot(identical(tools::md5sum(files), hash))
cat("Independent endpoint checked; its range remains NA and primary outcomes are unchanged.\n")
