#!/usr/bin/env Rscript
# Numerical development on a previously generated case; not confirmation.
source("scripts/analysis/spatial_censored_likelihood_reference.R")
pilot <- "output/spatial_likelihood_reflected_2026-11-01"
out <- "output/spatial_likelihood_precision_2026-11-01"
if (dir.exists(out) && length(list.files(out))) stop("Use a fresh precision directory")
dir.create(out, recursive = TRUE)
x <- read.csv(file.path(pilot, "pilot_results.csv"), stringsAsFactors = FALSE)
row <- x[which.max(x$audit_loglik_mc_se), ]
case <- readRDS(file.path(pilot, "cases", paste0(row$key, ".rds")))
field <- readRDS(file.path(pilot, "fields", paste0(sub("_(moderate|heavy)$", "", row$key), ".rds")))
writeLines(c("Numerical precision refinement — November 1, 2026",
  "Selected the completed reflected-pilot case with the largest audit loglik Monte Carlo SE.",
  paste("Case:", row$key),
  "Refit with 262144 and 1048576 uniforms, seeds 360010001 and 360010002.",
  "Audit the final location and offsets -0.01/0/+0.01 using 16 independent 65536-draw batches.",
  "Audit seed 370010000+batch. All covariance remains known; no new response or confidence interval.",
  "This is post-hoc numerical development, not independent confirmation or a global integration certificate."),
  file.path(out, "PRECISION_PROTOCOL.txt"))
model <- spatial_censored_setup_reference(case$lower, case$upper, field$covariance)
fits <- list()
for (i in 1:2) {
  draws <- c(262144L, 1048576L)[i]
  u <- spatial_censored_uniforms_reference(model, draws, 360010000L + i)
  fit <- spatial_censored_fit_mu_reference(model, u)
  stopifnot(fit$status == "finite")
  fits[[i]] <- data.frame(draws = draws, seed = 360010000L + i,
    theta = fit$theta, loglik = fit$loglik)
  cat("Precision refit", draws, "draws, location", fit$theta, "\n"); flush.console()
  rm(u); invisible(gc())
}
fits <- do.call(rbind, fits)
center <- fits$theta[2L]
audit <- matrix(NA_real_, 16L, 3L)
for (b in seq_len(nrow(audit))) {
  u <- spatial_censored_uniforms_reference(model, 65536L, 370010000L + b)
  audit[b, ] <- vapply(center + c(-.01, 0, .01), function(mu)
    spatial_censored_loglik_reference(model, mu, u), 0.0)
}
delta <- sweep(audit[, c(1L, 3L), drop = FALSE], 1L, audit[, 2L])
summary <- data.frame(key = row$key, pilot_theta = row$spatial_theta,
  theta_262144 = fits$theta[1L], theta_1048576 = center,
  change_from_pilot = center - row$spatial_theta,
  audit_loglik_mc_se = sd(audit[, 2L]) / sqrt(nrow(audit)),
  lower_delta = mean(delta[, 1L]), lower_delta_mc_se = sd(delta[, 1L]) / sqrt(nrow(delta)),
  upper_delta = mean(delta[, 2L]), upper_delta_mc_se = sd(delta[, 2L]) / sqrt(nrow(delta)))
write.csv(fits, file.path(out, "refits.csv"), row.names = FALSE)
write.csv(audit, file.path(out, "audit_batches.csv"), row.names = FALSE)
write.csv(summary, file.path(out, "precision_summary.csv"), row.names = FALSE)
file.copy(c("scripts/analysis/refine_spatial_likelihood_precision.R",
  "scripts/analysis/spatial_censored_likelihood_reference.R"), out)
saveRDS(list(case = case, field = field, selection = "maximum prior integration MC SE"),
        file.path(out, "selected_input.rds"))
writeLines("Precision development case completed; no statistical interval claim.", file.path(out, "COMPLETED.txt"))
cat("Independent numerical refinement/audit complete.\n")
