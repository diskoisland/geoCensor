#!/usr/bin/env Rscript
# Gaussian numerical SCPC benchmark; no censored or real-USGS intervals.
source("scripts/analysis/scpc_envelope_reference.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Supply at most one fresh output directory")
out <- if (length(args)) args[1L] else "output/scpc_envelope_validation_2026-11-01"
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) {
  stop("Preserve existing results and use a fresh output directory")
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
sources <- c("scripts/analysis/scpc_gaussian_reference.R",
             "scripts/analysis/scpc_envelope_reference.R",
             "scripts/analysis/validate_scpc_envelope.R")
snapshot <- file.path(out, "source_snapshot")
dir.create(snapshot)
file.copy(sources, snapshot)
write.csv(data.frame(file = sources, md5 = unname(tools::md5sum(sources))),
          file.path(out, "source_manifest.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
designs <- list()
for (n in c(36L, 64L)) {
  set.seed(11021000L + n)
  side <- as.integer(sqrt(n))
  s <- as.matrix(expand.grid(x = seq(0, 1, length.out = side),
                            y = seq(0, 1, length.out = side)))
  s <- s + matrix(runif(2L * n, -0.015, 0.015), ncol = 2L)
  designs[[paste0("jittered_", n)]] <- s
  dense <- as.integer(3 * n / 4)
  designs[[paste0("clustered_", n)]] <- rbind(
    matrix(runif(2L * dense, 0, 0.25), ncol = 2L),
    matrix(runif(2L * (n - dense), 0.25, 1), ncol = 2L))
}
saveRDS(designs, file.path(out, "designs.rds"))
alpha <- 0.05
q_candidates <- 1:20
fits <- list(); candidates <- list(); calibrations <- list(); audits <- list()
for (design in names(designs)) {
  s <- designs[[design]]
  c0 <- scpc_decay_from_average_reference(s, 0.03)
  coarse <- scpc_select_reference(s, c0, q_candidates, grid_points = 31L, alpha = alpha)
  fine <- scpc_select_reference(s, c0, q_candidates, grid_points = 61L, alpha = alpha)
  audit <- scpc_envelope_audit_reference(fine$selected, grid_points = 241L)
  error <- max(abs(coarse$candidates$critical / fine$candidates$critical - 1))
  stable <- coarse$selected$basis$q == fine$selected$basis$q
  stopifnot(error < 5e-5, stable, !fine$upper_candidate_boundary,
            audit$max_rejection <= alpha + 1e-7, audit$tail_upper_bound <= alpha)
  b <- fine$selected$basis; cal <- fine$selected$calibration
  # Check coordinate scale/rotation invariance of the selected eigenspace.
  transformed <- 7 * (s %*% matrix(c(0, -1, 1, 0), 2)) + 100
  moved_c0 <- scpc_decay_from_average_reference(transformed, 0.03)
  moved <- scpc_basis_reference(transformed, moved_c0, b$q)
  stopifnot(abs(moved_c0 * 7 / c0 - 1) < 1e-8,
    max(abs(tcrossprod(b$W[, -1L, drop = FALSE]) -
            tcrossprod(moved$W[, -1L, drop = FALSE]))) < 1e-8)
  fits[[design]] <- fine$selected
  candidates[[design]] <- cbind(design = design, fine$candidates)
  calibrations[[design]] <- data.frame(design = design, n_sites = nrow(s),
    c0 = c0, benchmark_average_correlation = 0.03, selected_q = b$q,
    critical = cal$critical, raw_ratio_critical = cal$critical / sqrt(b$q),
    expected_iid_length = fine$selected$expected_length,
    q_search_upper_boundary = fine$upper_candidate_boundary,
    grid_relative_difference = error, selected_q_stable = stable,
    audit_max_rejection = audit$max_rejection, tail_rejection_bound = audit$tail_upper_bound,
    tail_tv = cal$tail_tv, worst_decay_over_c0 = cal$worst_decay / c0)
  audits[[design]] <- cbind(design = design, audit$grid)
  cat("Calibrated and audited", design, "q =", b$q, "critical =", cal$critical, "\n")
  flush.console()
}
write.csv(do.call(rbind, candidates), file.path(out, "component_candidates.csv"), row.names = FALSE)
write.csv(do.call(rbind, calibrations), file.path(out, "calibration_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, audits), file.path(out, "dense_envelope_audit.csv"), row.names = FALSE)
# Freeze before any development/confirmation response draws.
saveRDS(fits, file.path(out, "frozen_calibrations.rds"))
fits <- readRDS(file.path(out, "frozen_calibrations.rds"))

rows <- list(); row <- 0L
for (phase in c("development", "confirmation")) {
  draws <- if (phase == "development") 20000L else 100000L
  base_seed <- if (phase == "development") 11022000L else 11023000L
  for (id in seq_along(designs)) {
    design <- names(designs)[id]; s <- designs[[design]]
    fit <- fits[[design]]; b <- fit$basis; n <- b$n; q <- b$q
    # Full orthonormal basis for the iid mean-t comparator's analytic check.
    fullW <- qr.Q(qr(cbind(rep(1 / sqrt(n), n), diag(n)[, -1L])))
    for (j in seq_along(c(0.5, 1, 2, 10, Inf))) {
      multiple <- c(0.5, 1, 2, 10, Inf)[j]
      Sigma <- if (is.infinite(multiple)) diag(n) else exp(-multiple * b$c0 * b$distance)
      Omega <- crossprod(b$W, Sigma %*% b$W)
      values <- c(scpc_numerical_envelope = fit$calibration$critical,
                  projection_student = qt(0.975, q), iid_mean_t = qt(0.975, n - 1L))
      exact <- c(scpc_gaussian_rejection_reference(Omega, values[1L]),
                 scpc_gaussian_rejection_reference(Omega, values[2L]),
                 scpc_gaussian_rejection_reference(crossprod(fullW, Sigma %*% fullW), values[3L]))
      names(exact) <- names(values)
      hits <- setNames(numeric(3L), names(values))
      widths <- width_squares <- hits
      set.seed(base_seed + id * 100L + j)
      for (batch in seq_len(draws / 5000L)) {
        y <- matrix(rnorm(5000L * n), ncol = n) %*% chol(Sigma)
        projection <- y %*% b$W
        mean <- rowMeans(y)
        projection_se <- sqrt(rowSums(projection[, -1L, drop = FALSE]^2) / (q * n))
        mean_se <- sqrt(rowSums((y - mean)^2) / ((n - 1L) * n))
        ses <- cbind(projection_se, projection_se, mean_se)
        for (m in seq_along(values)) {
          hits[m] <- hits[m] + sum(abs(mean) > values[m] * ses[, m])
          w <- 2 * values[m] * ses[, m]
          widths[m] <- widths[m] + sum(w)
          width_squares[m] <- width_squares[m] + sum(w^2)
        }
      }
      for (method in names(values)) {
        p <- exact[method]; empirical <- hits[method] / draws
        mc_se <- sqrt(p * (1 - p) / draws)
        stopifnot(abs(empirical - p) < 5 * mc_se + 1e-4)
        if (method == "scpc_numerical_envelope" && multiple >= 1) stopifnot(p <= alpha + 1e-7)
        ci <- binom.test(draws - hits[method], draws)$conf.int
        width_mean <- widths[method] / draws
        width_se <- sqrt(max(0, width_squares[method] / draws - width_mean^2) / draws)
        if (method == "scpc_numerical_envelope" && is.infinite(multiple)) {
          stopifnot(abs(width_mean - fit$expected_length) < 5 * width_se + 1e-5)
        }
        row <- row + 1L
        rows[[row]] <- data.frame(phase = phase, design = design, n_sites = n,
          actual_decay_over_c0 = multiple, inside_assumed_envelope = multiple >= 1,
          method = method, q = if (method == "iid_mean_t") n - 1L else q,
          critical = values[method], draws = draws, analytic_coverage = 1 - p,
          empirical_coverage = 1 - empirical, mc_low = ci[1L], mc_high = ci[2L],
          mean_interval_width = width_mean, width_mc_standard_error = width_se)
      }
    }
    cat("Gaussian", phase, design, "complete\n"); flush.console()
  }
}
write.csv(do.call(rbind, rows), file.path(out, "gaussian_validation.csv"), row.names = FALSE)
writeLines(c("Gaussian numerical SCPC envelope study complete.",
  "Compact interval checked by grids and local searches; no global continuum certificate.",
  "Tail controlled by explicit Gaussian total-variation bound relative to independence.",
  "No censored-score, estimated-weight, preferential-sampling, or real-data interval claim."),
  file.path(out, "COMPLETED.txt"))
