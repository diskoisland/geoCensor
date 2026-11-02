#!/usr/bin/env Rscript
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
source("scripts/analysis/scpc_weighted_reference.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Supply an optional fresh output directory")
out <- if (length(args)) args[1L] else "output/scpc_weighted_gaussian_2026-11-01"
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Use a fresh output directory")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
sources <- c("scripts/analysis/scpc_gaussian_reference.R", "scripts/analysis/scpc_envelope_reference.R",
             "scripts/analysis/scpc_weighted_reference.R", "scripts/analysis/validate_scpc_weighted_gaussian.R",
             "R/declustering.R")
dir.create(file.path(out, "source_snapshot"))
invisible(file.copy(sources, file.path(out, "source_snapshot")))
write.csv(data.frame(file = sources, md5 = unname(tools::md5sum(sources))),
          file.path(out, "source_manifest.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
designs <- list(); weights <- list(); fits <- list(); tables <- list(); audit_rows <- list(); grids <- list()
for (n in c(36L, 64L)) {
  set.seed(11031000L + n)
  s <- as.matrix(expand.grid(x = seq(0, 1, length.out = sqrt(n)),
                            y = seq(0, 1, length.out = sqrt(n))))
  s <- s + matrix(runif(2L * n, -0.015, 0.015), ncol = 2L)
  designs[[paste0("jittered_", n)]] <- s
  dense <- as.integer(3 * n / 4)
  designs[[paste0("clustered_", n)]] <- rbind(matrix(runif(2L * dense, 0, 0.25), ncol = 2L),
    matrix(runif(2L * (n - dense), 0.25, 1), ncol = 2L))
}
for (design in names(designs)) {
  s <- designs[[design]]; n <- nrow(s)
  a <- cell_declustering_weights(data.frame(x = s[, 1L], y = s[, 2L]), "x", "y", 0.25, 10L)
  # Four colocated records must aggregate to the same site weights.
  id <- rep(seq_len(n), each = 4L)
  obs <- cell_declustering_weights(data.frame(x = s[id, 1L], y = s[id, 2L]), "x", "y", 0.25, 10L)
  stopifnot(max(abs(as.vector(rowsum(obs, id)) - a)) < 1e-12)
  c0 <- scpc_decay_from_average_reference(s, 0.03)
  coarse <- scpc_weighted_select_reference(s, a, c0, 31L)
  fine <- scpc_weighted_select_reference(s, a, c0, 61L)
  audit <- scpc_weighted_audit_reference(fine$selected)
  relative <- max(abs(coarse$candidates$critical / fine$candidates$critical - 1))
  stopifnot(relative < 5e-5, coarse$selected$basis$q == fine$selected$basis$q,
            audit$max_rejection <= 0.05 + 1e-7, audit$iid_tail_bound <= 0.05)
  naive <- scpc_select_reference(s, c0, 1:20, 61L)$selected
  fits[[design]] <- list(weighted = fine$selected, geometry_only = naive)
  weights[[design]] <- a
  tables[[design]] <- cbind(design = design, fine$candidates)
  audit_rows[[design]] <- data.frame(design = design, n_sites = n, c0 = c0,
    cell_size = 0.25, origins = 100L, kish_weight_n = 1 / sum(a^2), min_weight = min(a), max_weight = max(a),
    selected_q = fine$selected$basis$q, critical = fine$selected$calibration$critical,
    iid_point_critical = fine$selected$calibration$iid_critical,
    expected_iid_length = fine$selected$calibration$expected_iid_length,
    geometry_only_q = naive$basis$q, geometry_only_critical = naive$calibration$critical,
    grid_relative_difference = relative, max_audited_rejection = audit$max_rejection,
    iid_tail_bound = audit$iid_tail_bound)
  grids[[design]] <- cbind(design = design, audit$grid)
  cat("Weighted calibration audited:", design, "q =", fine$selected$basis$q,
      "Kish weight n =", 1 / sum(a^2), "\n"); flush.console()
}
saveRDS(list(designs = designs, weights = weights, fits = fits), file.path(out, "frozen_calibrations.rds"))
write.csv(do.call(rbind, tables), file.path(out, "component_candidates.csv"), row.names = FALSE)
write.csv(do.call(rbind, audit_rows), file.path(out, "calibration_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, grids), file.path(out, "dense_envelope_audit.csv"), row.names = FALSE)
rows <- list(); k <- 0L
for (phase in c("development", "confirmation")) {
  draws <- if (phase == "development") 20000L else 100000L
  seed <- if (phase == "development") 11032000L else 11033000L
  for (d in seq_along(designs)) {
    design <- names(designs)[d]; fit <- fits[[design]]; b <- fit$weighted$basis
    a <- b$weights; n <- b$n; naive <- fit$geometry_only
    # Alternative transform retains weighted estimation but ignores weights
    # when choosing/calibrating the spatial components.
    altW <- cbind(sqrt(n) * a, t(b$B) %*% naive$basis$W[, -1L, drop = FALSE])
    fullV <- qr.Q(qr(cbind(rep(1 / sqrt(n), n), diag(n)[, -1L])))[, -1L, drop = FALSE]
    clusterW <- cbind(sqrt(n) * a, t(b$B) %*% fullV)
    for (j in seq_along(c(0.5, 1, 2, Inf))) {
      multiple <- c(0.5, 1, 2, Inf)[j]
      S <- if (is.infinite(multiple)) diag(n) else exp(-multiple * b$c0 * b$distance)
      values <- c(weight_aware = fit$weighted$calibration$critical,
                  geometry_only = naive$calibration$critical, site_cluster_t = qt(0.975, n - 1L))
      exact <- c(scpc_gaussian_rejection_reference(crossprod(b$W, S %*% b$W), values[1L]),
                 scpc_gaussian_rejection_reference(crossprod(altW, S %*% altW), values[2L]),
                 scpc_gaussian_rejection_reference(crossprod(clusterW, S %*% clusterW), values[3L]))
      names(exact) <- names(values)
      hits <- width_sum <- width_square <- setNames(numeric(3L), names(values))
      set.seed(seed + d * 100L + j)
      for (batch in seq_len(draws / 5000L)) {
        y <- matrix(rnorm(5000L * n), ncol = n) %*% chol(S)
        theta <- as.vector(y %*% a)
        psi <- (y - theta) %*% diag(n * a, n)
        stopifnot(max(abs(rowSums(psi))) < 1e-9)
        ses <- cbind(sqrt(rowSums((psi %*% b$V)^2) / (b$q * n)),
          sqrt(rowSums((psi %*% naive$basis$W[, -1L, drop = FALSE])^2) / (naive$basis$q * n)),
          sqrt(rowSums(psi^2) / (n * (n - 1L))))
        for (m in seq_along(values)) {
          hits[m] <- hits[m] + sum(abs(theta) > values[m] * ses[, m])
          widths <- 2 * values[m] * ses[, m]
          width_sum[m] <- width_sum[m] + sum(widths)
          width_square[m] <- width_square[m] + sum(widths^2)
        }
      }
      for (method in names(values)) {
        p <- exact[method]; observed <- hits[method] / draws
        stopifnot(abs(observed - p) < 5 * sqrt(p * (1 - p) / draws) + 1e-4)
        if (method == "weight_aware" && multiple >= 1) stopifnot(p <= 0.05 + 1e-7)
        width <- width_sum[method] / draws
        width_se <- sqrt(max(0, width_square[method] / draws - width^2) / draws)
        if (method == "weight_aware" && is.infinite(multiple)) stopifnot(
          abs(width - fit$weighted$calibration$expected_iid_length) < 5 * width_se + 1e-5)
        ci <- binom.test(draws - hits[method], draws)$conf.int
        k <- k + 1L
        rows[[k]] <- data.frame(phase = phase, design = design, n_sites = n,
          actual_decay_over_c0 = multiple, inside_assumed_envelope = multiple >= 1,
          method = method, critical = values[method], draws = draws,
          analytic_coverage = 1 - p, coverage = 1 - observed, mc_low = ci[1L], mc_high = ci[2L],
          mean_width = width, width_mc_standard_error = width_se)
      }
    }
    cat("Weighted Gaussian", phase, design, "complete\n"); flush.console()
  }
}
write.csv(do.call(rbind, rows), file.path(out, "gaussian_validation.csv"), row.names = FALSE)
writeLines(c("Conditional weighted Gaussian validation complete.",
  "Weights are derived from each fixed location pattern, not concentrations.",
  "Compact envelope is numerically searched, not globally certified.",
  "No claim for random-location weight uncertainty or preferential sampling."), file.path(out, "COMPLETED.txt"))
