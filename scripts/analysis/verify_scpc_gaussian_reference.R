# Run from the repository root. This verifies fixed-q formulas only; no
# nonlinear/censored inference, estimated weights, or production CIs are tested.
source("scripts/analysis/scpc_gaussian_reference.R")
out <- commandArgs(trailingOnly = TRUE)
if (length(out) > 1L) stop("Supply at most one output directory")
out <- if (length(out)) out else "output/scpc_gaussian_formula_checks_2026-11-01"
if (dir.exists(out) && length(list.files(out, all.files = FALSE))) stop("Use a fresh output directory")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

iid <- do.call(rbind, lapply(c(1L, 2L, 3L, 5L, 10L), function(q) {
  critical <- scpc_gaussian_critical_reference(diag(q + 1L))
  error <- abs(critical - qt(0.975, q))
  stopifnot(error < 1e-7)
  data.frame(q = q, reference_critical = critical, student_critical = qt(0.975, q),
             archive_raw_ratio_critical = critical / sqrt(q), absolute_error = error)
}))
write.csv(iid, file.path(out, "iid_checks.csv"), row.names = FALSE)

set.seed(110126)
regular <- as.matrix(expand.grid(x = seq(0, 1, length.out = 5L),
                                y = seq(0, 1, length.out = 5L)))
regular <- regular + matrix(runif(50L, -0.02, 0.02), ncol = 2L)
clustered <- rbind(matrix(runif(40L, 0, 0.25), ncol = 2L),
                   matrix(runif(10L, 0.25, 1), ncol = 2L))
settings <- list(jittered_grid = regular, exogenous_cluster = clustered)
results <- list()
row <- 0L
draws <- 100000L

for (design in names(settings)) {
  s <- settings[[design]]
  b <- scpc_basis_reference(s, c0 = 5, q = 3L)
  stopifnot(max(abs(crossprod(b$W) - diag(4L))) < 1e-10)
  moved <- scpc_basis_reference((s %*% matrix(c(0, -1, 1, 0), 2)) + 10,
                                c0 = 5, q = 3L)
  stopifnot(max(abs(tcrossprod(b$W[, -1L]) - tcrossprod(moved$W[, -1L]))) < 1e-9)
  for (decay in c(5, 20, Inf)) {
    Sigma <- if (is.infinite(decay)) diag(nrow(s)) else exp(-decay * b$distance)
    Omega <- crossprod(b$W, Sigma %*% b$W)
    critical <- scpc_gaussian_critical_reference(Omega)
    exact <- scpc_gaussian_rejection_reference(Omega, critical)
    stopifnot(abs(exact - 0.05) < 1e-8)
    # Numerical equality with the covariance eigenvalues used in the archive:
    # there the statistic is X0/sqrt(sum Xj^2), with critical/sqrt(q).
    raw_D <- diag(c(1, rep(-(critical / sqrt(b$q))^2, b$q)))
    archive_e <- sort(Re(eigen(Omega %*% raw_D, only.values = TRUE)$values))
    L <- t(chol(Omega))
    symmetric_e <- sort(eigen(crossprod(L, raw_D %*% L),
                             symmetric = TRUE, only.values = TRUE)$values)
    stopifnot(max(abs(archive_e - symmetric_e)) < 1e-9)

    hits <- 0L
    numerator_sum <- 0
    numerator_squares <- 0
    for (batch in seq_len(draws / 5000L)) {
      # Simulate the full observation vector, independently of the projected
      # covariance/eigenvalue integral used to calibrate the test.
      y <- matrix(rnorm(5000L * nrow(s)), ncol = nrow(s)) %*% chol(Sigma)
      projections <- y %*% b$W
      T <- sqrt(b$q) * projections[, 1L] / sqrt(rowSums(projections[, -1L]^2))
      se <- sqrt(rowSums(projections[, -1L]^2) / (b$q * b$n))
      stopifnot(max(abs(T - rowMeans(y) / se)) < 1e-8)
      hits <- hits + sum(abs(T) > critical)
      numerator_sum <- numerator_sum + sum(projections[, 1L])
      numerator_squares <- numerator_squares + sum(projections[, 1L]^2)
    }
    empirical <- hits / draws
    mc_se <- sqrt(exact * (1 - exact) / draws)
    stopifnot(abs(empirical - exact) < 4 * mc_se + 1e-4)
    row <- row + 1L
    results[[row]] <- data.frame(design = design, n_sites = b$n, q = b$q,
      actual_decay = decay, critical = critical,
      raw_ratio_critical = critical / sqrt(b$q), analytic_rejection = exact,
      simulated_rejection = empirical, draws = draws, mc_standard_error = mc_se)
  }
}
write.csv(do.call(rbind, results), file.path(out, "gaussian_checks.csv"), row.names = FALSE)
scripts <- c("scripts/analysis/scpc_gaussian_reference.R",
             "scripts/analysis/verify_scpc_gaussian_reference.R")
write.csv(data.frame(file = scripts, md5 = unname(tools::md5sum(scripts))),
          file.path(out, "source_manifest.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("Verified iid t quantiles, normalization, quadratic-form eigenvalues,\n",
    "coordinate invariance, and six independent 100,000-draw Gaussian checks.\n")
