# Compare with the previously checked implementation retained in its snapshot.
reference <- "output/scpc_weighted_gaussian_2026-11-01/source_snapshot/scpc_gaussian_reference.R"
if (!file.exists(reference)) stop("Original checked source snapshot is required")
old <- new.env(parent = globalenv())
sys.source(reference, envir = old)
source("scripts/analysis/scpc_gaussian_reference.R")
set.seed(200010001L)
error <- critical_error <- 0
for (q in c(1L, 2L, 3L, 5L, 10L, 20L)) for (i in 1:5) {
  z <- matrix(rnorm((q + 1L)^2), q + 1L)
  Omega <- crossprod(z) + diag(0.1, q + 1L)
  for (cv in c(1, 2, 4, 10)) {
    before <- old$scpc_gaussian_rejection_reference(Omega, cv)
    after <- scpc_gaussian_rejection_reference(Omega, cv)
    error <- max(error, abs(before - after))
  }
  critical_error <- max(critical_error, abs(old$scpc_gaussian_critical_reference(Omega) -
    scpc_gaussian_critical_reference(Omega)))
}
stopifnot(error < 1e-10, critical_error < 1e-7)
cat("Checked 120 probabilities and 30 critical values against the original source snapshot.\n",
    "Maximum probability difference:", error, "; critical-value difference:", critical_error, "\n")
