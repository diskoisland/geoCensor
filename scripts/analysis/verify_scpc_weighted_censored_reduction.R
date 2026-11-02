suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
source("scripts/analysis/scpc_weighted_reference.R")
source("scripts/analysis/conley_survreg_reference.R")
set.seed(190010001L)
s <- matrix(runif(32), ncol = 2L)
n <- nrow(s); id <- rep(seq_len(n), each = 4L)
a <- cell_declustering_weights(data.frame(x = s[, 1L], y = s[, 2L]), "x", "y", 0.3, 10L)
weights <- a[id] / 4
logvalue <- log(50) + rnorm(length(id), sd = 0.6)
x <- data.frame(detected_val = exp(logvalue), left_cen = FALSE, left_limit = NA_real_,
                right_cen = FALSE, right_limit = NA_real_)
fit <- fit_censored_score_reference(x, weights)
theta <- sum(weights * logvalue)
site_mean <- as.vector(rowsum(logvalue, id, reorder = FALSE)) / 4
expected <- n * a * (site_mean - theta)
influence <- as.vector(n * rowsum(fit$scores, id, reorder = FALSE) %*% fit$bread[1L, ])
stopifnot(abs(fit$theta[1L] - theta) < 1e-8, max(abs(influence - expected)) < 1e-7)
scaled_x <- x; scaled_x$detected_val <- x$detected_val * 1000
scaled <- fit_censored_score_reference(scaled_x, weights * 123)
scaled_influence <- as.vector(n * rowsum(scaled$scores, id, reorder = FALSE) %*% scaled$bread[1L, ])
stopifnot(abs(scaled$theta[1L] - theta - log(1000)) < 1e-8,
          max(abs(scaled_influence - influence)) < 1e-7)
cat("Uncensored weighted survreg location/influences agree with exact weighted-mean algebra.\n",
    "Concentration-unit shifts and common weight scaling preserve influences.\n")
