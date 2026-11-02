# Research reference: covariance shape fixed, one common scale estimated.
source("scripts/analysis/spatial_censored_likelihood_reference.R")

spatial_censored_scaled_loglik_reference <- function(model, location_scaled,
    precision_multiplier, uniforms = NULL, engine = "simulation") {
  if (length(location_scaled) != 1L || !is.finite(location_scaled) ||
      length(precision_multiplier) != 1L || !is.finite(precision_multiplier) ||
      precision_multiplier <= 0) stop("Invalid scaled location/precision")
  # If the covariance is S/kappa^2, delta=mu*kappa. This reparameterization
  # leaves the likelihood unchanged and avoids optimizing mu with enormous SD.
  kappa <- precision_multiplier
  scaled <- model
  scaled$constant <- model$constant + model$no * log(kappa)
  scaled$inv <- model$inv * kappa^2
  scaled$conditional <- model$conditional / kappa^2
  scaled$S <- model$S / kappa^2
  spatial_censored_loglik_reference(scaled, location_scaled / kappa, uniforms, engine)
}

spatial_censored_fit_scale_reference <- function(model, uniforms = NULL,
    starts = NULL, engine = "simulation", maxit = 100L) {
  if (!model$no) {
    all_left <- all(model$lower == -Inf)
    all_right <- all(model$upper == Inf)
    left_limits <- unique(model$upper[is.infinite(model$lower)])
    right_limits <- unique(model$lower[is.infinite(model$upper)])
    common_window <- all(is.infinite(model$lower) | is.infinite(model$upper)) &&
      length(left_limits) == 1L && length(right_limits) == 1L && left_limits < right_limits
    if (all_left || all_right || common_window)
      return(list(status = "no_finite_common_window_scale_mle",
        theta = NA_real_, multiplier = NA_real_, candidates = data.frame()))
    stop("No-exact scale classification implemented only for common left/right windows")
  }
  if (is.null(starts)) {
    known <- spatial_censored_fit_mu_reference(model, uniforms, engine)
    if (known$status != "finite") stop("Starting location unavailable")
    starts <- rbind(c(known$theta, 1), c(known$theta / 2, .5), c(known$theta * 2, 2))
  }
  starts <- as.matrix(starts)
  if (ncol(starts) != 2L || any(!is.finite(starts)) || any(starts[, 2L] <= 0))
    stop("Invalid optimization starts")
  candidates <- vector("list", nrow(starts))
  for (i in seq_len(nrow(starts))) {
    failures <- 0L
    objective <- function(par) {
      value <- tryCatch(spatial_censored_scaled_loglik_reference(model,
        par[1L], par[2L], uniforms, engine), error = function(e) {
          failures <<- failures + 1L; -Inf
        })
      -value
    }
    result <- tryCatch(optim(starts[i, ], objective, method = "L-BFGS-B",
      lower = c(-Inf, 1e-5), upper = c(Inf, 1e3),
      control = list(maxit = maxit, factr = 1e7, pgtol = 1e-6,
                     ndeps = c(1e-5, 1e-5))), error = function(e) e)
    if (inherits(result, "error")) {
      candidates[[i]] <- data.frame(start = i, status = "numerical_failure", delta = NA_real_,
        precision = NA_real_, theta = NA_real_, multiplier = NA_real_,
        loglik = NA_real_, convergence = NA_integer_, invalid_evaluations = failures,
        message = conditionMessage(result))
    } else {
      boundary <- result$par[2L] <= 1.01e-5 || result$par[2L] >= .99e3
      good <- result$convergence == 0L && !boundary && is.finite(result$value)
      candidates[[i]] <- data.frame(start = i, status = if (good) "finite" else "not_accepted",
        delta = result$par[1L], precision = result$par[2L],
        theta = result$par[1L] / result$par[2L], multiplier = 1 / result$par[2L],
        loglik = -result$value, convergence = result$convergence,
        invalid_evaluations = failures, message = if (is.null(result$message)) "" else result$message)
    }
  }
  table <- do.call(rbind, candidates)
  ok <- which(table$status == "finite")
  if (!length(ok)) return(list(status = "no_accepted_fit", theta = NA_real_,
    multiplier = NA_real_, candidates = table))
  best <- ok[which.max(table$loglik[ok])]
  list(status = "finite", theta = table$theta[best], multiplier = table$multiplier[best],
    loglik = table$loglik[best], candidates = table)
}
