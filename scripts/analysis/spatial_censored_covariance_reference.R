# Research only: exponential spatial covariance plus independent record noise.
# Four estimated parameters; no prediction at unsampled sites or CI procedure.
source("scripts/analysis/spatial_censored_likelihood_reference.R")

spatial_covariance_matrix_reference <- function(coordinates, sigma, range,
    nugget_fraction) {
  xy <- as.matrix(coordinates)
  if (!is.numeric(xy) || ncol(xy) != 2L || !nrow(xy) || nrow(xy) > 96L ||
      any(!is.finite(xy))) stop("Provide 1 to 96 finite two-dimensional record coordinates")
  if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0 ||
      length(range) != 1L || !is.finite(range) || range <= 0 ||
      length(nugget_fraction) != 1L || !is.finite(nugget_fraction) ||
      nugget_fraction < 0 || nugget_fraction > 1) stop("Invalid covariance parameters")
  covariance <- sigma^2 * ((1 - nugget_fraction) *
    exp(-as.matrix(dist(xy)) / range) + diag(nugget_fraction, nrow(xy)))
  chol(covariance) # No silent ridge, including for coincident noiseless records.
  covariance
}

spatial_covariance_loglik_reference <- function(lower, upper, coordinates, mu,
    sigma, range, nugget_fraction, uniforms = NULL, engine = "simulation") {
  S <- spatial_covariance_matrix_reference(coordinates, sigma, range, nugget_fraction)
  model <- spatial_censored_setup_reference(lower, upper, S)
  value <- spatial_censored_loglik_reference(model, mu, uniforms, engine)
  if (!is.finite(value)) stop("Nonfinite covariance likelihood")
  value
}

spatial_covariance_select_reference <- function(candidates, lower_bounds,
    upper_bounds, boundary_fraction = .001) {
  # The highest finite candidate controls status, even if it is on a boundary
  # or did not converge. Never replace it with a worse interior/converged fit.
  finite <- which(is.finite(candidates$loglik))
  missing <- list(theta = NA_real_, sigma = NA_real_, range = NA_real_,
    nugget_fraction = NA_real_, loglik = NA_real_, candidate = NULL)
  if (!length(finite)) return(c(list(status = "no_finite_candidate"), missing))
  winner <- finite[which.max(candidates$loglik[finite])]
  row <- candidates[winner, , drop = FALSE]
  par <- unlist(row[c("mu", "log_sigma", "log_range", "logit_nugget")], use.names = FALSE)
  distance <- pmin(par - lower_bounds, upper_bounds - par) / (upper_bounds - lower_bounds)
  at_boundary <- any(distance <= boundary_fraction)
  status <- if (is.na(row$convergence) || row$convergence != 0L) "optimizer_unresolved" else
    if (at_boundary) "search_boundary" else "finite_interior"
  candidate <- list(theta = row$mu, sigma = exp(row$log_sigma), range = exp(row$log_range),
    nugget_fraction = plogis(row$logit_nugget), loglik = row$loglik,
    par = par, start = row$start, at_boundary = at_boundary)
  if (status != "finite_interior") {
    missing$candidate <- candidate
    return(c(list(status = status), missing))
  }
  list(status = status, theta = candidate$theta, sigma = candidate$sigma,
    range = candidate$range, nugget_fraction = candidate$nugget_fraction,
    loglik = candidate$loglik, candidate = candidate)
}

spatial_censored_fit_covariance_reference <- function(lower, upper, coordinates,
    uniforms = NULL, starts = NULL, maxit = 1200L, reltol = 1e-6,
    engine = "simulation") {
  xy <- as.matrix(coordinates)
  empty <- list(theta = NA_real_, sigma = NA_real_, range = NA_real_,
    nugget_fraction = NA_real_, loglik = NA_real_, candidate = NULL,
    candidates = data.frame(), evaluations = data.frame())
  if (!length(lower) && !length(upper) && nrow(xy) == 0L)
    return(c(list(status = "no_data"), empty))
  # Validate response bounds with a unit covariance, without supplying a
  # generating covariance or estimating anything from missing concentrations.
  model <- spatial_censored_setup_reference(lower, upper, diag(length(lower)))
  if (!is.numeric(xy) || nrow(xy) != model$n || ncol(xy) != 2L ||
      any(!is.finite(xy))) stop("Invalid record coordinates")
  if (!model$no) {
    all_left <- all(lower == -Inf); all_right <- all(upper == Inf)
    a <- unique(upper[is.infinite(lower)]); b <- unique(lower[is.infinite(upper)])
    common_window <- all(is.infinite(lower) | is.infinite(upper)) &&
      length(a) == 1L && length(b) == 1L && a < b
    if (all_left || all_right || common_window)
      return(c(list(status = "no_exact_scale_identification"), empty))
    stop("No-exact classification only implemented for common reporting windows")
  }
  endpoints <- c(lower[is.finite(lower)], upper[is.finite(upper)])
  center <- median(endpoints)
  initial_mu <- median(lower[model$exact])
  lower_bounds <- c(center - 8, log(.03), log(.01), qlogis(.02))
  upper_bounds <- c(center + 8, log(10), log(2), qlogis(.98))
  if (is.null(starts)) starts <- cbind(rep(initial_mu, 3L),
    log(c(.3, 1, 2)), log(c(.03, .2, .8)), qlogis(c(.05, .3, .8)))
  starts <- as.matrix(starts)
  if (!is.numeric(starts) || !nrow(starts) || ncol(starts) != 4L || any(!is.finite(starts)) ||
      any(sweep(starts, 2L, lower_bounds, "-") < 0) ||
      any(sweep(starts, 2L, upper_bounds, "-") > 0)) stop("Invalid covariance starts")
  records <- list(); candidates <- list()
  for (i in seq_len(nrow(starts))) {
    objective <- function(par) {
      outside <- any(par < lower_bounds | par > upper_bounds)
      message <- ""; value <- NA_real_
      if (!outside) value <- tryCatch(spatial_covariance_loglik_reference(lower, upper, xy,
        par[1L], exp(par[2L]), exp(par[3L]), plogis(par[4L]), uniforms, engine),
        error = function(e) {message <<- conditionMessage(e); NA_real_})
      records[[length(records) + 1L]] <<- data.frame(start = i,
        mu = par[1L], log_sigma = par[2L], log_range = par[3L], logit_nugget = par[4L],
        loglik = value, status = if (outside) "outside_search" else
          if (is.finite(value)) "finite" else "likelihood_unavailable", message = message)
      # Inf is an explicit rejected trial, never a finite likelihood floor.
      if (is.finite(value)) -value else Inf
    }
    fit <- tryCatch(optim(starts[i, ], objective, method = "Nelder-Mead",
      control = list(maxit = maxit, reltol = reltol)), error = function(e) e)
    if (inherits(fit, "error")) {
      candidates[[i]] <- data.frame(start = i, mu = NA_real_, log_sigma = NA_real_,
        log_range = NA_real_, logit_nugget = NA_real_, loglik = NA_real_,
        convergence = NA_integer_, message = conditionMessage(fit))
    } else candidates[[i]] <- data.frame(start = i, mu = fit$par[1L], log_sigma = fit$par[2L],
      log_range = fit$par[3L], logit_nugget = fit$par[4L], loglik = -fit$value,
      convergence = fit$convergence, message = if (is.null(fit$message)) "" else fit$message)
  }
  candidates <- do.call(rbind, candidates)
  result <- spatial_covariance_select_reference(candidates, lower_bounds, upper_bounds)
  c(result, list(candidates = candidates, evaluations = do.call(rbind, records),
    search_lower = lower_bounds, search_upper = upper_bounds,
    maxit = maxit, reltol = reltol))
}
