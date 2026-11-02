# Research diagnostic: profile the common covariance scale at a fixed location.
# Numerical search limits are reported, never treated as statistical constraints.
source("scripts/analysis/spatial_censored_scale_reference.R")

spatial_censored_profile_scale_reference <- function(model, mu, uniforms = NULL,
    engine = "simulation", log_multiplier_grid = seq(-4, 4, length.out = 17L),
    tolerance = 1e-5) {
  if (length(mu) != 1L || !is.finite(mu)) stop("Provide one finite location")
  if (!is.numeric(log_multiplier_grid) || length(log_multiplier_grid) < 3L ||
      any(!is.finite(log_multiplier_grid)) || any(diff(log_multiplier_grid) <= 0))
    stop("Provide an increasing finite scale-search grid")
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0)
    stop("Provide a positive finite tolerance")
  evaluations <- list()
  evaluate <- function(eta, stage) {
    message <- ""
    value <- tryCatch({
      multiplier <- exp(eta)
      value <- spatial_censored_scaled_loglik_reference(model,
        mu / multiplier, 1 / multiplier, uniforms, engine)
      if (!is.finite(value)) stop("Nonfinite profile likelihood")
      value
    }, error = function(e) {message <<- conditionMessage(e); NA_real_})
    evaluations[[length(evaluations) + 1L]] <<- data.frame(stage = stage,
      log_multiplier = eta, loglik = value,
      status = if (is.finite(value)) "finite" else "unavailable", message = message)
    value
  }
  grid_ll <- vapply(log_multiplier_grid, evaluate, 0.0, stage = "grid")
  finite <- which(is.finite(grid_ll))
  if (!length(finite)) return(list(status = "no_finite_grid_value", mu = mu,
    multiplier = NA_real_, loglik = NA_real_, evaluations = do.call(rbind, evaluations)))
  best <- finite[which.max(grid_ll[finite])]
  optimized <- NULL; optimizer_message <- ""
  boundary <- best %in% c(1L, length(log_multiplier_grid))
  if (!boundary) {
    optimized <- tryCatch(withCallingHandlers(optimize(function(eta) {
      value <- evaluate(eta, "refinement")
      if (is.finite(value)) value else -Inf
    }, interval = log_multiplier_grid[c(best - 1L, best + 1L)],
      maximum = TRUE, tol = tolerance), warning = function(w) {
        optimizer_message <<- paste(optimizer_message, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) {
        optimizer_message <<- conditionMessage(e); NULL
      })
  }
  table <- do.call(rbind, evaluations)
  available <- which(is.finite(table$loglik))
  winner <- available[which.max(table$loglik[available])]
  eta <- table$log_multiplier[winner]
  at_limit <- min(abs(eta - range(log_multiplier_grid))) <= tolerance
  status <- if (at_limit) "profile_search_boundary" else if (is.null(optimized) ||
    !is.finite(optimized$objective)) "finite_grid_only" else "finite"
  list(status = status, mu = mu, multiplier = exp(eta), loglik = table$loglik[winner],
    log_multiplier = eta, invalid_evaluations = sum(!is.finite(table$loglik)),
    optimizer_message = optimizer_message, evaluations = table)
}
