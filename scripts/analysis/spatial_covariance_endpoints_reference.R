# Research revision: explicit independent and site-only limits. No prediction/CI.
source("scripts/analysis/spatial_censored_covariance_reference.R")

spatial_endpoint_setup_reference <- function(lower, upper, coordinates) {
  xy <- as.matrix(coordinates)
  if (!length(lower) && !length(upper) && nrow(xy) == 0L)
    return(list(status = "no_data"))
  unit <- spatial_censored_setup_reference(lower, upper, diag(length(lower)))
  if (!is.numeric(xy) || nrow(xy) != unit$n || ncol(xy) != 2L || any(!is.finite(xy)))
    stop("Invalid observed coordinates")
  xy[xy == 0] <- 0
  key <- paste(sprintf("%.17g", xy[, 1L]), sprintf("%.17g", xy[, 2L]), sep = ":")
  groups <- unname(split(seq_along(lower), factor(key, levels = unique(key))))
  if (any(lengths(groups) > 2L)) stop("Site-only reference supports at most two records per exact coordinate")
  center <- median(c(lower[is.finite(lower)], upper[is.finite(upper)]))
  status <- "usable"
  if (!unit$no) {
    a <- unique(upper[is.infinite(lower)]); b <- unique(lower[is.infinite(upper)])
    common <- all(is.infinite(lower) | is.infinite(upper)) &&
      length(a) == 1L && length(b) == 1L && a <= b
    if (all(lower == -Inf) || all(upper == Inf) || common) status <- "no_exact_scale_identification" else
      stop("No-exact classification only implemented for common reporting windows")
  }
  list(status = status, lower = lower, upper = upper, coordinates = xy, groups = groups,
    unit = unit, center = center, initial_mu = if (unit$no) median(lower[unit$exact]) else NA_real_)
}

spatial_endpoint_univariate_logmass_reference <- function(lower, upper, mu, sigma) {
  if (lower == upper) return(dnorm(lower, mu, sigma, log = TRUE))
  a <- (lower - mu) / sigma; b <- (upper - mu) / sigma
  if (a > 0) {
    x <- pnorm(a, lower.tail = FALSE, log.p = TRUE)
    y <- pnorm(b, lower.tail = FALSE, log.p = TRUE)
  } else {
    x <- pnorm(b, log.p = TRUE); y <- pnorm(a, log.p = TRUE)
  }
  x + log(-expm1(y - x))
}

spatial_endpoint_independent_loglik_reference <- function(setup, mu, sigma) {
  if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0) stop("Invalid location/SD")
  values <- vapply(seq_along(setup$lower), function(i)
    spatial_endpoint_univariate_logmass_reference(setup$lower[i], setup$upper[i], mu, sigma), 0.0)
  value <- sum(values)
  if (!is.finite(value)) stop("Independent endpoint likelihood unavailable")
  value
}

spatial_endpoint_site_loglik_reference <- function(setup, mu, sigma, nugget_fraction,
    engine = c("tvpack", "miwa"), miwa_steps = 4096L) {
  engine <- match.arg(engine)
  if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0 ||
      !is.finite(nugget_fraction) || nugget_fraction <= 0 || nugget_fraction > 1)
    stop("Invalid site-only parameters; zero-noise coincident records are singular")
  if (nugget_fraction == 1) return(spatial_endpoint_independent_loglik_reference(setup, mu, sigma))
  rho <- 1 - nugget_fraction; variance_factor <- 1 - rho^2
  cache <- new.env(parent = emptyenv()); total <- 0
  for (index in setup$groups) {
    lower <- setup$lower[index]; upper <- setup$upper[index]; exact <- lower == upper
    if (length(index) == 1L) {
      value <- spatial_endpoint_univariate_logmass_reference(lower, upper, mu, sigma)
    } else if (all(exact)) {
      z <- (lower - mu) / sigma
      value <- -log(2 * pi) - 2 * log(sigma) - .5 * log(variance_factor) -
        (z[1L]^2 - 2 * rho * z[1L] * z[2L] + z[2L]^2) / (2 * variance_factor)
    } else if (sum(exact) == 1L) {
      obs <- which(exact); cen <- which(!exact)
      conditional_mean <- mu + rho * (lower[obs] - mu)
      value <- dnorm(lower[obs], mu, sigma, log = TRUE) +
        spatial_endpoint_univariate_logmass_reference(lower[cen], upper[cen],
          conditional_mean, sigma * sqrt(variance_factor))
    } else {
      # Many sites share the same two reporting intervals; reuse the probability
      # only within this evaluation, without caching across parameter values.
      key <- paste(sprintf("%.17g", c(lower, upper)), collapse = ":")
      if (exists(key, envir = cache, inherits = FALSE)) value <- get(key, envir = cache) else {
        if (engine == "miwa") {
          covariance <- sigma^2 * matrix(c(1, rho, rho, 1), 2L, 2L)
          model <- spatial_censored_setup_reference(lower, upper, covariance)
          value <- spatial_censored_loglik_reference(model, mu, engine = "miwa", miwa_steps = miwa_steps)
        } else {
          lo <- (lower - mu) / sigma; up <- (upper - mu) / sigma
          sign <- ifelse(is.infinite(up) & is.finite(lo), -1, 1)
          lo_reflected <- ifelse(sign > 0, lo, -up)
          up_reflected <- ifelse(sign > 0, up, -lo)
          if (any(is.finite(lo_reflected))) stop("TVPACK site branch requires one-sided reporting bounds")
          probability <- mvtnorm::pmvnorm(lower = rep(-Inf, 2L), upper = up_reflected,
            corr = matrix(c(1, rho, rho, 1), 2L, 2L) * outer(sign, sign),
            algorithm = mvtnorm::TVPACK(abseps = 1e-12))
          if (!is.finite(probability) || probability <= 0) stop("TVPACK pair probability unavailable")
          value <- log(as.numeric(probability))
        }
        assign(key, value, envir = cache)
      }
    }
    if (!is.finite(value)) stop("Site-only block likelihood unavailable")
    total <- total + value
  }
  total
}

spatial_endpoint_optimize_reference <- function(objective, starts, lower, upper,
    maxit, reltol = 1e-7) {
  starts <- as.matrix(starts)
  evaluations <- list(); candidates <- list()
  for (i in seq_len(nrow(starts))) {
    f <- function(par) {
      outside <- any(!is.finite(par)) || any(par < lower | par > upper)
      message <- ""; value <- NA_real_
      if (!outside) value <- tryCatch(objective(par), error = function(e) {
        message <<- conditionMessage(e); NA_real_
      })
      evaluations[[length(evaluations) + 1L]] <<- data.frame(start = i,
        parameter = paste(sprintf("%.17g", par), collapse = ";"), loglik = value,
        status = if (outside) "outside_search" else if (is.finite(value)) "finite" else "unavailable",
        message = message)
      if (is.finite(value)) -value else Inf
    }
    fit <- tryCatch(optim(starts[i, ], f, method = "Nelder-Mead",
      control = list(maxit = maxit, reltol = reltol)), error = function(e) e)
    if (inherits(fit, "error")) candidates[[i]] <- list(start = i, par = rep(NA_real_, ncol(starts)),
      loglik = NA_real_, convergence = NA_integer_, message = conditionMessage(fit)) else
      candidates[[i]] <- list(start = i, par = fit$par, loglik = -fit$value,
        convergence = fit$convergence, message = if (is.null(fit$message)) "" else fit$message)
  }
  finite <- which(vapply(candidates, function(x) is.finite(x$loglik), TRUE))
  if (!length(finite)) return(list(status = "no_finite_candidate", candidate = NULL,
    candidates = candidates, evaluations = do.call(rbind, evaluations)))
  best <- candidates[[finite[which.max(vapply(candidates[finite], `[[`, 0.0, "loglik"))]]]
  distance <- pmin(best$par - lower, upper - best$par) / (upper - lower)
  status <- if (best$convergence != 0L) "optimizer_unresolved" else if (any(distance <= .001))
    "artificial_search_boundary" else "finite"
  list(status = status, candidate = best, candidates = candidates,
    evaluations = do.call(rbind, evaluations), lower = lower, upper = upper)
}

spatial_endpoint_fit_independent_reference <- function(setup) {
  fit <- tryCatch(withCallingHandlers(survival::survreg(
    survival::Surv(setup$lower, setup$upper, type = "interval2") ~ 1, dist = "gaussian"),
    warning = function(w) stop(conditionMessage(w))), error = function(e) e)
  if (inherits(fit, "error")) return(list(status = "optimizer_unresolved", candidate = NULL,
    message = conditionMessage(fit)))
  mu <- unname(coef(fit)); sigma <- fit$scale
  residual <- residuals(fit, type = "matrix", weighted = FALSE)
  info <- -matrix(c(sum(residual[, "ddg"]), sum(residual[, "dsg"]),
    sum(residual[, "dsg"]), sum(residual[, "dds"])), 2L, 2L)
  if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0 || any(!is.finite(info)) ||
      min(eigen(info, symmetric = TRUE, only.values = TRUE)$values) <= 0 ||
      max(abs(colSums(residual[, c("dg", "ds"), drop = FALSE]))) / length(setup$lower) > 1e-5)
    return(list(status = "optimizer_unresolved", candidate = NULL))
  lower <- c(setup$center - 8, log(.03)); upper <- c(setup$center + 8, log(10))
  par <- c(mu, log(sigma)); distance <- pmin(par - lower, upper - par) / (upper - lower)
  status <- if (any(distance <= .001)) "artificial_search_boundary" else "finite"
  list(status = status, candidate = list(par = par, mu = mu, sigma = sigma, range = NA_real_,
    nugget = 1, loglik = spatial_endpoint_independent_loglik_reference(setup, mu, sigma)))
}

spatial_endpoint_compare_reference <- function(branches, audit, deterministic_error,
    numerical_floor = 1e-6) {
  empty <- list(structure = NA_character_, theta = NA_real_, sigma = NA_real_,
    range = NA_real_, nugget_fraction = NA_real_)
  if (any(vapply(branches, function(x) is.null(x$candidate), TRUE)))
    return(c(list(status = "branch_search_unresolved", candidate_structure = NA_character_), empty))
  if (any(!is.finite(audit))) return(c(list(status = "audit_unavailable",
    candidate_structure = NA_character_), empty))
  mean_ll <- colMeans(audit); winner <- which.max(mean_ll)
  others <- setdiff(seq_along(branches), winner)
  gap <- sweep(audit[, others, drop = FALSE], 1L, audit[, winner], function(other, best) best - other)
  se <- apply(gap, 2L, sd) / sqrt(nrow(audit))
  margin <- 3 * se + 5 * (deterministic_error[winner] + deterministic_error[others]) + numerical_floor
  comparison <- data.frame(winner = names(branches)[winner], other = names(branches)[others],
    mean_gap = colMeans(gap), paired_mc_se = se, numerical_margin = margin)
  if (any(comparison$mean_gap <= comparison$numerical_margin)) status <- "comparison_unresolved" else
    if (branches[[winner]]$status != "finite") status <- "winning_branch_unresolved" else status <- "finite_numerical_winner"
  result <- c(list(status = status, candidate_structure = names(branches)[winner], comparison = comparison), empty)
  if (status == "finite_numerical_winner") {
    best <- branches[[winner]]$candidate
    result$structure <- names(branches)[winner]; result$theta <- best$mu; result$sigma <- best$sigma
    result$range <- if (result$structure == "spatial") best$range else NA_real_
    result$nugget_fraction <- best$nugget
  }
  result
}

spatial_censored_fit_endpoints_reference <- function(lower, upper, coordinates,
    coarse_seed, fine_seed, audit_seed, coarse_draws = 2048L, fine_draws = 8192L,
    audit_draws = 32768L, audit_batches = 8L) {
  setup <- spatial_endpoint_setup_reference(lower, upper, coordinates)
  empty <- list(structure = NA_character_, theta = NA_real_, sigma = NA_real_,
    range = NA_real_, nugget_fraction = NA_real_, branches = list())
  if (setup$status != "usable") return(c(list(status = setup$status), empty))
  center <- setup$center; initial <- setup$initial_mu
  independent <- spatial_endpoint_fit_independent_reference(setup)
  site_objective <- function(par) spatial_endpoint_site_loglik_reference(setup,
    par[1L], exp(par[2L]), plogis(par[3L]))
  site <- spatial_endpoint_optimize_reference(site_objective,
    rbind(c(initial, log(.5), qlogis(.2)), c(initial, log(1), qlogis(.8))),
    c(center - 8, log(.03), qlogis(.0001)), c(center + 8, log(10), qlogis(.9999)),
    maxit = 1800L, reltol = 1e-8)
  if (!is.null(site$candidate)) {
    par <- site$candidate$par
    site$candidate$mu <- par[1L]; site$candidate$sigma <- exp(par[2L])
    site$candidate$range <- NA_real_; site$candidate$nugget <- plogis(par[3L])
  }
  # Without repeats the site-only nugget drops out; do not invent an estimate.
  if (all(lengths(setup$groups) == 1L)) site <- list(status = "unidentified_site_nugget", candidate = NULL)
  lo <- c(center - 8, log(.03), log(.001), qlogis(.0001))
  hi <- c(center + 8, log(10), log(10), qlogis(.9999))
  coarse_u <- spatial_censored_uniforms_reference(setup$unit, coarse_draws, coarse_seed)
  fine_u <- spatial_censored_uniforms_reference(setup$unit, fine_draws, fine_seed)
  objective <- function(par, u) spatial_covariance_loglik_reference(lower, upper, coordinates,
    par[1L], exp(par[2L]), exp(par[3L]), plogis(par[4L]), u)
  coarse <- spatial_endpoint_optimize_reference(function(par) objective(par, coarse_u),
    cbind(rep(initial, 3L), log(c(.3, 1, 2)), log(c(.03, .2, .8)), qlogis(c(.05, .3, .8))),
    lo, hi, maxit = 800L, reltol = 1e-7)
  if (is.null(coarse$candidate)) spatial <- coarse else {
    spatial <- spatial_endpoint_optimize_reference(function(par) objective(par, fine_u),
      rbind(coarse$candidate$par, c(initial, log(.3), log(.03), qlogis(.05))),
      lo, hi, maxit = 1600L, reltol = 1e-7)
  }
  if (!is.null(spatial$candidate)) {
    par <- spatial$candidate$par
    spatial$candidate$mu <- par[1L]; spatial$candidate$sigma <- exp(par[2L])
    spatial$candidate$range <- exp(par[3L]); spatial$candidate$nugget <- plogis(par[4L])
  }
  branches <- list(spatial = spatial, site_only = site, independent = independent)
  audit <- matrix(NA_real_, audit_batches, 3L, dimnames = list(NULL, names(branches)))
  deterministic_error <- setNames(rep(0, 3L), names(branches))
  for (b in seq_len(audit_batches)) {
    u <- spatial_censored_uniforms_reference(setup$unit, audit_draws, audit_seed + b)
    for (j in seq_along(branches)) {
      candidate <- branches[[j]]$candidate
      if (is.null(candidate)) next
      audit[b, j] <- tryCatch(switch(names(branches)[j],
        spatial = spatial_covariance_loglik_reference(lower, upper, coordinates, candidate$mu,
          candidate$sigma, candidate$range, candidate$nugget, u),
        site_only = spatial_endpoint_site_loglik_reference(setup, candidate$mu, candidate$sigma,
          candidate$nugget),
        independent = spatial_endpoint_independent_loglik_reference(setup, candidate$mu, candidate$sigma)),
        error = function(e) NA_real_)
    }
  }
  if (!is.null(site$candidate)) deterministic_error["site_only"] <- tryCatch(abs(
    spatial_endpoint_site_loglik_reference(setup, site$candidate$mu, site$candidate$sigma,
      site$candidate$nugget, engine = "miwa", miwa_steps = 4096L) - audit[1L, "site_only"]),
    error = function(e) Inf)
  result <- spatial_endpoint_compare_reference(branches, audit, deterministic_error)
  c(result, list(branches = branches, coarse = coarse, audit = audit,
    deterministic_error = deterministic_error, groups = setup$groups))
}
