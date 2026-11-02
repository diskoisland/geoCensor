# Small, known-covariance Gaussian reference on log concentration.
# Research only: no production API, fitted covariance or interval calibration.

spatial_censored_setup_reference <- function(lower, upper, covariance,
    max_records = 96L) {
  n <- length(lower)
  if (!is.numeric(lower) || !is.numeric(upper) || n < 1L || n > max_records ||
      length(upper) != n || anyNA(lower) || anyNA(upper) || any(lower > upper) ||
      any(lower == Inf) || any(upper == -Inf) ||
      any(lower == -Inf & upper == Inf)) stop("Invalid log-response bounds")
  S <- as.matrix(covariance)
  if (!is.numeric(S) || !identical(dim(S), c(n, n)) || any(!is.finite(S)) ||
      max(abs(S - t(S))) > 1e-12) stop("Invalid covariance")
  chol(S) # Require positive definiteness; never add a silent ridge.
  exact <- which(lower == upper)
  censored <- which(lower != upper)
  no <- length(exact); nc <- length(censored)
  if (no) {
    C <- chol(S[exact, exact, drop = FALSE])
    inv <- chol2inv(C)
    constant <- -no / 2 * log(2 * pi) - sum(log(diag(C)))
  } else {
    inv <- matrix(numeric(), 0L, 0L); constant <- 0
  }
  if (nc) {
    A <- if (no) S[censored, exact, drop = FALSE] %*% inv else matrix(0, nc, 0L)
    conditional <- S[censored, censored, drop = FALSE] -
      if (no) A %*% S[exact, censored, drop = FALSE] else matrix(0, nc, nc)
    conditional <- (conditional + t(conditional)) / 2
    L <- t(chol(conditional))
    packed <- mvtnorm::ltMatrices(L[lower.tri(L, diag = TRUE)], diag = TRUE)
    intercept <- if (no) as.vector(A %*% lower[exact]) else rep(0, nc)
    slope <- 1 - rowSums(A)
  } else {
    conditional <- A <- packed <- NULL
    intercept <- slope <- numeric()
  }
  structure(list(lower = lower, upper = upper, S = S, exact = exact,
    censored = censored, inv = inv, constant = constant, A = A,
    conditional = conditional, packed = packed, intercept = intercept,
    slope = slope, n = n, no = no, nc = nc), class = "spatial_censored_reference")
}

spatial_censored_loglik_reference <- function(model, mu, uniforms = NULL,
    engine = c("simulation", "miwa"), miwa_steps = 128L) {
  engine <- match.arg(engine)
  if (length(mu) != 1L || !is.finite(mu)) stop("Invalid location")
  d <- model$lower[model$exact] - mu
  exact_ll <- model$constant - if (model$no) drop(crossprod(d, model$inv %*% d)) / 2 else 0
  if (!model$nc) return(exact_ll)
  index <- model$censored
  mean <- model$intercept + model$slope * mu
  lower <- model$lower[index] - mean
  upper <- model$upper[index] - mean
  if (model$nc == 1L) {
    a <- lower / sqrt(model$conditional[1L, 1L])
    b <- upper / sqrt(model$conditional[1L, 1L])
    # Tail-aware log difference avoids subtracting two CDFs rounded to one.
    if (a > 0) {
      x <- pnorm(a, lower.tail = FALSE, log.p = TRUE)
      y <- pnorm(b, lower.tail = FALSE, log.p = TRUE)
    } else {
      x <- pnorm(b, log.p = TRUE); y <- pnorm(a, log.p = TRUE)
    }
    rectangle <- x + log(-expm1(y - x))
  } else if (engine == "miwa") {
    if (model$nc > 8L) stop("Deterministic check limited to eight censored dimensions")
    # Mixed left/right orthants are converted to all-upper orthants. General
    # finite boxes use pmvnorm's documented Miwa handling (warnings propagate).
    sign <- ifelse(is.infinite(upper) & is.finite(lower), -1, 1)
    transformed_lower <- ifelse(sign > 0, lower, -upper)
    transformed_upper <- ifelse(sign > 0, upper, -lower)
    probability <- mvtnorm::pmvnorm(lower = transformed_lower,
      upper = transformed_upper,
      sigma = model$conditional * outer(sign, sign),
      algorithm = mvtnorm::Miwa(steps = miwa_steps))
    if (!is.finite(probability) || probability <= 0) stop("Rectangle probability unavailable")
    rectangle <- log(as.numeric(probability))
  } else {
    if (!is.matrix(uniforms) || nrow(uniforms) != model$nc - 1L ||
        ncol(uniforms) < 2L || any(!is.finite(uniforms)) ||
        any(uniforms <= 0 | uniforms >= 1)) stop("Provide fixed integration uniforms")
    # Reflect right tails (and wholly positive finite intervals) so the small
    # probability is evaluated in the lower tail. CDF subtraction near one can
    # otherwise lose it, including when the routine returns a finite floor.
    sign <- ifelse((is.infinite(upper) & is.finite(lower)) |
                     (is.finite(lower) & is.finite(upper) & lower > 0), -1, 1)
    reflected_lower <- ifelse(sign > 0, lower, -upper)
    reflected_upper <- ifelse(sign > 0, upper, -lower)
    L <- t(chol(model$conditional * outer(sign, sign)))
    packed <- mvtnorm::ltMatrices(L[lower.tri(L, diag = TRUE)], diag = TRUE)
    had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_rng) saved_rng <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(if (had_rng) assign(".Random.seed", saved_rng, envir = .GlobalEnv) else
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv))
    rectangle <- as.numeric(mvtnorm::lpmvnorm(lower = matrix(reflected_lower, ncol = 1L),
      upper = matrix(reflected_upper, ncol = 1L), chol = packed,
      M = ncol(uniforms), w = uniforms, fast = FALSE,
      tol = 0)) # Underflow must be unavailable, not a finite tolerance floor.
  }
  value <- exact_ll + rectangle
  if (!is.finite(value)) stop("Nonfinite likelihood; do not clip silently")
  value
}

spatial_censored_uniforms_reference <- function(model, draws, seed) {
  # Preserve the caller's RNG state. Uniforms are reused across locations in
  # each optimization, and are ordinary Monte Carlo, not a claimed QMC lattice.
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had) saved <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had) assign(".Random.seed", saved, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv))
  set.seed(seed)
  if (model$nc < 2L) return(NULL)
  matrix(runif((model$nc - 1L) * draws), nrow = model$nc - 1L)
}

spatial_censored_fit_mu_reference <- function(model, uniforms = NULL,
    engine = "simulation", tolerance = 1e-6) {
  # With only upper (or only lower) censoring and no exact observations the
  # likelihood has a supremum at an infinite location, even with covariance known.
  if (!model$no && all(model$lower == -Inf))
    return(list(status = "no_finite_location_mle", theta = -Inf, loglik = 0))
  if (!model$no && all(model$upper == Inf))
    return(list(status = "no_finite_location_mle", theta = Inf, loglik = 0))
  observed <- c(model$lower[is.finite(model$lower)], model$upper[is.finite(model$upper)])
  center <- median(observed); step <- max(sqrt(diag(model$S)))
  loglik <- function(mu) spatial_censored_loglik_reference(model, mu, uniforms, engine)
  for (k in seq_len(6L)) {
    interval <- center + c(-1, 1) * step * 2^(k + 1L)
    fit <- optimize(loglik, interval = interval, maximum = TRUE, tol = tolerance)
    if (fit$maximum > interval[1L] + 0.01 * diff(interval) &&
        fit$maximum < interval[2L] - 0.01 * diff(interval)) {
      return(list(status = "finite", theta = fit$maximum, loglik = fit$objective,
                  interval = interval))
    }
  }
  stop("No interior maximum within expanding numerical bracket")
}
