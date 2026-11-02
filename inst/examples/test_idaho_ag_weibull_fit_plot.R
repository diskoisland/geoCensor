#!/usr/bin/env Rscript

# Idaho Ag censored distribution-fit comparison.
#
# Run from the USGS_sed_soil/ directory, where the default paths look like:
#   data/ngdbsoil/bestvalue.csv
#   data/ngdbsoil/main.csv
#
# This compares censored lognormal, Weibull, and log-logistic fits.

suppressPackageStartupMessages({
  library(geoCensor)
  library(dplyr)
  library(ggplot2)
  library(survival)
})

dir.create("output", showWarnings = FALSE, recursive = TRUE)

medium <- "soil"
element <- "Ag"
state_to_plot <- "ID"
chem_csv <- "data/ngdbsoil/bestvalue.csv"
main_csv <- "data/ngdbsoil/main.csv"
join_fun <- join_soil_main
country_codes <- c("US", "USA", "UNITED STATES", "U.S.")

build_intervals <- function(df) {
  left_cen <- as.logical(df$left_cen)
  right_cen <- as.logical(df$right_cen)
  left_cen[is.na(left_cen)] <- FALSE
  right_cen[is.na(right_cen)] <- FALSE

  is_left <- left_cen & is.finite(df$left_limit) & df$left_limit > 0
  is_right <- right_cen & is.finite(df$right_limit) & df$right_limit > 0
  is_exact <- is.finite(df$detected_val) & df$detected_val > 0 &
    !is_left & !is_right

  lower <- upper <- rep(NA_real_, nrow(df))
  lower[is_left] <- 0
  upper[is_left] <- df$left_limit[is_left]
  lower[is_right] <- df$right_limit[is_right]
  upper[is_right] <- Inf
  lower[is_exact] <- df$detected_val[is_exact]
  upper[is_exact] <- df$detected_val[is_exact]

  ok <- is_left | is_right | is_exact

  list(
    lower = lower,
    upper = upper,
    ok = ok,
    n_total = sum(ok),
    n_left = sum(is_left[ok]),
    n_right = sum(is_right[ok]),
    n_exact = sum(is_exact[ok])
  )
}

# 1. Build censoring-aware Ag chemistry table.
ag <- build_ag_from_csv(
  chem_csv,
  element = element,
  to_ugkg = TRUE
)

# 2. Read and clean NGDB sample metadata.
main <- read_clean_main(main_csv)

# 3. Join Ag chemistry to state/country/coordinate metadata.
ag_main <- join_fun(
  ag,
  main,
  extra_cols = c("state", "country", "lat_wgs84", "long_wgs84")
)

# 4. Keep U.S. Idaho records.
us_ag_main <- ag_main %>%
  filter(is.na(country) | toupper(country) %in% country_codes)

fit_df <- us_ag_main %>%
  filter(state == state_to_plot)

intervals <- build_intervals(fit_df)

lognormal_fit <- survreg(
  Surv(intervals$lower[intervals$ok], intervals$upper[intervals$ok], type = "interval2") ~ 1,
  dist = "lognormal"
)

weibull_fit <- survreg(
  Surv(intervals$lower[intervals$ok], intervals$upper[intervals$ok], type = "interval2") ~ 1,
  dist = "weibull"
)

loglogistic_lower <- intervals$lower
positive_floor <- min(c(
  intervals$lower[intervals$ok & is.finite(intervals$lower) & intervals$lower > 0],
  intervals$upper[intervals$ok & is.finite(intervals$upper) & intervals$upper > 0]
)) * 1e-6
loglogistic_lower[intervals$ok & loglogistic_lower <= 0] <- positive_floor

loglogistic_fit <- survreg(
  Surv(loglogistic_lower[intervals$ok], intervals$upper[intervals$ok], type = "interval2") ~ 1,
  dist = "loglogistic"
)

lognormal_meanlog <- unname(coef(lognormal_fit)[1])
lognormal_sdlog <- lognormal_fit$scale
lognormal_gm <- exp(lognormal_meanlog)

weibull_scale <- exp(unname(coef(weibull_fit)[1]))
weibull_shape <- 1 / weibull_fit$scale
weibull_median <- weibull_scale * (log(2))^(1 / weibull_shape)
euler_gamma <- -digamma(1)
weibull_geomean <- weibull_scale * exp(-euler_gamma / weibull_shape)

loglogistic_location <- unname(coef(loglogistic_fit)[1])
loglogistic_scale <- loglogistic_fit$scale
loglogistic_median <- exp(loglogistic_location)
loglogistic_geomean <- exp(loglogistic_location)

exact_vals <- intervals$lower[
  intervals$ok &
    is.finite(intervals$lower) &
    intervals$lower > 0 &
    intervals$lower == intervals$upper
]
left_limits <- intervals$upper[
  intervals$ok &
    intervals$lower == 0 &
    is.finite(intervals$upper) &
    intervals$upper > 0
]

plot_vals <- c(
  exact_vals,
  left_limits,
  lognormal_gm,
  weibull_median,
  weibull_geomean,
  loglogistic_median,
  loglogistic_geomean
)
plot_vals <- plot_vals[is.finite(plot_vals) & plot_vals > 0]

curve_x <- exp(seq(log(min(plot_vals) / 1.5), log(max(plot_vals) * 1.5), length.out = 500))

curve_df <- dplyr::bind_rows(
  tibble::tibble(
    value = curve_x,
    log10_value = log10(curve_x),
    cumulative = plnorm(curve_x, meanlog = lognormal_meanlog, sdlog = lognormal_sdlog),
    fit = "lognormal"
  ),
  tibble::tibble(
    value = curve_x,
    log10_value = log10(curve_x),
    cumulative = pweibull(curve_x, shape = weibull_shape, scale = weibull_scale),
    fit = "weibull"
  ),
  tibble::tibble(
    value = curve_x,
    log10_value = log10(curve_x),
    cumulative = stats::plogis(
      (log(curve_x) - loglogistic_location) / loglogistic_scale
    ),
    fit = "loglogistic"
  )
)

exact_cdf_df <- tibble::tibble(value = exact_vals) %>%
  mutate(log10_value = log10(value)) %>%
  arrange(log10_value) %>%
  mutate(detected_only_cdf = ppoints(n()))

rug_df <- tibble::tibble(
  value = left_limits,
  censoring = "left-censored LOD"
) %>%
  mutate(log10_value = log10(value))

lognormal_loglik <- stats::logLik(lognormal_fit)
weibull_loglik <- stats::logLik(weibull_fit)
loglogistic_loglik <- stats::logLik(loglogistic_fit)

aic_tbl <- tibble::tibble(
  distribution = c("lognormal", "weibull", "loglogistic"),
  n_total = intervals$n_total,
  n_left = intervals$n_left,
  n_right = intervals$n_right,
  n_exact = intervals$n_exact,
  logLik = c(
    as.numeric(lognormal_loglik),
    as.numeric(weibull_loglik),
    as.numeric(loglogistic_loglik)
  ),
  df = c(
    as.numeric(attr(lognormal_loglik, "df")),
    as.numeric(attr(weibull_loglik, "df")),
    as.numeric(attr(loglogistic_loglik, "df"))
  ),
  AIC = c(
    stats::AIC(lognormal_fit),
    stats::AIC(weibull_fit),
    stats::AIC(loglogistic_fit)
  ),
  note = c(
    NA_character_,
    NA_character_,
    paste0("left lower bound set to ", signif(positive_floor, 3), " for fit")
  )
) %>%
  arrange(AIC) %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    akaike_weight = {
      w <- exp(-0.5 * delta_AIC)
      w / sum(w)
    },
    best = delta_AIC == 0
  )

subtitle <- sprintf(
  "n = %d; exact = %d; left-censored = %d; right-censored = %d",
  intervals$n_total,
  intervals$n_exact,
  intervals$n_left,
  intervals$n_right
)

p <- ggplot() +
  geom_point(
    data = exact_cdf_df,
    aes(x = log10_value, y = detected_only_cdf),
    color = "grey35",
    alpha = 0.12,
    size = 0.65
  ) +
  geom_line(
    data = curve_df,
    aes(x = log10_value, y = cumulative, color = fit),
    linewidth = 1
  ) +
  geom_hline(yintercept = 0.5, color = "grey55", linetype = "dotted") +
  geom_vline(
    xintercept = log10(lognormal_gm),
    color = "#0072B2",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = log10(weibull_median),
    color = "#CC79A7",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = log10(loglogistic_median),
    color = "#E69F00",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_rug(
    data = rug_df,
    aes(x = log10_value),
    sides = "b",
    color = "#009E73",
    alpha = 0.45
  ) +
  annotate(
    "label",
    x = log10(lognormal_gm),
    y = 0.10,
    label = sprintf("Lognormal GM\n%.0f ug/kg", lognormal_gm),
    size = 3,
    fill = "white",
    color = "grey20"
  ) +
  annotate(
    "label",
    x = log10(weibull_median),
    y = 0.22,
    label = sprintf("Weibull median\n%.0f ug/kg", weibull_median),
    size = 3,
    fill = "white",
    color = "grey20"
  ) +
  annotate(
    "label",
    x = log10(loglogistic_median),
    y = 0.34,
    label = sprintf("Log-logistic median\n%.0f ug/kg", loglogistic_median),
    size = 3,
    fill = "white",
    color = "grey20"
  ) +
  scale_color_manual(
    values = c(
      lognormal = "#0072B2",
      weibull = "#CC79A7",
      loglogistic = "#E69F00"
    ),
    name = "Fitted distribution"
  ) +
  scale_x_continuous(
    name = "Concentration (ug kg^-1, log10 scale)",
    labels = function(x) format(signif(10^x, 3), scientific = FALSE)
  ) +
  scale_y_continuous(
    name = "Fitted cumulative probability",
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = function(x) paste0(round(100 * x), "%")
  ) +
  labs(
    title = "Ag in soil - Idaho censored distribution fits",
    subtitle = subtitle,
    caption = "Gray points are exact detections only; green rug marks are left-censoring limits."
  ) +
  theme_minimal(base_size = 11)

outfile <- file.path(
  "output",
  paste0(medium, "_", element, "_ID_lognormal_weibull_loglogistic_cdf.png")
)

ggsave(
  filename = outfile,
  plot = p,
  width = 8,
  height = 5,
  dpi = 300
)

cat(sprintf(
  "\nBuilt Idaho %s %s lognormal vs Weibull vs log-logistic CDF comparison\n",
  medium,
  element
))
cat(sprintf("Idaho fit_df rows: %d\n", nrow(fit_df)))
cat(sprintf("Idaho exact detections: %d\n", intervals$n_exact))
cat(sprintf("Idaho left-censored: %d\n", intervals$n_left))
cat(sprintf("Lognormal GM = %.3f ug/kg\n", lognormal_gm))
cat(sprintf("Weibull median = %.3f ug/kg\n", weibull_median))
cat(sprintf("Weibull geometric mean = %.3f ug/kg\n", weibull_geomean))
cat(sprintf("Log-logistic median/GM = %.3f ug/kg\n", loglogistic_median))
cat(sprintf("Log-logistic left lower bound for fit = %.6g ug/kg\n", positive_floor))
cat(sprintf("wrote: %s\n\n", outfile))

print(as.data.frame(aic_tbl), row.names = FALSE)
