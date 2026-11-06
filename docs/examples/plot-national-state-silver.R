#!/usr/bin/env Rscript
# Redraw the public figures from aggregate tables; no source database needed.
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 first")
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- if (length(arg)) sub("^--file=", "", arg[1]) else NA_character_
if (length(arg) && !file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
here <- if (length(arg)) dirname(normalizePath(script, mustWork = TRUE)) else "docs/examples"
figures <- file.path(here, "..", "figures")
dir.create(figures, showWarnings = FALSE)
national <- read.csv(file.path(here, "data", "silver-national-grid.csv"), na.strings = c("", "NA"))
states <- read.csv(file.path(here, "data", "silver-state-summary.csv"), na.strings = c("", "NA"))
national$medium <- factor(national$medium, levels = c("soil", "sed"), labels = c("Soil", "Sediment"))
national$domain <- factor(national$domain, levels = c("US50_DC", "CONUS"),
  labels = c("50 states + DC", "Contiguous 48 states + DC"))
colours <- c(Soil = "#087F8C", Sediment = "#9B542D")
theme_silver <- function() ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(
  panel.grid.minor = ggplot2::element_blank(), panel.grid.major = ggplot2::element_line(colour = "#E5EAF0"),
  plot.title = ggplot2::element_text(face = "bold", size = 19, colour = "#142D43"),
  plot.subtitle = ggplot2::element_text(colour = "#41576A", margin = ggplot2::margin(b = 12)),
  strip.text = ggplot2::element_text(face = "bold", colour = "#142D43"),
  plot.caption = ggplot2::element_text(hjust = 0, colour = "#41576A", size = 9),
  legend.position = "bottom", plot.margin = ggplot2::margin(16, 20, 16, 16))
baseline <- unique(national[, c("medium", "domain", "baseline_gm")])
p <- ggplot2::ggplot(national, ggplot2::aes(size_km, gm, colour = medium)) +
  ggplot2::geom_hline(data = baseline, ggplot2::aes(yintercept = baseline_gm, colour = medium),
    linetype = "dashed", linewidth = 0.65, alpha = 0.8) +
  ggplot2::geom_line(linewidth = 1.1) + ggplot2::geom_point(size = 3) +
  ggplot2::facet_wrap(~domain, nrow = 1) + ggplot2::scale_colour_manual(values = colours) +
  ggplot2::scale_x_continuous(breaks = c(20, 40, 80)) +
  ggplot2::scale_y_continuous(limits = c(0, 0.10), breaks = seq(0, 0.10, 0.02)) +
  ggplot2::labs(title = "Silver: how much does spatial weighting change the estimate?",
    subtitle = "Censored lognormal geometric means | same usable records within each domain and medium",
    x = "Declustering cell width (km)", y = "Silver geometric mean (mg/kg)", colour = NULL,
    caption = "Solid: cell-declustered fit across 100 placements. Dashed: unweighted fit.\nOccupied sampled support only; grid sensitivity is not a confidence interval.") + theme_silver()
ggplot2::ggsave(file.path(figures, "silver-national-grid.png"), p, width = 11, height = 5.3, dpi = 160, bg = "white")

states$medium <- factor(states$medium, levels = c("soil", "sed"), labels = c("Soil", "Sediment"))
if (anyNA(states$medium)) stop("Unexpected medium code")
state_order <- rev(sort(c(state.abb, "DC")))
states$state <- factor(states$state, levels = state_order)
states$display_gm <- states$geomean_20km_mgkg
finite <- is.finite(states$display_gm) & states$display_gm > 0
states$missing_label <- ifelse(states$n_usable == 0, "No usable Ag records", "Fit unavailable")
states$review <- ifelse(states$information_review_20km == "review_information", "Review flag", "Finite estimate")
p <- ggplot2::ggplot(states, ggplot2::aes(x = display_gm, y = state, colour = medium)) +
  ggplot2::geom_point(data = states[finite, ], ggplot2::aes(shape = review), size = 2.5, stroke = 0.9) +
  ggplot2::geom_text(data = states[!finite, ], ggplot2::aes(x = 0.00012, label = missing_label),
    inherit.aes = TRUE, colour = "#687787", hjust = 0, size = 2.8) +
  ggplot2::facet_wrap(~medium, nrow = 1, drop = FALSE) +
  ggplot2::scale_y_discrete(limits = state_order, drop = FALSE) +
  ggplot2::scale_x_log10(limits = c(0.0001, 1), breaks = c(0.0001, 0.001, 0.01, 0.1, 1),
    labels = c("0.0001", "0.001", "0.01", "0.1", "1")) +
  ggplot2::scale_colour_manual(values = colours, guide = "none") +
  ggplot2::scale_shape_manual(values = c("Finite estimate" = 16, "Review flag" = 1), name = NULL) +
  ggplot2::labs(title = "State silver estimates after 20 km cell declustering",
    subtitle = "Pooled recorded-method fits | unavailable estimates stay missing",
    x = "Silver geometric mean (mg/kg; logarithmic scale)", y = NULL,
    caption = "Open circles retain count-based information-review flags; filled circles do not certify reliability.\nOccupied sampled support only. No spatial confidence limits; see tables for 5/10/20 km sensitivity.") +
  theme_silver() + ggplot2::theme(panel.grid.major.y = ggplot2::element_line(linewidth = 0.25),
    axis.text.y = ggplot2::element_text(size = 10))
if (any(finite & (states$display_gm < 0.0001 | states$display_gm > 1))) stop("Finite state estimate outside plot range")
ggplot2::ggsave(file.path(figures, "silver-states.png"), p, width = 10.5, height = 13.5, dpi = 160, bg = "white")
cat("Redrew two silver figures from public aggregate tables.\n")
