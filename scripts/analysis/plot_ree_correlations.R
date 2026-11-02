args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1L) args[[1]] else {
  "/Users/rossedwards/.codex/attachments/a1b1afb5-a193-4984-9471-789f711f44db/pasted-text.txt"
}
output_dir <- if (length(args) >= 2L) args[[2]] else "results/ree_correlations"

if (!file.exists(input_path)) stop("Input file not found: ", input_path, call. = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2.", call. = FALSE)
if (!requireNamespace("scales", quietly = TRUE)) stop("Install scales.", call. = FALSE)

raw <- read.delim(
  input_path, header = TRUE, sep = "\t", quote = "", check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = c("N.D.", "N.D", "ND", "")
)

parse_result <- function(x) {
  x <- trimws(as.character(x))
  censored <- grepl("^<", x)
  value <- suppressWarnings(as.numeric(sub("^<", "", x)))
  # Half the reporting limit is used only to rank the small number of censored
  # observations. Results are rank correlations, not regressions on these values.
  value[censored & is.finite(value)] <- value[censored & is.finite(value)] / 2
  list(value = value, censored = censored & is.finite(value))
}

parsed <- lapply(raw, parse_result)
values <- as.data.frame(lapply(parsed, `[[`, "value"), check.names = FALSE)
censored <- as.data.frame(lapply(parsed, `[[`, "censored"), check.names = FALSE)

ree_fields <- c("A_Ce", "A_La", "A_Y")
missing_ree <- setdiff(ree_fields, names(values))
if (length(missing_ree)) stop("Missing REE columns: ", paste(missing_ree, collapse = ", "), call. = FALSE)

element_name <- function(x) sub("^A_", "", x)
bootstrap_replicates <- 2000L
set.seed(20261021)

correlate_pair <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  if (n < 10L || length(unique(x)) < 3L || length(unique(y)) < 3L) {
    return(c(n = n, rho = NA_real_, p_value = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  boot <- replicate(bootstrap_replicates, {
    index <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(x[index], y[index], method = "spearman"))
  })
  boot <- boot[is.finite(boot)]
  ci <- if (length(boot) >= 0.9 * bootstrap_replicates) {
    quantile(boot, c(0.025, 0.975), names = FALSE, type = 7)
  } else c(NA_real_, NA_real_)
  c(n = n, rho = unname(test$estimate), p_value = test$p.value,
    ci_low = ci[1], ci_high = ci[2])
}

results <- do.call(rbind, lapply(ree_fields, function(ree) {
  comparison_fields <- setdiff(names(values), ree)
  do.call(rbind, lapply(comparison_fields, function(field) {
    estimate <- correlate_pair(values[[ree]], values[[field]])
    data.frame(
      rare_earth_element = element_name(ree),
      correlated_element = element_name(field),
      n_pairs = as.integer(estimate["n"]),
      spearman_rho = unname(estimate["rho"]),
      bootstrap_ci_low = unname(estimate["ci_low"]),
      bootstrap_ci_high = unname(estimate["ci_high"]),
      p_value = unname(estimate["p_value"]),
      stringsAsFactors = FALSE
    )
  }))
}))
rownames(results) <- NULL
results$fdr_q_value <- NA_real_
testable <- is.finite(results$p_value)
results$fdr_q_value[testable] <- p.adjust(results$p_value[testable], method = "BH")
results$significant_fdr_0_05 <- results$fdr_q_value < 0.05
results <- results[order(results$rare_earth_element, -abs(results$spearman_rho)), ]

top_n <- 12L
plot_data <- do.call(rbind, lapply(split(results, results$rare_earth_element), function(z) {
  z <- z[is.finite(z$spearman_rho), , drop = FALSE]
  z[seq_len(min(top_n, nrow(z))), , drop = FALSE]
}))
rownames(plot_data) <- NULL
plot_data$display_label <- paste0(
  plot_data$correlated_element,
  ifelse(plot_data$significant_fdr_0_05, " *", "")
)

# Prefix the category with its facet so each panel can have its own ordering.
plot_data$panel_label <- paste(plot_data$rare_earth_element, plot_data$display_label, sep = "___")
ordered_levels <- unlist(lapply(split(plot_data, plot_data$rare_earth_element), function(z) {
  z$panel_label[order(z$spearman_rho)]
}), use.names = FALSE)
plot_data$panel_label <- factor(plot_data$panel_label, levels = unique(ordered_levels))
plot_data$rare_earth_element <- factor(plot_data$rare_earth_element, levels = c("Ce", "La", "Y"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(results, file.path(output_dir, "rare_earth_element_correlations.csv"), row.names = FALSE, na = "")
write.csv(
  data.frame(
    element = names(values),
    n = vapply(values, function(x) sum(is.finite(x)), integer(1)),
    n_below_limit = vapply(censored, sum, integer(1)),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "rare_earth_element_data_summary.csv"), row.names = FALSE
)

correlation_plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = spearman_rho, y = panel_label)
) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.45) +
  ggplot2::geom_segment(
    ggplot2::aes(x = bootstrap_ci_low, xend = bootstrap_ci_high, yend = panel_label),
    colour = "grey60", linewidth = 0.8
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = spearman_rho, yend = panel_label),
    colour = "#2E8B57", linewidth = 0.75
  ) +
  ggplot2::geom_point(
    ggplot2::aes(shape = significant_fdr_0_05),
    colour = "#2E8B57", fill = "white", size = 2.9, stroke = 0.9
  ) +
  ggplot2::facet_wrap(~rare_earth_element, nrow = 1, scales = "free_y") +
  ggplot2::scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1),
    labels = c(`TRUE` = "FDR q < 0.05", `FALSE` = "Not FDR-significant"),
    name = "Correlation status"
  ) +
  ggplot2::scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  ggplot2::scale_x_continuous(
    limits = c(-1, 1), breaks = seq(-1, 1, 0.5),
    labels = scales::label_number(accuracy = 0.1),
    expand = ggplot2::expansion(mult = c(0.02, 0.02))
  ) +
  ggplot2::labs(
    title = "Elements most strongly correlated with rare-earth elements",
    subtitle = paste0(
      "Top ", top_n, " associations per REE; points are Spearman correlations and grey whiskers are 95% bootstrap intervals."
    ),
    x = "Spearman rank correlation (ρ)", y = NULL,
    caption = paste0(
      "* Benjamini–Hochberg FDR q < 0.05 across all REE–element tests. ",
      "Nondetects were ranked at one-half the reporting limit; N.D. values were omitted pairwise."
    )
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_line(linewidth = 0.25, colour = "grey86"),
    strip.background = ggplot2::element_rect(fill = "#EAF1F7", colour = "grey65"),
    strip.text = ggplot2::element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = ggplot2::element_text(face = "bold"),
    plot.caption = ggplot2::element_text(hjust = 0)
  )

ggplot2::ggsave(
  file.path(output_dir, "rare_earth_element_correlations.png"), correlation_plot,
  width = 12.5, height = 7.8, units = "in", dpi = 240, bg = "white"
)
ggplot2::ggsave(
  file.path(output_dir, "rare_earth_element_correlations.pdf"), correlation_plot,
  width = 12.5, height = 7.8, units = "in", device = grDevices::cairo_pdf
)

cat("Wrote", nrow(results), "correlations and", nrow(plot_data), "plotted associations to", output_dir, "\n")
