#!/usr/bin/env Rscript

# Plot the censored lognormal MLE fit for Ag.
#
# Run from the USGS_sed_soil/ directory, where the default data paths look like:
#   data/ngdbsoil/bestvalue.csv
#   data/ngdbsoil/main.csv
#
# By default this plots pooled soil Ag for the 10 seeding states. Set
# state_to_plot to a two-letter state code to make a state-specific plot.

suppressPackageStartupMessages({
  library(geoCensor)
  library(dplyr)
  library(ggplot2)
})

dir.create("output", showWarnings = FALSE, recursive = TRUE)

medium <- "soil"
element <- "Ag"
chem_csv <- "data/ngdbsoil/bestvalue.csv"
main_csv <- "data/ngdbsoil/main.csv"
join_fun <- join_soil_main

# Sediment version:
# medium <- "sed"
# chem_csv <- "data/ngdbsed/bestvalue.csv"
# main_csv <- "data/ngdbsed/main.csv"
# join_fun <- join_sed_main

seed_states <- c("CA", "CO", "ID", "KS", "ND", "NM", "NV", "TX", "UT", "WY")
state_to_plot <- NULL
country_codes <- c("US", "USA", "UNITED STATES", "U.S.")

ag <- build_ag_from_csv(
  chem_csv,
  element = element,
  to_ugkg = TRUE
)

main <- read_clean_main(main_csv)

ag_main <- join_fun(
  ag,
  main,
  extra_cols = c("state", "country", "lat_wgs84", "long_wgs84")
)

plot_df <- ag_main %>%
  filter(is.na(country) | toupper(country) %in% country_codes)

if (is.null(state_to_plot)) {
  plot_df <- plot_df %>%
    filter(state %in% seed_states)
  plot_label <- "10-state pooled"
  outfile <- file.path(
    "output",
    paste0(medium, "_", element, "_seedstates_lognormal_fit.png")
  )
} else {
  plot_df <- plot_df %>%
    filter(state == state_to_plot)
  plot_label <- state_to_plot
  outfile <- file.path(
    "output",
    paste0(medium, "_", element, "_", state_to_plot, "_lognormal_fit.png")
  )
}

fit_summary <- fit_cen_lognorm(plot_df)

p <- plot_lognormal_fit(
  plot_df,
  unit_label = "ug kg^-1",
  title = paste0(
    element,
    " in ",
    medium,
    " - ",
    plot_label,
    " censored lognormal fit"
  )
)

ggsave(
  filename = outfile,
  plot = p,
  width = 8,
  height = 5,
  dpi = 300
)

cat(sprintf(
  "\n%s %s %s lognormal MLE geomean = %.3f ug/kg; 95%% CI %.3f-%.3f ug/kg; n = %d\n",
  toupper(medium),
  element,
  plot_label,
  fit_summary$geomean,
  fit_summary$ci_low,
  fit_summary$ci_high,
  fit_summary$n_total
))
cat(sprintf("wrote: %s\n", outfile))
