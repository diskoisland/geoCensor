#!/usr/bin/env Rscript

# Idaho Ag censored lognormal fit test script.
#
# Run from the USGS_sed_soil/ directory, where the default paths look like:
#   data/ngdbsoil/bestvalue.csv
#   data/ngdbsoil/main.csv
#
# This keeps the intermediate objects visible:
#   ag, main, ag_main, us_ag_main, fit_df, fit_summary, p

suppressPackageStartupMessages({
  library(geoCensor)
  library(dplyr)
  library(ggplot2)
})

dir.create("output", showWarnings = FALSE, recursive = TRUE)

medium <- "soil"
element <- "Ag"
state_to_plot <- "ID"
chem_csv <- "data/ngdbsoil/bestvalue.csv"
main_csv <- "data/ngdbsoil/main.csv"
join_fun <- join_soil_main

# Sediment version:
# medium <- "sed"
# chem_csv <- "data/ngdbsed/bestvalue.csv"
# main_csv <- "data/ngdbsed/main.csv"
# join_fun <- join_sed_main

country_codes <- c("US", "USA", "UNITED STATES", "U.S.")

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

# 4. Keep U.S. records.
us_ag_main <- ag_main %>%
  filter(is.na(country) | toupper(country) %in% country_codes)

# 5. Build the Idaho data frame used for the fit and plot.
fit_df <- us_ag_main %>%
  filter(state == state_to_plot)

# 6. Fit the same censored lognormal MLE used by the plot.
fit_summary <- fit_cen_lognorm(fit_df)

# 7. Make and save the default CDF plot.
p <- plot_lognormal_fit(
  fit_df,
  unit_label = "ug kg^-1",
  title = paste0(
    element,
    " in ",
    medium,
    " - Idaho censored lognormal fit"
  )
)

outfile <- file.path(
  "output",
  paste0(medium, "_", element, "_ID_lognormal_fit.png")
)

ggsave(
  filename = outfile,
  plot = p,
  width = 8,
  height = 5,
  dpi = 300
)

cat(sprintf(
  "\nBuilt Idaho %s %s lognormal fit data frame\n",
  medium,
  element
))
cat(sprintf("ag rows: %d\n", nrow(ag)))
cat(sprintf("main rows: %d\n", nrow(main)))
cat(sprintf("ag_main rows: %d\n", nrow(ag_main)))
cat(sprintf("us_ag_main rows: %d\n", nrow(us_ag_main)))
cat(sprintf("Idaho fit_df rows: %d\n", nrow(fit_df)))
cat(sprintf("Idaho exact detections: %d\n", sum(!is.na(fit_df$detected_val) & fit_df$detected_val > 0)))
cat(sprintf("Idaho left-censored: %d\n", sum(fit_df$left_cen, na.rm = TRUE)))
cat(sprintf("Idaho right-censored: %d\n", sum(fit_df$right_cen, na.rm = TRUE)))
cat(sprintf(
  "Idaho geomean = %.3f ug/kg (%.3f mg/kg); 95%% CI %.3f-%.3f ug/kg\n",
  fit_summary$geomean,
  fit_summary$geomean / 1000,
  fit_summary$ci_low,
  fit_summary$ci_high
))
cat(sprintf("wrote: %s\n\n", outfile))

print(as.data.frame(fit_summary), row.names = FALSE)
