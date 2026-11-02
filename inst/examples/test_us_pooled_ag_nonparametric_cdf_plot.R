#!/usr/bin/env Rscript

# U.S. pooled Ag lognormal CDF with nonparametric censored-CDF overlay.
#
# Run from the USGS_sed_soil/ directory, where the default paths look like:
#   data/ngdbsoil/bestvalue.csv
#   data/ngdbsoil/main.csv
#
# This keeps the intermediate objects visible:
#   ag, main, ag_main, fit_df, fit_summary, p

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
country_codes <- c("US", "USA", "UNITED STATES", "U.S.")

# Sediment version:
# medium <- "sed"
# chem_csv <- "data/ngdbsed/bestvalue.csv"
# main_csv <- "data/ngdbsed/main.csv"
# join_fun <- join_sed_main

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

fit_df <- ag_main %>%
  filter(is.na(country) | toupper(country) %in% country_codes)

fit_summary <- fit_cen_lognorm(fit_df)

p <- plot_lognormal_fit(
  fit_df,
  unit_label = "ug kg^-1",
  title = paste0(
    element,
    " in ",
    medium,
    " - U.S. pooled censored CDF diagnostic"
  ),
  show_nonparametric = TRUE
)

outfile <- file.path(
  "output",
  paste0(medium, "_", element, "_US_pooled_lognormal_nonparametric_cdf.png")
)

ggsave(
  filename = outfile,
  plot = p,
  width = 8,
  height = 5,
  dpi = 300
)

cat(sprintf(
  "\nBuilt U.S. pooled %s %s lognormal CDF with nonparametric overlay\n",
  medium,
  element
))
cat(sprintf("ag rows: %d\n", nrow(ag)))
cat(sprintf("main rows: %d\n", nrow(main)))
cat(sprintf("ag_main rows: %d\n", nrow(ag_main)))
cat(sprintf("fit_df rows: %d\n", nrow(fit_df)))
cat(sprintf("exact detections: %d\n", sum(!is.na(fit_df$detected_val) & fit_df$detected_val > 0)))
cat(sprintf("left-censored: %d\n", sum(fit_df$left_cen, na.rm = TRUE)))
cat(sprintf("right-censored: %d\n", sum(fit_df$right_cen, na.rm = TRUE)))
cat(sprintf(
  "lognormal geomean = %.3f ug/kg (%.3f mg/kg); 95%% CI %.3f-%.3f ug/kg\n",
  fit_summary$geomean,
  fit_summary$geomean / 1000,
  fit_summary$ci_low,
  fit_summary$ci_high
))
cat(sprintf("wrote: %s\n\n", outfile))

print(as.data.frame(fit_summary), row.names = FALSE)
