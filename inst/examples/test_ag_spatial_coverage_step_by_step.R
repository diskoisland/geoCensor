#!/usr/bin/env Rscript

# Step-by-step Ag spatial coverage test script.
#
# Run from the USGS_sed_soil/ directory, where the default paths look like:
#   data/ngdbsoil/bestvalue.csv
#   data/ngdbsoil/main.csv
#
# This intentionally keeps the intermediate objects visible:
#   ag, main, ag_main, us_ag_main, mle_state_df, coverage_df, mle_coverage_df

suppressPackageStartupMessages({
  library(geoCensor)
  library(readr)
  library(dplyr)
})

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# Change these lines to test sediment instead of soil.
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

grid_deg <- 0.5
min_uncensored <- 5
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

us_ag_main <- ag_main %>%
  filter(is.na(country) | toupper(country) %in% country_codes)

mle_state_df <- mle_by_state(
  ag_main,
  country_codes = country_codes
) %>%
  mutate(n_uncensored = n_total - n_left - n_right)

too_few <- is.na(mle_state_df$n_uncensored) |
  mle_state_df$n_uncensored < min_uncensored
mle_state_df$geomean[too_few] <- NA_real_
mle_state_df$ci_low[too_few] <- NA_real_
mle_state_df$ci_high[too_few] <- NA_real_
mle_state_df$method[too_few] <- paste0(
  mle_state_df$method[too_few],
  "_min_uncensored_not_met"
)

coverage_df <- spatial_coverage_by_state(
  us_ag_main,
  grid_deg = grid_deg
)

mle_coverage_df <- mle_state_df %>%
  left_join(coverage_df, by = "state") %>%
  mutate(
    geomean_mgkg = geomean / 1000,
    ci_low_mgkg = ci_low / 1000,
    ci_high_mgkg = ci_high / 1000,
    pct_censored = if_else(
      n_total > 0,
      100 * (n_left + n_right) / n_total,
      NA_real_
    ),
    grid_cells_per_1000_samples = if_else(
      n_spatial > 0,
      1000 * occupied_grid_cells / n_spatial,
      NA_real_
    )
  ) %>%
  arrange(desc(geomean_mgkg))

outfile <- file.path(
  "output",
  paste0(medium, "_", element, "_mle_spatial_coverage_step_by_step.csv")
)
write_csv(mle_coverage_df, outfile)

cat(sprintf(
  "\nBuilt %s %s state MLE + spatial coverage table\n",
  toupper(medium),
  element
))
cat(sprintf("ag rows: %d\n", nrow(ag)))
cat(sprintf("main rows: %d\n", nrow(main)))
cat(sprintf("ag_main rows: %d\n", nrow(ag_main)))
cat(sprintf("US ag_main rows: %d\n", nrow(us_ag_main)))
cat(sprintf("states in mle_state_df: %d\n", nrow(mle_state_df)))
cat(sprintf("states in coverage_df: %d\n", nrow(coverage_df)))
cat(sprintf("wrote: %s\n\n", outfile))

print(
  as.data.frame(
    mle_coverage_df %>%
      select(
        state,
        geomean_mgkg,
        ci_low_mgkg,
        ci_high_mgkg,
        n_total,
        n_uncensored,
        n_spatial,
        occupied_grid_cells,
        samples_per_grid_cell,
        grid_cells_per_1000_samples
      ) %>%
      head(12)
  ),
  row.names = FALSE
)
