#!/usr/bin/env Rscript

# Build Ag state MLE tables with spatial coverage diagnostics.
#
# Run from the USGS_sed_soil/ directory, where the default data paths look like:
#   data/ngdbsoil/bestvalue.csv
#   data/ngdbsoil/main.csv
#   data/ngdbsed/bestvalue.csv
#   data/ngdbsed/main.csv
#
# Outputs include the censoring-aware MLE fields plus valid coordinate counts,
# occupied lat/long grid cells, coordinate ranges, and simple coverage ratios.

suppressPackageStartupMessages({
  library(geoCensor)
  library(readr)
  library(dplyr)
})

dir.create("output", showWarnings = FALSE, recursive = TRUE)

element <- "Ag"
grid_deg <- 0.5
min_uncensored <- 5
country_codes <- c("US", "USA", "UNITED STATES", "U.S.")
seed_states <- c("CA", "CO", "ID", "KS", "ND", "NM", "NV", "TX", "UT", "WY")

apply_min_uncensored <- function(x, min_uncensored) {
  x <- x %>%
    mutate(n_uncensored = n_total - n_left - n_right)

  too_few <- is.na(x$n_uncensored) | x$n_uncensored < min_uncensored
  x$geomean[too_few] <- NA_real_
  x$ci_low[too_few] <- NA_real_
  x$ci_high[too_few] <- NA_real_
  x$method[too_few] <- paste0(x$method[too_few], "_min_uncensored_not_met")
  x
}

build_medium <- function(medium, chem_csv, main_csv, join_fun) {
  element_data <- build_ag_from_csv(
    chem_csv,
    element = element,
    to_ugkg = TRUE
  )

  main_data <- read_clean_main(main_csv)

  element_main <- join_fun(
    element_data,
    main_data,
    extra_cols = c("state", "country", "lat_wgs84", "long_wgs84")
  )

  us_element_main <- element_main %>%
    filter(is.na(country) | toupper(country) %in% country_codes)

  mle_state_df <- mle_by_state(
    element_main,
    country_codes = country_codes
  ) %>%
    apply_min_uncensored(min_uncensored)

  coverage_df <- spatial_coverage_by_state(
    us_element_main,
    grid_deg = grid_deg
  )

  out <- mle_state_df %>%
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
    )

  outfile <- file.path(
    "output",
    paste0(medium, "_", element, "_mle_state_df_with_spatial_coverage.csv")
  )
  write_csv(out, outfile)

  cat(sprintf(
    "[%s] %d states, %d fitted, grid %.2f degrees, wrote %s\n",
    toupper(medium),
    nrow(out),
    sum(!is.na(out$geomean)),
    grid_deg,
    outfile
  ))

  out
}

soil <- build_medium(
  medium = "soil",
  chem_csv = "data/ngdbsoil/bestvalue.csv",
  main_csv = "data/ngdbsoil/main.csv",
  join_fun = join_soil_main
)

sed <- build_medium(
  medium = "sed",
  chem_csv = "data/ngdbsed/bestvalue.csv",
  main_csv = "data/ngdbsed/main.csv",
  join_fun = join_sed_main
)

seed_coverage <- soil %>%
  select(
    state,
    gm_soil = geomean_mgkg,
    lo_soil = ci_low_mgkg,
    hi_soil = ci_high_mgkg,
    n_total_soil = n_total,
    cen_soil = pct_censored,
    n_spatial_soil = n_spatial,
    grid_cells_soil = occupied_grid_cells,
    cells_per_1000_soil = grid_cells_per_1000_samples
  ) %>%
  full_join(
    sed %>%
      select(
        state,
        gm_sed = geomean_mgkg,
        lo_sed = ci_low_mgkg,
        hi_sed = ci_high_mgkg,
        n_total_sed = n_total,
        cen_sed = pct_censored,
        n_spatial_sed = n_spatial,
        grid_cells_sed = occupied_grid_cells,
        cells_per_1000_sed = grid_cells_per_1000_samples
      ),
    by = "state"
  ) %>%
  filter(state %in% seed_states) %>%
  mutate(
    across(starts_with("gm_"), ~ round(.x, 3)),
    across(starts_with("lo_"), ~ round(.x, 3)),
    across(starts_with("hi_"), ~ round(.x, 3)),
    across(starts_with("cen_"), ~ round(.x, 1)),
    across(starts_with("cells_per_1000"), ~ round(.x, 1))
  ) %>%
  arrange(desc(gm_soil))

write_csv(
  seed_coverage,
  "output/seedstate_soil_vs_sed_Ag_mle_spatial_coverage.csv"
)

cat("\n=== seeding-state Ag MLE and spatial coverage ===\n")
print(as.data.frame(seed_coverage), row.names = FALSE)
