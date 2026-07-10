#!/usr/bin/env Rscript

# Build per-state and pooled Ag MLE geomeans for USGS NGDB soil and sediment.
#
# Run from the USGS_sed_soil/ directory, where the data paths look like:
#   data/ngdbsoil/bestvalue.csv
#   data/ngdbsoil/main.csv
#   data/ngdbsed/bestvalue.csv
#   data/ngdbsed/main.csv

suppressPackageStartupMessages({
  library(geoCensor)
  library(readr)
  library(dplyr)
})

dir.create("output", showWarnings = FALSE, recursive = TRUE)

country_codes <- c("US", "USA", "UNITED STATES", "U.S.")
seed_states <- c("CA", "CO", "ID", "KS", "ND", "NM", "NV", "TX", "UT", "WY")
min_uncensored <- 5

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

run_medium <- function(chem_csv, main_csv, label, outfile, join_fun) {
  ag <- build_ag_from_csv(
    chem_csv,
    element = "Ag",
    to_ugkg = TRUE
  )

  main <- read_clean_main(main_csv)

  ag_main <- join_fun(
    ag,
    main,
    extra_cols = c("state", "country")
  )

  us_data <- ag_main %>%
    filter(is.na(country) | toupper(country) %in% country_codes)

  state_res <- mle_by_state(
    ag_main,
    country_codes = country_codes
  ) %>%
    apply_min_uncensored(min_uncensored)

  write_csv(state_res, outfile)

  national <- fit_cen_lognorm(us_data) %>%
    apply_min_uncensored(min_uncensored)

  cat(sprintf(
    "[%s] %d states, national pooled GM = %.1f ug/kg (%.4f mg/kg), n=%d, min uncensored=%d\n",
    label,
    sum(!is.na(state_res$geomean)),
    national$geomean,
    national$geomean / 1000,
    national$n_total,
    min_uncensored
  ))

  list(
    state = state_res,
    us_data = us_data,
    national = national
  )
}

percent_censored <- function(x) {
  dplyr::if_else(
    x$n_total > 0,
    round(100 * (x$n_left + x$n_right) / x$n_total),
    NA_real_
  )
}

state_meta <- function(x, suffix) {
  x %>%
    transmute(
      state,
      gm = round(geomean / 1000, 3),
      lo = round(ci_low / 1000, 3),
      hi = round(ci_high / 1000, 3),
      n_total = n_total,
      cen = dplyr::if_else(
        n_total > 0,
        round(100 * (n_left + n_right) / n_total),
        NA_real_
      )
    ) %>%
    rename_with(~ paste0(.x, "_", suffix), -state)
}

pooled_meta <- function(x, suffix) {
  tibble(
    gm = round(x$geomean / 1000, 3),
    lo = round(x$ci_low / 1000, 3),
    hi = round(x$ci_high / 1000, 3),
    n_total = x$n_total,
    cen = percent_censored(x)
  ) %>%
    rename_with(~ paste0(.x, "_", suffix))
}

soil_run <- run_medium(
  chem_csv = "data/ngdbsoil/bestvalue.csv",
  main_csv = "data/ngdbsoil/main.csv",
  label = "SOIL",
  outfile = "output/soil_Ag_geomean_by_state_MLE_v2.csv",
  join_fun = join_soil_main
)

sed_run <- run_medium(
  chem_csv = "data/ngdbsed/bestvalue.csv",
  main_csv = "data/ngdbsed/main.csv",
  label = "SEDIMENT",
  outfile = "output/sed_Ag_geomean_by_state_MLE_v2.csv",
  join_fun = join_sed_main
)

cmp_states <- state_meta(soil_run$state, "soil") %>%
  full_join(state_meta(sed_run$state, "sed"), by = "state") %>%
  filter(state %in% seed_states) %>%
  arrange(desc(gm_soil))

pooled_seed_states <- bind_cols(
  tibble(state = "10-state pooled"),
  soil_run$us_data %>%
    filter(state %in% seed_states) %>%
    fit_cen_lognorm() %>%
    apply_min_uncensored(min_uncensored) %>%
    pooled_meta("soil"),
  sed_run$us_data %>%
    filter(state %in% seed_states) %>%
    fit_cen_lognorm() %>%
    apply_min_uncensored(min_uncensored) %>%
    pooled_meta("sed")
)

pooled_us <- bind_cols(
  tibble(state = "US pooled"),
  pooled_meta(soil_run$national, "soil"),
  pooled_meta(sed_run$national, "sed")
)

cmp <- bind_rows(
  cmp_states,
  pooled_seed_states,
  pooled_us
)

cat("\n=== seeding-state SOIL vs SEDIMENT MLE Ag (mg/kg) ===\n")
print(as.data.frame(cmp), row.names = FALSE)

write_csv(cmp, "output/seedstate_soil_vs_sed_Ag_mgkg.csv")
