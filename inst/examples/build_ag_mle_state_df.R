#!/usr/bin/env Rscript

# Build an mle_state_df data frame for Ag from USGS NGDB bestvalue + main files.
#
# Run from the USGS_sed_soil/ directory, where the default data paths look like:
#   data/ngdbsoil/bestvalue.csv
#   data/ngdbsoil/main.csv
#
# The output object is named mle_state_df and has the same structure returned by
# geoCensor::mle_by_state(). It can be used directly with plot_ranked_bar() and
# plot_us_map().

suppressPackageStartupMessages({
  library(geoCensor)
  library(readr)
  library(dplyr)
})

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# Change these four lines to build a different medium or element.
medium <- "soil"
element <- "Ag"
chem_csv <- "data/ngdbsoil/bestvalue.csv"
main_csv <- "data/ngdbsoil/main.csv"
min_uncensored <- 5

# For sediment, use:
# medium <- "sed"
# chem_csv <- "data/ngdbsed/bestvalue.csv"
# main_csv <- "data/ngdbsed/main.csv"

country_codes <- c("US", "USA", "UNITED STATES", "U.S.")
join_fun <- if (identical(medium, "sed")) join_sed_main else join_soil_main

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

element_data <- build_ag_from_csv(
  chem_csv,
  element = element,
  to_ugkg = TRUE
)

main_data <- read_clean_main(main_csv)

element_main <- join_fun(
  element_data,
  main_data,
  extra_cols = c("state", "country")
)

mle_state_df <- mle_by_state(
  element_main,
  country_codes = country_codes
) %>%
  apply_min_uncensored(min_uncensored)

mle_state_df_mgkg <- mle_state_df %>%
  mutate(
    geomean_mgkg = geomean / 1000,
    ci_low_mgkg = ci_low / 1000,
    ci_high_mgkg = ci_high / 1000,
    pct_censored = dplyr::if_else(
      n_total > 0,
      100 * (n_left + n_right) / n_total,
      NA_real_
    )
  )

base_name <- paste0(medium, "_", element, "_mle_state_df")

write_csv(mle_state_df, file.path("output", paste0(base_name, "_ugkg.csv")))
write_csv(mle_state_df_mgkg, file.path("output", paste0(base_name, "_mgkg.csv")))
saveRDS(mle_state_df, file.path("output", paste0(base_name, ".rds")))

cat(sprintf(
  "Built %s for %s %s: %d states (%d fitted, min uncensored = %d), output/%s_ugkg.csv\n",
  "mle_state_df",
  toupper(medium),
  element,
  nrow(mle_state_df),
  sum(!is.na(mle_state_df$geomean)),
  min_uncensored,
  base_name
))

print(head(as.data.frame(mle_state_df)))
