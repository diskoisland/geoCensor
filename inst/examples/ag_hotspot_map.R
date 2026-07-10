#!/usr/bin/env Rscript

# Build an Ag hotspot table and state map from USGS NGDB bestvalue + main files.
#
# Run from the USGS_sed_soil/ directory, where the default data paths look like:
#   data/ngdbsed/bestvalue.csv
#   data/ngdbsed/main.csv
#
# The default example maps sediment Ag hotspots in Idaho at >= 2 mg/kg.

suppressPackageStartupMessages({
  library(geoCensor)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# Change these settings for a different map.
medium <- "sed"
element <- "Ag"
state <- "ID"
threshold_mgkg <- 2
chem_csv <- "data/ngdbsed/bestvalue.csv"
main_csv <- "data/ngdbsed/main.csv"

# For soil, use:
# medium <- "soil"
# chem_csv <- "data/ngdbsoil/bestvalue.csv"
# main_csv <- "data/ngdbsoil/main.csv"

join_fun <- if (identical(medium, "sed")) join_sed_main else join_soil_main
value_col <- paste0(element, "_mgkg")
state_code <- norm_state(state)

element_data <- build_ag_from_csv(
  chem_csv,
  element = element,
  to_ugkg = TRUE
) %>%
  add_substitution() %>%
  transmute(
    lab_id,
    !!value_col := val_sub_mgkg,
    detected_val,
    left_cen,
    left_limit,
    right_cen,
    right_limit,
    value_unit
  )

main_data <- read_clean_main(main_csv)

element_main <- join_fun(
  element_data,
  main_data,
  extra_cols = c("state", "country", "lat_wgs84", "long_wgs84")
)

hotspot_df <- element_main %>%
  filter(
    state == state_code,
    is.finite(lat_wgs84),
    is.finite(long_wgs84),
    .data[[value_col]] >= threshold_mgkg
  ) %>%
  arrange(desc(.data[[value_col]]))

base_name <- paste0(
  medium, "_", element, "_", state_code, "_hotspots_ge",
  gsub("\\.", "p", threshold_mgkg), "mgkg"
)

write_csv(hotspot_df, file.path("output", paste0(base_name, ".csv")))

hotspot_map <- state_hotspot_map(
  element_main,
  state = state_code,
  value_col = value_col,
  thr_min = threshold_mgkg,
  title = paste0(
    state_code, " ", toupper(medium), " ", element,
    " hotspots >= ", threshold_mgkg, " mg/kg"
  )
)

ggsave(
  filename = file.path("output", paste0(base_name, ".png")),
  plot = hotspot_map,
  width = 7,
  height = 5,
  dpi = 300
)

cat(sprintf(
  "Wrote %d hotspots and map for %s %s in %s: output/%s.csv / .png\n",
  nrow(hotspot_df),
  toupper(medium),
  element,
  state_code,
  base_name
))

print(head(as.data.frame(hotspot_df)))
