#!/usr/bin/env Rscript

# Convert the Great Salt Lake / Wasatch Front dust release to NGDB-like
# chemistry + main tables, then run geoCensor's censored lognormal MLE for Ag.
#
# Run from the geoCensor repository root:
#   Rscript inst/examples/utah_dust_ag_mle.R
#
# The default input path assumes the ScienceBase files are stored at:
#   ../../manuscript/Utahdust/SLVDust_DataRelease.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

load_geocensor <- function() {
  if (requireNamespace("geoCensor", quietly = TRUE)) {
    suppressPackageStartupMessages(library(geoCensor))
    return(invisible(TRUE))
  }

  source_files <- c(
    "R/parse_usgs_q.r",
    "R/censoring.R",
    "R/build_ag_from_chem.r",
    "R/build_element.r"
  )

  missing <- source_files[!file.exists(source_files)]
  if (length(missing) > 0) {
    stop(
      "geoCensor is not installed and source files are not available: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  lapply(source_files, source)
  invisible(TRUE)
}

make_lab_id <- function(x) {
  id <- ifelse(!is.na(x$Geochemical_ID) & x$Geochemical_ID != "",
               x$Geochemical_ID,
               x$Physical_Properties_ID)
  id <- ifelse(!is.na(id) & id != "",
               id,
               paste(x$Sample_ID, x$Size_Max, x$Size_Min, sep = "_"))
  make.unique(as.character(id), sep = "_rep")
}

size_fraction <- function(size_max, size_min) {
  max_chr <- as.character(size_max)
  min_chr <- as.character(size_min)

  dplyr::case_when(
    max_chr == "63" & min_chr == "0" ~ "<63 um",
    max_chr == "125" & min_chr == "63" ~ "63-125 um",
    max_chr == "250" & min_chr == "125" ~ "125-250 um",
    max_chr == "500" & min_chr == "250" ~ "250-500 um",
    max_chr == "500" & min_chr == "0" ~ "<500 um",
    max_chr == "bulk" & min_chr == "0" ~ "bulk",
    is.na(size_max) & min_chr == "63" ~ ">63 um",
    is.na(size_max) & min_chr == "250" ~ ">250 um",
    TRUE ~ paste0(min_chr, "-", max_chr, " um")
  )
}

screen_to_qualifier <- function(x) {
  dplyr::case_when(
    str_to_lower(str_trim(as.character(x))) == "bdl" ~ "L",
    TRUE ~ NA_character_
  )
}

convert_utah_dust_to_ngdb <- function(raw) {
  raw$lab_id <- make_lab_id(raw)

  raw <- raw %>%
    mutate(
      sample_id = Sample_ID,
      state = "UT",
      country = "US",
      latitude = suppressWarnings(as.numeric(Latitude)),
      longitude = suppressWarnings(as.numeric(Longitude)),
      sample_type = str_to_lower(as.character(Sample_Type)),
      size_fraction = size_fraction(Size_Max, Size_Min),
      source_dataset = "USGS ScienceBase 10.5066/P996NOES"
    )

  metadata_cols <- c(
    "lab_id", "sample_id", "Site_Name", "Site_ID", "state", "country",
    "latitude", "longitude", "sample_type", "Sample_Height",
    "Start_Date", "End_Date", "Physical_Properties_ID", "Size_Max",
    "Size_Min", "size_fraction", "Mass", "Geochemical_ID", "Element_ID",
    "Isotope_ID", "87Sr_86Sr", "87Sr_86Sr_error", "Sample_Comments",
    "source_dataset"
  )

  main <- raw %>%
    select(any_of(metadata_cols))

  skip_cols <- c(
    names(main),
    "Site_Name", "Site_ID", "Sample_ID", "Sample_Type", "Sample_Height",
    "Start_Date", "End_Date", "Physical_Properties_ID", "Size_Max",
    "Size_Min", "Mass", "Geochemical_ID", "Element_ID", "Isotope_ID",
    "87Sr_86Sr", "87Sr_86Sr_error", "Sample_Comments",
    "Latitude", "Longitude",
    "Li_SD", "Be_SD", "B_SD", "Na_SD", "Mg_SD", "Al_SD", "K_SD",
    "Ca_SD", "V_SD", "Cr_SD", "Mn_SD", "Fe_SD", "Co_SD", "Ni_SD",
    "Cu_SD", "Zn_SD", "As_SD", "Se_SD", "Rb_SD", "Sr_SD", "Y_SD",
    "Mo_SD", "Ag_SD", "Cd_SD", "Sb_SD", "Cs_SD", "Ba_SD", "La_SD",
    "Tl_SD", "Pb_SD", "Ce_SD", "Pr_SD", "Nd_SD", "Sm_SD", "Eu_SD",
    "Gd_SD", "Tb_SD", "Dy_SD", "Ho_SD", "Er_SD", "Yb_SD", "Lu_SD",
    "Th_SD", "U_SD",
    grep("_screen$", names(raw), value = TRUE)
  )

  element_cols <- setdiff(names(raw), skip_cols)
  element_cols <- element_cols[vapply(raw[element_cols], function(x) {
    any(!is.na(suppressWarnings(parse_number(as.character(x)))))
  }, logical(1))]

  chem <- raw %>%
    select(lab_id, all_of(element_cols), ends_with("_screen")) %>%
    pivot_longer(
      cols = all_of(element_cols),
      names_to = "species",
      values_to = "qvalue",
      values_transform = list(qvalue = as.character)
    ) %>%
    mutate(
      qvalue = as.character(qvalue),
      screen_col = paste0(species, "_screen")
    ) %>%
    rowwise() %>%
    mutate(
      qualifier = if (screen_col %in% names(raw)) {
        screen_to_qualifier(raw[[screen_col]][match(lab_id, raw$lab_id)])
      } else {
        NA_character_
      }
    ) %>%
    ungroup() %>%
    mutate(
      unit = "ppm",
      lab_id = as.character(lab_id),
      species = as.character(species)
    ) %>%
    filter(!is.na(parse_number(qvalue))) %>%
    select(lab_id, species, qvalue, qualifier, unit)

  list(main = main, chem = chem)
}

fit_by_group <- function(df, group_cols) {
  df %>%
    group_by(across(all_of(group_cols))) %>%
    group_modify(~ fit_cen_lognorm(.x)) %>%
    ungroup() %>%
    mutate(
      geomean_mgkg = geomean,
      ci_low_mgkg = ci_low,
      ci_high_mgkg = ci_high,
      pct_censored = if_else(n_total > 0, 100 * n_left / n_total, NA_real_)
    ) %>%
    select(all_of(group_cols), n_total, n_left, n_right, pct_censored,
           geomean_mgkg, ci_low_mgkg, ci_high_mgkg, method)
}

load_geocensor()

args <- commandArgs(trailingOnly = TRUE)
input_csv <- if (length(args) >= 1) {
  args[[1]]
} else {
  "../../manuscript/Utahdust/SLVDust_DataRelease.csv"
}
output_dir <- if (length(args) >= 2) args[[2]] else "output/utahdust_ngdb"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

raw <- read_csv(input_csv, show_col_types = FALSE, na = c("", "NA", "NaN"))
converted <- convert_utah_dust_to_ngdb(raw)

main_out <- file.path(output_dir, "utahdust_main_ngdb_like.csv")
chem_out <- file.path(output_dir, "utahdust_bestvalue_ngdb_like.csv")
ag_out <- file.path(output_dir, "utahdust_ag_censoring_mgkg.csv")
pooled_out <- file.path(output_dir, "utahdust_ag_mle_pooled_mgkg.csv")
group_out <- file.path(output_dir, "utahdust_ag_mle_by_sample_type_fraction_mgkg.csv")

write_csv(converted$main, main_out, na = "")
write_csv(converted$chem, chem_out, na = "")

ag <- build_element(converted$chem, element = "Ag", to_ugkg = FALSE) %>%
  left_join(
    converted$main %>%
      select(lab_id, state, country, latitude, longitude, sample_type,
             size_fraction, Site_Name, sample_id),
    by = "lab_id"
  )

write_csv(ag, ag_out, na = "")

pooled <- fit_cen_lognorm(ag) %>%
  mutate(
    geomean_mgkg = geomean,
    ci_low_mgkg = ci_low,
    ci_high_mgkg = ci_high,
    pct_censored = if_else(n_total > 0, 100 * n_left / n_total, NA_real_)
  ) %>%
  select(n_total, n_left, n_right, pct_censored,
         geomean_mgkg, ci_low_mgkg, ci_high_mgkg, method)

by_sample_fraction <- fit_by_group(
  ag %>% filter(!is.na(sample_type), !is.na(size_fraction)),
  c("sample_type", "size_fraction")
)

write_csv(pooled, pooled_out, na = "")
write_csv(by_sample_fraction, group_out, na = "")

cat("\nWrote NGDB-like Utah dust files:\n")
cat("  main:      ", main_out, "\n", sep = "")
cat("  bestvalue: ", chem_out, "\n", sep = "")
cat("  Ag table:  ", ag_out, "\n", sep = "")
cat("  pooled MLE: ", pooled_out, "\n", sep = "")
cat("  group MLE: ", group_out, "\n\n", sep = "")

cat("Pooled Ag censored lognormal MLE, mg/kg:\n")
print(as.data.frame(pooled), row.names = FALSE)

cat("\nAg MLE by sample type and fraction, mg/kg:\n")
print(as.data.frame(by_sample_fraction), row.names = FALSE)
