#!/usr/bin/env Rscript
# Synthetic software smoke check only; no real geochemical observations are used.
# Run from any working directory with Rscript /path/to/this/script.R.
args <- commandArgs(FALSE)
script <- sub("^--file=", "", args[startsWith(args, "--file=")])
if (length(script) != 1L) stop("Run this check with Rscript")
# Rscript encodes spaces as ~+~ in --file on some R builds. Prefer the literal
# path when it exists so a genuine ~+~ filename is not changed accidentally.
if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
workflow <- file.path(dirname(normalizePath(script, mustWork = TRUE)), "national-state-silver.R")
source(workflow)
for (pkg in c("geoCensor", "sf")) if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg, " first")
if (as.character(utils::packageVersion("geoCensor")) != "1.2.0") stop("Requires geoCensor 1.2.0")
expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  if (!inherits(error, "error") || !grepl(pattern, conditionMessage(error), ignore.case = TRUE)) {
    stop("Expected rejection matching: ", pattern)
  }
}
temporary_root <- Sys.getenv("TMPDIR")
if (!nzchar(temporary_root)) temporary_root <- "/tmp"
root <- tempfile("geocensor-silver-synthetic-", tmpdir = temporary_root)
stopifnot(dir.create(root))

# Hand-counted at a SINGLE fixed origin: A has a double-occupied negative
# cell plus one positive cell; B has one cell. Three cells share total mass.
geometry <- data.frame(component = c("A", "A", "A", "B"),
  x_m = c(-1, -2, 11, -1), y_m = c(-1, -2, 1, -1))
got <- silver_weights(geometry, 10, n_offsets = 1L, anchor = "fixed0")
stopifnot(max(abs(got$weights - c(1/6, 1/6, 1/3, 1/3))) < 1e-14, got$cells_mean == 3)
geometry$component <- "A"
got <- silver_weights(geometry, 10, n_offsets = 10L, anchor = "minimum")$weights
reference <- geoCensor::cell_declustering_weights(geometry, "x_m", "y_m", 10, 10)
stopifnot(max(abs(got - reference)) < 1e-13)

# Eight CA Ag records include one unlocated record. Six AK records are
# uncensored; six TX records are all left-censored. A blank-country Ag record
# must be excluded, while a CA metadata-only record must influence the CRS.
n <- 22L
metadata <- data.frame(lab_id = sprintf("synthetic_%02d", seq_len(n)),
  country = c(rep("US", 20), "", "US"),
  state = c(rep("CA", 8), rep("AK", 6), rep("TX", 6), "CA", "CA"),
  primary_class = "soil",
  lat_wgs84 = c(rep(36, 4), 36.1, 36.2, 36.3, NA, 64 + (0:5)/100, 31 + (0:5)/100, 36, 50),
  long_wgs84 = c(rep(-120, 4), -119.9, -119.8, -119.7, NA, -150 + (0:5)/100, -98 + (0:5)/100, -120, -100))
chemistry <- data.frame(lab_id = metadata$lab_id[1:21], species = "Ag", unit = "ppm",
  amethod = "MS_HF", qvalue = c("-0.01", "-0.01", "0.02", "0.03", "0.04", "0.05", "1.1111", "2.11",
    "0.01", "0.02", "0.03", "0.04", "0.05", "0.06", rep("-0.1", 6), "9999"), summary = "")
chemistry$amethod[3] <- "ES_SQ"
chemistry$summary[8] <- "2.1111,MS_HF"
chemistry$summary[6] <- "0.051111,UNKNOWN"
paths <- character()
for (medium in c("soil", "sediment")) {
  metadata$primary_class <- medium
  for (kind in c("main", "bestvalue")) {
    path <- file.path(root, paste0(medium, "_", kind, ".csv"))
    utils::write.csv(if (kind == "main") metadata else chemistry, path, row.names = FALSE, na = "")
    paths <- c(paths, path)
  }
}
before <- tools::md5sum(paths)
output <- file.path(root, "result")
# No small-offset shortcut: run the published CLI entry function unchanged,
# with its default 100 placements for every cell width.
suppressMessages(silver_run(c(paths, output)))
stopifnot(identical(before, tools::md5sum(paths)), file.exists(file.path(output, "COMPLETED.txt")))
for (medium in c("soil", "sediment")) {
  folder <- file.path(output, medium)
  d <- readRDS(file.path(folder, "analysis.rds"))
  stopifnot(nrow(d) == 20L, identical(d$chemistry_source_record, as.double(1:20)),
    identical(d$metadata_source_record, as.double(1:20)), identical(d$qvalue, chemistry$qvalue[1:20]),
    d$left_cen[1], d$right_cen[7], d$right_cen[8], d$censoring_source[8] == "summary_upper_limit",
    d$summary_recovery_status[8] == "recovered_summary_upper_limit", !d$right_cen[6],
    d$censoring_source[6] == "qvalue", d$summary_recovery_status[6] == "blocked_method_information_missing")
  projection <- utils::read.csv(file.path(folder, "projections.csv"))
  ca_crs <- projection$crs[projection$role == "state" & projection$domain == "CA"]
  stopifnot(length(ca_crs) == 1L, grepl("lat_0=36.0500000000", ca_crs, fixed = TRUE),
    grepl("lon_0=-119.9500000000", ca_crs, fixed = TRUE))
  z <- utils::read.csv(file.path(folder, "estimates.csv"))
  pooled <- subset(z, domain_type == "state" & analysis_role == "ALL_METHODS_REFERENCE")
  stopifnot(nrow(pooled) == 255L, length(unique(pooled$domain)) == 51L,
    all(is.na(z$ci_low)), all(is.na(z$ci_high)),
    all(is.na(pooled$geomean_mgkg[pooled$domain == "HI"])),
    all(pooled$data_status[pooled$domain == "HI"] == "no_reported_Ag"),
    all(pooled$n_metadata[pooled$domain == "HI"] == 0),
    all(is.na(pooled$geomean_mgkg[pooled$domain == "TX"])),
    all(pooled$data_status[pooled$domain == "TX"] == "observed_Ag_records"),
    all(pooled$method[pooled$domain == "TX"] == "insufficient_data"))
  baseline <- subset(z, estimate_basis == "unweighted_located")
  weighted <- subset(z, estimate_basis == "cell_declustered")
  stopifnot(all(weighted$n_total == baseline$n_total[match(weighted$case_id, baseline$case_id)]),
    all(weighted$n_origins == 100L),
    pooled$n_total[pooled$domain == "CA" & pooled$estimate_basis == "unweighted_all"] == 8L,
    pooled$n_total[pooled$domain == "CA" & pooled$estimate_basis == "unweighted_located"] == 7L)
  weights <- utils::read.csv(file.path(folder, "row_weights.csv"))
  totals <- stats::aggregate(weight ~ case_id + cell_size_m, weights, sum)
  stopifnot(all(abs(totals$weight - 1) < 1e-12), all(weights$chemistry_source_record != 8L))
}
expect_error(silver_run(c(paths, output)), "already exists")
bad <- file.path(root, "malformed_unselected.csv")
writeLines(c("lab_id,species,qvalue,unit,amethod,summary",
  '"keep","Ag","0.1","ppm","MS_HF",""',
  '"unselected","Cu","0.1","ppm","MS_HF","",EXTRA'), bad)
expect_error(silver_scan(bad, c("lab_id", "species"), function(x) x$species == "Ag"), "parsing|CSV")
stopifnot(identical(before, tools::md5sum(paths)))
cat("PASS: synthetic geometry, full workflow, source/censoring retention, missing-data distinctions,\n",
  "same-subset comparisons, 100 placements, CRS derivation, immutable inputs and rejection checks.\n",
  "Synthetic check outputs: ", root, "\n", sep = "")
