#!/usr/bin/env Rscript
# Export text fields and separate source-record ordinals; no statistical analysis.
args <- commandArgs(trailingOnly = TRUE)
elements <- NULL
option <- grepl("^--elements=", args)
if (sum(option) > 1L) stop("Supply --elements only once", call. = FALSE)
if (any(option)) {
  value <- sub("^--elements=", "", args[option])
  if (!nzchar(value) || grepl("(^,|,,|,$)", value)) stop("--elements requires nonblank comma-separated labels", call. = FALSE)
  elements <- strsplit(value, ",", fixed = TRUE)[[1]]
  args <- args[!option]
}
chunk_size <- 10000
option <- grepl("^--chunk-size=", args)
if (sum(option) > 1L) stop("Supply --chunk-size only once", call. = FALSE)
if (any(option)) {
  value <- sub("^--chunk-size=", "", args[option])
  if (!grepl("^[0-9]+$", value)) stop("--chunk-size requires a positive integer", call. = FALSE)
  chunk_size <- as.numeric(value)
  args <- args[!option]
}
if (length(args) != 2L) stop(
  "Usage: Rscript export_chemistry.R chemistry.csv new_output_dir [--elements=Ag,Pb] [--chunk-size=10000]",
  call. = FALSE)
suppressPackageStartupMessages(library(geoCensor))
helper <- system.file("examples", "workflow_requirements.R", package = "geoCensor")
if (!nzchar(helper)) stop("Install the final geoCensor 1.2.0 release and restart R; this installation lacks workflow support.", call. = FALSE)
sys.source(helper, envir = environment())
geocensor_require_workflow("export")
x <- export_usgs_chemistry(args[1], args[2], elements = elements, chunk_size = chunk_size)
message("Exported ", x$n_selected, " of ", x$n_input, " data records to ", x$output_dir)
