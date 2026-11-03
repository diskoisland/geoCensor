#!/usr/bin/env Rscript
# Audit metadata separately for each intended dataset/domain.
args <- commandArgs(trailingOnly = TRUE)
usage <- paste("Usage: Rscript audit_sampling_support.R metadata.csv output_dir",
  "[--group-cols=state,...] [--support-cols=depth,...] [--id-col=lab_id]",
  "[--coords=x_col,y_col] [--memberships]", sep = "\n")
if (identical(args, "--help")) { cat(usage, "\n"); quit(status = 0L) }
positional <- character(); options <- list(); seen <- character(); memberships <- FALSE
for (arg in args) {
  if (!startsWith(arg, "--")) { positional <- c(positional, arg); next }
  if (arg == "--memberships") { key <- "memberships"; value <- NULL } else {
    token <- regmatches(arg, regexec("^--([a-z-]+)=(.*)$", arg))[[1L]]
    if (length(token) != 3L || !token[2L] %in% c("group-cols", "support-cols", "id-col", "coords")) {
      stop("Unknown or malformed option: ", arg, call. = FALSE)
    }
    key <- token[2L]; value <- token[3L]
    if (!nzchar(trimws(value))) stop("Empty option: --", key, call. = FALSE)
  }
  if (key %in% seen) stop("Supply --", key, " only once", call. = FALSE)
  seen <- c(seen, key)
  if (key == "memberships") memberships <- TRUE else options[[key]] <- value
}
if (length(positional) != 2L) stop(usage, call. = FALSE)
parse_names <- function(value) {
  if (grepl("(^,|,,|,$)", value)) stop("Empty column name in option", call. = FALSE)
  out <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  if (any(!nzchar(out)) || anyDuplicated(out)) stop("Column names must be nonempty and distinct", call. = FALSE)
  out
}
groups <- if (is.null(options[["group-cols"]])) character() else parse_names(options[["group-cols"]])
support <- if (is.null(options[["support-cols"]])) NULL else parse_names(options[["support-cols"]])
coords <- if (is.null(options$coords)) NULL else parse_names(options$coords)
if (!is.null(coords) && length(coords) != 2L) stop("--coords requires exactly two columns", call. = FALSE)
if (memberships && is.null(coords)) stop("--memberships requires --coords", call. = FALSE)
suppressPackageStartupMessages(library(geoCensor))
if (utils::packageVersion("geoCensor") < "1.2.0.9005") stop("Install geoCensor 1.2.0.9005 or newer for this workflow", call. = FALSE)
audit_sampling_file(positional[1L], positional[2L], group_cols = groups, support_cols = support,
  id_col = if (is.null(options[["id-col"]])) "lab_id" else options[["id-col"]],
  x_col = if (is.null(coords)) NULL else coords[1L],
  y_col = if (is.null(coords)) NULL else coords[2L], include_memberships = memberships)
message("Sampling-support audit saved to ", normalizePath(positional[2L]))
