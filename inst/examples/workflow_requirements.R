# Shared by the installed command examples, not an exported package API.
geocensor_require_workflow <- function(workflow,
    version = utils::packageVersion("geoCensor"), namespace = asNamespace("geoCensor")) {
  requirements <- list(
    all = list(read_usgs_chemistry = "chem_csv", analyze_elements = c(
      "chem_df", "main_df", "cell_sizes", "stratify_methods",
      "audit_method_overlap", "audit_source_precision")),
    single = list(read_usgs_chemistry = "chem_csv",
      read_usgs_chemistry_select = c("chem_csv", "elements", "chunk_size"),
      analyze_elements = c("chem_df", "stratify_methods", "audit_source_precision")),
    audit = list(audit_sampling_file = c("input_file", "output_dir", "include_memberships")),
    export = list(export_usgs_chemistry = c("chem_csv", "output_dir", "elements", "chunk_size")),
    illinois = list(read_usgs_chemistry = "chem_csv", analyze_elements = c(
      "chem_df", "main_df", "cell_sizes", "stratify_methods",
      "audit_method_overlap", "audit_source_precision"),
      sampling_support_audit = "include_memberships")
  )
  if (!workflow %in% names(requirements)) stop("Unknown example workflow", call. = FALSE)
  needed <- requirements[[workflow]]
  missing <- names(needed)[!vapply(names(needed), function(name) {
    if (!exists(name, envir = namespace, inherits = FALSE)) return(FALSE)
    fun <- get(name, envir = namespace, inherits = FALSE)
    exported <- !isNamespace(namespace) || name %in% getNamespaceExports(namespace)
    is.function(fun) && exported && all(needed[[name]] %in% names(formals(fun)))
  }, logical(1))]
  if (utils::compareVersion(as.character(version), "1.2.0") < 0 || length(missing)) {
    stop("This command needs the final geoCensor 1.2.0 release or a compatible newer build. ",
      "An earlier local 1.2.0 candidate can have the same version but lack these features. ",
      if (length(missing)) paste0("Missing workflow capabilities: ", paste(missing, collapse = ", "), ". "),
      "Install the checked release archive and restart R.", call. = FALSE)
  }
  invisible(TRUE)
}
