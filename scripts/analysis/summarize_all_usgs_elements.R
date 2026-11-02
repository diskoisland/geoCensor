#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
root <- file.path('output','usgs_all_elements_2026-10-31')
stopifnot(file.exists(file.path(root,'COMPLETED.txt')))
combined <- list()
summary <- list()
for (medium in c('soil','sed')) {
  folder <- file.path(root,medium)
  pooled <- fread(file.path(folder,'pooled_estimates.csv'))
  state <- fread(file.path(folder,'state_estimates.csv'))
  spatial <- fread(file.path(folder,'state_declustering.csv'))
  validation <- fread(file.path(folder,'pooled_validation.csv'))
  status <- fread(file.path(folder,'run_status.csv'))
  inventory <- fread(file.path(folder,'input_inventory.csv'))
  stopifnot(nrow(pooled)==inventory$n_species, nrow(status)==inventory$n_species)
  stopifnot(all(pooled$n_input==pooled$n_total+pooled$n_excluded))
  totals <- state[,.(n_input=sum(n_input)),by=element]
  checked <- merge(pooled[,.(element,n_input)],totals,by='element')
  stopifnot(all(checked$n_input.x==checked$n_input.y))
  counts <- validation[,.(n_input=sum(n_observations)),by=element]
  checked <- merge(pooled[,.(element,n_input)],counts,by='element')
  stopifnot(all(checked$n_input.x==checked$n_input.y))
  stopifnot(nrow(spatial)==3*nrow(state), all(spatial$n_input==spatial$n_total+spatial$n_excluded))
  stopifnot(all(is.na(spatial$ci_low)),all(is.na(spatial$ci_high)))
  spatial[,medium:=medium]
  spatial[,ratio_to_matched_unweighted:=fifelse(is.finite(geomean)&is.finite(unweighted_geomean)&unweighted_geomean>0,geomean/unweighted_geomean,NA_real_)]
  spatial[,n_exact:=n_total-n_left-n_right]
  combined[[medium]] <- spatial
  summary[[medium]] <- data.table(
    medium=medium, input_records=inventory$input_chemistry_rows,
    us_records=inventory$us_chemistry_rows, us_metadata_samples=inventory$us_metadata_rows,
    reported_species=nrow(pooled), pooled_estimates=sum(is.finite(pooled$geomean)),
    state_species_rows=nrow(state), state_estimates=sum(is.finite(state$geomean)),
    spatial_fits_attempted=nrow(spatial), spatial_estimates=sum(is.finite(spatial$geomean)),
    recovered_upper_limits=sum(status$n_upper_limits_recovered),
    invalid_chemistry_records=sum(pooled$n_excluded)
  )
}
summary <- rbindlist(summary)
fwrite(summary,file.path(root,'run_overview.csv'))
fwrite(rbindlist(combined),file.path(root,'declustering_comparison.csv'))
format_number <- function(x) format(x,big.mark=',',scientific=FALSE,trim=TRUE)
lines <- c('# Full USGS all-species analysis', '',
  'Completed using the local geoCensor source. Soil and sediment were processed separately.',
  'The analysis includes records explicitly marked United States, including reported state/territory groups.', '',
  '| Medium | U.S. chemistry records | Reported species | Pooled estimates | State/territory estimates | Weighted spatial estimates |',
  '| --- | ---: | ---: | ---: | ---: | ---: |')
for (i in seq_len(nrow(summary))) {
  x<-summary[i]
  lines<-c(lines,sprintf('| %s | %s | %s | %s | %s | %s |',x$medium,format_number(x$us_records),x$reported_species,x$pooled_estimates,format_number(x$state_estimates),format_number(x$spatial_estimates)))
}
lines<-c(lines,'',
  'Spatial sensitivity uses 5, 10, and 20 km square cells with 100 grid placements per size, separately for each species and state/territory. Local equal-area projection definitions are saved with each medium.',
  'These are exploratory sensitivity scales, not selected optimal geological spacings. Weighted estimates use the same usable chemistry/coordinate subset as their unweighted comparison. They do not represent unsampled terrain.', '',
  paste('Recovered upper-censoring markers from rounded source exports:',format_number(sum(summary$recovered_upper_limits)),'. Original files were not changed.'),
  'The first summary value was used only when it contained an explicit upper-censoring suffix and agreed with rounded qvalue.', '',
  'All input species are retained, including compounds and ancillary species. Records reported in grams rather than concentration units cannot be modeled as concentration and remain excluded with explicit counts.',
  'Insufficient-data and failed-fit rows remain in exported results. Weighted confidence intervals are omitted pending spatial uncertainty validation.',
  'Pooled confidence intervals are model-based and do not account for spatial dependence. Distinct samples from different horizons or size fractions remain as supplied.', '',
  'Files: run_overview.csv; declustering_comparison.csv; soil/ and sed/ contain pooled and state estimates, distribution comparisons, validation summaries, spatial sensitivity, join audits, projections, checkpoints, and serialized results.',
  'Reproduce with scripts/analysis/run_all_usgs_elements.R; rerunning resumes completed checkpoints.',
  'Analysis checks passed: input/validation/grouped counts reconcile; every state/species has three spatial scales; weighted confidence intervals are missing as intended.')
writeLines(lines,file.path(root,'README.md'))
print(summary)
cat('Output validation passed.\n')
