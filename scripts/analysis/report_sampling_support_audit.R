#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
f<-'output/usgs_sampling_support_audit_2026-10-31'
stopifnot(file.exists(file.path(f,'COMPLETED.txt')))
x<-fread(file.path(f,'overview.csv'))
fmt<-function(z)format(z,big.mark=',',trim=TRUE,scientific=FALSE)
lines<-c('# USGS sampling-support and analytical-method audit','',
 'This descriptive audit examines U.S.-identified soil and sediment samples before defining comparable populations. No chemistry estimate was recalculated and no samples were deleted or reclassified.','',
 '| Medium | Metadata samples | Exact coordinate pairs | Pairs reused by multiple samples | Records at reused pairs | Pairs with differing recorded support |',
 '| --- | ---: | ---: | ---: | ---: | ---: |')
for(i in seq_len(nrow(x))){z<-x[i];lines<-c(lines,sprintf('| %s | %s | %s | %s | %s | %s |',z$medium,fmt(z$n_metadata_samples),fmt(z$n_exact_locations),fmt(z$n_repeated_locations),fmt(z$n_records_at_repeated_locations),fmt(z$n_mixed_support_locations)))}
lines<-c(lines,'','Exact coordinate reuse may reflect true co-location, coordinate rounding, or reused locations. It does not prove sample duplication. Differing support means more than one nonmissing recorded value in at least one inspected field; spelling or unit differences may contribute.','',
 '## Metadata coverage','',
 '| Medium | Field | Records with a value | Percent recorded |','| --- | --- | ---: | ---: |')
for(m in c('soil','sed')){
 z<-fread(file.path(f,m,'support_field_coverage.csv'))
 stopifnot(all(z$n_samples==x[medium==m]$n_metadata_samples),all(z$n_recorded+z$n_missing==z$n_samples))
 for(i in seq_len(nrow(z)))lines<-c(lines,sprintf('| %s | %s | %s | %.2f%% |',m,z$field[i],fmt(z$n_recorded[i]),100*z$fraction_recorded[i]))
 locations<-fread(file.path(f,m,'exact_locations.csv'))
 stopifnot(sum(locations$n_samples)==x[medium==m]$n_valid_coordinate_records)
 combinations<-fread(file.path(f,m,'support_combinations_by_state.csv'))
 stopifnot(sum(combinations$n_samples)==x[medium==m]$n_metadata_samples)
}
lines<-c(lines,'','Horizon is available only in the soil metadata layout. Missing metadata means comparability is unknown; an identical or missing label does not establish the same sampling support. The broad support-gap flag is intentionally sensitive, including sparsely populated protocol fields; use the per-field coverage table to see which information is actually missing.','',
 '## Analytical method mixtures','',
 '| Medium | Species/state groups | Groups with multiple recorded method codes | Percent |','| --- | ---: | ---: | ---: |')
for(i in seq_len(nrow(x))){z<-x[i];lines<-c(lines,sprintf('| %s | %s | %s | %.1f%% |',z$medium,fmt(z$n_species_state_groups),fmt(z$n_groups_multiple_methods),100*z$n_groups_multiple_methods/z$n_species_state_groups))}
lines<-c(lines,'',
 'A best-value table has one selected result per sample/species, but the selected analytical method can vary between samples. Method codes were counted without assuming that different digestions, detection limits, or measurement methods are equivalent. A rare method was not automatically treated as lower quality.','',
 '## Consequences for the next analysis','',
 '1. Define the intended medium, horizon/depth, size fraction, and geological domain before a reference population is estimated.',
 '2. For a well-described subset, compare results within recorded methods or method families whose comparability has been verified from laboratory definitions.',
 '3. Report known-support and unknown-support populations separately where scientifically justified; do not assign unknown samples a convenient support category.',
 '4. Treat co-located records as potentially different supports. Declustering can reduce repeated-location influence but does not resolve differences in what was sampled.',
 '5. Keep previous all-sample state estimates explicitly exploratory until these population decisions are reviewed. The audit alone does not justify an automatic exclusion rule.','',
 '## Outputs and reproduction','',
 'Each medium contains support_fields_by_state.csv, support_field_coverage.csv, support_combinations_by_state.csv, exact_locations.csv, location_summary_by_state.csv, methods_by_species_state.csv, method_mixtures.csv, and the serialized audit. overview.csv provides the totals above.',
 'Reproduce with scripts/analysis/audit_usgs_sampling_support.R followed by scripts/analysis/report_sampling_support_audit.R. The public sampling_support_audit() function provides the descriptive metadata audit on user-supplied data.','',
 'Checks passed: metadata combination counts and coordinate-location counts reconcile; recorded plus missing field counts equal the metadata totals.')
writeLines(lines,file.path(f,'README.md'))
cat('Sampling-support output reconciliation passed.\n')
