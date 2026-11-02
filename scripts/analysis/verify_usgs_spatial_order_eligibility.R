#!/usr/bin/env Rscript
# Independent record-count reconciliation and descriptive support review only.
suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1L] else 'output/usgs_spatial_order_eligibility_2026-11-01'
stopifnot(file.exists(file.path(out,'COMPLETED.txt')))
manifest <- fread(file.path(out,'source_input_manifest.csv'))
# File digests are verified by the final SHA-256 receipt. Only audited-readable
# metadata is opened; a timeout is never turned into an empty data table.
availability <- fread(file.path(out,'input_availability.csv'))
summary <- fread(file.path(out,'metadata_summary.csv'),na.strings='NA')
fields <- fread(file.path(out,'metadata_field_coverage.csv'))
multiplicity <- fread(file.path(out,'coordinate_multiplicity.csv'))
checks <- list(); support <- list(); labels <- list(); k <- 0L
for(audit_medium in summary[status=='available',medium]) {
 path <- file.path('..','data',paste0('ngdb',audit_medium),'main.csv')
 field_set <- c('lab_id','field_id','job_id','date_collect','date_submitted','project_name',
   'coord_prec_est','depth','mesh_size','prep','collection_protocol','previous_lab_id')
 if(audit_medium=='soil') field_set <- c(field_set,'horizon','sample_zone')
 d <- fread(path,select=c('country','lat_wgs84','long_wgs84',field_set),
   colClasses=setNames(rep('character',length(field_set)),field_set),na.strings=NULL)
 original_n <- nrow(d)
 d <- d[toupper(trimws(country))%in%c('US','USA','UNITED STATES','U.S.')]
 old <- summary[summary$medium==audit_medium]
 stopifnot(nrow(old)==1L,original_n==old$input_rows,nrow(d)==old$us_rows,
   !anyDuplicated(d$lab_id))
 valid <- is.finite(d$lat_wgs84)&abs(d$lat_wgs84)<=90&
   is.finite(d$long_wgs84)&abs(d$long_wgs84)<=180
 stopifnot(sum(valid)==old$valid_coordinate_rows)
 g <- d[valid,.(records=.N),by=.(lat_wgs84,long_wgs84)]
 stopifnot(nrow(g)==old$exact_coordinate_groups,
   sum(g$records>1)==old$repeated_coordinate_groups,
   sum(g$records[g$records>1])==old$records_at_repeated_coordinates)
 observed <- g[,.N,by=records][order(records)]
 expected <- multiplicity[multiplicity$medium==audit_medium][order(records_at_coordinate)]
 stopifnot(identical(observed$records,expected$records_at_coordinate),
   identical(observed$N,expected$n_coordinate_groups),sum(observed$records*observed$N)==sum(valid))
 for(audit_field in field_set) {
   values <- trimws(d[[audit_field]])
   recorded <- !is.na(values)&nzchar(values)
   old_field <- fields[fields$medium==audit_medium&fields$field==audit_field]
   if(nrow(old_field)) stopifnot(sum(recorded)==old_field$n_nonblank,
     abs(mean(recorded)-old_field$fraction_nonblank)<1e-14)
   tab <- data.table(value=values)[,.N,by=value][order(-N,value)]
   labels[[paste(audit_medium,audit_field)]] <- data.table(medium=audit_medium,field=audit_field,
     raw_trimmed_label=tab$value,n_records=tab$N)
   z <- data.table(lat=d$lat_wgs84[valid],lon=d$long_wgs84[valid],value=values[valid])
   at <- z[,.(n=.N,n_nonblank=sum(!is.na(value)&nzchar(value)),
     n_distinct_nonblank=uniqueN(value[!is.na(value)&nzchar(value)])),by=.(lat,lon)]
   repeated <- at[n>1]
   support[[paste(audit_medium,audit_field)]] <- data.table(medium=audit_medium,field=audit_field,
     n_us_records=nrow(d),n_nonblank=sum(recorded),fraction_nonblank=mean(recorded),
     n_repeated_coordinate_groups=nrow(repeated),
     repeated_groups_multiple_recorded_labels=sum(repeated$n_distinct_nonblank>1),
     repeated_groups_any_blank=sum(repeated$n_nonblank<repeated$n),
     repeated_groups_all_blank=sum(repeated$n_nonblank==0),
     interpretation='Exact coordinate groups; raw trimmed labels; missing is not agreement; no physical-site claim')
 }
 k <- k+1L
 checks[[k]] <- data.table(medium=audit_medium,input_rows=original_n,us_rows=nrow(d),
   exact_coordinate_groups=nrow(g),metadata_counts_reconciled=TRUE,
   all_reported_field_counts_reconciled=TRUE,coordinate_multiplicities_reconciled=TRUE,
   lab_id_unique_in_us_metadata=TRUE)
}
fwrite(rbindlist(checks),file.path(out,'independent_metadata_verification.csv'))
fwrite(rbindlist(support),file.path(out,'repeated_coordinate_support_review.csv'))
fwrite(rbindlist(labels),file.path(out,'metadata_raw_labels.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'verification_sessionInfo.txt'))
invisible(file.copy('scripts/analysis/verify_usgs_spatial_order_eligibility.R',
  file.path(out,'source_snapshot'),overwrite=FALSE))
writeLines(c('Independent metadata counts, field coverage and every coordinate multiplicity reconciled.',
 'No sample IDs duplicated in the audited U.S. metadata.',
 'Support contrasts and blanks are descriptive; no physical-site map or eligibility claim.',
 'Unreadable media and all p-values remain NA.'),file.path(out,'METADATA_VERIFIED.txt'))
cat('Independent metadata verification passed for',k,'readable medium/media.\n')
