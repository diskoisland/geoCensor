#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);pkgload::load_all('.',quiet=TRUE)})
setDTthreads(2)
folder<-'output/usgs_sampling_support_audit_2026-10-31'
dir.create(folder,recursive=TRUE,showWarnings=FALSE)
writeLines(capture.output(sessionInfo()),file.path(folder,'session_info.txt'))
overview<-list()
for(medium in c('soil','sed')) {
 cat('Auditing',medium,'metadata\n');flush.console()
 source<-file.path('..','data',paste0('ngdb',medium))
 header<-names(fread(file.path(source,'main.csv'),nrows=0))
 support<-intersect(c('depth','horizon','mesh_size','primary_class','sample_source','method_collected','collection_protocol'),header)
 metadata<-fread(file.path(source,'main.csv'),select=c('lab_id','state','country','lat_wgs84','long_wgs84',support),colClasses=c(lab_id='character'))
 metadata<-metadata[toupper(trimws(country))%in%c('UNITED STATES','US','USA','U.S.')]
 metadata[,state:=toupper(trimws(state))];metadata[is.na(state)|state=='',state:='UNKNOWN_STATE']
 metadata[!is.finite(lat_wgs84)|abs(lat_wgs84)>90,lat_wgs84:=NA_real_]
 metadata[!is.finite(long_wgs84)|abs(long_wgs84)>180,long_wgs84:=NA_real_]
 out<-sampling_support_audit(as.data.frame(metadata),group_cols='state',support_cols=support,x_col='long_wgs84',y_col='lat_wgs84')
 target<-file.path(folder,medium);dir.create(target,showWarnings=FALSE)
 fwrite(as.data.table(out$fields),file.path(target,'support_fields_by_state.csv'))
 fwrite(as.data.table(out$combinations),file.path(target,'support_combinations_by_state.csv'))
 locations<-as.data.table(out$locations)
 fwrite(locations,file.path(target,'exact_locations.csv'))
 location_summary<-locations[,.(n_locations=.N,n_repeated_locations=sum(n_samples>1),
  n_sample_records_at_repeated_locations=sum(n_samples[n_samples>1]),
  n_mixed_recorded_support_locations=sum(mixed_recorded_support),
  n_support_gap_locations=sum(has_support_gaps)),by=state]
 fwrite(location_summary,file.path(target,'location_summary_by_state.csv'))
 # Method codes are described; total/partial digestion equivalence is not assumed.
 cat('Auditing',medium,'analytical method mixtures\n');flush.console()
 chemistry<-fread(file.path(source,'bestvalue.csv'),select=c('lab_id','species','amethod','unit'),colClasses=c(lab_id='character'))
 chemistry<-chemistry[lab_id%in%metadata$lab_id]
 chemistry[,state:=metadata$state[match(lab_id,metadata$lab_id)]]
 chemistry[,amethod:=trimws(amethod)]
 methods<-chemistry[,.N,by=.(species,state,amethod,unit)]
 setnames(methods,'N','n_records')
 fwrite(methods,file.path(target,'methods_by_species_state.csv'))
 method_summary<-methods[,.(n_records=sum(n_records),
  n_recorded_methods=uniqueN(amethod[!is.na(amethod)&amethod!='']),
  n_missing_method_records=sum(n_records[is.na(amethod)|amethod=='']),
  largest_method_unit_share=max(n_records)/sum(n_records)),by=.(species,state)]
 fwrite(method_summary,file.path(target,'method_mixtures.csv'))
 field_totals<-as.data.table(out$fields)[,.(n_samples=sum(n_samples),n_recorded=sum(n_recorded),n_missing=sum(n_missing)),by=field]
 field_totals[,fraction_recorded:=n_recorded/n_samples]
 fwrite(field_totals,file.path(target,'support_field_coverage.csv'))
 overview[[medium]]<-data.table(medium=medium,n_metadata_samples=nrow(metadata),n_valid_coordinate_records=sum(locations$n_samples),
  n_exact_locations=nrow(locations),n_repeated_locations=sum(locations$n_samples>1),
  n_records_at_repeated_locations=sum(locations$n_samples[locations$n_samples>1]),
  n_mixed_support_locations=sum(locations$mixed_recorded_support),n_locations_with_support_gaps=sum(locations$has_support_gaps),
  n_species_state_groups=nrow(method_summary),n_groups_multiple_methods=sum(method_summary$n_recorded_methods>1))
 saveRDS(out,file.path(target,'support_audit.rds'))
 cat('Completed',medium,'audit\n');flush.console()
 rm(chemistry,metadata,out,methods);gc()
}
fwrite(rbindlist(overview),file.path(folder,'overview.csv'))
writeLines(c(
 'Descriptive audit of U.S.-identified metadata and best-value chemistry; original files read only.',
 'Coordinates are compared exactly, without rounding or proximity clustering. Reused coordinates are not necessarily the same sample support.',
 'Mixed recorded support means more than one nonmissing value for at least one inspected metadata field at the location.',
 'Missing metadata independently flags support gaps; identical or missing labels do not establish comparability.',
 'Raw labels are trimmed but not harmonized; apparent differences can reflect spelling, units, protocols, or real support differences.',
 'Method codes are counted, not assumed interchangeable. Different digestion/measurement methods need domain-specific review.',
 'Largest method/unit share is descriptive and is not a quality score. No observations were dropped for having rare methods.',
 'No pooled estimate or uncertainty interval was recalculated. See per-medium outputs for audit counts.'
),file.path(folder,'analysis_notes.txt'))
writeLines('Sampling-support and method audit completed.',file.path(folder,'COMPLETED.txt'))
print(rbindlist(overview))
