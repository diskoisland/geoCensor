#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(dplyr);library(readr);library(ggplot2)})
f<-'output/spatial_bootstrap_diagnostics_2026-10-31'
stopifnot(file.exists(file.path(f,'COMPLETED.txt')))
s<-read_csv(file.path(f,'summary.csv'),show_col_types=FALSE)
x<-read_csv(file.path(f,'results.csv'),show_col_types=FALSE)
stopifnot(nrow(s)==44,nrow(x)==36400)
stopifnot(all(s$n_intervals<=s$n_datasets),all(s$coverage>=0 & s$coverage<=1))
s<-s%>%mutate(se_ratio=mean_boot_se/empirical_sd_log)
write_csv(s,file.path(f,'summary.csv'))
labels<-c(uncensored_regular='Uncensored; regular density',censored_regular='Censored; regular density',
 censored_uneven='Censored; uneven density',correlated_blocks_1='Correlated; smaller blocks',
 correlated_blocks_2='Correlated; larger blocks')
a<-s%>%filter(phase=='confirmation')
comparison<-a%>%filter(method%in%c('percentile_log','basic_log','normal_log_centered','known_variance_z','complete_data_oracle'))%>%
 mutate(label=factor(labels[scenario],levels=unname(labels)),
 method_label=recode(method,percentile_log='Percentile',basic_log='Basic (log scale)',normal_log_centered='Normal (log scale)',
 known_variance_z='Complete-data oracle',complete_data_oracle='Complete-data oracle'))
p<-ggplot(comparison,aes(x=label,y=coverage,color=method_label))+
 geom_hline(yintercept=.95,linetype='dashed',color='grey40')+
 geom_errorbar(aes(ymin=mc_low,ymax=mc_high),width=.15,position=position_dodge(.65))+
 geom_point(position=position_dodge(.65),size=2)+coord_flip()+
 scale_y_continuous(limits=c(0,1),labels=function(z)paste0(round(z*100),'%'))+
 labs(title='Independent confirmation: coverage of nominal 95% intervals',
 subtitle='Oracle intervals use uncensored measurements and known covariance; they are diagnostic references',
 x=NULL,y='Coverage of the superpopulation geometric mean',color=NULL,
 caption='Uncensored: 3,000 datasets; censored: 200 per setting. Bars: 95% binomial Monte Carlo intervals.\nAll 64 sites observed. Results do not validate irregular or partly unsampled NGDB domains.')+
 theme_minimal(base_size=11)+theme(legend.position='bottom',panel.grid.major.y=element_blank())
ggsave(file.path(f,'confirmation_coverage.png'),p,width=12,height=6,dpi=160)
lines<-c('# Bootstrap benchmark and diagnostic study','',
 'The independent confirmation phase was run with disjoint seeds after all scenarios and candidate intervals were fixed. No candidate was selected or retuned using development results.','',
 '## Estimation target','',
 'The primary target is exp(E(log concentration)) = 50 for a repeatedly generated population. A second, explicitly different target is exp(mu + mean(realized site effects)), the realized equal-site domain mean. Secondary coverage is saved separately; the target was not changed to improve apparent coverage. Neither target establishes concentrations in unobserved terrain.','',
 '## Exact benchmark and implementation checks','',
 'A new regression test confirms that the package uncensored point estimate and every bootstrap replicate match resampling of site log means. Thus the basic weighting/multiplicity algorithm agrees with a mathematically transparent reference in this setting.','',
 'The uncensored benchmark uses independent Gaussian site log means with known variance. Analytic t and known-variance normal intervals have exact 95% coverage under that model. There were 2,000 development and 3,000 confirmation datasets, with 999 resamples each.','',
 '| Uncensored confirmation method | Coverage | 95% Monte Carlo interval |',
 '| --- | ---: | --- |')
for(i in which(a$scenario=='uncensored_regular')){
 z<-a[i,];lines<-c(lines,sprintf('| %s | %.2f%% | %.2f–%.2f%% |',z$method,100*z$coverage,100*z$mc_low,100*z$mc_high))
}
lines<-c(lines,'','## Censoring, density, and correlation diagnostics','',
 'Each of four settings used 200 development and 200 confirmation datasets, with 199 actual package bootstrap refits per dataset. Matched scenarios share underlying measurements where possible. The larger-block correlated setting differs only in block width.','',
 '| Confirmation setting | Method | Coverage | Bootstrap/oracle SE divided by empirical SD |',
 '| --- | --- | ---: | ---: |')
for(i in which(a$scenario!='uncensored_regular')){
 z<-a[i,];lines<-c(lines,sprintf('| %s | %s | %.1f%% | %.2f |',labels[z$scenario],z$method,100*z$coverage,z$se_ratio))
}
lines<-c(lines,'','![Confirmation coverage](confirmation_coverage.png)','',
 '## Interpretation and decision','',
 'Confirmation findings: the analytic t interval covered 95.13% and the studentized bootstrap 95.00% of 3,000 uncensored datasets. The current percentile method covered 94.57%. The original 88% result did not persist in this larger simple benchmark.', '',
 'For censored independent sites, percentile coverage was 96.0% with regular density and 96.5% with uneven density (200 datasets each). Censoring and uneven density alone did not reproduce the major shortfall in these settings.', '',
 'With spatial correlation, percentile coverage fell to 77% and 85% for smaller and larger blocks. Basic intervals gave 79% and 83.5%; centered normal intervals gave 79.5% and 84.5%. Changing the interval formula did not resolve the spatial design problem.', '',
 'The correlated point-estimate log bias was about 0.001, while bootstrap standard errors were only 0.63 and 0.76 times the empirical standard deviation. Known-covariance complete-data oracle intervals covered 97% in both correlated settings. The main unresolved problem here is underestimated dependence-related variation, not a large point-estimate bias.', '',
 'Compare coverage and Monte Carlo uncertainty rather than fit convergence alone. A standard-error ratio below one indicates underestimated variation across generated datasets. Complete-data oracles are diagnostic references only: their hidden measurements and known covariance are unavailable in real censored data.','',
 'The original 100-dataset result of 88% independent-site coverage is not by itself evidence of a coding defect or solely a spatial block-size problem. This larger uncensored benchmark isolates finite-sample/tail effects and simulation variation. Censored cases are evaluated separately.','',
 'Keep the production spatial bootstrap experimental. Changing percentile limits to basic or normal limits alone does not establish a design-appropriate uncertainty method. Next steps should address dependence across blocks and verify performance under irregular locations and genuinely preferential omission of sites. Studentized uncensored results do not validate studentization of censored spatial models.','',
 'All sites were observed on a regular 8x8 grid. Uneven density changes sample counts at those sites. This study does not test missing-area bias or automatic block-size selection. Prior full-data USGS estimates were not recalculated or assigned these intervals.','',
 '## Reproduction','',
 'Run scripts/analysis/diagnose_spatial_bootstrap.R and scripts/analysis/report_spatial_bootstrap_diagnostics.R from the repository root. Use a new output folder if the implementation or study design changes; current runs resume their checkpoints. The software suite and local package check were run separately.','',
 'Files: design_notes.txt, settings.csv, results.csv, summary.csv, confirmation_coverage.png, session_info.txt, and simulation checkpoints.')
writeLines(lines,file.path(f,'README.md'))
print(a%>%select(scenario,method,n_datasets,coverage,mc_low,mc_high,se_ratio))
cat('Diagnostic output checks passed.\n')
