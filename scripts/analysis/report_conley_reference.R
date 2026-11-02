#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(dplyr);library(readr);library(ggplot2)})
f<-'output/conley_hac_validation_2026-10-31'
stopifnot(file.exists(file.path(f,'COMPLETED.txt')))
s<-read_csv(file.path(f,'summary.csv'),show_col_types=FALSE)
x<-read_csv(file.path(f,'results.csv'),show_col_types=FALSE)
stopifnot(nrow(x)==28800,nrow(s)==144,all(s$n_datasets==200),all(s$n_intervals+s$n_invalid==200))
a<-s%>%filter(phase=='confirmation')
a$kernel<-factor(a$kernel,levels=c('bartlett_product','bartlett_radial','epanechnikov'))
base_theme<-theme_minimal(base_size=11)+theme(legend.position='bottom')
p<-ggplot(a,aes(x=bandwidth,y=coverage_counting_invalid_as_miss,color=kernel))+
 geom_hline(yintercept=.95,linetype='dashed',color='grey45')+geom_line()+geom_point()+
 scale_x_continuous(trans='log2',breaks=c(.5,1,2,4,8,16))+
 scale_y_continuous(limits=c(0,1),labels=function(z)paste0(round(z*100),'%'))+
 facet_wrap(~scenario,ncol=2)+labs(title='Conley-style score HAC: independent confirmation coverage',
 x='Fixed bandwidth (synthetic distance units)',y='Coverage; unavailable intervals count as misses',
 caption='200 confirmation datasets per setting; nominal 95% normal intervals. No automatic bandwidth selection.\nEpanechnikov covariances that were not positive semidefinite were withheld, not clipped.')+base_theme
q<-ggplot(a,aes(x=bandwidth,y=mean_se,color=kernel))+
 geom_hline(data=distinct(a,scenario,empirical_sd),aes(yintercept=empirical_sd),linetype='dashed',color='grey35')+
 geom_line()+geom_point()+scale_x_continuous(trans='log2',breaks=c(.5,1,2,4,8,16))+
 facet_wrap(~scenario,ncol=2,scales='free_y')+
 labs(title='Estimated versus actual variation of log geometric means',x='Fixed bandwidth (synthetic distance units)',
 y='Mean estimated standard error',caption='Dashed line: empirical standard deviation across confirmation datasets.\nMean estimated SE uses valid covariance results only; invalid counts are saved separately.')+base_theme
ggsave(file.path(f,'confirmation_coverage.png'),p,width=10,height=7,dpi=160)
ggsave(file.path(f,'bandwidth_standard_errors.png'),q,width=10,height=7,dpi=160)
lines<-c('# Conley-style spatial HAC for weighted censored estimates','',
 'This is a research reference, not a production geoCensor interval API. It applies the weights to the joint location and log-scale likelihood scores, uses inverse weighted observed information, sums co-located scores, and forms the spatial sandwich. It retains the original censored point estimator.','',
 '## Source and compatibility checks','',
 'The user supplied the [SpatialInference bandwidth vignette](https://cran.r-universe.dev/articles/SpatialInference/spatial_HAC.html). Its high-level wrappers target felm linear models, so those wrappers are not a direct replacement for weighted interval-censored survreg. Its kernel and bandwidth diagnostics are useful reference points. This independent reference does not copy external implementation code.',
 'The six checks in output/conley_reference_checks_2026-10-31/checks.csv passed: numerical joint scores; bread agreement with survreg; isolated-location HC0; invariance to the overall weight scale; the uncensored location-mean identity; and wide-bandwidth shrinkage.',
 'A local compatibility check with sandwich 3.1.3 and survival 3.8.6 reproduced a scale-score weighting discrepancy in sandwich::estfun() for a weighted fit. Direct survival matrix derivatives multiplied by case weights for both parameters matched numerical differentiation to about 1e-10. This is a specific, reproducible version/case observation; the verification script saves the comparison.','',
 '## Prespecified study','',
 'Each setting used 200 development plus 200 independently seeded confirmation datasets. Six fixed bandwidths and three fixed kernels were assessed, with no selection or retuning based on coverage. There were 1,600 fitted datasets and 28,800 covariance evaluations. Both censoring directions and unequal within-site sample counts were included.','',
 'The product Bartlett kernel supplies a PSD reference. The radial Bartlett and Epanechnikov alternatives assess the bandwidth/kernel sensitivity discussed in the supplied vignette. Invalid covariance estimates are retained as failures rather than repaired by eigenvalue clipping. Dense calculations are guarded at 2,000 unique sites and are not suitable for a national run.','',
 '## Product Bartlett confirmation results','',
 '| Setting | Bandwidth | Coverage | Mean SE | Empirical SD | Log-scale bias |','| --- | ---: | ---: | ---: | ---: | ---: |')
for(i in which(a$kernel=='bartlett_product')){z<-a[i,];lines<-c(lines,sprintf('| %s | %.1f | %.1f%% | %.4f | %.4f | %.4f |',z$scenario,z$bandwidth,100*z$coverage,z$mean_se,z$empirical_sd,z$bias_log))}
lines<-c(lines,'','![Confirmation coverage](confirmation_coverage.png)','',
 '![Bandwidth versus standard errors](bandwidth_standard_errors.png)','',
 '## Interpretation and decision','',
 'The research calculation explicitly handles censoring, unequal case weights, nuisance scale estimation, and exact-location score aggregation. A passing score identity does not establish finite-sample interval coverage.',
 'In the correlated 64-site confirmation case, the product Bartlett candidate reached about 80% coverage at bandwidth 4, then fell to 57.5% at bandwidth 16. With 256 correlated sites, the same kernel reached 89.5% at bandwidths 4 and 8. Epanechnikov at bandwidth 4 gave 92.5% in the 256-site case, but this is one prespecified candidate result, not a selected or validated production rule. Some wider Epanechnikov covariances were non-PSD.',
 'Excessively wide cutoffs can reduce SEs substantially. A nearly constant kernel approaches the outer product of summed fitted scores, which are near zero. Therefore the widest bandwidth must not be assumed conservative.',
 'The preferential-omission negative control retained substantial point-estimate bias and zero coverage across these candidates. A covariance correction does not by itself repair a biased estimating equation under informative site selection.',
 'Keep spatial HAC experimental. Next assess bandwidth diagnostics on the estimator-specific influence field, domain size and stationarity assumptions, finite-sample calibration, and range-estimation uncertainty. Do not choose the largest SE or the candidate with coverage closest to 95% from this study and call it validated. A separate prespecified follow-up is needed for any data-driven rule.','',
 '## Reproduction','',
 'Run scripts/analysis/verify_conley_score_reference.R, scripts/analysis/validate_conley_reference.R, and scripts/analysis/report_conley_reference.R. The calculation is in scripts/analysis/conley_survreg_reference.R. Use fresh output folders for changed code/design. No full-data USGS estimate or routine confidence interval was replaced.',
 'Saved files: design_notes.txt, settings.csv, results.csv, summary.csv, figures, session_info.txt, and phase checkpoints. Output checks reconcile all dataset/covariance rows and withheld intervals.')
writeLines(lines,file.path(f,'README.md'))
cat('Conley reference report and output checks complete.\n')
