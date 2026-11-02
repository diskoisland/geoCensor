#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(dplyr);library(readr);library(ggplot2)})
args<-commandArgs(trailingOnly=TRUE)
f<-if(length(args))args[1]else 'output/gaussian_hac_calibration_2026-10-31'
stopifnot(file.exists(file.path(f,'COMPLETED.txt')))
c<-read_csv(file.path(f,'critical_values.csv'),show_col_types=FALSE)
g<-read_csv(file.path(f,'gaussian_validation.csv'),show_col_types=FALSE)
s<-read_csv(file.path(f,'censored_summary.csv'),show_col_types=FALSE)
x<-read_csv(file.path(f,'censored_results.csv'),show_col_types=FALSE)
range_file<-if(file.exists(file.path(f,'influence_ranges_bins_checked.csv')))'influence_ranges_bins_checked.csv'else 'influence_ranges.csv'
d<-read_csv(file.path(f,range_file),show_col_types=FALSE)
stopifnot(nrow(c)==32,nrow(g)==128,nrow(s)==192,nrow(x)==76800,nrow(d)==4800)
stopifnot(all(is.finite(x$se)),all(x$se>0),all(s$n_datasets%in%c(300,500)))
a<-s%>%filter(phase=='confirmation')%>%mutate(mean_log_width=2*critical*mean_se)
labels<-c(normal='Normal 1.96',gaussian_iid='Gaussian iid calibration',gaussian_reference_0.8='Gaussian reference range 0.8',gaussian_range_grid_max='Maximum over fixed Gaussian range grid')
a$method_label<-factor(labels[a$method],levels=unname(labels))
p<-ggplot(a,aes(x=bandwidth,y=coverage,color=method_label))+
 geom_hline(yintercept=.95,linetype='dashed',color='grey40')+geom_line()+geom_point()+
 scale_x_continuous(trans='log2',breaks=c(2,4,8,16))+
 scale_y_continuous(limits=c(0,1),labels=function(z)paste0(round(100*z),'%'))+
 facet_grid(n_sites~true_range,labeller=label_both)+
 labs(title='Gaussian-reference calibration: independent censored confirmation',
 subtitle='Same HAC variances; different prespecified critical-value assumptions',
 x='Fixed bandwidth',y='Coverage of nominal 95% intervals',color=NULL,
 caption='500 confirmation datasets per setting. Correlation range is the exponential decay parameter.\nCalibration assumes a specified Gaussian covariance class; it is not estimated from these censored datasets.')+
 theme_minimal(base_size=11)+theme(legend.position='bottom')+guides(color=guide_legend(nrow=2))
ggsave(file.path(f,'censored_confirmation_coverage.png'),p,width=12,height=8,dpi=160)
q<-ggplot(c,aes(x=bandwidth,y=critical,color=factor(assumed_range)))+
 geom_hline(yintercept=qnorm(.975),linetype='dashed',color='grey40')+geom_line()+geom_point()+
 scale_x_continuous(trans='log2',breaks=c(2,4,8,16))+facet_wrap(~n_sites,scales='free_y')+
 labs(title='Critical values conditional on the assumed Gaussian covariance',x='Fixed bandwidth',y='Simulated 95% critical value',
 color='Assumed range',caption='10,000 calibration draws per design/range; range 0 denotes independent sites.\nDashed line: 1.96. A larger critical value compensates for a reference distribution, not a selection-biased point estimate.')+
 theme_minimal(base_size=11)+theme(legend.position='bottom')
ggsave(file.path(f,'critical_values.png'),q,width=10,height=5.5,dpi=160)
range_summary<-d%>%group_by(phase,n_sites,true_range,range_status)%>%summarise(n=n(),
 median_zero_crossing=if(all(is.na(estimated_zero_crossing)))NA_real_ else median(estimated_zero_crossing,na.rm=TRUE),.groups='drop')
write_csv(range_summary,file.path(f,'influence_range_summary.csv'))
lines<-c('# Gaussian-reference calibration and influence-range diagnostics','',
 'This is a deliberately assumption-dependent research experiment, not a production geoCensor interval method and not a claimed implementation of published fixed-b or Müller–Watson procedures. It asks whether the critical value, in addition to the covariance estimate, matters for the observed undercoverage.','',
 '## Design and safeguards','',
 'The product Bartlett kernel and bandwidths 2, 4, 8, and 16 were fixed in advance. Designs have 64 or 256 observed sites on a regular grid. No site was preferentially omitted in this study. The Gaussian reference assumes latent variance 0.55^2, four observations per site with noise variance 0.2^2, and exponential covariance with assumed decay parameters 0, 0.25, 0.8, or 2. Range 0 is an independent-site reference.',
 'Each design/range used 10,000 Gaussian calibration draws, then 5,000 independently seeded Gaussian validation draws. Critical values were frozen before the censored experiment: 300 development and 500 independently seeded confirmation datasets for each of six site-count/range settings. All three actual correlation ranges and all four critical-value methods were prespecified.',
 'The range-grid maximum takes the largest calibrated critical value over the finite assumed grid. It is not a distribution-free bound for arbitrary correlation, heterogeneity, censoring, or sampling selection. A reference with range 0.8 is correctly matched only when that covariance assumption describes the generating model.','',
 '## Exact Gaussian benchmark, bandwidth 4, actual range 0.8','',
 '| Sites | Method | Independent validation coverage | Critical value |','| --- | --- | ---: | ---: |')
for(i in which(g$true_range==.8 & g$bandwidth==4)){
 z<-g[i,];lines<-c(lines,sprintf('| %d | %s | %.2f%% | %.3f |',z$n_sites,labels[z$method],100*z$coverage,z$critical))}
lines<-c(lines,'','## Censored confirmation, bandwidth 4','',
 'This one prespecified bandwidth is displayed for readability; summary.csv and the figure retain every candidate. It is not a selected production bandwidth.','',
 '| Sites | Actual range | Method | Coverage | 95% Monte Carlo interval | Mean log interval width |',
 '| --- | ---: | --- | ---: | --- | ---: |')
for(i in which(a$bandwidth==4)){
 z<-a[i,];lines<-c(lines,sprintf('| %d | %.2f | %s | %.1f%% | %.1f–%.1f%% | %.3f |',
 z$n_sites,z$true_range,labels[z$method],100*z$coverage,100*z$mc_low,100*z$mc_high,z$mean_log_width))}
lines<-c(lines,'','![Confirmation coverage](censored_confirmation_coverage.png)','',
 '![Conditional critical values](critical_values.png)','',
 '## Interpretation','',
 'For actual range 0.8 at bandwidth 4, the censored normal-critical-value coverage was 83.8% and 86.0% for 64 and 256 sites. The matched Gaussian-reference critical values raised coverage to 96.2% and 96.4%. This is encouraging evidence for investigating calibration under the stated reference assumptions; it is not evidence that those assumptions are known for NGDB.',
 'When the actual range was 2 but the assumed reference range remained 0.8, coverage fell to 83.8% and 84.8%. Taking the maximum critical over the fixed Gaussian range grid gave 96.4% and 95.2% in those settings, with conservative or even 100% coverage under weaker correlation. Robustness to a broader assumed class has a cost in interval width. The study does not establish protection outside the finite assumed class.',
 'The Gaussian benchmark isolates finite-design calibration of the HAC studentized statistic. Transfer to nonlinear censored scores was evaluated empirically rather than assumed. The method retains the original point estimates and has no mechanism for correcting omitted-site selection bias.','',
 '## Influence-based range diagnostic','',
 'The location influence contribution uses the location row of the full joint bread times each aggregated site score, including nuisance-scale effects. An identity check verifies that its kernel quadratic form reproduces the joint HAC location variance. Diagnostics use 0.5-wide lag bins to lag 6, at least 20 pairs per bin, and the first positive-to-nonpositive crossing. Unavailable ranges remain missing; there is no silent 1 km or other fallback.',
 'The exponential decay parameter is not a population zero-crossing distance: exponential covariance remains positive at all finite lags. Empirical zero crossings depend on sampling variation, estimated centering, domain size, and binning. They must not be interpreted as unbiased estimates of that decay parameter.',
 'Lag bins were verified for translation invariance and corrected for floating-point boundary assignment. The completed interval calculations did not use those bins. The corrected diagnostics were regenerated from identical seeds in a separate bins-checked study folder, leaving the original files intact. The interval results were checked against the original study. A source/runtime manifest now prevents silent reuse of changed-code checkpoints; reproduction should specify a fresh output directory.', '',
 'No interval in this experiment uses a bandwidth selected from its estimated influence range. That rule needs a separate prespecified validation. Weak correlation can produce no initial positive bin, and stronger correlation can have no zero crossing within the chosen maximum lag. Their exact counts are retained in influence_range_summary.csv rather than replaced with fallback ranges.',
 '## Decision and next work','',
 'Keep these calculations research-only. The experiment identifies a plausible calibration route, not a ready default. Next read the fixed-b and spatial correlation-robust papers in PAPER_FETCH_LIST.md; compare theory-supported calibration and influence-based bandwidth rules under estimated covariance, irregular locations, support heterogeneity, and separately modeled preferential site selection. Do not tune a correlation class or critical value to obtain 95% on this confirmation set.','',
 'Reproduce with scripts/analysis/validate_gaussian_hac_calibration.R and scripts/analysis/report_gaussian_hac_calibration.R. Reference functions are in conley_survreg_reference.R and influence_covariogram_reference.R. The verification script verify_influence_covariogram.R checks the joint influence/HAC identity and translation invariance.','',
 'Output checks reconcile 76,800 censored interval evaluations, 4,800 fitted censored datasets, 192 summary rows, 32 calibrated critical values, and 128 Gaussian validation summaries. No routine package API, full-data USGS estimate, or production spatial interval was changed.')
writeLines(lines,file.path(f,'README.md'))
cat('Gaussian-reference calibration report and reconciliation passed.\n')
