#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(readr);library(dplyr);library(ggplot2)})
f<-'output/irregular_spatial_validation_2026-10-31'
stopifnot(file.exists(file.path(f,'COMPLETED.txt')))
s<-read_csv(file.path(f,'summary.csv'),show_col_types=FALSE)
x<-read_csv(file.path(f,'results.csv'),show_col_types=FALSE)
stopifnot(nrow(x)==1000,nrow(s)==10,all(s$n_datasets==100),all(s$n_intervals<=100))
labels<-c(irregular_independent_uniform='Independent; uniformly selected sites',
 irregular_correlated_uniform='Correlated; uniformly selected sites',
 irregular_correlated_uniform_shifted='Correlated; uniform, shifted blocks',
 irregular_correlated_preferential='Correlated; preferentially selected sites',
 irregular_correlated_preferential_shifted='Correlated; preferential, shifted blocks')
a<-s%>%filter(phase=='confirmation')%>%mutate(label=factor(labels[scenario],levels=rev(unname(labels))))
p<-ggplot(a,aes(y=label,x=coverage))+
 geom_vline(xintercept=.95,linetype='dashed',color='grey45')+
 geom_segment(aes(x=mc_low,xend=mc_high,yend=label),linewidth=.8,color='#276b8e')+
 geom_point(size=3,color='#276b8e')+
 scale_x_continuous(limits=c(0,1),labels=function(z)paste0(round(100*z),'%'))+
 labs(title='Irregular-site confirmation: nominal 95% bootstrap coverage',
 subtitle='64 observed sites selected from a 144-site domain; censored data and unequal records per site',
 x='Coverage of the superpopulation geometric mean',y=NULL,
 caption='100 confirmation datasets per setting; 199 resamples each. Bars: exact 95% Monte Carlo intervals.\nPreferential selection favors high latent concentrations. No inclusion-probability correction was supplied.')+
 theme_minimal(base_size=11)+theme(panel.grid.major.y=element_blank())
ggsave(file.path(f,'confirmation_coverage.png'),p,width=11,height=5.5,dpi=160)
lines<-c('# Irregular sampling and preferential site-omission study','',
 'Prespecified development and independently seeded confirmation phases used 100 datasets per setting, with 199 actual package bootstrap refits each. Candidate locations were jittered from a 12x12 grid; only 64 of the 144 candidate sites were observed. All recorded observations at a site shared the same coordinate and had unequal counts independent of concentration.','',
 'Uniform site selection was independent of concentration. Preferential selection sampled without replacement with probabilities proportional to exp(2 * latent effect / 0.55). This mechanism omits many low-concentration sites rather than merely adding more measurements at already observed high-concentration sites.','',
 'The primary target is the superpopulation geometric mean of 50. The secondary target is the realized equal-site geometric mean over all 144 candidate sites, including unobserved ones. Neither was replaced by the observed-site mean to improve reported coverage.','',
 '| Confirmation setting | Intervals | Population coverage | Finite-domain coverage | Mean point estimate | Mean log error for finite domain |',
 '| --- | ---: | ---: | ---: | ---: | ---: |')
for(i in seq_len(nrow(a))){z<-a[i,];lines<-c(lines,sprintf('| %s | %d/100 | %.1f%% | %.1f%% | %.1f | %.3f |',
 labels[z$scenario],z$n_intervals,100*z$coverage,100*z$finite_domain_coverage,z$mean_geomean,z$mean_log_error_finite))}
lines<-c(lines,'','![Confirmation coverage](confirmation_coverage.png)','',
 '## Interpretation','',
 'Coverage, point-estimate error, and block-origin sensitivity are separate diagnostics. A change of block origin leaves the original declustered point estimate unchanged but can change its bootstrap distribution. Successful refitting is not evidence of valid interval coverage.',
 'The hidden complete-data observed-site mean is included in summary.csv only to diagnose selection bias that persists even without censoring. It is not an available correction for real data and is not given a purported oracle confidence interval under informative selection.',
 'Declustering repeated records within observed locations cannot by itself reconstruct the unsampled low-concentration sites. Enlarging intervals is not a substitute for defining the target population and addressing informative site selection.','',
 '## Decision and next work','',
 'Keep the bootstrap explicitly experimental and retain missing routine spatial intervals in the full NGDB outputs. Review support and method mixtures before defining comparable populations. For geographic-domain inference, investigate sampling-design information or an explicitly stated spatial/selection model and validate it with separate simulations. A generalized bootstrap for irregular but exogenous sampling must not be assumed to solve concentration-dependent site selection.',
 'The present study is one synthetic design with a fixed correlation scale and selection strength. It does not identify or quantify the true selection mechanism in NGDB data. Its uncertainty summaries include binomial Monte Carlo intervals; no automatic block size or origin was chosen from these results.','',
 '## Reproduction and checks','',
 'Run scripts/analysis/validate_irregular_spatial_sampling.R and scripts/analysis/report_irregular_spatial_sampling.R from the repository root. Use a new output directory for changed design or code. Source data were not modified and previous full-data estimates were not recalculated.',
 'Files include design_notes.txt, settings.csv, results.csv, summary.csv, confirmation_coverage.png, session_info.txt, and phase checkpoints. Output checks reconcile 1,000 dataset/setting rows and 100 datasets per phase/setting.')
writeLines(lines,file.path(f,'README.md'))
print(a%>%select(scenario,n_intervals,coverage,finite_domain_coverage,mean_geomean,mean_log_error_finite))
cat('Irregular-study output checks passed.\n')
