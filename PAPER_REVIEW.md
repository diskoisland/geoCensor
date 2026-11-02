# Spatial methods: review for geoCensor

Reviewed November 1, 2026. This is a methods-focused first pass through the nine
original PDFs and three additional papers that arrived during review: twelve
distinct papers/versions across fourteen PDFs, including second copies of Ordoñez
et al. and Qiang & Craigmile. Selected assumptions, methods, simulation results,
and discussion were examined. This is not a line-by-line verification of every
proof. PDF page numbers
below count from the first page of each supplied file, including any publisher
cover. Paper findings and our proposed adaptations are distinguished explicitly.

The source inventory, SHA-256 fingerprints, page counts, and pages consulted are
saved in [the review record](output/literature_review_2026-11-01/README.md).
The supplied PDFs have not been renamed or modified.

**Supplement follow-up, November 1:** both requested supplements are now received
and reviewed. See [the follow-up and numerical formula checks](output/supplement_review_2026-11-01/README.md).
The Qiang supplement qualifies the prior/parameter-recovery evidence. The
Müller–Watson archive supplies Fortran replication sources and clarifies the
statistic's scaling. A bounded fixed-component Gaussian formula check passed;
the complete SCPC procedure and censored adaptation remain future work.

**Development follow-up, November 1:** completed a bounded numerical
[component-selection/envelope benchmark](output/scpc_envelope_validation_2026-11-01/README.md)
and an independent
[equal-weight censored transfer study](output/scpc_censored_transfer_site_control_2026-11-01/README.md).
The compact envelope remains a numerical grid/local-search approximation;
weighted and preferential designs are not validated. These are subsequent
experiments, not additional claims made by the papers.

**Weighted follow-up, November 1:** completed an original
[conditional weighted-Gaussian extension](output/scpc_weighted_gaussian_2026-11-01/README.md)
and [weighted censored transfer](output/scpc_weighted_censored_2026-11-01/README.md).
The fitted residual transformation and independent reference distribution are
recalculated for the weights. Results remain conditional on four fixed exogenous
site patterns; random-location weight uncertainty and preferential sampling are
not established. These studies do not change the papers' assumptions or promote
a production interval.

**Additional literature, November 1:** the published CensSpatial SoftwareX
article (2024) has now arrived and all four pages were reviewed. See
[the version-specific assessment and supplement request](output/literature_additions_2026-11-01/README.md).
Its prediction example is not evidence for common-parameter interval coverage
under heavy censoring. Additional Lahiri–Zhu and Ibragimov–Müller copies were
identified and preserved.

**CensSpatial supplement follow-up, November 1:** all nine pages of the newly
supplied 2024 supplement are now reviewed, with equations/figure on pages 4, 5
and 7 visually checked. It supplies the conditional censored-Gaussian likelihood
and SAEM steps; S.2 is rainfall trend/variogram diagnostics. It does not supply
heavy-censoring common-parameter interval calibration. A dimensional discrepancy
in printed prediction equation S.1.9 is recorded for independent checking, not
asserted to be an implementation defect. See
[the full assessment and benchmark implications](output/censspatial_supplement_review_2026-11-01/README.md).
Ross reports no 2018 Ordoñez supplement/code is available; the fetch request is
closed. Lehner's preprint is now received; its full methods review is pending.

**Other additions reviewed, November 1:** completed a methods-focused assessment
of Lehner's 34-page bandwidth preprint, Fouedjio's censored-random-forest paper,
the historical geoR article, and the newly identified Conley–Kelly online
appendix. Page-specific scope, practical implications, reported failures and
original deterministic counterchecks are in
[the additional methods review](output/additional_methods_review_2026-11-01/README.md).
Lehner's county-centroid table retains nominal 5% rejection rates up to 12.8%
for Epanechnikov; it does not establish universal calibration. A fixed smooth
taper can attenuate nonzero covariance even at an oracle support cutoff.
Fouedjio's geochemical prediction example is useful for a future mapping
comparison, but constraint-compatible ensembles need probability calibration:
random quadratic projection is not generally truncated-normal sampling.
The geoR article supplies uncensored Gaussian groundwork. The appendix adds
trend/basis/placebo diagnostics; its association with Conley–Kelly is now
supported by publisher metadata and the matching author-hosted Appendix C.
The spatial censored-likelihood benchmark remains the immediate priority.

**Geochemical database additions, November 1:** reviewed Zhang et al. (2005)
on the USGS NGS database and sample-size effects, Olea (2021) on stratified
declustering, and the methods/validation portions of Lawley et al. (2022) on
Zn–Pb prospectivity. See
[the focused review](output/geochemical_database_papers_2026-11-01/README.md)
for page scopes and implications. Zhang's 27-element analysis excluded the
hardest censoring examples, including Au/Ag, and used half-LOD values for
graphical display. Small-sample non-rejection is not model validation;
large-sample rejection is a power/model-adequacy issue distinct from spatial
independence and our correctly generated Gaussian censoring stress test.
Olea offers a useful resampling comparator, with geometry and empty-cell limits
that need checking. Exact enumeration verifies its linear mean/CDF connection
to occupied-cell weights and shows that rank-median recombination need not
preserve one observation per cell. Lawley's spatial predictive validation is
useful for a later mapping extension, not a concentration-mean confidence method.
No public estimator, filter, dataset or confidence limit changed in this review.

## What this changes for our work

The papers support continuing three separate lines of investigation:

1. **Spatial uncertainty:** verify a published covariance-and-critical-value
   procedure under its own assumptions before adapting it to censored scores.
2. **Uneven sampling independent of chemistry:** investigate design-aware
   resampling and the effect of the sampling density on calibration.
3. **Concentration-dependent sampling:** investigate a model for locations and
   responses together. A variance adjustment does not itself repair a biased
   estimating equation or reconstruct unsampled terrain.

These are development decisions, not claims of validated NGDB inference. Cell
declustering remains a descriptive sensitivity analysis. Its occupied-cell
weights are not known survey inclusion probabilities. No paper reviewed here
establishes coverage for our combination of censored likelihood, estimated
declustering weights, repeated locations, and mixed historical sample support.

## Paper-by-paper assessment

### Conley (1999): spatial GMM/HAC foundations

[Supplied paper](Papers/1-s2.0-S0304407698000840-main.pdf), especially PDF pp. 6-7,
10-12, and 22; [DOI](https://doi.org/10.1016/S0304-4076(98)00084-0).

**Paper:** spatial dependence enters through the covariance of moment conditions.
The directing process that selects observed locations is independent of the
underlying field in the stated model. The paper distinguishes increasing-domain
covariance consistency from merely computing a feasible matrix. The product of
one-dimensional Bartlett windows gives a positive-semidefinite construction
(equation 3.14). The conclusion explicitly identifies endogenous locations and
minimum separation as restrictions requiring further work.

**For geoCensor:** using joint location/log-scale scores in a sandwich is a sound
algebraic starting point. Product Bartlett is a useful reference kernel. This
does not validate normal critical values at our sample sizes, preferential site
selection, or many nearly coincident observations. Exact-coordinate aggregation
does not establish the minimum-separation assumption for distinct nearby sites.
Keep checking the full covariance matrix; a nonnegative diagonal alone is
insufficient.

### Bester, Conley, Hansen & Vogelsang (2016): fixed-b spatial HAC

[Supplied paper](Papers/fixed-b-asymptotics-for-spatially-dependent-robust-nonparametric-covariance-matrix-estimators.pdf),
PDF pp. 4-6, 10-13, 18, and Appendix B, pp. 31-33;
[DOI](https://doi.org/10.1017/S0266466614000814).

**Paper:** the bandwidths grow proportionally to the linear extent of the
expanding region: `h1 = b1*n`, `h2 = b2*n`, with sample size of order `n^2`.
The limiting covariance estimate remains random, and the test distribution
depends on the kernel, bandwidth fractions, and region shape. The stated main
results impose stationary scores, moment/mixing conditions, homogeneous limiting
derivatives, and a particular boundary geometry. Theorem 1's smooth-kernel
condition does not directly cover Bartlett; the authors discuss that extension
separately. Appendix B permits certain exogenous missing-location processes,
subject to the transformed scores satisfying the assumptions.

**For geoCensor:** this gives a theoretical reason to investigate nonstandard
critical values. It does not identify our earlier Gaussian-reference experiment
as fixed-b inference. State polygons, nonuniform density, and estimated weights
need explicit treatment. The iid-bootstrap argument for simulating critical
values on PDF p. 13 is presented as a conjecture, not a theorem to inherit.

### Müller & Watson (2022): spatial correlation principal components (SCPC)

[Supplied paper](<Papers/Econometrica - 2022 - Müller - Spatial Correlation Robust Inference.pdf>),
PDF pp. 4-9, 11, 15-16, and 21-22;
[DOI](https://doi.org/10.3982/ECTA19465).

**Paper:** construct leading eigenvectors of a demeaned benchmark covariance
matrix, estimate variance from projections onto those eigenvectors, and choose a
critical value controlling rejection over exponential benchmarks with
`c >= c0`. Choose the component count by expected interval length in the iid
reference model. Here `exp(-c*distance)` uses inverse range: smaller `c` means
more persistent dependence. Geometry and sampling density matter. Section 6
describes a GMM/influence extension, conditional on appropriate assumptions for
the resulting score process. Finite-sample guarantees concern specified Gaussian
covariance classes; the weak-correlation asymptotic argument is distinct.

**For geoCensor:** this is the strongest next published benchmark to reproduce
on small uncensored Gaussian datasets. Our existing product-Bartlett Gaussian
calibration is not SCPC: its variance construction and finite range grid differ.
Maximizing over four ranges is not the paper's supremum over a continuum. A
subsequent censored-score adaptation needs its own validation, including the
effects of fitted nuisance scale and estimated declustering weights. Do not
transfer finite-Gaussian guarantees to nonlinear estimated influences by name.

### Bester, Conley & Hansen (2011): a fixed number of large groups

[Supplied paper](Papers/1-s2.0-S0304407611000431-main.pdf), PDF pp. 3-6 and
Appendix C on p. 14; [DOI](https://doi.org/10.1016/j.jeconom.2011.01.007).

**Paper:** group-score sums can become asymptotically independent when contiguous
groups grow and boundary contributions become negligible under weak dependence.
The approach uses the pooled estimator with a cluster covariance. The familiar
scalar t reference requires the appropriate covariance scaling and homogeneity
conditions; more general heterogeneity produces a nonpivotal limit. The paper
also treats unequal groups and sketches nonlinear M-estimator extensions.

**For geoCensor:** worth a controlled comparator, not a reason to replace 1.96
with a t quantile in any existing HAC calculation. A raw cluster covariance
needs the specified group correction. Adjacent small blocks, uneven support,
and groups with different limiting information can defeat the simple reference.
Appendix C is a useful derivation route for censored scores, not a completed
proof for our weights and sampling history.

### Ibragimov & Müller (2010): t tests of group estimates

[Supplied paper](<Papers/t-Statistic Based Correlation and Heterogeneity Robust Inference.pdf>),
PDF pp. 4-5 and 11; [DOI](https://doi.org/10.1198/jbes.2009.08046).

**Paper:** estimate the same scalar parameter separately in a fixed number of
groups, then use a t statistic of those group estimates. Independent normal
group estimates can have unequal variances; their common mean is essential.
The large-sample version requires suitable joint convergence, including
asymptotic independence. Spatial groups need to grow enough that boundaries do
not dominate. The procedure uses the equally weighted mean of group estimates,
which is different from the pooled-estimator procedure above.

**For geoCensor:** sparse or heavily censored groups may not yield approximately
normal estimates. Separately fitted geological domains may estimate different
parameters, so their differences cannot automatically be treated as sampling
noise around one mean. Specify groups and target before confirmation; do not
search groupings until a preferred interval appears.

### Lahiri & Zhu (2006): bootstrap for irregular stochastic designs

[Supplied paper, currently Download.pdf](Papers/Download.pdf), PDF pp. 4-6 and
10-19; [DOI](https://doi.org/10.1214/009053606000000551).

**Paper:** nonuniform random locations are independent of the underlying field.
The data-site-shifted block bootstrap can fail under this nonuniform design.
The proposed grid-based block bootstrap (GBBB) draws blocks from a spatial grid,
handles boundary shapes, permits random bootstrap sample size, and centers the
bootstrap estimating equations. Section 5 gives the centering terms explicitly,
as well as cubic-block and nonoverlapping variants. These details matter along
with block size and domain growth.

**For geoCensor:** our occupied-block bootstrap is not GBBB. The published
inconsistency example is for a particular data-site-shifted method; it does not
by itself prove every other block method inconsistent. Equally, our use of square
blocks does not inherit GBBB validity. A reproduction should first cover an
uncensored scalar mean with exogenous nonuniform locations, then establish how
score centering, censoring, and estimated weights would be handled. Empty areas
and eligible source blocks must be defined by the observation window rather
than inferred solely from occupied blocks.

### Diggle, Menezes & Su (2010): joint modeling of locations and concentrations

[Supplied paper, including discussion and rejoinder](<Papers/Royal Stata Society Series C - 2010 - Diggle - Geostatistical inference under preferential sampling.pdf>),
PDF pp. 2-4, 8-9, 15-16, and discussion/rejoinder pp. 39-40;
[DOI](https://doi.org/10.1111/j.1467-9876.2009.00701.x).

**Paper:** preferential sampling is dependence between the location process and
the latent field, not simply an uneven map. A shared latent Gaussian process
drives a log-Gaussian Cox intensity and responses. Joint likelihood includes the
location pattern over the observation region. The discussion is material:
sampling based on earlier observed measurements can be ignorable when the
relevant conditioning history is retained. A shared dependence on measured
explanatory variables is also distinguished from residual preferential sampling.
Reported plug-in predictive distributions do not include parameter uncertainty.

**For geoCensor:** a high-grade cluster is a reason to investigate the sampling
mechanism, not proof of preferential bias. Project boundaries, sampling campaigns,
retained earlier measurements, geological covariates, and access constraints may
matter. A state boundary is an analysis domain, not evidence that the historical
survey could have sampled everywhere inside it. Repeated coordinates and varied
sample support complicate a simple point-process model.

### Qiang & Craigmile (2026): flexible preferential-sampling models

[Supplied paper](<Papers/Environmetrics - 2026 - Qiang - Accounting for Preferential Sampling in Geostatistical Inference.pdf>),
PDF pp. 2-7 and 11-12; [DOI](https://doi.org/10.1002/env.70081).

**Paper:** the framework links a transformation of point intensity to the latent
response process, permits non-Gaussian observations, and explores location
processes including LGCP and Thomas clustering. Section 4 offers a simpler
predictor using estimated log intensity as a covariate in universal kriging.
The simulations compare prediction error and overprediction, not just nominal
interval coverage. In the Gaussian/LGCP example (200 simulations), the plug-in
method improves on ignoring preference but retains more directional prediction
error than joint Bayesian fitting. The Thomas implementation conditions on
known or pre-estimated parent locations, a limitation acknowledged in discussion.

**For geoCensor:** the plug-in intensity approach is a tractable sensitivity
candidate, not an automatic correction or a diagnostic proving preference.
Intensity-estimation error, boundaries, coordinate reuse, and chemical covariates
need controls. Replacing conditional observation densities with interval
probabilities could accommodate censoring in principle; this is our proposed
extension, not a censored implementation established by this paper. No numerical
zero-imputation should be introduced to make the model fit.

### Ordoñez et al. (2018): a spatial likelihood retaining censoring

[Supplied paper](Papers/1-s2.0-S221167531730043X-main.pdf), PDF pp. 3-7 and 9-13;
[DOI](https://doi.org/10.1016/j.spasta.2017.12.001).

**Paper:** a Gaussian spatial linear model combines mean structure, nugget,
and spatial covariance. Observed likelihood multiplies the density of exact
responses by the conditional probability of censored responses lying inside
their intervals (equations 5 and 7). SAEM uses truncated multivariate-normal
simulation rather than detection-limit substitution. The application has 1,500
arsenic observations, with 392 left-censored; the simulations use 500 datasets
per setting, sample sizes 50/200, and 25%/45% left censoring. The general model
allows left, right, and interval censoring, but these simulations do not validate
every direction. The paper explicitly notes limitations for very large datasets.

**For geoCensor:** a useful model-based comparator on small, comparable domains.
For positive concentrations modeled on the log scale, transform finite bounds
consistently and preserve infinite bounds. This differs from adding HAC to a
working-independence fit and does not itself model preferential locations.
Prediction error improvements are not proof of coverage for our geometric-mean
interval. Validate covariance estimation, convergence, and parameter uncertainty
separately; distinguish prediction of a latent field from a future measurement.

### Ordoñez, Galarza & Lachos (2021): CensSpatial software preprint

[Supplied preprint](Papers/2110.05570v1.pdf), PDF pp. 2-3, 6, and 9-11;
[arXiv record](https://arxiv.org/abs/2110.05570).

**Paper:** documents estimation, prediction, and perturbation/case-influence tools
for censored spatial Gaussian models. The rainfall illustration artificially
censors 25% of the first 100 observations and predicts the remaining 43. The
diagnostic construction uses the expected complete-data log-likelihood (the
EM Q function), its Hessian, and perturbations of responses, covariance, and
covariates. The paper explicitly treats diagnostic thresholds as context-dependent.

**For geoCensor:** an independent small-domain implementation and influence
comparison would be useful. These local perturbation diagnostics are different
from our estimator influence contributions used for HAC variance. A Q-function
Hessian must not be assumed to equal the observed-data likelihood information;
verify derivatives, signs, and uncertainty calculations before adapting them.
Review flags should trigger inspection, not automatic removal of a high-grade
sample. Paper-level interval support still needs a separate audit of the actual
software version we choose to run.

### Zeileis (2006), supplied modified article: sandwich software structure

[Supplied modified article](Papers/sandwich-OOP.pdf), PDF pp. 2, 4-7, and 13;
[original article DOI](https://doi.org/10.18637/jss.v016.i09).

**Paper:** separates extraction of empirical estimating functions, derivative
information, covariance construction, and inference. The supplied document
identifies itself as a slightly modified version of the published article.
Its normalized bread and meat use averages, and the parameter covariance adds
the sample-size normalization. Its time-series HAC example uses observation
ordering, not two-dimensional geography. The discussion warns that robust
covariance is useful only when the target parameter remains consistently estimated.

**For geoCensor:** preserve the independently checked weighted score extractor,
including nuisance scale. Our research code uses inverse summed information and
unnormalized score cross-products: retain this convention consistently rather
than adding the library formula's sample-size factors again. Ordinary ordered
HAC is not a substitute for a spatial distance matrix. This conceptual article
does not verify any particular current `estfun.survreg` implementation; our earlier
version-specific numerical checks remain separate evidence.

### Kelly (2019): an earlier persistence working paper

[Supplied working paper](Papers/WP2019_13.pdf), PDF pp. 1-3, 16-18, and 34.
The cover identifies UCD WP19/13, April 2019, by Morgan Kelly alone. This is not
the requested Conley & Kelly (2025) published version.

**Paper:** uses spatially correlated artificial variables to investigate spurious
regression significance. Its examples show that cutoff choice and clustered
locations affect the performance of conventional Conley intervals. It advocates
residual spatial diagnostics and spatial-noise placebo experiments; ordinary
random permutations destroy the dependence structure the exercise needs to retain.

**For geoCensor:** this is useful when testing future element-element or
covariate associations. Such simulations should preserve geography and the
relevant censoring mechanism. A Moran statistic is a dependence diagnostic, not
a correction for preferential sampling, a calibrated global significance test,
or proof that an association is spurious. The older working paper's empirical
results are not assumed to be the final published version's results.

## Supplements and replication material to fetch

The first two were received and reviewed later on November 1, 2026. Their detailed
assessment is in the supplement follow-up linked above. Other entries remain
requests unless indicated otherwise.

| Priority | Item | Where / why |
| --- | --- | --- |
| 1 | RECEIVED / REVIEWED: Müller & Watson (2022), `ecta200510-sup-0001-dataandprograms.zip` | Fortran sources, data, and Intel/IMSL build requirements inspected; original program not run. [Publisher](https://onlinelibrary.wiley.com/doi/10.3982/ECTA19465). |
| 2 | RECEIVED / REVIEWED: Qiang & Craigmile (2026), Data S1 | Local filename `env70081-sup-0001-supplementalmaterial.pdf`; all eight pages reviewed. [Publisher](https://onlinelibrary.wiley.com/doi/10.1002/env.70081). The separate [PrefSampling code/data](https://github.com/petercraigmile/PrefSampling/) has not been audited here. |
| 3 | Bester et al. (2016), any author simulation/critical-value code | [Article DOI](https://doi.org/10.1017/S0266466614000814). Useful if supplied; a separate supplement was not confirmed. Appendices A/B are already in our PDF. |
| 4 | CLOSED: Ordoñez et al. (2018), application/simulation scripts | Ross reports none available. Main paper is already reviewed; no further fetch requested. The 2024 CensSpatial supplement is received and reviewed separately. |
| 5 | Ibragimov & Müller (2010), supplementary materials | [Article DOI](https://doi.org/10.1198/jbes.2009.08046); supplementary calculations are referenced in the paper. Exact download not checked. |

The Diggle PDF already contains discussion and rejoinder; these need not be
fetched again. Bester et al. (2011) already includes the unequal-group and
nonlinear appendices.

## Next experiment and promotion criteria

**First priority: an independently checked SCPC benchmark.** Reproduce the
published scalar Gaussian construction on small fixed designs with uniform and
uneven locations. Specify the dependence envelope and component selection before
examining confirmation data. Compare projection variance, critical values, and
rejection probabilities with author material. Include iid and strong-dependence
controls, location-scale/rotation checks, and nonuniform exogenous sampling.

Only then try the existing full-joint-score influence representation for a
censored lognormal location parameter. First use known/fixed weights, then
estimated declustering weights; report failures and coverage separately. Retain
independent calibration/development/confirmation seeds. The previous confirmation
set must not become a tuning set followed by a claim of new confirmation.

**Second priority: an exogenous-nonuniform GBBB reproduction.** Verify source
block eligibility, centering, random sample counts, boundaries, and the scalar
mean before proposing a censored estimating-equation version.

**Separate selection experiment:** generate clustered exogenous locations and
response-dependent locations as different controls; compare no-preference,
intensity-covariate, and joint-model candidates on a declared observation window.
Score prediction and domain-average bias separately from interval coverage. The
realized-domain spatial average and a superpopulation location parameter are
different targets. Retain a no-preference control to detect damage caused by an
unnecessary correction. This simulation cannot establish the amount of bias in
the actual USGS compilation.

Production promotion requires a credible target and assumptions, numerical
agreement with an independent reference where available, and satisfactory
prespecified validation including irregular locations and censoring. Successful
optimization, a smooth map, or one favorable covariance setting is insufficient.

## Incoming files and remaining original requests

The subsequent weighted SCPC experiments have now reached independently generated
exogenous survey layouts, with weights and calibration recalculated for each
layout. Across 200 layouts and 30,000 censored analyses, confirmation interval
coverage was 96.25–99.05% under moderate censoring but 85.50–93.75% under heavy
censoring. Available-only heavy-censoring coverage was still 89.39–93.94%, with
downward point-estimate bias. The Gaussian benchmark met the 95% target within
the tested dependence envelope. These results identify a failure of the censored
transfer; they do not establish production intervals or robustness to preferential
sampling. The next priority is a separately validated spatial-likelihood or
bias-aware benchmark for heavy censoring, retaining these confirmation results
as completed evidence rather than using them for further tuning. See
[the changing-layout study](output/scpc_random_designs_2026-11-01/README.md).

The CensSpatial preprint, modified Zeileis article, and Kelly (2019) working paper
arrived during review and are assessed above. The extra
`YGd5HK-1-s2.0-S221167531730043X-main.pdf` is a second copy of the Ordoñez article,
not its supplement. A second Qiang PDF ending in `-2.pdf` is also a copy of the
main article, not Data S1. Data S1 subsequently arrived as its own PDF and is
reviewed in the follow-up. All copies remain in place. The review inventory records their file
fingerprints and full extracted-text comparisons (allowing for download stamps).

Lehner's bandwidth preprint is now received and methods-reviewed. The Conley &
Kelly online appendix has been identified and its method reviewed; the published
main PDF is still not present locally. The fetch list distinguishes receipt
from review; new arrivals can be added without renaming existing files.


## Site-block permutation methods — November 1, 2026

Online method excerpts consulted (not a full local-PDF review): Hemerik & Goeman
(2018), [Exact testing with random permutations](https://doi.org/10.1007/s11749-017-0571-1),
especially the transformation-group/identity-inclusion result; Phipson & Smyth
(2010), [Permutation P-values Should Never Be Zero, author preprint](https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf),
for nonzero randomly drawn permutation p-values. No new supplied PDF or supplement
is required for the present calculation.

These support the randomization machinery under null invariance, not our
censor-aware order score's power or real-USGS eligibility. The new research
procedure preserves whole repeated-site vectors and reporting bounds, calibrates
the maximum over fixed distance scales, requires an exchangeable site-vector null
and common reporting design, and leaves absent information/eligibility as NA.
This differs from the Kelly warning above: permutation here intentionally tests
a specified exchangeable null; it is not used to mimic a spatial sampling
distribution for geochemical confidence limits.

The new score and exact6-site checks are original reference calculations. Frozen
development/confirmation experiments and their assumptions are documented in
[the study report](output/spatial_order_signal_2026-11-01/README.md). No general
NGDB significance test or concentration confidence interval follows.
