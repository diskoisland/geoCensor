# Papers to fetch for geoCensor's spatial-inference work

Prepared October 31, 2026. Ordered by the implementation decisions they help
resolve. This is a screened reading/fetch list, not a claim that every paper has
been read in full. Abstracts, publisher/author metadata, and selected software
reference material were inspected. The relevance notes are our assessment.

## Receipt and review update — November 1, 2026

The nine original PDFs plus three papers that arrived during review have received
a methods-focused first review (twelve distinct papers/versions; fourteen PDFs
including second copies of Ordoñez et al. and Qiang & Craigmile). See [PAPER_REVIEW.md](PAPER_REVIEW.md) for findings, assumptions, page
references, and specific supplement requests. This does not claim a complete
proof audit or review of supplements that have not yet arrived. The original
access notes below describe the October 31 search, not today's local holdings.

| Original item | Local receipt / review |
| --- | --- |
| Lehner (2026) | Received as `2603.03997v1.pdf`; methods and simulation tables reviewed |
| Bester et al. (2016), fixed-b | Received; methods reviewed; appendices included |
| Müller & Watson (2022) | Received; methods reviewed; replication ZIP also received and source-inspected |
| Bester et al. (2011) | Received; methods reviewed; nonlinear appendix included |
| Ibragimov & Müller (2010) | Received; methods reviewed |
| Conley (1999) | Received; methods reviewed |
| Lahiri & Zhu (2006) | Received as `Download.pdf`; methods reviewed |
| Diggle et al. (2010) | Received with discussion/rejoinder; methods and selected discussion reviewed |
| Qiang & Craigmile (2026) | Received; methods and Data S1 reviewed |
| Ordoñez et al. (2018) | Received; methods reviewed |
| CensSpatial software preprint (2021) | Received and reviewed; the published SoftwareX (2024) article also arrived and was reviewed separately |
| Zeileis (2006) | Received during review as `sandwich-OOP.pdf`, a modified article version; methods reviewed |
| Conley & Kelly (2025) | Online appendix identified and method reviewed; published main PDF/code not present. Earlier Kelly (2019) working paper reviewed separately |
| Ribeiro & Diggle (2001), geoR | Newsletter article received and reviewed as historical Gaussian-model groundwork |
| Fouedjio, censored random forest | Received; methods and prediction examples reviewed; no calibrated interval claim inferred |
| Zhang et al. (2005), geochemical database/sample size | Received; sampling, censoring conventions, sample-size results and limitations reviewed |
| Olea (2021), preferential-sampling declustering | Received; method/examples reviewed; USGS software-release availability verified, code not run |
| Lawley et al. (2022), Zn–Pb prospectivity | Received; training/preprocessing/validation methods screened; supplements/data not reviewed |

**Supplement receipt update:** Müller & Watson's
`ecta200510-sup-0001-dataandprograms.zip` and Qiang & Craigmile's
Data S1 (`env70081-sup-0001-supplementalmaterial.pdf`) have now arrived and been
reviewed. See [the supplement assessment](output/supplement_review_2026-11-01/README.md).
They no longer need fetching. Other requests and paper-version distinctions
remain in the review. New files may arrive after these inventory snapshots.

**New published-paper follow-up:** the CensSpatial SoftwareX article (2024), DOI
[10.1016/j.softx.2024.101762](https://doi.org/10.1016/j.softx.2024.101762), is now
received and reviewed. Its **Supplementary Material** has also arrived as
`1-s2.0-S235271102400133X-mmc1.pdf`; all nine pages, including S.2, are reviewed.
See [the supplement review](output/censspatial_supplement_review_2026-11-01/README.md).
Ross reports no supplements or code are available for Ordoñez et al. (2018);
that request is closed. This is the user's availability report, not an
independent claim that no such material could exist anywhere.

## Remaining fetch priorities — November 1, 2026

The main papers in the original four-paper shortlist are already received.
After the random-layout study exposed heavy-censoring bias and undercoverage,
the most useful outstanding material is:

1. **Conley & Kelly (2025): published article and replication files.** Lower
   priority for later spatial-trend diagnostics. The Kelly (2019) working paper
   is already received; it is a different version.
   [Published article](https://doi.org/10.1016/j.jinteco.2024.104027) ·
   [Replication files](https://doi.org/10.17632/p2yt92jnnr.1).
   The 18-page online appendix is now received as
   `qxKeoa-1-s2.0-S0022199624001545-mmc1.pdf`. Its association was checked
   against publisher metadata and the matching author-hosted working-paper
   Appendix C. Its method was reviewed; it does not replace the published
   main PDF or constitute a code archive.

**Now received and methods-reviewed:** Alexander Lehner's *Bandwidth
Selection for Spatial HAC Standard Errors* (2026),
[arXiv:2603.03997](https://arxiv.org/abs/2603.03997), as `2603.03997v1.pdf`.
The method, assumptions, simulation tables and appendix were reviewed, with
definition/calibration cautions recorded in
[the additional methods assessment](output/additional_methods_review_2026-11-01/README.md).
This is useful supporting work on HAC bandwidths, rather than a correction
for the newly observed point-estimation bias.

**No need to fetch again:** Müller–Watson's replication ZIP, Qiang–Craigmile's
Data S1, the 2024 CensSpatial supplement, Lehner's preprint, or the sandwich
paper/package. The three supplements and the sandwich paper
are already reviewed; sandwich 3.1.3 is installed. Optional author code for
Bester et al. (2016) and supplementary calculations for Ibragimov–Müller (2010)
remain useful, but are below the likelihood benchmark in priority.

The sections below retain the original reading sequence and historical access
notes. Use the receipt table and remaining priorities above for current fetch
status. Do not purchase anything solely for this list.

No additional paper is required to start the small spatial censored-likelihood
benchmark. Optional author replication scripts for Lehner would help clarify
the difference between the printed near-zero rule and the described sign
crossing, plus binning, empty-bin and covariance-definiteness handling.

The three newly supplied geochemical references are assessed in
[the focused database/sample-size review](output/geochemical_database_papers_2026-11-01/README.md).
Olea's Fortran decluster release is publicly listed by
[USGS](https://www.usgs.gov/software/computer-programs-assessment-coal-resources-ver-20-april-2021)
under DOI 10.5066/P971RL9L; it has not been downloaded or compiled here.
Lawley's paper references supplementary tables/data, but these are not needed
for the current concentration/uncertainty gate and are not an urgent fetch
request. None of these references calls for refetching the main PDFs.

## First reading sequence: spatial uncertainty and finite-sample calibration

### 1. Alexander Lehner (2026)

**Bandwidth Selection for Spatial HAC Standard Errors.** arXiv:2603.03997.
[Paper/preprint](https://arxiv.org/abs/2603.03997) ·
[PDF](https://arxiv.org/pdf/2603.03997) ·
[SpatialInference vignette](https://cran.r-universe.dev/articles/SpatialInference/spatial_HAC.html).

- **Access:** open preprint; identified as a preprint rather than a journal publication.
- **Why:** directly motivates the user's bandwidth/kernel diagnostic suggestion.
  Read the covariogram rule, simulation design, and failure cases. Its linear-model
  validation does not establish validity for our weighted censored estimator.
- **Question:** can an estimator-specific influence covariogram support a useful
  bandwidth rule, and how much uncertainty is introduced by choosing that range?

### 2. C. Alan Bester, Timothy G. Conley, Christian B. Hansen & Timothy J. Vogelsang (2016)

**Fixed-b Asymptotics for Spatially Dependent Robust Nonparametric Covariance
Matrix Estimators.** *Econometric Theory* 32(1), 154–186.
[DOI/publisher](https://doi.org/10.1017/S0266466614000814).

- **Access:** fetch published full text and supplementary critical-value/code material
  if available. Journal year is 2016; the article appeared online in 2014.
- **Why:** our normal-critical-value HAC intervals undercovered even with correct
  score algebra. This addresses the distribution of test statistics when smoothing
  bandwidth is not negligible relative to the spatial domain.
- **Question:** which nonstandard critical values, domain-shape assumptions, and
  bandwidth fractions are applicable to our score/influence representation?

### 3. Ulrich K. Müller & Mark W. Watson (2022)

**Spatial Correlation Robust Inference.** *Econometrica* 90(6), 2901–2935.
[DOI/publisher](https://doi.org/10.3982/ECTA19465) ·
[Author PDF](https://www.princeton.edu/~mwatson/papers/ECTA19465.pdf).

- **Access:** author full-text PDF and publisher replication material are linked.
- **Why:** treats both standard-error construction and critical-value calibration;
  useful as a principled comparator to simply using 1.96 with an HAC variance.
- **Question:** can its scalar-mean/GMM arguments be adapted to estimated censored
  influence contributions, and which restrictions on spatial covariance are needed?
- **Fetch if convenient:** the published supplement and replication files.

### 4. C. Alan Bester, Timothy G. Conley & Christian B. Hansen (2011)

**Inference with Dependent Data Using Cluster Covariance Estimators.**
*Journal of Econometrics* 165(2), 137–151.
[DOI/publisher](https://doi.org/10.1016/j.jeconom.2011.01.007) ·
[Author PDF route](https://economics.uwo.ca/people/conley_docs/inferencewithdependentdata.pdf).

- **Access:** author PDF route was found but redirected unsuccessfully in the tool;
  fetch the published PDF if that route does not work tomorrow.
- **Why:** candidate small-number-of-large-spatial-groups alternative, including
  the role of t/F reference distributions.
- **Question:** what must be true of the groups and their parameter estimates?
  Merely replacing a normal critical value with a t critical value is not enough.

### 5. Rustam Ibragimov & Ulrich K. Müller (2010)

**t-Statistic Based Correlation and Heterogeneity Robust Inference.**
*Journal of Business & Economic Statistics* 28(4), 453–468.
[DOI](https://doi.org/10.1198/jbes.2009.08046) ·
[Author PDF](https://www.princeton.edu/~umueller/tstat.pdf).

- **Access:** open author PDF.
- **Why:** independently estimated group parameters provide another scalar-inference
  benchmark, with important independence, normality, and common-target assumptions.
- **Question:** could sufficiently separated, comparable domains support this
  method without changing the target parameter or assuming independent adjacent blocks?

### 6. Timothy G. Conley (1999)

**GMM Estimation with Cross Sectional Dependence.** *Journal of Econometrics*
92(1), 1–45.
[DOI/publisher](https://doi.org/10.1016/S0304-4076(98)00084-0).

- **Access:** abstract/metadata inspected; fetch published full text and original
  implementation notes if available.
- **Why:** foundational estimating-equation, distance, kernel, asymptotic-domain,
  and positive-semidefinite covariance requirements.
- **Question:** which assumptions remain credible after declustering and censoring,
  and which do not survive concentration-dependent sampling of sites?

## Second reading sequence: irregular and preferential sampling

### 7. S. N. Lahiri & Jun Zhu (2006)

**Resampling Methods for Spatial Regression Models under a Class of Stochastic
Designs.** *Annals of Statistics* 34(4), 1774–1813.
[DOI](https://doi.org/10.1214/009053606000000551) ·
[Open preprint](https://arxiv.org/abs/math/0611261) ·
[Publisher PDF](https://projecteuclid.org/journals/annals-of-statistics/volume-34/issue-4/Resampling-methods-for-spatial-regression-models-under-a-class-of/10.1214/009053606000000551.pdf).

- **Access:** open versions linked.
- **Why:** explains why a straightforward transfer of grid block bootstrapping can
  fail for irregular nonuniform sampling; proposes a different blocking mechanism.
- **Question:** can it support our nonlinear censored M-estimator under exogenous
  site selection? Nonuniform sampling and concentration-dependent preferential
  sampling must not be treated as the same assumption.

### 8. Peter J. Diggle, Raquel Menezes & Ting-li Su (2010)

**Geostatistical Inference under Preferential Sampling.** *JRSS Series C*
59(2), 191–232.
[DOI/publisher](https://doi.org/10.1111/j.1467-9876.2009.00701.x) ·
[University-hosted PDF](https://www.stat.ubc.ca/~jim/Diggle2010.pdf).

- **Access:** open university-hosted full-text route found.
- **Why:** a model for the sampling locations is relevant to the severe omitted-site
  bias demonstrated in our simulations; a variance correction alone is insufficient.
- **Question:** what sampling-location information and joint likelihood assumptions
  would be needed for a historical geochemical compilation?
- **Fetch if convenient:** the published discussion/rejoinder and supplements as well
  as the main article, since the page range includes discussion material.

### 9. Rui Qiang & Peter F. Craigmile (2026)

**Accounting for Preferential Sampling in Geostatistical Inference.**
*Environmetrics* 37(2), e70081.
[DOI/publisher](https://doi.org/10.1002/env.70081).

- **Access:** publisher abstract/metadata inspected; full text and supplements requested.
- **Why:** recent joint modeling of locations and responses offers a modern
  comparator to the foundational preferential-sampling model.
- **Question:** which extensions or practical estimation tools are useful for
  uneven historical sampling, and can detection intervals be accommodated?

## Third reading sequence: spatial models with detection limits

### 10. José A. Ordoñez, Dipankar Bandyopadhyay, Victor H. Lachos & Celso R. B. Cabral (2018)

**Geostatistical Estimation and Prediction for Censored Responses.**
*Spatial Statistics* 23, 109–123.
[DOI/publisher](https://doi.org/10.1016/j.spasta.2017.12.001).

- **Access:** publisher abstract/metadata inspected; fetch full text, supplement,
  and reproducible code if available.
- **Why:** a spatial likelihood model preserves censoring while estimating covariance,
  instead of adjusting a working-independence fit after estimation.
- **Question:** how are both mean and covariance uncertainty estimated under SAEM,
  and what is the practical cost at the sample sizes needed for domain analyses?

### 11. Jose A. Ordonez, Christian E. Galarza & Victor H. Lachos (2021)

**Spatial Censored Regression Models in R: The CensSpatial Package.**
arXiv:2110.05570.
[Open preprint](https://arxiv.org/abs/2110.05570) ·
[PDF](https://arxiv.org/pdf/2110.05570).

- **Access:** open software preprint; listed as the verified arXiv version, without
  assuming a later journal publication.
- **Why:** implementation, estimation algorithms, simulation utilities, and influence
  diagnostics could provide an independent model-based comparator.
- **Question:** does its present implementation support our censoring directions,
  observation supports, uncertainty target, and feasible domain sizes?

## Supporting computational and diagnostic references

### 12. Achim Zeileis (2006)

**Object-oriented Computation of Sandwich Estimators.** *Journal of Statistical
Software* 16(9), 1–16.
[DOI/publisher](https://doi.org/10.18637/jss.v016.i09) ·
[Author/package PDF](https://sandwich.r-forge.r-project.org/articles/pdf/sandwich-OOP.pdf).

- **Access:** open article; already used for the joint-score/bread framework.
- **Why:** supports the nonlinear score-sandwich adaptation and independent numerical
  verification of weighting for location and nuisance-scale parameters.

### 13. Timothy G. Conley & Morgan Kelly (2025)

**The Standard Errors of Persistence.** *Journal of International Economics*
153, 104027.
[DOI/publisher](https://doi.org/10.1016/j.jinteco.2024.104027) ·
[2024 working-paper PDF](https://www.ucd.ie/economics/t4media/WP2024_17.pdf) ·
[Replication files](https://doi.org/10.17632/p2yt92jnnr.1).

- **Access:** open working-paper and replication routes found; published article
  is dated 2025 (not the working paper's 2024 date).
- **Why:** spatial trends, spurious relationships, and simulation/placebo diagnostics
  matter if the package expands into element–element or covariate inference.
- **Question:** can useful diagnostics be transferred without confusing a spatial
  trend with stationary local dependence? This is supporting reading, not a
  direct estimator for our censored concentration mean.

## What to return

PDFs plus supplements/replication files are ideal. Use the DOI/title to avoid
mixing a working-paper version with a different final article. Filenames such as
`Bester_etal_2016_fixed_b.pdf` and `Ordonez_etal_2018_censored_spatial.pdf` will
be easy to identify. The access notes are observations from this search, not
promises about access from a different institution or browser.

The list can guide reading while development continues. No paper's abstract or
package example is being treated as a validation of geoCensor's weighted,
censored, preferentially sampled setting.


### Optional online permutation references — November 1, 2026

No blocking paper or supplement request for the current site-block evidence
study. Online methods excerpts were sufficient to implement/check the restricted
randomization machinery; local archival copies would be optional:

- Hemerik & Goeman (2018), Exact testing with random permutations,
  https://doi.org/10.1007/s11749-017-0571-1 (open publisher text).
- Phipson & Smyth (2010), Permutation P-values Should Never Be Zero,
  https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf (author preprint).

These are not reported as received local Papers files or full-PDF reviews.
