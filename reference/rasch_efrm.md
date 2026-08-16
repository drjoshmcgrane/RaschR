# Fit the extended frame of reference model

Fits Humphry's extended frame of reference model, in which the unit can
differ across item-set by person-group frames. For item \\i\\ in set
\\s\\ and person \\n\\ in group \\g\\,
\$\$P(X\_{ni}=x)=\frac{\exp\\\rho\_{sg}\[x\theta_n-
\sum\_{k=1}^{x}\delta\_{ik}\]\\}
{\sum\_{y=0}^{m_i}\exp\\\rho\_{sg}\[y\theta_n-
\sum\_{k=1}^{y}\delta\_{ik}\]\\},\qquad \rho\_{sg}=\alpha_s\phi_g.\$\$

## Usage

``` r
rasch_efrm(
  data,
  item_sets,
  groups,
  id = NULL,
  factors = NULL,
  items = NULL,
  n_groups = NULL,
  adjust_N = NA,
  na_codes = -1,
  maxit = 50,
  tol = 1e-07,
  min_link_persons = 30,
  se_method = c("hybrid", "bootstrap"),
  boot_reps = NULL
)
```

## Arguments

- data:

  Persons-by-items data (matrix or data frame, like
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)),
  plus a person-group column.

- item_sets:

  A named list mapping set names to item-column names, or a named
  character vector mapping item names to set names. Items not mentioned
  form their own set `"(rest)"` when a list is given.

- groups:

  Name of the person-group column in `data`, or a vector with one entry
  per person. Several columns define crossed group cells. Their units
  are returned in `phi_table`; `phi_factorial` and `phi_factorial_tests`
  contain the GLS factorial decomposition and omnibus Wald tests.
  Structurally unidentified units are refused. Very imprecise but
  identified units are retained with a warning.

- id, factors, items, n_groups, adjust_N, na_codes:

  As in
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- maxit, tol:

  Outer iteration cap and convergence tolerance of the bilinear pairwise
  stage.

- min_link_persons:

  Minimum number of common persons required for a set pair to contribute
  to the unit linking.

- se_method:

  `"hybrid"` (sandwich + linking bootstrap + delta propagation; fast,
  default) or `"bootstrap"` (full person bootstrap of all stages).

- boot_reps:

  Bootstrap replicates; defaults to 300 for the linking bootstrap and
  200 for the full bootstrap.

## Value

An object of classes `"rasch_efrm"` and `"rasch"`. Model-specific
components include `frames`, `phi_table`, `alpha_table`, `set_table`,
common-unit item and threshold tables, group-specific `score_curves`,
`efrm_vs_rasch`, and `linking`. See the extended frame of reference
vignette for their interpretation.

## Details

The partial credit model holds within each frame in its natural unit.
Person-group units \\\phi_g\\ and centred set thresholds are estimated
by within-frame pairwise conditional maximum likelihood. Item-set units
\\\alpha_s\\ and set locations are then estimated from persons common to
linked sets by true-score variance ratios, with the true-score variance
recovered by a truncated-score-moment correction: the mean and variance
of the weighted likelihood score map over the non-extreme scores are
exact functions of the person location given the fitted thresholds, and
their person-distribution expectations are estimated through score
weights that are unbiased for any person distribution because the raw
score is sufficient. The naive \\var(\hat u) - mean(SE^2)\\ correction
is badly calibrated on short tests (the reported error variance
overstates the actual one and weighted likelihood shrinkage makes the
errors covary negatively with the locations); its residual distortion
biased the log unit ratio upward by about 0.05 at eight dichotomous
items per set, confirmed against an external TAM 2PL slope-group anchor,
while the corrected estimator is unbiased there. The linking graph must
connect all sets to a common scale.

The default hybrid standard errors combine the pairwise Godambe
covariance, a person bootstrap for set linking, and delta-method
propagation. Each linking replicate also draws the within-frame
parameters from their joint stage-one covariance: without that redraw
the set-unit standard errors understate by about 20% and the unit tests
reject a true null at 9-10%; with it they reject at 4.9% over 1,200
simulated replicates, stable across item counts, sample sizes,
imbalance, and weak linking, and matching a full-bootstrap benchmark.
The corrected set-unit estimator holds this calibration across designs:
null size 3–5% with 93–99% coverage over 5–15 items per set, unit ratios
1–2, partial credit items, booklet missingness, pairwise-only person
overlap, and person skewness to 2.8 – where it stays unbiased while a
normal-population MML anchor drifts. With `se_method = "bootstrap"`, the
complete model is refitted to each person resample and all reported
covariance comes from the bootstrap distribution.

The dichotomous model and the theory of frame-dependent units follow
Humphry (2005) and Humphry and Andrich (2008). The item-set linking step
is an error-corrected method-of-moments implementation of the
variance-ratio argument in Humphry (2005), rather than the likelihood
proposed in section 5.3 of that thesis. The polytomous, multigroup, and
crossed-frame forms are extensions implemented in this package. For
these designs, the full bootstrap gives the least conditional account of
uncertainty.

## References

Andrich, D. (1982). An extension of the Rasch model for ratings
providing both location and dispersion parameters. Psychometrika, 47(1),
105–113.

Andrich, D. and Luo, G. (2003). Conditional pairwise estimation in the
Rasch model for ordered response categories using principal components.
Journal of Applied Measurement, 4(3), 205–221.

Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
Measuring in the Educational, Social and Health Sciences. Springer.

Humphry, S. M. (2005). Maintaining a Common Arbitrary Unit in Social
Measurement. PhD thesis, Murdoch University.

Humphry, S. M. (2010). Modeling the effects of person group factors on
discrimination. Educational and Psychological Measurement, 70(2),
215–231.

Humphry, S. M. (2012). Item set discrimination and the unit in the Rasch
model. Journal of Applied Measurement, 13(2), 165–180.

Montuoro, P. and Humphry, S. M. (2024). Modeling the effect of reading
item clarity on item discrimination. Journal of Applied Measurement,
24(3/4), 121–132.

Humphry, S. M. and Andrich, D. (2008). Understanding the unit in the
Rasch model. Journal of Applied Measurement, 9(3), 249–264.

## See also

[`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md),
[`rasch_mfrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md),
[`test_information`](https://drjoshmcgrane.github.io/rasch/reference/test_information.md),
and
[`simulate_efrm`](https://drjoshmcgrane.github.io/rasch/reference/simulate_efrm.md).

## Examples

``` r
# \donttest{
set.seed(1); Np <- 400
simP <- function(th, tau, r) { x <- 0:length(tau)
  p <- exp(r * (x * th - c(0, cumsum(tau)))); p / sum(p) }
grp <- rep(c("A", "B"), each = Np / 2)
phi <- c(A = 0.8, B = 1.25)
d <- seq(-1.5, 1.5, length.out = 10)
theta <- rnorm(Np)
X <- sapply(seq_along(d), function(i) sapply(seq_len(Np), function(n)
  sample(0:1, 1, prob = simP(theta[n], d[i], phi[grp[n]]))))
colnames(X) <- sprintf("I%02d", seq_along(d))
fit <- rasch_efrm(data.frame(X, grp = grp), item_sets = list(core = colnames(X)),
                  groups = "grp")
fit$phi_table
#>   group       phi se_log_phi
#> 1     A 0.8422988 0.04222911
#> 2     B 1.1872271 0.04222911
# }
```
