# Fit the extended frame of reference model

Estimates Humphry's extended frame of reference model, in which the unit
of the latent scale differs across frames (item-set by person-group
cells). For item \\i\\ in set \\s\\ and person \\n\\ in group \\g\\, the
response model is \$\$P(X\_{ni}=x)=\frac{\exp\\\rho\_{sg}\[x\theta_n-
\sum\_{k=1}^{x}\delta\_{ik}\]\\}
{\sum\_{y=0}^{m_i}\exp\\\rho\_{sg}\[y\theta_n-
\sum\_{k=1}^{y}\delta\_{ik}\]\\},\qquad \rho\_{sg}=\alpha_s\phi_g.\$\$
Within frames the partial credit model holds in the frame's natural
unit, so item thresholds and the person group units `phi` are estimated
by within-frame pairwise conditional maximum likelihood (the person
parameter cancels; Andrich and Luo 2003), jointly across frames through
the sets shared by several groups. Item-set units `alpha` and set
locations are then estimated from persons common to pairs of sets, using
error-corrected true-score variances, and reconciled over the linking
graph by weighted least squares. Everything is reported in a common
arbitrary unit, and the returned object is also a full
[`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md) fit
at the item-by-group level, so the package's diagnostic tables and plots
apply.

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
  per person. Several column names may be given: the frames are then
  their crossed cells, per-cell units appear in `phi_table`, and a
  factorial decomposition of the cell units (sum-coded main effects, and
  the interaction when every cell is observed) is returned in
  `phi_factorial`. The decomposition is a generalised least-squares fit
  of the cell log-units using their joint covariance (bootstrap
  replicates when available, otherwise the analytic centred covariance,
  inverted spectrally along its identified directions); coefficient rows
  are descriptive, and inference is carried by the
  multi-degree-of-freedom Wald test per term in `phi_factorial_tests`.
  Group units are checked for identification on the joint information: a
  flat direction along a unit (structural non-identification) is refused
  with an error naming the group, since every common-unit quantity would
  silently depend on it. A unit whose analytic standard error exceeds 5
  log-units (uncertain beyond a factor of about 150) is practically
  uninformative but not structurally unidentified: its estimate is kept
  for sensitivity work, with a warning and a note. Weakly identified
  units with real threshold spread are kept, with standard errors that
  say how weak they are.

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

An object of classes `"rasch_efrm"` and `"rasch"`. In addition to the
standard components (computed over item-by-group virtual columns with
the frame units carried in `disc`), it has `frames` (one row per frame:
units, origin, pooled fit; under `se_method = "bootstrap"` the
frame-level `se_log_rho` comes from the joint replicate draws of
`log(alpha) + log(phi)`, capturing cross-stage dependence, while the
hybrid fallback combines the stagewise errors as if uncorrelated),
`phi_table`, `alpha_table`, `set_table`, `item_arbitrary` and
`thresholds_arbitrary` (the structural parameters in the common unit),
`score_curves` (per-group score-to-measure curves, replacing the
raw-score table), `efrm_vs_rasch` (fit comparison against the equal-unit
model on the same conditional information, omnibus Wald tests for the
unit families, and Holm-adjusted exploratory unit contrasts), and
`linking` (the linking evidence).

## Details

Humphry (2005) states the model for dichotomous responses and names it:
the thesis defines an extended frame of reference (EFR) as the union of
two or more compatible frames and calls the model over it the extended
frame of reference model. The journal statements use different
terminology for the same machinery – Humphry and Andrich (2008) present
the frame unit as a scale parameter, and Humphry (2009), which develops
and applies the person-group side, calls the model the logistic
measurement function: with a single item set, the model here is exactly
that article's model, with `phi` its person-group discrimination
parameter under the same product-one identification. The polytomous form
fitted here, with the frame unit multiplying the whole exponent over the
item's partial-credit thresholds, is this package's extension of that
statement. It is the form characterised by preserving the two properties
the model's logic rests on: the partial credit model holds within every
frame in the frame's natural unit (so the pairwise conditional
cancellation remains valid), and the weighted score remains sufficient
for the person parameter. It reduces exactly to the dichotomous model
when items are scored 0/1 and to the ordinary partial credit model when
all units equal one. One interpretive consequence: category widths in
natural units scale with the frame unit, so a high-unit frame makes
proportionally sharper category distinctions; frame-level fit and the
per-frame category curves are where a violation of this would appear.

Estimation order: the within-frame pairwise stage establishes the
centred set thresholds and the person-group units `phi`; the person-side
linking stage then establishes the item-set units `alpha` and set
locations. The reported item parameters and all person measures are
computed only after every unit is established: item thresholds are
mapped into the common arbitrary unit using `alpha` and the set
locations, and person measures are weighted-score weighted likelihood
estimates evaluated under the final units `rho = alpha * phi`. The
per-frame person estimates used inside the linking stage are interim
quantities for the unit ratios only and are discarded. The within-frame
stage needs no re-estimation once `alpha` is known, because the pairwise
likelihood is invariant to the within-set rescaling that `alpha`
represents; the units' own uncertainty is reported in `alpha_table` and
`phi_table` and folded into the common-unit standard errors as described
below.

Standard errors: under `se_method = "hybrid"` (default) the group units
carry sandwich standard errors from the pairwise stage; the set units
carry standard errors from a linking-stage bootstrap in which each
replicate resamples persons and also redraws the within-frame thresholds
and group units jointly from their estimated stage-1 covariance before
rebuilding the person estimates. The redraw matters because the set unit
is a scale: error in the estimated threshold spread moves every person
estimate's variance coherently, which person resampling alone cannot see
(without it the log-alpha standard error understates by about 20 in
simulation and the unit tests reject at 9-10 they reject at 4.9 counts,
sample sizes, imbalance, and weak linking). The unit uncertainty is
propagated into the common-unit threshold and item standard errors by
the delta method, treating the stage-1 and person-side linking
information as independent – the one remaining approximation of the
hybrid method. Under `se_method = "bootstrap"` all stages are
re-estimated on `boot_reps` person resamples and every standard error
and the threshold covariance come from the replicate spread; slower, but
captures all cross-dependencies jointly.

Relation to Humphry (2005): the within-frame stage follows the thesis's
conditional separation logic. Humphry (2012) is the published statement
of the item-set discrimination side, with the score vector across item
sets as the sufficient statistic; Montuoro and Humphry (2024) apply it
with sets formed a priori by qualitative item review and set units
estimated from ratios of person-location standard deviations across
sets. The linking stage implemented here is that estimator in
error-corrected method-of-moments form, based on the true-score variance
ratios in equations 2.28–2.29 of the thesis (after Andrich 1982) – the
correction removes the attenuation that raw estimate-standard-deviation
ratios carry. It is not the distinct likelihood equation proposed in
section 5.3. The multigroup, polytomous, and crossed-frame
implementation is therefore an experimental package extension whose
sampling performance should be checked for the intended design,
preferably with the full person bootstrap. The standard errors go
further than the thesis's section 5.4, which inverts each diagonal
element of the joint-likelihood information separately and therefore
conditions on the remaining parameters, including the person locations,
being treated as known. Here full covariance matrices are used
throughout; the item-side covariance carries the Godambe sandwich
correction required for a pairwise composite likelihood; the unit
uncertainty that the thesis's transformation treats as fixed is
propagated into the common-unit parameters; and resampling replaces
analytic plug-in variances for the person-side linking stage.

Measurement-theoretic status: within every frame the model is strictly
Rasch, with person-free item comparisons by conditioning. Across frames
it is an argued extension of the theory of the unit (Humphry 2005;
Humphry and Andrich 2008): on this account the unit was always a
frame-dependent empirical property that the ordinary model leaves
implicit, and the extension makes it explicit; the orthodox reading of
Rasch measurement contests this, and applied reports should present it
as an extension rather than settled doctrine. Two concessions are
intrinsic to the model rather than to this implementation: the item-set
units are identified only from the person side (their conditional
identification is impossible, as documented above), so that step uses
distributional information; and person measures rest on weighted-score
sufficiency with estimated weights, whose uncertainty is propagated
rather than ignored.

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

Humphry, S. M. (2009). Modeling the effects of person group factors on
discrimination. Educational and Psychological Measurement, 70(2),
215–231.

Humphry, S. M. (2012). Item set discrimination and the unit in the Rasch
model. Journal of Applied Measurement, 13(2), 165–180.

Montuoro, P. and Humphry, S. M. (2024). Modeling the effect of reading
item clarity on item discrimination. Journal of Applied Measurement,
24(3/4), 121–132.

Humphry, S. M. and Andrich, D. (2008). Understanding the unit in the
Rasch model. Journal of Applied Measurement, 9(3), 249–264.

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
