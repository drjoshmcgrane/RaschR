# Bootstrap null distribution for the item fit statistics

Every item fit statistic
[`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
reports is computed at estimated person locations and referred to a
distribution derived as though those locations were known, and each is
miscalibrated by an amount that grows with the sample. The item-trait
chi-square, whose class intervals are formed by selecting on the
estimates, grows about linearly in N against a fixed reference. The fit
residual's null SD runs from about 0.71 at 250 persons to 1.00 at 4,000,
so the conventional \\\pm\\2.5 cut means different things at different
sizes; infit z beyond 1.96 flags 12 of correctly fitting items at 250
persons and 69 replaces the reference distribution for all of them at
once.

## Usage

``` r
fit_bootstrap(
  fit,
  B = 200,
  theta = c("resample", "fixed", "normal"),
  workers = 4L,
  seed = NULL
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).
  Extended-frame, many-facet and explanatory fits are not supported.

- B:

  Number of bootstrap replicates.

- theta:

  How to generate person locations for a replicate: `"resample"` (the
  default) resamples the person estimates, `"fixed"` reuses them as they
  stand, and `"normal"` draws from a normal with the error-corrected
  variance the fitted model implies. The three calibrate the tests
  within Monte Carlo error of each other. The first two generate at the
  spread the data showed; `"normal"` generates at a narrower one, which
  leaves an extreme category unvisited often enough on polytomous items
  to lose a substantial share of replicates, so it is not the default.

- workers:

  Number of parallel bootstrap workers. The default is four, reduced
  when fewer physical cores are available or the R process has a lower
  system limit. Starting them costs about half a second, which even the
  smallest useful run earns back; a larger count pays at large samples
  or large `B`, where eight workers run about five times faster than
  one. Person locations and the per-replicate seeds are drawn before
  distribution, so a fixed seed gives the same result for any worker
  count.

- seed:

  Optional seed. The caller's random stream is restored on exit.

## Value

A list with `items`, one row per item carrying each observed statistic
beside its bootstrap probability and that probability's Holm adjustment
across items; `total` for the whole-test chi-square and for the mean and
SD of the item fit residuals, each against its own bootstrap null;
`replicates`, the raw replicated statistics as one matrix per statistic;
`B` requested and `B_used` (replicates that estimated); and the `theta`
scheme.

## Details

Each replicate generates responses from the fitted thresholds at person
locations drawn under `theta`, then re-estimates the item side,
re-estimates person locations, re-forms class intervals and re-takes
residuals exactly as the observed fit did, so the same bias enters the
null as entered the observed statistics. One set of replicates serves
every statistic, since a replicate must refit the whole model in any
case. The class-interval count is held at the fitted value, and observed
missing data are carried into every replicate.

The chi-square is read in its upper tail alone. The fit residual, the
mean squares and the standardised statistics depart in both directions –
above for an item flatter than the model predicts, below for a steeper
one – and both are misfit, so their p-values are equal-tailed. A
bootstrap p-value is `(1 + r) / (1 + B)`, so it is never zero and its
resolution is `1 / (1 + B)`, halved again for the two-sided statistics:
`B = 200` cannot distinguish anything below about .01 two-sided, and
testing at .01 or with a familywise adjustment across many items wants
more.

Calibration is not power. A flatter-than-Rasch item carries less of the
class-interval selection bias than a fitting item does, so its
chi-square comes out *smaller* than a fitting item's and no reference
distribution can make it significant; the fit residual detects that
departure readily. Calibrating both is what covers the two blind spots.
A null estimated from data that contain a misfitting item is also mildly
contaminated by it, which leaves the remaining items flagging somewhat
above nominal.

## References

Andrich, D. and Marais, I. (2019) *A Course in Rasch Measurement
Theory*. Springer.

## See also

[`chisq_detail`](https://drjoshmcgrane.github.io/rasch/reference/chisq_detail.md)
for the class-interval breakdown behind one item's chi-square.

## Examples

``` r
set.seed(1)
d <- seq(-1.5, 1.5, length.out = 6)
X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
colnames(X) <- paste0("I", 1:6)
bs <- fit_bootstrap(rasch(X), B = 49, seed = 1)
bs$items[c("item", "chisq", "chisq_p_boot", "fit_resid", "fit_resid_p_boot")]
#>  item chisq chisq_p_boot fit_resid fit_resid_p_boot
#>    I1 2.839        0.860    -0.090            0.840
#>    I2 3.867        0.660     0.679            0.440
#>    I3 5.085        0.720     0.313            1.000
#>    I4 6.088        0.500     1.172            0.400
#>    I5 3.258        0.900     0.023            0.880
#>    I6 7.682        0.180    -0.200            0.680
```
