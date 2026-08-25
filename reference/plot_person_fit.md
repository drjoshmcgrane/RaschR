# Plot person fit

Person locations against a person fit statistic with the +/- 2.5 band;
persons beyond the band respond erratically (positive) or too
deterministically (negative). The default statistic is the
log-of-mean-square fit residual; `"infit"` and `"outfit"` display the
Wilson–Hilferty standardised mean squares, to which the same band
convention applies.

## Usage

``` r
plot_person_fit(fit, statistic = c("residual", "infit", "outfit"), band = 2.5)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- statistic:

  `"residual"` (the default fit residual), `"infit"`, or `"outfit"`.

- band:

  Acceptance band for the standardised statistic.

## Value

Called for its plotting side effect; invisibly `NULL`.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
colnames(X) <- paste0("I", 1:6)
plot_person_fit(rasch(X))
```
