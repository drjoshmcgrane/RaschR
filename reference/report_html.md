# Write a self-contained HTML report of a Rasch analysis

Writes one HTML file containing the summary statistics, diagnostic
tables, and test-level plots. Images and styles are embedded in the
file.

## Usage

``` r
report_html(
  fit,
  file,
  title = "Rasch measurement analysis",
  dpi = 150,
  dif = NULL,
  bootstrap = NULL,
  dif_bootstrap = NULL,
  dimensionality = NULL,
  invariance = NULL,
  subtest = NULL
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- file:

  Path of the HTML file to write.

- title:

  Report title.

- dpi:

  Resolution of the embedded plots.

- dif, bootstrap:

  Optional computed
  [`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  and
  [`fit_bootstrap`](https://drjoshmcgrane.github.io/rasch/reference/fit_bootstrap.md)
  results from this fit, exported as run; the DIF table is otherwise
  recomputed at defaults when the fit carries person factors.

- dif_bootstrap:

  Optional
  [`dif_bootstrap`](https://drjoshmcgrane.github.io/rasch/reference/dif_bootstrap.md)
  sensitivity analysis from this fit and DIF specification.

- dimensionality:

  Optional computed
  [`plot_scree`](https://drjoshmcgrane.github.io/rasch/reference/plot_scree.md)
  result from this fit. Supplying it keeps the table and figure
  identical to the analysis already run.

- invariance:

  Optional computed
  [`frame_invariance`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
  result from an EFRM fit.

- subtest:

  Optional computed
  [`dimensionality_test`](https://drjoshmcgrane.github.io/rasch/reference/dimensionality_test.md)
  result from this fit. Supplying it keeps the item split and bootstrap
  calibration used in the analysis.

## Value

Invisibly, `file`.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(150 * 6, 1, plogis(outer(rnorm(150), d, "-"))), 150, 6)
colnames(X) <- paste0("I", 1:6)
out <- file.path(tempdir(), "report.html")
report_html(rasch(X), out)
```
