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
  bootstrap = NULL
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
