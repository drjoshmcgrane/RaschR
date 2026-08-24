# Write a self-contained HTML report of a Rasch analysis

Writes one HTML file containing the summary statistics, diagnostic
tables, and test-level plots. Images and styles are embedded in the
file.

## Usage

``` r
report_html(fit, file, title = "Rasch measurement analysis", dpi = 150)
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
