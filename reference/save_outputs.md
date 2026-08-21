# Save the outputs of a Rasch analysis

Writes the summary, estimates, diagnostic tables, person measures, and
model-specific results as CSV. Plots are written as PNG and, optionally,
PDF, together with a plain-text analysis summary. For MFRM and EFRM
fits, item estimates and response-cell diagnostics are saved separately.

## Usage

``` r
save_outputs(
  fit,
  dir,
  formats = c("png", "pdf"),
  width = 9,
  height = 6,
  dpi = 300,
  item_plots = TRUE
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- dir:

  Output directory; created if absent.

- formats:

  Plot formats, any of `"png"` and `"pdf"`.

- width, height:

  Plot size in inches.

- dpi:

  PNG resolution.

- item_plots:

  Also write the per-item plot set (one ICC, category curve, threshold
  curve, and frequency chart per item).

## Value

Invisibly, the vector of files written.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
colnames(X) <- paste0("I", 1:6)
out <- file.path(tempdir(), "rasch-out")
save_outputs(rasch(X), out, formats = "png", item_plots = FALSE)
```
