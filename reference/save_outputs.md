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
  item_plots = TRUE,
  dif = NULL,
  bootstrap = NULL
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

- dif:

  Optional
  [`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  result to export as computed — an application analysis carries the DIF
  model the analyst chose, which a default recomputation would silently
  replace. `NULL` computes the default when the fit carries person
  factors.

- bootstrap:

  Optional
  [`fit_bootstrap`](https://drjoshmcgrane.github.io/rasch/reference/fit_bootstrap.md)
  result; its item table and whole-test readings join the exported
  tables.

## Value

Invisibly, the vector of files written.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(150 * 6, 1, plogis(outer(rnorm(150), d, "-"))), 150, 6)
colnames(X) <- paste0("I", 1:6)
out <- file.path(tempdir(), "rasch-out")
save_outputs(rasch(X), out, formats = "png", item_plots = FALSE, dpi = 96)
```
