# Write an editable or print-ready analysis report

Renders the active Rasch or paired-comparison fit as a self-contained
HTML document, an editable Word document, or a PDF. The report contains
the principal estimates, model-specific tables, diagnostic figures, and
software provenance. Complete machine-readable results remain available
from
[`save_outputs`](https://drjoshmcgrane.github.io/rasch/reference/save_outputs.md).

## Usage

``` r
report_document(
  fit,
  file,
  format = c("auto", "html", "docx", "pdf"),
  title = "Rasch measurement analysis",
  dif = NULL,
  bootstrap = NULL,
  dif_bootstrap = NULL
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md),
  [`rasch_mfrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md),
  [`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md),
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md), or
  [`btl_efrm`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md).

- file:

  Output path ending in `.html`, `.docx`, or `.pdf`.

- format:

  Output format. By default it is inferred from `file`.

- title:

  Report title.

- dif, bootstrap:

  Optional computed
  [`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  and
  [`fit_bootstrap`](https://drjoshmcgrane.github.io/rasch/reference/fit_bootstrap.md)
  results from this fit, rendered as run rather than recomputed at
  defaults.

- dif_bootstrap:

  Optional
  [`dif_bootstrap`](https://drjoshmcgrane.github.io/rasch/reference/dif_bootstrap.md)
  sensitivity analysis from this fit and DIF specification.

## Value

Invisibly, the output path.

## Details

Word and HTML output require Pandoc, supplied with RStudio and available
through rmarkdown. PDF output also requires a LaTeX installation such as
TinyTeX.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- rasch(matrix(rbinom(3000, 1, .5), 300, 10))
report_document(fit, file.path(tempdir(), "analysis.docx"))
} # }
```
