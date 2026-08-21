# Differential item functioning by residual analysis of variance

Tests uniform and non-uniform DIF by analysing each item's standardised
residuals over person factors and trait class intervals (Andrich and
Marais 2019, ch. 16). Several person factors are fitted jointly. The
function also supports designs containing both between-person and
within-person factors.

## Usage

``` r
dif_anova(
  fit,
  factors = NULL,
  n_groups = NULL,
  p_adjust = "BH",
  alpha = 0.05,
  effects = c("main", "factorial"),
  sizes = FALSE,
  id = NULL,
  within = NULL,
  pool_facets = TRUE
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md),
  [`rasch_mfrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md),
  or
  [`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md).

- factors:

  A vector (one factor), a data frame of person factors, or a character
  vector naming factor columns nominated in the fit. Defaults to every
  factor stored in the fit.

- n_groups:

  Number of trait class intervals. The default uses the smallest joint
  factor cell to retain about 30 expected responses per interval and
  cell, with between 2 and 10 intervals. The selected value is returned
  in `n_groups`.

- p_adjust:

  Multiplicity adjustment over all item-by-term tests; default `"BH"`.

- alpha:

  Significance level applied to the adjusted probabilities.

- effects:

  `"main"` (default) models several factors additively (each factor's
  main effect and its class-interval interaction, but no
  factor-by-factor terms); `"factorial"` also crosses the factors with
  each other. Immaterial with a single factor.

- sizes:

  If `TRUE`, refit each flagged item-term and calculate pairwise DIF
  differences in logits using
  [`dif_size`](https://drjoshmcgrane.github.io/rasch/reference/dif_size.md).

- id:

  Person identifier for stacked or repeated-measures data. It may be a
  column name stored in the fit or a vector with one value per row; by
  default the identifier carried by the fit is used.

- within:

  Names of within-person factors. With repeated identifiers, varying
  factors are detected automatically when this is omitted. See Details
  for the mixed-design analysis.

- pool_facets:

  For MFRM fits: pool residuals to the underlying items (the default),
  so DIF is tested per item rather than per item-by-facet cell; `FALSE`
  tests each cell as its own item. EFRM response cells are always pooled
  by item, so this argument does not alter EFRM fits. Ignored for
  ordinary fits.

## Value

A list with:

- `summary`:

  One row per item and group term, containing the uniform and
  non-uniform tests, partial eta-squared, adjusted probabilities, DIF
  flags, and supersession flag.

- `terms`:

  The complete item-wise analysis-of-variance tables.

- `tukey`:

  Residual-mean Tukey comparisons retained for compatibility in
  between-person designs. Use `posthoc` for logit-scale follow-ups.

- `sizes`:

  When requested, pairwise logit differences for the significant,
  non-superseded item-terms.

- `posthoc`:

  When `sizes = TRUE`, marginal pairwise differences for main effects
  and difference-in-differences magnitudes for interactions, calculated
  by
  [`dif_posthoc`](https://drjoshmcgrane.github.io/rasch/reference/dif_posthoc.md).

The remaining components record the factors, class intervals,
adjustment, significance level, and design settings.

## Details

With one factor \\G\\ and class interval \\C\\, the residual model is
\$\$z=\mu+G+C+G\mathbin{:}C+\varepsilon.\$\$ The factor term tests
uniform DIF and its interaction with class interval tests non-uniform
DIF. With several factors, `effects = "main"` fits
`(f1 + f2 + ...) * ci`; `effects = "factorial"` also includes
factor-by-factor interactions. Type II sums of squares are used. The
multiplicity adjustment covers all item-by-DIF-term tests, including
both uniform and non-uniform DIF; the class-interval main effect is a
nuisance term and is not included.

When identifiers repeat, the person is the unit of analysis.
Between-person terms use person means and the between-person error
stratum. Within-person terms use orthonormal contrasts of person-by-cell
means. A Greenhouse–Geisser correction is applied to within-person
factors with more than two levels. Persons missing a required cell are
excluded from the corresponding within-person test. In incomplete mixed
designs, within-cell effects are removed before the between-person
analysis.

A significant higher-order factor term supersedes its component terms in
the summary. For EFRM fits, frame-defining factors are excluded because
they define the model rather than a separate DIF contrast; testing such
a factor means stepping outside the model, which is what
[`frame_invariance`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
does. MFRM residuals are pooled to underlying items unless
`pool_facets = FALSE`. EFRM response cells are always pooled by item;
the frame-defining factors remain excluded. Inference is available only
from a converged calibration.

## References

Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery
rate: a practical and powerful approach to multiple testing. Journal of
the Royal Statistical Society: Series B, 57(1), 289–300.

Hagquist, C. and Andrich, D. (2017). Recent advances in analysis of
differential item functioning in health research using the Rasch model.
Health and Quality of Life Outcomes, 15, 181.

Maxwell, S. E. and Delaney, H. D. (2004). Designing Experiments and
Analyzing Data: A Model Comparison Perspective (2nd ed.). Lawrence
Erlbaum.

## See also

[`dif_size`](https://drjoshmcgrane.github.io/rasch/reference/dif_size.md),
[`dif_contrasts`](https://drjoshmcgrane.github.io/rasch/reference/dif_contrasts.md),
and
[`resolve_dif`](https://drjoshmcgrane.github.io/rasch/reference/resolve_dif.md);
and
[`frame_invariance`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
for the frame-defining factor this function excludes.

## Examples

``` r
set.seed(1); n <- 800
d <- seq(-1.5, 1.5, length.out = 6)
g1 <- rep(c("a", "b"), each = n / 2)
g2 <- rep(c("x", "y"), times = n / 2)
sh <- matrix(0, n, 6); sh[g1 == "b", 2] <- 0.8
X <- matrix(rbinom(n * 6, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 6)
colnames(X) <- paste0("I", 1:6)
fit <- rasch(data.frame(X, g1 = g1, g2 = g2), factors = c("g1", "g2"))
dif_anova(fit)$summary
#>  item term F_uniform p_uniform p_uniform_adj eta2_uniform uniform_DIF
#>    I1   g1     0.612     0.434         0.651        0.001            
#>    I1   g2     1.341     0.247         0.494        0.002            
#>    I2   g1    22.600   < 0.001       < 0.001        0.031           *
#>    I2   g2     2.931     0.087         0.419        0.004            
#>    I3   g1     2.623     0.106         0.423        0.004            
#>    I3   g2     0.324     0.570         0.719        0.000            
#>    I4   g1     2.279     0.132         0.451        0.003            
#>    I4   g2     1.428     0.232         0.494        0.002            
#>    I5   g1     0.833     0.362         0.620        0.001            
#>    I5   g2     1.384     0.240         0.494        0.002            
#>    I6   g1     0.680     0.410         0.651        0.001            
#>    I6   g2     3.954     0.047         0.283        0.006            
#>  F_nonuniform p_nonuniform p_nonuniform_adj eta2_nonuniform nonuniform_DIF
#>         0.471        0.757            0.865           0.003               
#>         0.117        0.977            0.977           0.001               
#>         2.571        0.037            0.283           0.014               
#>         0.224        0.925            0.965           0.001               
#>         1.216        0.303            0.559           0.007               
#>         0.582        0.676            0.811           0.003               
#>         0.841        0.499            0.705           0.005               
#>         1.582        0.177            0.494           0.009               
#>         0.754        0.556            0.719           0.004               
#>         0.386        0.819            0.893           0.002               
#>         3.544        0.007            0.085           0.020               
#>         1.406        0.230            0.494           0.008               
#>  superseded
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            

# \donttest{
# Mixed design: group is between persons and occasion is within persons.
N <- 320; theta <- rnorm(N); group <- rep(c("A", "B"), each = N / 2)
make_wave <- function(occasion_shift) {
  shift <- matrix(0, N, 6)
  shift[group == "B", 2] <- 0.9
  shift[, 5] <- occasion_shift
  matrix(rbinom(N * 6, 1,
         plogis(outer(theta, d, "-") - shift)), N, 6)
}
Xm <- rbind(make_wave(0), make_wave(1.0))
colnames(Xm) <- paste0("I", 1:6)
repeated <- data.frame(Xm, group = rep(group, 2),
                       occasion = rep(c("T1", "T2"), each = N))
mixed_fit <- rasch(repeated, id = rep(seq_len(N), 2),
                   factors = c("group", "occasion"))
mixed_dif <- dif_anova(mixed_fit, within = "occasion")
subset(mixed_dif$summary, uniform_DIF | nonuniform_DIF)
#>  item     term F_uniform p_uniform p_uniform_adj eta2_uniform uniform_DIF
#>    I2    group    24.646   < 0.001       < 0.001        0.075           *
#>    I5 occasion    19.702   < 0.001       < 0.001        0.061           *
#>  F_nonuniform p_nonuniform p_nonuniform_adj eta2_nonuniform nonuniform_DIF
#>         1.368        0.245            0.653           0.018               
#>         0.196        0.940            0.950           0.003               
#>  superseded
#>            
#>            
# }
```
