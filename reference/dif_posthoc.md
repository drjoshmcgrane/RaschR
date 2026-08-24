# Pairwise follow-up comparisons for a DIF term

Resolves one item's locations over the complete person-factor design and
follows up a selected main effect or interaction. Main effects are
pairwise marginal differences. Interactions are differences between
those differences, providing a logit-scale magnitude for the interaction
itself.

## Usage

``` r
dif_posthoc(
  fit,
  item,
  term,
  factors = NULL,
  within = NULL,
  id = NULL,
  p_adjust = "holm",
  alpha = 0.05,
  flag_logits = 0.5,
  min_n = 20
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md) or
  [`rasch_mfrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md).
  EFRM fits are excluded because resolved comparisons would discard
  their frame units.

- item:

  Item name or index.

- term:

  A factor name for a main effect, or a character vector of factor names
  for an interaction. A single colon-separated string is also accepted
  when the factor names themselves contain no colon.

- factors:

  The complete person-factor design, specified as for
  [`dif_contrasts`](https://drjoshmcgrane.github.io/rasch/reference/dif_contrasts.md).
  Other factors are retained when calculating marginal comparisons.

- within:

  Within-person factor names, specified as for
  [`dif_contrasts`](https://drjoshmcgrane.github.io/rasch/reference/dif_contrasts.md).

- id:

  Person identifiers, specified as for
  [`dif_contrasts`](https://drjoshmcgrane.github.io/rasch/reference/dif_contrasts.md).

- p_adjust:

  Familywise adjustment over this post-hoc family; default `"holm"`.

- alpha:

  Significance level for adjusted probabilities.

- flag_logits:

  Absolute logit magnitude flagged as practically important.

- min_n:

  Minimum responders required in a resolved design cell.

## Value

An object of class `"rasch_dif_posthoc"`, extending the
[`dif_contrasts`](https://drjoshmcgrane.github.io/rasch/reference/dif_contrasts.md)
result. Its `table` contains the pairwise marginal differences or
interaction contrasts, with logit estimates, standard errors where
available, confidence intervals, raw and adjusted probabilities, and
statistical and practical flags.

## Details

For levels \\a,b\\ of one factor, the comparison is
\$\$\Delta\_{ba}=\bar\delta_b-\bar\delta_a,\$\$ where the bars average
equally over complete cells of the other nominated factors. For a
two-factor interaction, levels \\a,b\\ and \\c,d\\ give
\$\$\Delta\_{ba\mathbin{:}dc}=
(\delta\_{bd}-\delta\_{ad})-(\delta\_{bc}-\delta\_{ac}).\$\$
Higher-order interactions use the corresponding tensor-product contrast.
Standard errors use the full covariance of the resolved locations.

This is the follow-up to a significant DIF term with more than two
levels. It reports effects in Rasch logits, adjusts the chosen family of
comparisons, and uses person-level scores with the same equal-cell
marginal weights in repeated-measures designs.

## References

Holm, S. (1979). A simple sequentially rejective multiple test
procedure. Scandinavian Journal of Statistics, 6(2), 65–70.

## See also

[`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md),
[`dif_size`](https://drjoshmcgrane.github.io/rasch/reference/dif_size.md),
and
[`dif_contrasts`](https://drjoshmcgrane.github.io/rasch/reference/dif_contrasts.md).

## Examples

``` r
set.seed(1); n <- 800
g <- factor(rep(c("A", "B", "C", "D"), each = n / 4))
sex <- factor(rep(c("female", "male"), length.out = n))
d <- seq(-1.5, 1.5, length.out = 6)
sh <- matrix(0, n, 6); sh[g == "D", 2] <- 0.8
X <- matrix(rbinom(n * 6, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 6)
colnames(X) <- paste0("I", 1:6)
fit <- rasch(data.frame(X, group = g, sex = sex),
             factors = c("group", "sex"))
dif_posthoc(fit, "I2", term = "group")
#> DIF follow-up for group (pairwise marginal difference; holm)
#>  item contrast estimate    se statistic p_adj significant practical
#>    I2    B - A   -0.314 0.270    -1.164 0.619                      
#>    I2    C - A    0.028 0.268     0.104 0.917                      
#>    I2    D - A    0.596 0.262     2.276 0.114                     *
#>    I2    C - B    0.342 0.271     1.264 0.619                      
#>    I2    D - B    0.910 0.265     3.437 0.004           *         *
#>    I2    D - C    0.568 0.264     2.155 0.125                     *
```
