# Fit an explanatory comparative judgement model

Constrains Bradley–Terry–Luce object locations to linear functions of
observed object characteristics. The formulation applies to dichotomous
and ordered comparative judgements; ordered-response thresholds retain
the structure selected in `thresholds`.

## Usage

``` r
btl_explanatory(
  data,
  predictors,
  formula,
  object_a,
  object_b,
  winner = NULL,
  response = NULL,
  margin = NULL,
  judge = NULL,
  count = NULL,
  order = NULL,
  position = FALSE,
  ties = c("drop", "half", "error"),
  thresholds = c("free", "pc"),
  maxit = 60,
  tol = 1e-08
)
```

## Arguments

- data:

  A data frame with one comparison per row.

- predictors:

  Data frame with one row per object, an `object` column, and the
  predictors named in `formula`. Column names must be unique.

- formula:

  One-sided explanatory formula, including selected interactions if
  required.

- object_a, object_b:

  Names of the columns holding the two objects compared. Columns used
  for the comparison roles must be distinct.

- winner:

  Name of the column holding the winner of each row: its value must
  equal one of the two objects. `"tie"` and `"draw"` mark ties. Do not
  supply both `winner` and `response`.

- response:

  Optional ordered response favouring `object_a` over `object_b`: an
  ordered factor from least to greatest preference for `object_a`, or
  integer scores `0..m`.

- margin:

  Optional ordered margin-of-victory column, combined with `winner` to
  construct an orientation-invariant response. Use an ordered factor
  (levels from the smallest to largest margin) or a positive numeric
  magnitude. Margins on ties and rows excluded from the analysis are
  ignored.

- judge:

  Optional name of a judge column; enables the judge fit table and
  clusters the sandwich standard errors by judge.

- count:

  Optional name of a column of replication counts (a row standing for
  several identical comparisons). Counts greater than one cannot be
  combined with `order`, because a compressed row does not retain the
  sequence of the comparisons it represents.

- order:

  Optional column giving each judge's comparison sequence; requires
  `judge`. See Details. Incompatible with `ties = "half"`.

- position:

  If `TRUE`, estimate a first-presentation advantage, treating
  `object_a` as the first object in each comparison.

- ties:

  How to treat ties in the dichotomous analysis: `"drop"` (default,
  removed with a note), `"half"` (half a win each way, a common
  pragmatic device; the two halves remain one sampling unit in the
  sandwich because they are not independent Bernoulli trials), or
  `"error"`. With polytomous responses, code ties as a middle category
  instead.

- thresholds:

  `"free"` (default) estimates every symmetric threshold; `"pc"` retains
  only the symmetric spread component.

- maxit, tol:

  Newton-Raphson iteration cap and convergence tolerance.

## Value

An object of class `"rasch_btl_explanatory"`, inheriting from
`"rasch_btl"`.

## Details

For objects \\a\\ and \\b\\, \$\$\log\\P(a \succ b)/P(b \succ
a)\\=\beta_a-\beta_b,\qquad \beta_i=\mathbf
z_i^{T}\boldsymbol\gamma.\$\$ The scale origin is fixed at mean object
location zero. Numeric predictors are continuous, unordered factors are
categorical, and ordered factors use successive contrasts between
adjacent levels. Character predictors are converted to unordered
factors. Selected interactions may be included in `formula`. A free
calibration is retained for
[`explanatory_test()`](https://drjoshmcgrane.github.io/rasch/reference/explanatory_test.md).
Standard errors use the same sandwich covariance as
[`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md); when
judges are identified, coefficient tests use the judge-clustered
covariance and a \\t\\ reference with judge-cluster degrees of freedom.
Holm adjustment covers the coefficient family.

## References

Bradley, R. A. and Terry, M. E. (1952). Rank analysis of incomplete
block designs: I. The method of paired comparisons. Biometrika, 39,
324–345.

Fischer, G. H. (1973). The linear logistic test model as an instrument
in educational research. Acta Psychologica, 37, 359–374.

## See also

[`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md),
[`explanatory_test`](https://drjoshmcgrane.github.io/rasch/reference/explanatory_test.md),
[`explanatory_diagnostics`](https://drjoshmcgrane.github.io/rasch/reference/explanatory_diagnostics.md),
and
[`relax_btl_explanatory`](https://drjoshmcgrane.github.io/rasch/reference/relax_btl_explanatory.md).

## Examples

``` r
set.seed(1)
q <- data.frame(object = LETTERS[1:6],
                domain = rep(0:1, each = 3))
beta <- setNames(0.8 * q$domain, q$object)
pr <- t(combn(q$object, 2))
d <- data.frame(a = rep(pr[, 1], each = 20),
                b = rep(pr[, 2], each = 20))
p <- plogis(beta[d$a] - beta[d$b])
d$winner <- ifelse(runif(nrow(d)) < p, d$a, d$b)
fit <- btl_explanatory(d, q, ~ domain, "a", "b", winner = "winner")
fit$object_coefficients
#>    term estimate    se     t  df       p   p_adj
#>  domain    1.040 0.170 6.130 Inf < 0.001 < 0.001
explanatory_test(fit)
#>           model parameters free_parameters r_squared r_squared_adj
#>  Explanatory CJ          1               5     0.920         0.900
#>            r2_basis chisq df p_naive chisq_kent     p p_kent
#>  object calibration 3.715  4   0.446      3.730 0.444  0.444
```
