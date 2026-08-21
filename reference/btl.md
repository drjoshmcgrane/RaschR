# Fit comparative judgement models to paired comparisons

Fits the Bradley–Terry–Luce model to dichotomous comparisons or its
ordered-response extension (Tutz 1986). Object, judge, and pair fit are
reported with an object separation index and design diagnostics.

## Usage

``` r
btl(
  data,
  object_a,
  object_b,
  winner = NULL,
  response = NULL,
  margin = NULL,
  judge = NULL,
  count = NULL,
  order = NULL,
  position = FALSE,
  anchors = NULL,
  ties = c("drop", "half", "error"),
  thresholds = c("free", "pc"),
  maxit = 60,
  tol = 1e-08
)
```

## Arguments

- data:

  A data frame with one comparison per row.

- object_a, object_b:

  Names of the columns holding the two objects compared.

- winner:

  Name of the column holding the winner of each row: its value must
  equal one of the two objects. `"tie"` and `"draw"` mark ties. Ignored
  when `response` is supplied.

- response:

  Optional ordered response favouring `object_a` over `object_b`: an
  ordered factor from least to greatest preference for `object_a`, or
  integer scores `0..m`.

- margin:

  Optional ordered margin-of-victory column, combined with `winner` to
  construct an orientation-invariant response.

- judge:

  Optional name of a judge column; enables the judge fit table and
  clusters the sandwich standard errors by judge.

- count:

  Optional name of a column of replication counts (a row standing for
  several identical comparisons).

- order:

  Optional column giving each judge's comparison sequence; requires
  `judge`. See Details. Incompatible with `ties = "half"`.

- position:

  If `TRUE`, estimate a first-presentation advantage, treating
  `object_a` as the first object in each comparison.

- anchors:

  Optional named numeric vector of fixed object locations. Anchored
  objects have standard error zero and must not be boundary objects.

- ties:

  How to treat ties in the dichotomous analysis: `"drop"` (default,
  removed with a note), `"half"` (half a win each way, a common
  pragmatic device – flagged in the notes because the halves are not
  independent Bernoulli trials), or `"error"`. With polytomous
  responses, code ties as a middle category instead.

- thresholds:

  `"free"` (default) estimates every symmetric threshold; `"pc"` retains
  only the symmetric spread component.

- maxit, tol:

  Newton-Raphson iteration cap and convergence tolerance.

## Value

A `"rasch_btl"` object. Principal components are `objects`, `pairs`,
`judges`, the total pair-fit test, `osi`, `loglik`, composite-likelihood
information `cl`, convergence details, and `notes`. Ordered-response
fits also contain `thresholds`, `m`, and `categories`. Fits using
`order` contain `dependence` and `dependence_data`.

## Details

For objects \\a\\ and \\b\\, the dichotomous model is \$\$P(a\succ
b)=\frac{\exp(\beta_a)} {\exp(\beta_a)+\exp(\beta_b)}.\$\$ This is the
conditional form of the dichotomous Rasch model (Andrich 1978). For an
ordered response \\Y=0,\ldots,m\\,
\$\$\log\\P(Y=r)/P(Y=r-1)\\=\beta_a-\beta_b-\tau_r,\$\$ with thresholds
constrained to be symmetric under reversal of presentation order. Two
categories reproduce the dichotomous model.

Locations are identified by a sum-zero constraint unless anchors are
supplied. The comparison graph must be connected, and the directed win
graph must be strongly connected for all free locations to be finite
(Ford 1957). Boundary objects are removed when this leaves an identified
model; otherwise fitting stops.

Standard errors use the Godambe sandwich covariance. When `judge` is
supplied, the covariance is clustered by judge. Clustered inference is
withheld when there are fewer than ten judges, fewer than eight
effective judges, or no residual cluster degrees of freedom. A caution
is attached when the effective count is below 9.5 or one judge supplies
more than 20 per cent of the comparisons.

Dichotomous data may be supplied as a winner, with ties dropped or
divided equally between the two outcomes. Ordered data may instead be
supplied directly as scores from 0 to \\m\\, or assembled from the
winner and an ordered margin of victory. Plain factors are refused
because alphabetical ordering can reverse the response scale. The `"pc"`
threshold option retains the symmetric spread component, which can
stabilise thin categories.

If comparison order is supplied, exposure and carry-over effects are
estimated from each judge's preceding comparisons. The `position` term
estimates a first-presentation effect. These coefficients enter the
model jointly with the object locations and are reported in logits. The
carry-over estimate and clustered SE remain descriptive below 30 judges;
its probability is withheld because null calibration is mildly
anti-conservative at smaller judge counts. Anchors fix nominated object
locations and replace the sum-zero origin.

## References

Bradley, R. A. and Terry, M. E. (1952). Rank analysis of incomplete
block designs: I. The method of paired comparisons. Biometrika, 39,
324–345.

Luce, R. D. (1959). Individual Choice Behavior. Wiley.

Andrich, D. (1978). Relationships between the Thurstone and Rasch
approaches to item scaling. Applied Psychological Measurement, 2,
451–462.

Tutz, G. (1986). Bradley-Terry-Luce models with an ordered response.
Journal of Mathematical Psychology, 30(3), 306–316.

Agresti, A. (1992). Analysis of ordinal paired comparison data. Journal
of the Royal Statistical Society C, 41(2), 287–297.

Davidson, R. R. (1970). On extending the Bradley-Terry model to
accommodate ties in paired comparison experiments. Journal of the
American Statistical Association, 65(329), 317–328.

Ford, L. R. (1957). Solution of a ranking problem from binary
comparisons. American Mathematical Monthly, 64(8), 28–33.

Davidson, R. R. and Beaver, R. J. (1977). On extending the Bradley-Terry
model to incorporate within-pair order effects. Biometrics, 33(4),
693–702.

## See also

[`btl_dif`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md),
[`btl_efrm`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md),
[`btl_information`](https://drjoshmcgrane.github.io/rasch/reference/btl_information.md),
[`btl_transitivity`](https://drjoshmcgrane.github.io/rasch/reference/btl_transitivity.md),
and
[`simulate_btl`](https://drjoshmcgrane.github.io/rasch/reference/simulate_btl.md).

## Examples

``` r
set.seed(1)
beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
pairs <- t(combn(names(beta), 2))
d <- data.frame(a = rep(pairs[, 1], each = 30),
                b = rep(pairs[, 2], each = 30))
p <- plogis(beta[d$a] - beta[d$b])
d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
btl(d, object_a = "a", object_b = "b", winner = "win")
#> Bradley-Terry-Luce analysis: 4 objects, 180 comparisons
#> Conditional ML: converged in 6 iterations; sandwich SEs
#> Object separation index 0.963; pairwise chi-square 1.07 on 3 df, p = 0.783
#>  object location    se comparisons wins fit_resid
#>       A   -1.238 0.214          90   16    -0.111
#>       B   -0.354 0.186          90   36     0.526
#>       C    0.448 0.180          90   56    -0.056
#>       D    1.144 0.209          90   72     0.037
```
