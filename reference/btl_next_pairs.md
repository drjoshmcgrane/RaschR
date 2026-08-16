# Recommend the next informative comparisons (adaptive step)

Ranks candidate object pairs by the information expected from one
additional comparison at the current estimates (Pollitt 2012). By
default, priority is the one-step reduction in total location variance
from a rank-one covariance update. This favours close pairs and objects
measured with less precision.

## Usage

``` r
btl_next_pairs(fit, n = 10, weight_se = TRUE)
```

## Arguments

- fit:

  A paired-comparison fit from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md).

- n:

  Number of pairs to return.

- weight_se:

  If `TRUE` (the default), rank pairs by their one-step reduction in
  total location variance. When the fit has no covariance, the fallback
  priority is expected information multiplied by the sum of the two
  squared standard errors. If `FALSE`, rank pairs by expected
  information alone.

## Value

A data frame of the top `n` candidate pairs, each oriented to its
stronger object: `object_a`, `object_b`, the location `gap`,
`n_existing` (replications already observed for the pair),
`expected_information` (of one new comparison), and `priority`. Sorted
by `priority` (or by `expected_information` when `weight_se = FALSE`).

## Details

The procedure is a greedy, one-step ranking rather than a jointly
optimal design. Applied to a sandwich covariance, the update ranks pairs
but does not give an exact variance reduction. Adaptive selection can
also inflate a separation reliability calculated from the same
comparisons (Bramley 2015).

## References

Pollitt, A. (2012). The method of adaptive comparative judgement.
*Assessment in Education*, 19(3), 281-300. Bramley, T. (2015).
Investigating the reliability of Adaptive Comparative Judgment.
*Cambridge Assessment Research Report*.

## See also

[`btl_information`](https://drjoshmcgrane.github.io/rasch/reference/btl_information.md),
[`plot_btl_targeting`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_targeting.md)

## Examples

``` r
set.seed(1)
beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
pr <- t(combn(names(beta), 2))
d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
btl_next_pairs(btl(d, "a", "b", "win"), n = 5)
#>   object_a object_b       gap n_existing expected_information     priority
#> 1        B        A 0.8837971         30            0.2068987 0.0010120610
#> 2        D        C 0.6961084         30            0.2220026 0.0009330875
#> 3        D        B 1.4978719         30            0.1493481 0.0009279925
#> 4        C        A 1.6855607         30            0.1319119 0.0008269402
#> 5        C        B 0.8017636         30            0.2137663 0.0007881426
```
