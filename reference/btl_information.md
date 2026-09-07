# Information and targeting of a paired-comparison design

Calculates the Fisher information supplied by the observed comparison
design. For location difference \\d=\beta_a-\beta_b\\, one dichotomous
comparison contributes \$\$I(d)=P(a\succ b)\\1-P(a\succ b)\\.\$\$ For an
ordered comparison, the contribution is the variance of the response
score. Information is summed over the comparisons involving each object,
including replication counts. Observed-design information retains the
fitted position and dependence effects; it is conditional on the
recorded comparison history.

## Usage

``` r
btl_information(fit)
```

## Arguments

- fit:

  A paired-comparison fit from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md).

## Value

A list of class `"rasch_btl_info"`: `objects` (per object: `location`,
the fit's `se`, `n_comparisons`, the design `information`, and
`se_naive`); `pairs` (per observed pair: `n`, the mean location `gap`,
and the pair's `information`); `comparisons` (per comparison: the signed
`gap`, `weight`, and the single-comparison `information`); the scalar
`total` information; `m`; the `clustered` flag; and `notes`.

## Details

`se_naive = 1/sqrt(information)` treats each object's comparisons in
isolation. It is a description of the design, not the fitted standard
error or a bound on it. The fitted standard error also reflects joint
estimation, the identifying constraint, and judge clustering. The fitted
model must have converged.

## References

Pollitt, A. (2012). The method of adaptive comparative judgement.
*Assessment in Education*, 19(3), 281-300.

## See also

[`plot_btl_targeting`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_targeting.md),
[`btl_next_pairs`](https://drjoshmcgrane.github.io/rasch/reference/btl_next_pairs.md)

## Examples

``` r
set.seed(1)
beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
pr <- t(combn(names(beta), 2))
d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
btl_information(btl(d, "a", "b", "win"))
#> Paired-comparison design information: 4 objects, total 30.04
#> One-comparison Fisher information about the location gap (dichotomous: P(1 - P))
#>  object location    se n_comparisons information se_naive
#>       A   -1.238 0.214            90      12.487    0.283
#>       B   -0.354 0.186            90      17.100    0.242
#>       C    0.448 0.180            90      17.030    0.242
#>       D    1.144 0.209            90      13.463    0.273
#> Note: se is the Godambe sandwich standard error; se_naive = 1/sqrt(information) is a design-only yardstick (the object's comparisons treated in isolation), not a bound -- the fitted se can sit below or above it
```
