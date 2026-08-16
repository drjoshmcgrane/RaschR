# Transitivity of paired comparisons

Summarises circular triads in the observed paired comparisons. A triad
is circular when A is preferred to B, B to C, and C to A. For a complete
tournament, the function reports Kendall's coefficient of consistency
(Kendall and Babington Smith 1940). Judge-specific summaries are
returned when judges are available.

## Usage

``` r
btl_transitivity(fit, min_triples = 5L)
```

## Arguments

- fit:

  A paired-comparison fit from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md).

- min_triples:

  A judge is reported only if this many complete triples (all three
  pairs judged) are available.

## Value

A list of class `"rasch_btl_transitivity"`: `summary` (one row: objects,
pairs compared, complete triples, circular triads, the circular rate,
the chance rate 0.25, the consistency index `1 - rate/0.25`, and
Kendall's `zeta` when the design is a complete round-robin with no
exactly-tied pair – `NA` otherwise); `objects` (each object's
circular-triad involvement); `judges` (per-judge consistency, when
judges exist); and `notes`.

## Details

A circular-triad rate of one quarter is the benchmark for a random
tournament. It is not the expected rate under a fitted BTL model with
unequal object locations, so this function is a descriptive consistency
measure rather than a calibrated goodness-of-fit test.

## References

Kendall, M. G., & Babington Smith, B. (1940). On the method of paired
comparisons. *Biometrika*, 31(3/4), 324-345.

## Examples

``` r
set.seed(1); objs <- LETTERS[1:6]; beta <- setNames(seq(-1.5, 1.5, len = 6), objs)
pr <- t(utils::combn(objs, 2))
d <- data.frame(a = rep(pr[, 1], each = 20), b = rep(pr[, 2], each = 20))
d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
btl_transitivity(btl(d, "a", "b", "win"))
#> Paired-comparison transitivity: 6 objects, 20 complete triples
#> Circular triads: 0 (0.0% of triples; random-tournament benchmark 25%) -> consistency 1.00
#> Kendall coefficient of consistency (complete design): 1.000
#> Note: the 0.25 chance rate is a random-tournament benchmark, not the fitted BTL expected circular rate; transitivity is descriptive 
```
