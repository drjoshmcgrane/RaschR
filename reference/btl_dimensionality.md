# Experimental residual dimensionality of paired comparisons

Decomposes the skew-symmetric matrix of observed-minus-expected pair
log-odds into Gower's (1977) rotational planes, or bimensions. A large
leading bimension indicates a structured cycle in the residual
comparisons. Its strength is compared with simulations from the fitted
one-dimensional model using the observed comparison counts.

## Usage

``` r
btl_dimensionality(fit, reps = 200L)
```

## Arguments

- fit:

  A paired-comparison fit from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md).

- reps:

  Model-simulated replicates for the noise reference.

## Value

A list of class `"rasch_btl_dim"`: `bimensions` (per bimension: strength
and share of residual size; the reference mean, 95th percentile, and the
clears-the-reference flag are reported for the leading bimension and
`NA` for the rest); `coords` (each object's position in the leading
bimension plane, for the residual map); `leading_structured` (whether
bimension 1 clears its reference); `residual_matrix`; and `notes`.

## Details

This is an experimental diagnostic. The reference is conditional on the
fitted point estimates because the model is not re-estimated in each
replicate. Ordered-response fits use the same points-proportion residual
in the data and simulations. Fits with exposure or carry-over effects
simulate those effects through each judge's observed sequence. The
fitted model must have converged.

A categorical result is withheld if any object pair is unobserved, or if
an ordered analysis contains count-weighted rows whose within-row
sequence is unavailable. It is also withheld when every judge receives
essentially the same comparison sequence and an order effect is fitted,
because order and residual structure are then confounded.

## References

Gower, J. C. (1977). The analysis of asymmetry and orthogonality. In J.
R. Barra et al. (Eds.), *Recent Developments in Statistics* (pp.
109-123). North-Holland.

## Examples

``` r
set.seed(1); objs <- LETTERS[1:6]; beta <- setNames(seq(-1.5, 1.5, len = 6), objs)
pr <- t(utils::combn(objs, 2))
d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
btl_dimensionality(btl(d, "a", "b", "win"), reps = 20)
#> Paired-comparison residual dimensionality: 3 bimension(s)
#> Leading bimension strength 1.370 (91% of residual; reference 95%: 2.440) -> within the conditional reference
```
