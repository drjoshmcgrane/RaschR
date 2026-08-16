# DIF analysis for paired comparisons

Tests whether object locations differ across groups of judges. Several
judge factors can be fitted jointly, with optional factor-by-factor
interactions. Uniform DIF is a judge-factor effect; non-uniform DIF is
its interaction with opponent-strength band.

## Usage

``` r
btl_dif(
  fit,
  factors,
  objects = NULL,
  effects = c("main", "factorial"),
  p_adjust = "BH",
  alpha = 0.05,
  flag_logits = 0.5,
  min_n = 20,
  maxit = 60,
  tol = 1e-08
)
```

## Arguments

- fit:

  An object from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md).

- factors:

  One judge factor, or a named list containing several. Each factor may
  have one value per comparison row or be a vector named by judge.

- objects:

  Objects to test; all by default.

- effects:

  `"main"` (default) models several factors additively (each factor's
  main effect and its band interaction); `"factorial"` also crosses the
  factors with one another.

- p_adjust:

  Multiplicity adjustment across objects within each term; the
  resolved-size probabilities are adjusted in one pool over all objects,
  terms, and cell pairs.

- alpha:

  Significance level for adjusted probabilities.

- flag_logits:

  Absolute resolved difference flagged as practically significant.

- min_n:

  Term cells with fewer comparisons involving the object are dropped
  from its resolution, with a note.

- maxit, tol:

  Newton controls for the resolution refits.

## Value

A list of class `"rasch_btl_dif"`: `summary` (one row per object and
group term with the uniform F, adjusted p and partial eta-squared – the
term itself – the non-uniform ones – the term crossed with the opponent
band – plus `uniform_DIF`, `nonuniform_DIF` and `superseded` flags);
`terms` (the full per-object analysis-of-variance table); `levels`
(resolved location and SE per object, term and cell); `sizes` (per
object, term and cell pair: difference in logits, SE, t, degrees of
freedom, adjusted p, significance and practical flags); `effects`,
`factors`, and `notes`.

## Details

Judges are the independent units. For each object, oriented residuals
are aggregated to one weighted mean per judge and opponent band. A
split-plot analysis then tests judge factors between judges and band
effects within judges. Each factor level requires at least two judges.
Confirmatory Wald tests are available only when the base fit supplies a
valid judge-clustered covariance.

A significant uniform term is followed by a joint refit in which the
object has one location per factor cell. Differences between these
locations are reported in logits with clustered Wald tests. Higher-order
terms supersede their component terms. Models fitted with `order` retain
the exposure and carry-over effects in both the residual analysis and
refit.

Objects are resolved one at a time against the common locations of the
remaining objects. With DIF in several objects, this can induce
compensating apparent DIF in invariant objects (Andrich and Hagquist
2012, 2015).

## References

Andrich, D., & Hagquist, C. (2012). Real and artificial differential
item functioning. *Journal of Educational and Behavioral Statistics*,
37(3), 387-416.

Dittrich, R., Hatzinger, R., & Katzenbeisser, W. (1998). Modelling the
effect of subject-specific covariates in paired comparison studies with
an application to university rankings. *Journal of the Royal Statistical
Society C*, 47(4), 511-525.

## Examples

``` r
set.seed(1)
beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
pr <- t(combn(names(beta), 2))
d <- data.frame(a = rep(pr[, 1], each = 60), b = rep(pr[, 2], each = 60),
                judge = sample(sprintf("J%02d", 1:12), 360, TRUE))
shift <- ifelse(d$judge %in% sprintf("J%02d", 1:6) & d$a == "C", 0.9,
         ifelse(d$judge %in% sprintf("J%02d", 1:6) & d$b == "C", -0.9, 0))
p <- plogis(beta[d$a] - beta[d$b] + shift)
d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
f <- btl(d, "a", "b", winner = "win", judge = "judge")
grp <- setNames(rep(c("g1", "g2"), each = 6), sprintf("J%02d", 1:12))
btl_dif(f, grp, objects = "C")
#> DIF for paired comparisons: 1 factor(s) [group], main effects
#> Residual ANOVA per object and term (uniform = term; non-uniform = term x opponent band)
#>  object  term F_uniform p_uniform_adj uniform_DIF F_nonuniform p_nonuniform_adj
#>       C group    10.576         0.009           *        0.685            0.514
#>  nonuniform_DIF
#>                
#> 
#> Resolved locations (logits; BH over 1 comparison(s); practical 0.50)
#>  object  term level_a level_b difference    se     t df   p_adj significant
#>       C group      g1      g2      1.404 0.283 4.965 11 < 0.001           *
#>  practical
#>          *
```
