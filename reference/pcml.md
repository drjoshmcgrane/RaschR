# Estimate Rasch thresholds by pairwise conditional maximum likelihood

Estimates PCM or RSM thresholds by Newton–Raphson maximisation of the
pairwise conditional likelihood (Andrich and Luo 2003; Zwinderman 1995).

## Usage

``` r
pcml(X, model = c("PCM", "RSM"), anchors = NULL, maxit = 60, tol = 1e-08)
```

## Arguments

- X:

  Persons-by-items integer score matrix (categories from 0). Missing
  values are handled by pairwise deletion, so linked booklet designs and
  random missingness estimate without imputation; the item-pair graph
  must be connected (some person answering items in both of any two
  blocks), otherwise relative locations between blocks are unidentified
  and the fit stops with an error naming the blocks – unless `anchors`
  fix an item in every block, the disjoint-form equating case.

- model:

  `"PCM"` or `"RSM"`.

- anchors:

  Optional anchor table for equating: a data frame with columns `item`
  (name or column index), `k`, and `tau` (the fixed value). A numeric
  `k` fixes that single threshold (individual anchoring); `k = NA` fixes
  the item's mean location at `tau` while its thresholds remain free
  (average anchoring). The remaining parameters are estimated on the
  anchored scale and no recentring is applied. PCM only.

- maxit, tol:

  Newton-Raphson iteration cap and convergence tolerance.

## Value

A list containing the threshold table `thr`, covariance matrix
`cov_tau`, pairwise conditional log-likelihood, iteration count,
convergence flag, notes, and maximum scores `m`. In `thr`, `weak` marks
all thresholds of an item with fewer than eight responses in any
category, or a threshold adjacent to a category with fewer than three
responses. Standard errors for weak thresholds are reported as `NA`.

## Details

For the PCM, the adjacent-category log odds are
\$\$\log\\P(X\_{ni}=k)/P(X\_{ni}=k-1)\\=\theta_n-\delta\_{ik}.\$\$
Conditioning on the score for an item pair removes \\\theta_n\\. The PCM
estimates each \\\delta\_{ik}\\; the RSM imposes
\\\delta\_{ik}=\beta_i+\tau_k\\ through a design matrix.

## References

Andrich, D. and Luo, G. (2003). Conditional pairwise estimation in the
Rasch model for ordered response categories using principal components.
Journal of Applied Measurement, 4(3), 205–221.

Zwinderman, A. H. (1995). Pairwise parameter estimation in Rasch models.
Applied Psychological Measurement, 19(4), 369–375.

## Examples

``` r
set.seed(1)
d <- seq(-1.5, 1.5, length.out = 6)
X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
colnames(X) <- paste0("I", 1:6)
pcml(X)$thr
#>   id item k        tau        se anchored  weak
#> 1  1    1 1 -1.5174835 0.1297013    FALSE FALSE
#> 2  2    2 1 -0.7730904 0.1121427    FALSE FALSE
#> 3  3    3 1 -0.2093836 0.1023242    FALSE FALSE
#> 4  4    4 1  0.3667512 0.1047318    FALSE FALSE
#> 5  5    5 1  0.7889434 0.1105680    FALSE FALSE
#> 6  6    6 1  1.3442630 0.1221786    FALSE FALSE
# anchor two items at fixed values (equating)
anchors <- data.frame(item = c("I1", "I6"), k = 1,
                      tau = c(-1.5, 1.5))
pcml(X, anchors = anchors)$thr
#>   id item k        tau        se anchored  weak
#> 1  1    1 1 -1.5000000 0.0000000     TRUE FALSE
#> 2  2    2 1 -0.6880959 0.1613982    FALSE FALSE
#> 3  3    3 1 -0.1211423 0.1463167    FALSE FALSE
#> 4  4    4 1  0.4586830 0.1532659    FALSE FALSE
#> 5  5    5 1  0.8835376 0.1568160    FALSE FALSE
#> 6  6    6 1  1.5000000 0.0000000     TRUE FALSE
```
