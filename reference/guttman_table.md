# Guttman-ordered response matrix and reproducibility

Orders persons by descending location and dichotomous items by ascending
location (the Guttman scalogram arrangement), computing the
reproducibility coefficient against the deterministic pattern implied by
each person's total score. PCM category steps can interleave across
items, so polytomous fits are rejected rather than forced into an
invalid whole-item deterministic order.

## Usage

``` r
guttman_table(fit)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  whose columns form one administered item set. Expanded EFRM and MFRM
  response-cell matrices are not accepted; one-cell-per-item reductions
  are.

## Value

A list with the ordered score matrix `matrix` (persons by items, row and
column names carrying the ID and item labels), the person and item
orderings, and the coefficient of reproducibility `CR`.

## References

Guttman, L. (1944). A basis for scaling qualitative data. American
Sociological Review, 9(2), 139–150.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(200 * 6, 1, plogis(outer(rnorm(200), d, "-"))), 200, 6)
colnames(X) <- paste0("I", 1:6)
guttman_table(rasch(X))$CR
#> [1] 0.8116667
```
