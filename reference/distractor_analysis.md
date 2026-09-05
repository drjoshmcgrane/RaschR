# Distractor analysis for multiple-choice items

For every keyed item and response option: the count and proportion
choosing it among respondents with a non-extreme rest measure, their
mean location, and the point-biserial correlation between choosing the
option and the person measure. These summaries use the rest measure (the
person estimate from the other items), so the analysed item cannot
credit its own takers. The keyed option should attract the ablest
persons and carry the only positive point-biserial; a distractor whose
takers are abler than the keyed option's (with at least `min_n` takers)
is flagged as a possible miskey.

## Usage

``` r
distractor_analysis(fit, items = NULL, min_n = 10)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  run with a `key`.

- items:

  Optional subset of item names; defaults to every keyed item.

- min_n:

  Minimum takers for an option to be eligible for the miskey flag.

## Value

A data frame with one row per item-option: `item`, `option`, its
assigned `score`, `keyed` (full credit), `n`, `prop`, `mean_location`,
`point_biserial`, and `flag`.

## Examples

``` r
set.seed(1); Np <- 400
th <- rnorm(Np)
raw <- sapply(seq(-1, 1, length.out = 6), function(d) {
  ok <- rbinom(Np, 1, plogis(th - d))
  ifelse(ok == 1, "A", sample(c("B", "C", "D"), Np, replace = TRUE))
})
colnames(raw) <- paste0("M", 1:6)
fit <- rasch(raw, key = setNames(rep("A", 6), colnames(raw)))
head(distractor_analysis(fit))
#>   item option score keyed   n       prop mean_location point_biserial  flag
#> 1   M1      A     1  TRUE 243 0.71260997     0.2559400      0.2368696 FALSE
#> 2   M1      B     0 FALSE  33 0.09677419    -0.1202694     -0.1070635 FALSE
#> 3   M1      C     0 FALSE  28 0.08211144    -0.2374100     -0.1422987 FALSE
#> 4   M1      D     0 FALSE  37 0.10850440    -0.1274202     -0.1172765 FALSE
#> 5   M2      A     1  TRUE 218 0.62643678     0.2520763      0.2455266 FALSE
#> 6   M2      B     0 FALSE  49 0.14080460    -0.1177247     -0.1061147 FALSE
```
