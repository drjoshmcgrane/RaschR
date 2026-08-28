# Draw a Wright map with WrightMap

Prepares person estimates and item locations from a fitted model and
passes them to
[`WrightMap::wrightMap()`](https://rdrr.io/pkg/WrightMap/man/wrightMap.html).
Person panels may be formed from variables retained in the fit or
supplied as a matrix. Item panels use the `item.groups` facility in
WrightMap 1.5.

## Usage

``` r
wright_map(
  fit,
  type = c("thresholds", "locations"),
  person_panels = NULL,
  item_panels = NULL,
  ...
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md),
  including explanatory, MFRM and EFRM fits. Comparative-judgement
  models do not estimate person locations and are not supported.

- type:

  Plot category `"thresholds"` or one `"locations"` estimate per item.

- person_panels:

  Optional person-panel specification. Supply the name of one or more
  variables retained in `fit$person`, a vector or factor of panel
  memberships with one value per person, or a numeric matrix or data
  frame whose columns contain the estimates for the panels. Missing
  estimates are permitted. For EFRM fits, `"groups"` uses the fitted
  person groups.

- item_panels:

  Optional vector or factor assigning each item row to a panel. A named
  vector is matched to item names; an unnamed vector is used in item
  order. A named list may instead map panel names to item names. For
  EFRM fits, use `"sets"`, `"groups"`, or `c("sets", "groups")` to
  arrange the calibrated response columns by item set, person group, or
  frame. This option requires WrightMap 1.5 or later.

- ...:

  Further arguments passed to
  [`WrightMap::wrightMap()`](https://rdrr.io/pkg/WrightMap/man/wrightMap.html),
  such as `person.side`, `item.side`, `main.title`, or graphical
  settings.

## Value

Invisibly, the threshold matrix returned by
[`WrightMap::wrightMap()`](https://rdrr.io/pkg/WrightMap/man/wrightMap.html)
when its `return.thresholds` argument is true; otherwise `NULL`.

## Details

For partial credit and rating scale models, `type = "thresholds"`
displays each item's estimated category thresholds, labelled `t1`, `t2`,
and so on. A wholly dichotomous scale omits the redundant `t1` labels;
in a mixed scale, dichotomous items retain `t1` to align them with the
polytomous items. `type = "locations"` displays one location per item.
In MFRM and EFRM fits, the rows are the calibrated item-by-facet or
item-by-frame response columns. For EFRM fits,
`person_panels = "groups"` and `item_panels = "sets"` use the fitted
frame design; both may be specified together. A single person panel has
no heading by default. When several person panels are requested, their
panel labels are printed above the distributions. The plot has no
overall title unless one is supplied through `main.title`.

## References

Torres Irribarra, D., and Freund, R. (2025). WrightMap: IRT item-person
map with ConQuest integration. R package version 1.5.

## See also

[`plot_wright`](https://drjoshmcgrane.github.io/rasch/reference/plot_wright.md)

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
colnames(X) <- paste0("I", 1:6)
fit <- rasch(X)
if (requireNamespace("WrightMap", quietly = TRUE)) {
  wright_map(fit)
}
```
