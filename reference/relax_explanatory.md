# Relax a nominated explanatory restriction

Adds either one fixed item-location departure or the part of an item's
threshold-structure block not already represented by the predictor
design, then repeats the complete conditional calibration and downstream
Rasch analysis. The departure is fixed rather than random; raw-score
sufficiency and the common discrimination remain. Earlier DIF splits and
superitem definitions are retained.

## Usage

``` r
relax_explanatory(fit, item, component = c("location", "thresholds"))
```

## Arguments

- fit:

  A fitted explanatory Rasch model.

- item:

  Item name.

- component:

  Either `"location"` or `"thresholds"`.

## Value

A partially relaxed `"rasch_explanatory"` fit.
