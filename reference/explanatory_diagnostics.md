# Diagnose fixed departures from an explanatory model

Fits each available item-location, polytomous threshold-structure or
comparative-judgement object departure separately from the active model.
Probabilities use Kent calibration and Holm adjustment over the complete
candidate family. A candidate with a withheld probability remains in
that family.

## Usage

``` r
explanatory_diagnostics(fit, p_adjust = "holm")
```

## Arguments

- fit:

  A fitted explanatory Rasch or comparative judgement model.

- p_adjust:

  Multiplicity adjustment over the candidate departures.

## Value

A data frame ordered by adjusted probability. For item fits, a `weak`
column marks items whose thresholds the calibration flags as weakly
identified; their probabilities are withheld, since the departure test
rests on the same sparse categories, and a note on the table records the
withholding.
