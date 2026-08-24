# Compare an explanatory model with its free calibration

Tests explanatory item, threshold or object restrictions against the
corresponding free calibration of the same responses. The inferential
result uses the first-order Kent calibration for the fitted likelihood
and sandwich covariance. The calibration coefficient of determination is
\$\$R^2\_{cal}=1-\frac{\sum_j(\hat\eta^{free}\_j-
\hat\eta^{expl}\_j-\bar d)^2}{\sum_j(\hat\eta^{free}\_j-
\bar\eta^{free})^2},\$\$ where \\\bar d\\ removes the arbitrary scale
origin. It describes the proportion of variation in the well-determined
free threshold calibration (Rasch models) or free object calibration
(comparative judgement) reproduced by the explanatory model. It is at
most one and may be negative. It is not adjusted for the number of
predictors, so with few calibrated parameters it reads above zero even
for an uninformative design. `r_squared_adj` divides the unexplained
proportion by its share of the degrees of freedom,
\\1-(1-R^2\_{cal})(n-1)/\mathit{df}\\, where \\n\\ counts the calibrated
parameters compared and \\\mathit{df}\\ is \\n\\ minus the rank of the
retained explanatory design with its origin, so exclusions that remove a
level's only support reduce it. The correction is exact for independent
homoskedastic estimates fitted by least squares, which these
calibrations are not, so read it as a descriptive optimism adjustment.
Read either beside the test rather than in place of it.

## Usage

``` r
explanatory_test(fit)
```

## Arguments

- fit:

  A fitted explanatory Rasch or comparative judgement model.

## Value

A one-row data frame containing the raw and Kent-calibrated statistics,
degrees of freedom and parameter counts. The primary `p` and the
retained `p_kent` are the Kent-calibrated probability. `p_naive` is the
unscaled composite-likelihood probability and is provided for
methodological inspection, not inference. `r_squared` is the calibration
coefficient of determination, `r_squared_adj` its
degrees-of-freedom-adjusted counterpart, and `r2_basis` names the
calibrated parameters used.
