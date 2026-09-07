# Compare the partial credit and rating scale models

Compares a fitted partial credit model with the rating scale
reparameterisation of the same data. Both raw and composite-likelihood
adjusted statistics are returned.

## Usage

``` r
lr_test(fit, maxit = 60, tol = 1e-08)
```

## Arguments

- fit:

  An unrestricted, unanchored `"PCM"` fit from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  with equal maximum scores across items (the rating parameterisation
  requires them).

- maxit, tol:

  Passed to the rating-scale refit.

## Value

A list of class `"rasch_lr"`: raw `chisq`, `df`, `p` (the conventional
display); adjusted `chisq_adj`, `p_adj`, and the eigenvalues `lambda`;
the two log-likelihoods; and the rating-scale refit (`fit_rsm`),
retaining keyed scoring, DIF-split records and superitem definitions for
subsequent analyses.

## Details

The pairwise conditional likelihood is a composite likelihood: each
response contributes to every item pair in which it appears.
Consequently, the raw statistic \\W=2(cl\_{PCM}-cl\_{RSM})\\ does not
have an ordinary chi-square reference distribution. Its limiting
distribution is \\\sum_j\lambda_j\chi^2_1\\ (Kent 1982; Varin, Reid and
Firth 2011), where the \\\lambda_j\\ are obtained from the sensitivity
matrix \\H\\, variability matrix \\J\\, and the constraints defining the
RSM. The mean-matched statistic is \$\$W\_{adj}=rW/\sum_j\lambda_j,\$\$
with \\r\\ degrees of freedom.

Use `p_adj` for inference. The unadjusted `p` is retained for
descriptive comparison with conventional displays. The adjustment is a
first-order approximation and can be mildly anti-conservative in small
samples with long polytomous tests. Interpret values near the nominal
level cautiously in such designs.

## References

Kent, J. T. (1982). Robust properties of likelihood ratio tests.
Biometrika, 69, 19-27. Varin, C., Reid, N. and Firth, D. (2011). An
overview of composite likelihood methods. Statistica Sinica, 21, 5-42.

## Examples

``` r
set.seed(1)
tau <- c(-0.7, 0.7)
X <- sapply(seq(-1, 1, length.out = 6), function(d) vapply(rnorm(300),
  function(b) sample(0:2, 1, prob = item_moments(b, tau + d)$P), 0L))
colnames(X) <- paste0("Q", 1:6)
lr_test(rasch(X, model = "PCM"))
#> Likelihood-ratio test: partial credit vs rating parameterisation
#>   Raw composite chi-square 8.741 on 5 df, p = 0.120 (conventional display; anticonservative)
#>   Adjusted chi-square 2.062 on 5 df, p = 0.840 (Kent 1982 first-order calibration)
#>   log-likelihood (pairwise composite): PCM -2814.271, RSM -2818.642
```
