# rasch: Models and Diagnostics for Rasch Measurement Theory <img src="man/figures/logo.png" align="right" height="139" alt="rasch hex logo" />

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/rasch)](https://CRAN.R-project.org/package=rasch)
[![R-CMD-check](https://github.com/drjoshmcgrane/rasch/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/drjoshmcgrane/rasch/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/drjoshmcgrane/rasch/actions/workflows/pkgdown.yaml/badge.svg)](https://drjoshmcgrane.github.io/rasch/)
<!-- badges: end -->

`rasch` fits and evaluates models within Rasch Measurement Theory. It includes
models for item responses, ratings, linked frames of reference, and paired
comparisons, with a common set of functions for examining fit, invariance,
targeting, dimensionality, and local dependence.

For an item with ordered scores \(x=0,\ldots,m_i\), the partial credit model
can be written as

$$
\log\frac{P(X_{ni}=x)}{P(X_{ni}=x-1)}=\theta_n-\delta_{ix},
$$

where \(\theta_n\) is the location of person \(n\) and \(\delta_{ix}\) is
the threshold between categories \(x-1\) and \(x\). Dichotomous items have
one threshold. The rating scale model imposes a common threshold structure
across items. The total score is sufficient for \(\theta_n\), allowing item
parameters to be estimated without assuming a population distribution for
the person locations.

The package treats fit to the model as an empirical question. Its diagnostics
examine whether comparisons remain invariant across persons, items, groups,
occasions, raters, and other parts of the measurement design.

## Models

| Function | Model |
|---|---|
| `rasch()` | Dichotomous Rasch, partial credit, and rating scale models |
| `rasch_mfrm()` | Many-facet Rasch model |
| `rasch_efrm()` | Extended frame of reference model |
| `btl()` | Comparative judgement models for dichotomous and polytomous paired comparisons |
| `btl_efrm()` | Extended frame of reference model for paired comparisons |

The item-response models use pairwise conditional maximum likelihood.
Person locations are estimated by Warm's weighted likelihood method. Anchored
estimation is available for equating, and incomplete linked designs can be
fitted when their observed response structure identifies a common scale.

The suite follows Rasch (1960) and Andrich and Marais (2019); the frame of
reference models follow Humphry's (2005) thesis, where the model is
introduced and named, and Humphry and Andrich (2008); the dichotomous
comparative judgement model is the conditional form of the dichotomous
Rasch model (Andrich, 1978), and the polytomous comparative judgement
model is its adjacent-categories extension (Tutz, 1986).

## Shiny application

The package includes a graphical interface for analysts who do not normally
work in R. Launch it after installation with:

```r
rasch::run_app()
```

The application imports data, assigns variables to their measurement roles,
fits the selected model, and displays the resulting tables and plots. Results
can be downloaded, and the R call for each analysis is shown in the interface.

<p align="center">
  <img src="man/figures/app-items.png" alt="Item statistics and an item characteristic curve in the rasch Shiny application" width="90%" />
</p>

## Installation

Install the CRAN release with:

```r
install.packages("rasch")
```

The development version is available from GitHub:

```r
# install.packages("remotes")
remotes::install_github("drjoshmcgrane/rasch")
```

## Example

```r
library(rasch)

d <- simulate_rasch(
  n_persons = 500,
  n_items = 10,
  n_groups = 2,
  seed = 1
)

fit <- rasch(d, model = "PCM", id = "id", factors = "group")

summary(fit)
fit_summary_table(fit)
targeting_table(fit)
dif_anova(fit)
residual_correlations(fit)
dimensionality_test(fit)

plot_pimap(fit)
plot_icc(fit, "I05", group = "group")
```

The [function reference](https://drjoshmcgrane.github.io/rasch/reference/index.html)
documents the data requirements and returned values for each analysis. The
[articles](https://drjoshmcgrane.github.io/rasch/articles/) give worked
examples for the main model families, DIF with repeated measures, and
simulation.

## References

Andrich, D. (1978). Relationships between the Thurstone and Rasch approaches
to item scaling. *Applied Psychological Measurement*, 2(3), 451--462.

Andrich, D., and Luo, G. (2003). Conditional pairwise estimation in the Rasch
model for ordered response categories using principal components. *Journal
of Applied Measurement*, 4(3), 205--221.

Andrich, D., and Marais, I. (2019). *A Course in Rasch Measurement Theory:
Measuring in the Educational, Social and Health Sciences*. Springer.

Bradley, R. A., and Terry, M. E. (1952). Rank analysis of incomplete block
designs: I. The method of paired comparisons. *Biometrika*, 39, 324--345.

Humphry, S. M. (2005). *Maintaining a Common Arbitrary Unit in Social
Measurement*. PhD thesis, Murdoch University.

Humphry, S. M., and Andrich, D. (2008). Understanding the unit in the Rasch
model. *Journal of Applied Measurement*, 9(3), 249--264.

Humphry, S. M. (2010). Modeling the effects of person group factors on
discrimination. *Educational and Psychological Measurement*, 70(2), 215--231.

Humphry, S. M. (2012). Item set discrimination and the unit in the Rasch
model. *Journal of Applied Measurement*, 13(2), 165--180.

Linacre, J. M. (1989). *Many-Facet Rasch Measurement*. MESA Press.

Montuoro, P., and Humphry, S. M. (2024). Modeling the effect of reading item
clarity on item discrimination. *Journal of Applied Measurement*, 24(3/4),
121--132.

Rasch, G. (1960). *Probabilistic Models for Some Intelligence and Attainment
Tests*. Danish Institute for Educational Research. Expanded edition,
University of Chicago Press, 1980.

Rasch, G. (1977). On specific objectivity: An attempt at formalizing the
request for generality and validity of scientific statements. *Danish
Yearbook of Philosophy*, 14, 58--94.

Tutz, G. (1986). Bradley-Terry-Luce models with an ordered response. *Journal
of Mathematical Psychology*, 30(3), 306--316.

Warm, T. A. (1989). Weighted likelihood estimation of ability in item
response theory. *Psychometrika*, 54(3), 427--450.

Zwinderman, A. H. (1995). Pairwise parameter estimation in Rasch models.
*Applied Psychological Measurement*, 19(4), 369--375.

## Citation

Use `citation("rasch")` to obtain the citation for the installed version.
