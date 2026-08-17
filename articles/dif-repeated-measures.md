# DIF with several person factors and repeated measures

``` r

library(rasch)
```

## Design the analysis around the person

Differential item functioning is a violation of the Rasch model’s
requirement of invariant comparison: item locations should not depend on
which persons respond (Rasch 1961). When several person factors are
relevant, they should enter one model. Repeated observations require a
further distinction: group is a between-person factor, whereas occasion
varies within person.

For a factor \\G\\, class interval \\C\\, and standardised residual
\\z\\, the single-factor model is

\\ z=\mu+G+C+G\mathbin{:}C+\varepsilon. \\

The \\G\\ term tests uniform DIF. The \\G\mathbin{:}C\\ term tests
non-uniform DIF. With several factors, `dif_anova` fits their terms
jointly and uses Type II sums of squares.

The following dataset has two observations per person. Item I03 has a
group shift and item I06 has an occasion shift.

``` r

set.seed(21)
N <- 240
difficulty <- seq(-1.5, 1.5, length.out = 8)
theta <- rnorm(N)
group <- rep(c("A", "B"), each = N / 2)

make_wave <- function(occasion_shift) {
  shift <- matrix(0, N, 8)
  shift[group == "B", 3] <- 0.8
  shift[, 6] <- occasion_shift
  matrix(rbinom(N * 8, 1,
                plogis(outer(theta, difficulty, "-") - shift)), N, 8)
}

X <- rbind(make_wave(0), make_wave(0.9))
colnames(X) <- sprintf("I%02d", 1:8)
dat <- data.frame(
  pid = rep(sprintf("P%03d", seq_len(N)), 2),
  X,
  group = rep(group, 2),
  occasion = rep(c("T1", "T2"), each = N)
)
```

## Fit once and test both factors

The repeated person identifier is carried into the fit. `dif_anova`
detects occasion as within-person; it can also be declared explicitly.
Persons, not stacked rows, are the units of analysis. Between-person
tests use Type II sums of squares. Within-person tests use person-level
contrasts, with a Greenhouse–Geisser correction when a factor has more
than two levels. This mixed-design analysis extends the single-factor
residual analysis of variance described by Andrich and Marais (2019).
Its F references are large-sample approximations.

``` r

fit <- rasch(dat, id = "pid", factors = c("group", "occasion"),
             items = sprintf("I%02d", 1:8))
da <- dif_anova(fit, within = "occasion")
da$summary
#>    item     term F_uniform p_uniform p_uniform_adj eta2_uniform uniform_DIF
#> 1   I01    group  0.752609 3.866e-01     0.8302540    3.247e-03       FALSE
#> 2   I01 occasion  0.002245 9.623e-01     0.9622515    9.718e-06       FALSE
#> 3   I02    group  0.098668 7.537e-01     0.8764326    4.269e-04       FALSE
#> 4   I02 occasion  0.748150 3.880e-01     0.6207361    3.228e-03       FALSE
#> 5   I03    group 15.829349 9.281e-05     0.0007425    6.413e-02        TRUE
#> 6   I03 occasion  0.520999 4.711e-01     0.6281947    2.250e-03       FALSE
#> 7   I04    group  0.059232 8.079e-01     0.8764326    2.563e-04       FALSE
#> 8   I04 occasion  0.354111 5.524e-01     0.6312890    1.531e-03       FALSE
#> 9   I05    group  3.973740 4.739e-02     0.1895652    1.691e-02       FALSE
#> 10  I05 occasion  6.852965 9.434e-03     0.0377370    2.881e-02        TRUE
#> 11  I06    group  0.359684 5.493e-01     0.8764326    1.555e-03       FALSE
#> 12  I06 occasion 19.380615 1.637e-05     0.0001310    7.740e-02        TRUE
#> 13  I07    group  0.024232 8.764e-01     0.8764326    1.049e-04       FALSE
#> 14  I07 occasion  1.524786 2.182e-01     0.4363025    6.558e-03       FALSE
#> 15  I08    group  0.666467 4.151e-01     0.8302540    2.877e-03       FALSE
#> 16  I08 occasion  3.943691 4.823e-02     0.1286127    1.679e-02       FALSE
#>    F_nonuniform p_nonuniform p_nonuniform_adj eta2_nonuniform nonuniform_DIF
#> 1        1.7981       0.1483           0.6567        0.022819          FALSE
#> 2        0.8939       0.4450           0.7722        0.011476          FALSE
#> 3        1.1534       0.3284           0.6567        0.014758          FALSE
#> 4        0.6573       0.5791           0.7722        0.008464          FALSE
#> 5        0.5694       0.6357           0.8476        0.007340          FALSE
#> 6        0.2067       0.8917           0.8917        0.002677          FALSE
#> 7        0.6183       0.6038           0.8476        0.007966          FALSE
#> 8        1.4413       0.2315           0.7722        0.018375          FALSE
#> 9        0.1963       0.8989           0.9402        0.002542          FALSE
#> 10       0.2284       0.8766           0.8917        0.002957          FALSE
#> 11       0.1332       0.9402           0.9402        0.001726          FALSE
#> 12       2.0156       0.1125           0.7722        0.025509          FALSE
#> 13       1.3740       0.2514           0.6567        0.017532          FALSE
#> 14       0.7690       0.5124           0.7722        0.009888          FALSE
#> 15       1.3748       0.2512           0.6567        0.017541          FALSE
#> 16       0.8186       0.4847           0.7722        0.010520          FALSE
#>    superseded
#> 1       FALSE
#> 2       FALSE
#> 3       FALSE
#> 4       FALSE
#> 5       FALSE
#> 6       FALSE
#> 7       FALSE
#> 8       FALSE
#> 9       FALSE
#> 10      FALSE
#> 11      FALSE
#> 12      FALSE
#> 13      FALSE
#> 14      FALSE
#> 15      FALSE
#> 16      FALSE
```

The planted effects are I03 by group and I06 by occasion, and both are
flagged. This run also flags I05 by occasion, where nothing was planted:
with eight items tested against every term, an occasional spurious flag
is expected at any fixed level, which is why adjusted probabilities and
effect sizes should be read together before an item is acted on. The
`dif_size` step below shows the planted effects are large where the
spurious flag is not.

The multiplicity adjustment is applied across items separately within
each term. Uniform DIF is a factor effect that is stable over the trait;
a factor-by-class-interval effect is non-uniform DIF. If person-factor
interactions are substantively required, use `effects = "factorial"` and
interpret a significant higher-order term before its component main
effects.

## Quantify the departure

ANOVA identifies evidence against invariance; it does not state the size
of the departure in logits. `dif_size` resolves an item by a
between-person factor. `dif_contrasts` is the corresponding
planned-contrast approach and uses person-level differencing for
within-person questions.

``` r

dif_size(fit, "I03", by = "group")
#> DIF size for I03 by group (resolved locations, logits)
#>  level location se weak   n
#>      A   -0.549 NA    0 240
#>      B    0.295 NA    0 240
#>  level_a level_b difference se  z  p p_adj lower upper significant practical
#>        A       B     -0.844 NA NA NA    NA    NA    NA          NA   >= 0.50
#> p adjusted by holm over 1 pairwise comparison(s); practical criterion 0.50 logits
#> notes: person identifiers repeat across response rows: resolved point differences remain descriptive, but sampling SEs, confidence intervals and Wald tests are withheld; use dif_contrasts for person-level inference or a whole-person bootstrap
dc <- dif_contrasts(fit, items = c("I03", "I06"), within = "occasion",
                    id = fit$person$id)
dc$table
#>   item                         contrast within estimate se statistic    df
#> 1  I03                     group: B - A  FALSE  0.84130 NA    4.0899 238.0
#> 2  I03                occasion: T2 - T1   TRUE -0.16638 NA   -0.7172 239.0
#> 3  I03 group(B - A) x occasion(T2 - T1)   TRUE -0.07149 NA   -0.2301 233.6
#> 4  I06                     group: B - A  FALSE  0.04635 NA    0.4839 237.6
#> 5  I06                occasion: T2 - T1   TRUE  1.04276 NA    4.6501 239.0
#> 6  I06 group(B - A) x occasion(T2 - T1)   TRUE -0.24577 NA   -1.1188 226.5
#>           p     p_adj lower upper significant practical
#> 1 5.907e-05 2.953e-04    NA    NA        TRUE      TRUE
#> 2 4.739e-01 1.000e+00    NA    NA       FALSE     FALSE
#> 3 8.182e-01 1.000e+00    NA    NA       FALSE     FALSE
#> 4 6.289e-01 1.000e+00    NA    NA       FALSE     FALSE
#> 5 5.491e-06 3.295e-05    NA    NA        TRUE      TRUE
#> 6 2.644e-01 1.000e+00    NA    NA       FALSE     FALSE
```

For repeated-person contrasts, significance comes from person-level
residual contrast scores. The resolved logit difference remains the
magnitude, but its row-independent calibration covariance is not a
repeated-measures standard error; the package therefore withholds the
logit SE and interval in this case.

A statistical flag should be considered with the logit magnitude,
targeting, item content, and the intended use of the scale. Resolving an
item changes the measurement model and should follow a substantive
account of why the item is not invariant.

## References

Andrich, D., and Marais, I. (2019). *A Course in Rasch Measurement
Theory: Measuring in the Educational, Social and Health Sciences*.
Springer.

Rasch, G. (1961). On general laws and the meaning of measurement in
psychology. In *Proceedings of the Fourth Berkeley Symposium on
Mathematical Statistics and Probability* (Vol. 4, pp. 321–333).
Berkeley: University of California Press.
