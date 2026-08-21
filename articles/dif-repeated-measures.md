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
da <- dif_anova(fit, within = "occasion", sizes = TRUE)
da$summary
#>  item     term F_uniform p_uniform p_uniform_adj eta2_uniform uniform_DIF
#>   I01    group     0.753     0.387         0.840        0.003            
#>   I01 occasion     0.002     0.962         0.962        0.000            
#>   I02    group     0.099     0.754         0.959        0.000            
#>   I02 occasion     0.748     0.388         0.840        0.003            
#>   I03    group    15.829   < 0.001         0.001        0.064           *
#>   I03 occasion     0.521     0.471         0.840        0.002            
#>   I04    group     0.059     0.808         0.959        0.000            
#>   I04 occasion     0.354     0.552         0.840        0.002            
#>   I05    group     3.974     0.047         0.309        0.017            
#>   I05 occasion     6.853     0.009         0.101        0.029            
#>   I06    group     0.360     0.549         0.840        0.002            
#>   I06 occasion    19.381   < 0.001       < 0.001        0.077           *
#>   I07    group     0.024     0.876         0.959        0.000            
#>   I07 occasion     1.525     0.218         0.731        0.007            
#>   I08    group     0.666     0.415         0.840        0.003            
#>   I08 occasion     3.944     0.048         0.309        0.017            
#>  F_nonuniform p_nonuniform p_nonuniform_adj eta2_nonuniform nonuniform_DIF
#>         1.798        0.148            0.678           0.023               
#>         0.894        0.445            0.840           0.011               
#>         1.153        0.328            0.840           0.015               
#>         0.657        0.579            0.840           0.008               
#>         0.569        0.636            0.848           0.007               
#>         0.207        0.892            0.959           0.003               
#>         0.618        0.604            0.840           0.008               
#>         1.441        0.231            0.731           0.018               
#>         0.196        0.899            0.959           0.003               
#>         0.228        0.877            0.959           0.003               
#>         0.133        0.940            0.962           0.002               
#>         2.016        0.112            0.600           0.026               
#>         1.374        0.251            0.731           0.018               
#>         0.769        0.512            0.840           0.010               
#>         1.375        0.251            0.731           0.018               
#>         0.819        0.485            0.840           0.011               
#>  superseded
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#> 
```

The multiplicity adjustment is applied across items separately within
each term. Uniform DIF is a factor effect that is stable over the trait;
a factor-by-class-interval effect is non-uniform DIF. If person-factor
interactions are substantively required, use `effects = "factorial"` and
interpret a significant higher-order term before its component main
effects. Read adjusted probabilities with effect sizes before changing
an item.

## Quantify the departure

ANOVA identifies evidence against invariance; it does not state the size
of the departure in logits. `dif_size` resolves an item by a
between-person factor. `dif_contrasts` provides planned contrasts and
uses person-level differencing for within-person questions.

``` r

dif_size(fit, "I03", by = "group")
#> DIF size for I03 by group (resolved locations, logits)
#>  level location se weak   n
#>      A   -0.549       0 240
#>      B    0.295       0 240
#>  level_a level_b difference se z p p_adj lower upper significant practical  ets
#>        A       B     -0.844                                        >= 0.50 <NA>
#>  signed_area
#>             
#> p adjusted by holm over 1 pairwise comparison(s); practical criterion 0.50 logits
#> notes: person identifiers repeat across response rows: resolved point differences remain descriptive, but sampling SEs, confidence intervals and Wald tests are withheld; use dif_contrasts for person-level inference or a whole-person bootstrap
dc <- dif_contrasts(fit, items = c("I03", "I06"), within = "occasion")
dc$table
#>  item                         contrast within estimate se statistic      df
#>   I03                     group: B - A           0.841        4.090 237.998
#>   I03                occasion: T2 - T1      *   -0.166       -0.717 239.000
#>   I03 group(B - A) x occasion(T2 - T1)      *   -0.071       -0.230 233.639
#>   I06                     group: B - A           0.046        0.484 237.564
#>   I06                occasion: T2 - T1      *    1.043        4.650 239.000
#>   I06 group(B - A) x occasion(T2 - T1)      *   -0.246       -1.119 226.506
#>        p   p_adj lower upper significant practical
#>  < 0.001 < 0.001                       *         *
#>    0.474   1.000                                  
#>    0.818   1.000                                  
#>    0.629   1.000                                  
#>  < 0.001 < 0.001                       *         *
#>    0.264   1.000
da$posthoc
#>  item     term item.1 contrast within estimate se statistic      df       p
#>   I03    group    I03    B - A           0.841        4.090 237.998 < 0.001
#>   I06 occasion    I06  T2 - T1      *    1.043        4.650 239.000 < 0.001
#>    p_adj lower upper significant practical
#>  < 0.001                       *         *
#>  < 0.001                       *         *
```

[`dif_posthoc()`](https://drjoshmcgrane.github.io/rasch/reference/dif_posthoc.md)
is the general follow-up for a significant term. A main effect with more
than two levels is reported as pairwise marginal differences over the
other fitted factors. An interaction is reported as a
difference-in-differences, or its higher-order counterpart. These
comparisons use the joint covariance of the resolved item locations and
Holm adjustment over the stated family. They are preferable to an
ordinary Tukey comparison of residual means because the result is on the
logit scale and respects the factor structure used in the DIF analysis.

For repeated-person contrasts, significance comes from person-level
residual contrast scores. The resolved logit difference remains the
magnitude, but its row-independent calibration covariance is not a
repeated-measures standard error; the package therefore withholds the
logit SE and interval in this case. The fitted person identifier is used
unless `id` is supplied explicitly.

For a many-facet fit, follow-ups may name the underlying item; its
virtual facet cells are pooled with common weights so facet severity
cancels from the group contrast. Extended-frame fits support the
residual ANOVA for factors outside the frame definition, but not an
ordinary resolved-item magnitude: that refit would discard the fitted
frame units.

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
