# Tailored analysis for guessing

Runs the four-step tailored procedure of Andrich, Marais and Humphry
(2012) on a dichotomous analysis. Step 1 is the supplied fit. Step 2
(tailored) sets to missing every observed response whose modelled
probability of success, at the step-1 person and item estimates, is
below `chance`, and re-estimates items and persons. Step 3
(origin-equated) re-analyses the *original* data with the mean location
of the anchor items fixed at its tailored value by average item
anchoring (see
[`pcml`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md)):
every item keeps its initial position relative to the others and the
calibration as a whole moves onto the tailored origin, so the two
calibrations can be compared item by item. Step 4 (all-anchored) fixes
every item at its tailored difficulty and re-estimates persons on the
original data. Guessing is indicated when difficult items are estimated
harder in the tailored analysis than in the origin-equated one; the
comparison table and
[`plot_equate`](https://drjoshmcgrane.github.io/rasch/reference/plot_equate.md)
on the two calibrations show it directly.

## Usage

``` r
tailored_analysis(
  fit,
  chance = 0.25,
  anchor_items = NULL,
  se_method = c("none", "bootstrap"),
  boot_reps = 999L
)
```

## Arguments

- fit:

  An unanchored, unconstrained dichotomous fit from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).
  The procedure estimates its own common origin.

- chance:

  The guessing floor: the probability of success by chance (1/number of
  options; default 0.25).

- anchor_items:

  Items whose mean location fixes the common origin in step 3. The
  default takes the third of the test (at least two items) least
  affected by tailoring – fewest responses removed, ties broken towards
  the easier tailored location – which are the easy items the procedure
  trusts.

- se_method:

  `"none"` (default) reports the item shifts descriptively.
  `"bootstrap"` resamples persons and repeats the complete four-step
  procedure, including automatic anchor selection, to obtain standard
  errors, percentile intervals, and Holm-adjusted tests. When a person
  identifier occurs on several rows, all of that person's rows are
  resampled together.

- boot_reps:

  Person-bootstrap replicates when `se_method = "bootstrap"`; at least
  50, default 999. The sign-count bootstrap p-value has resolution floor
  `2/(boot_reps + 1)`, so after the Holm adjustment across m items the
  smallest achievable adjusted p is `2m/(boot_reps + 1)`; a warning
  fires when that floor exceeds 0.05 (detection would be impossible).

## Value

A list of class `"rasch_tailored"`: `tailored`, `origin_equated`, and
`anchored` fits, the comparison `table` (initial, tailored,
origin-equated locations, the tailored-minus-equated `shift`; bootstrap
uncertainty columns when requested), the number of responses removed,
the anchor items used, `se_method`, and bootstrap accounting: requested,
usable, non-converged, other failures, and the minimum usable count.

## References

Waller, M. I. (1989). Modeling guessing behavior: A comparison of two
IRT models. Applied Psychological Measurement, 13, 233-243. Andrich, D.,
Marais, I. and Humphry, S. (2012). Using a theorem by Andersen and the
dichotomous Rasch model to assess the presence of random guessing in
multiple choice items. Journal of Educational and Behavioral Statistics,
37, 417-442.

## Examples

``` r
set.seed(1); N <- 800
d <- seq(-2, 2.5, length.out = 10); th <- rnorm(N)
P <- plogis(outer(th, d, "-"))
P <- 0.25 + 0.75 * P            # uniform guessing floor
X <- matrix(rbinom(N * 10, 1, P), N, 10)
colnames(X) <- paste0("I", 1:10)
ta <- tailored_analysis(rasch(X), chance = 0.25)
ta$table
#>  item initial tailored origin_equated removed  shift se ci_low ci_high p p_adj
#>    I1  -1.673   -1.770         -1.767       0 -0.003                          
#>    I2  -1.221   -1.360         -1.315       0 -0.044                          
#>    I3  -0.751   -0.826         -0.845       7  0.019                          
#>    I4  -0.340   -0.405         -0.434      33  0.029                          
#>    I5   0.112    0.048          0.018      33  0.031                          
#>    I6   0.272    0.265          0.177      90  0.088                          
#>    I7   0.542    0.513          0.448      90  0.065                          
#>    I8   0.901    0.883          0.806     206  0.077                          
#>    I9   1.037    1.284          0.943     206  0.341                          
#>   I10   1.120    1.369          1.026     206  0.343                          
#>  significant
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
