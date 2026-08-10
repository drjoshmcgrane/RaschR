# Apply a statistic across a simulation batch, resiliently

Applies `FUN` to each replicate of a
[`sim_replicate`](https://drjoshmcgrane.github.io/rasch/reference/sim_replicate.md)
batch, catching replicates on which `FUN` errors – for example a small
or disconnected draw the estimator refuses as unidentified – so a single
failure does not abort the whole Monte-Carlo run. Failed replicates
contribute `NA`; the number of failures and the distinct error messages
are attached as attributes.

## Usage

``` r
sim_apply(batch, FUN, ...)
```

## Arguments

- batch:

  A `"rasch_sim_batch"` from
  [`sim_replicate`](https://drjoshmcgrane.github.io/rasch/reference/sim_replicate.md)
  (or any list of datasets).

- FUN:

  A function of one dataset returning a scalar statistic.

- ...:

  Further arguments passed to `FUN`.

## Value

A vector of per-replicate statistics, with `NA` where the function
failed. Attribute `n_failed` gives the failure count; `failure_messages`
contains the distinct messages.

## Examples

``` r
batch <- sim_replicate(simulate_rasch, 10, n_persons = 300, n_items = 8,
                       seed = 1)
psi <- sim_apply(batch, function(d) rasch(d)$psi$PSI)
mean(psi, na.rm = TRUE)
#> [1] 0.4981672
```
