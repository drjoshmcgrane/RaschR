# STUDY: alpha-correction-limits
#
# Deliberate stress tests of the corrected EFRM item-set unit estimator:
# where does it break, and does it break loudly (refusals) or silently
# (bias without warning)? Cells target the construction's own weak
# points: differential mistargeting, extreme ratios, minimal score
# bases, tiny linking samples, heavy-tailed persons, and the two
# inherited-model-violation modes (within-set misfit from
# discrimination heterogeneity or guessing, and ability-dependent
# missingness), which the correction cannot see because it computes
# score distributions from the fitted model.
# Serial. Rscript tools/simval/studies/alpha-correction-limits.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "alpha-correction-limits"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()
tick <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M")),
                          sprintf(...), "\n")

# hand generator: dichotomous EFRM, full control over targeting, per-item
# discrimination jitter, guessing, and ability-dependent missingness
gen <- function(n, ips, ratio, seed, theta_fun = function(n) rnorm(n, 0, 1.3),
                offset = c(0, 0), disc_jitter = 0, guess2 = 0,
                inform_miss = 0) {
  set.seed(seed)
  al <- c(ratio^-0.5, ratio^0.5)
  th <- theta_fun(n)
  items <- paste0("S", rep(1:2, each = ips), "I",
                  sprintf("%02d", rep(seq_len(ips), 2)))
  X <- matrix(NA_integer_, n, 2L * ips, dimnames = list(NULL, items))
  for (s in 1:2) {
    delta <- seq(-1.5, 1.5, length.out = ips) + offset[s]
    dj <- if (disc_jitter > 0)
      exp(rnorm(ips, 0, disc_jitter) - disc_jitter^2 / 2) else rep(1, ips)
    dj <- dj / exp(mean(log(dj)))            # set geomean preserved = truth
    for (i in seq_len(ips)) {
      p <- plogis(al[s] * dj[i] * (th - delta[i]))
      if (s == 2 && guess2 > 0 && i <= ips / 2) p <- guess2 + (1 - guess2) * p
      x <- rbinom(n, 1, p)
      if (inform_miss > 0) {
        pm <- plogis(-2 + inform_miss * (delta[i] - th))  # low th skip hard items
        x[runif(n) < pm] <- NA
      }
      X[, (s - 1L) * ips + i] <- x
    }
  }
  list(d = data.frame(id = sprintf("P%05d", seq_len(n)), X, group = "g1",
                      check.names = FALSE),
       item_sets = list(set1 = items[seq_len(ips)],
                        set2 = items[ips + seq_len(ips)]))
}

cell <- function(scen, R, seed0, ..., ratio = 1.4, n = 500, ips = 8) {
  lr <- rep(NA_real_, R); n_ref <- 0L
  for (r in seq_len(R)) {
    g <- gen(n, ips, ratio, seed0 + r, ...)
    fe <- tryCatch(rasch_efrm(g$d, item_sets = g$item_sets, groups = "group",
                              id = "id", boot_reps = 0),
                   error = function(e) NULL)
    if (is.null(fe)) { n_ref <- n_ref + 1L; next }
    a <- fe$alpha_table$alpha
    lr[r] <- log(a[2] / a[1])
  }
  ok <- is.finite(lr)
  add(scen, "log unit ratio bias", sum(ok),
      bias = mean(lr[ok]) - log(ratio), emp_sd = sd(lr[ok]),
      n_attempted = R, n_refused = n_ref,
      refusal_rate = n_ref / R)
  tick("%s: bias %+.4f (sd %.4f) refusals %d/%d", scen,
       mean(lr[ok]) - log(ratio), sd(lr[ok]), n_ref, R)
}

## targeting
cell("common mistargeting +2 logits both sets", 60L, 90e3,
     offset = c(2, 2))
cell("differential mistargeting: set2 +1.5 logits", 60L, 91e3,
     offset = c(0, 1.5))
cell("differential mistargeting: set2 +2.5 logits", 60L, 92e3,
     offset = c(0, 2.5))
## extreme ratios
cell("ratio 2.5", 60L, 93e3, ratio = 2.5)
cell("ratio 3.5", 60L, 94e3, ratio = 3.5)
## minimal score basis
cell("4 items per set, N=1000", 80L, 95e3, ips = 4, n = 1000)
cell("3 items per set, N=1000", 80L, 96e3, ips = 3, n = 1000)
## tiny linking samples
cell("N=150", 100L, 97e3, n = 150)
cell("N=80", 100L, 98e3, n = 80)
## heavy tails
cell("t3-tailed persons", 60L, 99e3,
     theta_fun = function(n) rt(n, 3) / sqrt(3) * 1.3)
cell("wide bimodal persons (modes +/-2.2)", 60L, 100e3,
     theta_fun = function(n) sample(c(-2.2, 2.2), n, TRUE) + rnorm(n, 0, 0.6))
## inherited model violations
cell("within-set discrimination jitter sd 0.25", 60L, 101e3,
     disc_jitter = 0.25)
cell("within-set discrimination jitter sd 0.5", 60L, 102e3,
     disc_jitter = 0.5)
cell("guessing 0.15 on half of set 2", 60L, 103e3, guess2 = 0.15)
cell("ability-dependent missingness (low skip hard)", 60L, 104e3,
     inform_miss = 0.8)


## citation rows: the fixed-design floor curve (alpha-floor-curve.R,
## N = 20,000 x 25 replicates per length, isolating the probability-limit
## offset from sampling noise)
add("floor curve, normal persons (citation row)",
    "probability-limit offset by items/set", 25L,
    notes = paste("corrected: -0.009(4) +0.011(6) +0.0124(8) +0.0092(12)",
                  "+0.0119(16) +0.0085(24) +0.0064(32) ~ 1/sqrt(I);",
                  "raw: NaN at 4-6 items (negative variances),",
                  "+0.057(8) then -0.024..-0.035 persisting to 32 items --",
                  "the raw construction never converges in this range"))
add("floor curve, skewed persons chi-sq(3) (citation row)",
    "probability-limit offset by items/set", 25L,
    notes = paste("corrected: -0.003(4) +0.0156(8) +0.0017(16) -0.0067(32);",
                  "raw: NaN(4) +0.093(8) -0.025(16) -0.023(32).",
                  "At practical N the small-sample Jensen term offsets the",
                  "corrected floor, giving net ~0 from 8 items"))

sv_write(do.call(rbind, rows), "alpha-correction-limits")
cat(sprintf("TOTAL elapsed: %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
