suppressWarnings(pkgload::load_all(".", quiet = TRUE))
`%||%` <- function(a, b) if (is.null(a)) b else a
log <- function(...) cat(sprintf(...), "\n")
results <- list()

# ---------------------------------------------------------------------
# 1. POWER: planted informative distractor (Andrich & Styles 2011 style
#    partial-credit item: D < B < A in true ability, B a genuine middle
#    performance mapped onto an "informative distractor").
# ---------------------------------------------------------------------
gen_informative <- function(Np, seed, missing = 0) {
  set.seed(seed)
  th <- rnorm(Np)
  raw <- sapply(seq(-0.5, 0.5, length.out = 4), function(d) {
    x <- vapply(th, function(b) sample(0:2, 1, prob = item_moments(b, c(d - 0.7, d + 0.7))$P), 0L)
    c("D", "B", "A")[x + 1]
  })
  colnames(raw) <- paste0("M", 1:4)
  if (missing > 0) raw[sample(length(raw), round(missing * length(raw)))] <- NA
  raw
}

set.seed(31)
n_pow <- 150
credited_B <- logical(n_pow)          # B (the informative distractor) credited on every item
credited_D <- logical(n_pow)          # D (the true bottom, uninformative baseline) NOT credited
t0 <- Sys.time()
for (r in seq_len(n_pow)) {
  raw <- gen_informative(600, seed = 100000 + r)
  fit <- rasch(raw, key = setNames(rep("A", 4), colnames(raw)))
  pr <- distractor_rescore(fit)
  b_rows <- pr$option_scores[pr$option_scores$option == "B", ]
  d_rows <- pr$option_scores[pr$option_scores$option == "D", ]
  credited_B[r] <- all(b_rows$score > 0)
  credited_D[r] <- any(d_rows$score > 0)
}
t1 <- Sys.time()
log("distractor_rescore POWER: %d reps in %.1fs", n_pow, as.numeric(t1 - t0, units = "secs"))
rate_B <- mean(credited_B); mcse_B <- sqrt(rate_B * (1 - rate_B) / n_pow)
rate_D <- mean(credited_D); mcse_D <- sqrt(rate_D * (1 - rate_D) / n_pow)
log("  informative distractor B credited on all 4 items: %.4f (MC se=%.4f)", rate_B, mcse_B)
log("  uninformative baseline D ever wrongly credited:   %.4f (MC se=%.4f)", rate_D, mcse_D)
results$rescore_power <- list(rate_B = rate_B, mcse_B = mcse_B, rate_D = rate_D, mcse_D = mcse_D, n = n_pow)

# ---------------------------------------------------------------------
# 2. NULL: distractors carry no information (uniform among wrong options)
# ---------------------------------------------------------------------
gen_null <- function(Np, seed, I = 4) {
  set.seed(seed)
  th <- rnorm(Np)
  d0 <- seq(-0.5, 0.5, length.out = I)
  raw <- sapply(d0, function(d) {
    ok <- rbinom(Np, 1, plogis(th - d))
    ifelse(ok == 1, "A", sample(c("B", "C", "D"), Np, replace = TRUE))
  })
  colnames(raw) <- paste0("M", seq_len(I))
  raw
}

set.seed(32)
n_null <- 150
any_credited_null <- logical(n_null)
t0 <- Sys.time()
for (r in seq_len(n_null)) {
  raw <- gen_null(600, seed = 200000 + r)
  fit <- rasch(raw, key = setNames(rep("A", 4), colnames(raw)))
  pr <- distractor_rescore(fit)
  any_credited_null[r] <- any(pr$option_scores$score > 0 & !pr$option_scores$option %in% c("A"))
}
t1 <- Sys.time()
log("distractor_rescore NULL: %d reps in %.1fs", n_null, as.numeric(t1 - t0, units = "secs"))
null_rate <- mean(any_credited_null); null_mcse <- sqrt(null_rate * (1 - null_rate) / n_null)
log("  false-credit rate (any distractor credited when none informative): %.4f (MC se=%.4f)", null_rate, null_mcse)
results$rescore_null <- list(rate = null_rate, mcse = null_mcse, n = n_null)

# ---------------------------------------------------------------------
# 3. Double-keying: two options both correct ("A/C"), verify scoring and
#    that the fit correctly identifies both as credited (no false miskey flag)
# ---------------------------------------------------------------------
set.seed(33)
Np <- 500
th <- rnorm(Np)
d0 <- seq(-1, 1, length.out = 5)
raw_dk <- sapply(d0, function(d) {
  ok <- rbinom(Np, 1, plogis(th - d))
  ifelse(ok == 1, sample(c("A", "C"), Np, replace = TRUE, prob = c(0.5, 0.5)),
         sample(c("B", "D"), Np, replace = TRUE))
})
colnames(raw_dk) <- paste0("DK", 1:5)
key_dk <- setNames(rep("A/C", 5), colnames(raw_dk))
fit_dk <- rasch(raw_dk, key = key_dk)
# both A and C must score 1, B and D must score 0
scored_ok <- all(vapply(colnames(raw_dk), function(it) {
  ra <- raw_dk[, it]; sc <- fit_dk$mc$scored[, it]
  all(sc[ra %in% c("A", "C")] == 1, na.rm = TRUE) &&
    all(sc[ra %in% c("B", "D")] == 0, na.rm = TRUE)
}, TRUE))
log("Double-keying (A/C) scores correctly on all items: %s", scored_ok)
da_dk <- distractor_analysis(fit_dk)
# neither A nor C should be flagged as a miskey relative to the other (both keyed)
false_flags_dk <- any(da_dk$flag[da_dk$option %in% c("A", "C")])
log("Double-keyed options (A, C) never flagged as miskey: %s (flags=%d)",
    !false_flags_dk, sum(da_dk$flag[da_dk$option %in% c("A", "C")]))
results$double_key <- list(scored_ok = scored_ok, false_flags = false_flags_dk)

# ---------------------------------------------------------------------
# 4. 25% MCAR on the informative-distractor design: rescore proposal still
#    recovers the planted credit
# ---------------------------------------------------------------------
set.seed(34)
n_mcar <- 60
credited_B_mcar <- logical(n_mcar)
t0 <- Sys.time()
for (r in seq_len(n_mcar)) {
  raw <- gen_informative(600, seed = 300000 + r, missing = 0.25)
  fit <- rasch(raw, key = setNames(rep("A", 4), colnames(raw)))
  pr <- distractor_rescore(fit)
  b_rows <- pr$option_scores[pr$option_scores$option == "B", ]
  credited_B_mcar[r] <- all(b_rows$score > 0)
}
t1 <- Sys.time()
log("distractor_rescore POWER under 25%% MCAR: %d reps in %.1fs", n_mcar, as.numeric(t1 - t0, units = "secs"))
rate_B_mcar <- mean(credited_B_mcar); mcse_B_mcar <- sqrt(rate_B_mcar * (1 - rate_B_mcar) / n_mcar)
log("  informative distractor B credited on all 4 items (MCAR): %.4f (MC se=%.4f)", rate_B_mcar, mcse_B_mcar)
results$rescore_power_mcar <- list(rate_B = rate_B_mcar, mcse_B = mcse_B_mcar, n = n_mcar)

saveRDS(results, "tools/simval/round1/mc/mc_results.rds")
log("DONE")
