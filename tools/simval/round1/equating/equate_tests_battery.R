suppressWarnings(pkgload::load_all(".", quiet = TRUE))
`%||%` <- function(a, b) if (is.null(a)) b else a

alpha <- 0.05
I <- 15
N <- 500
delta0 <- seq(-2.5, 2.5, length.out = I)
drift_item_idx <- 8L
drift_item <- sprintf("I%02d", drift_item_idx)
drift_size <- 0.5

results <- list()

log <- function(...) cat(sprintf(...), "\n")

# ---------------------------------------------------------------------
# 1. NULL calibration, complete data: two independent calibrations of the
#    SAME items (same true difficulties), independent = TRUE.
# ---------------------------------------------------------------------
run_pair <- function(seedA, seedB, deltaA, deltaB, missing = 0) {
  dA <- simulate_rasch(N, I, difficulty = deltaA, missing = missing, seed = seedA)
  dB <- simulate_rasch(N, I, difficulty = deltaB, missing = missing, seed = seedB)
  fitA <- rasch(dA, id = "id")
  fitB <- rasch(dB, id = "id")
  equate_tests(fitA, fitB, independent = TRUE)
}

set.seed(1)
n_null <- 300
null_flags <- vector("list", n_null)
t0 <- Sys.time()
for (r in seq_len(n_null)) {
  eq <- run_pair(100000 + r, 200000 + r, delta0, delta0)
  null_flags[[r]] <- eq$table$drift
}
t1 <- Sys.time()
log("equate_tests NULL complete: %d reps in %.1fs", n_null, as.numeric(t1 - t0, units = "secs"))
nf <- unlist(null_flags)
nf <- nf[!is.na(nf)]
null_rate <- mean(nf)
null_mcse <- sqrt(null_rate * (1 - null_rate) / length(nf))
log("  item-level drift flag rate under null: %.4f (n_tests=%d, MC se=%.4f)", null_rate, length(nf), null_mcse)
results$equate_null <- list(rate = null_rate, n = length(nf), mcse = null_mcse)

# ---------------------------------------------------------------------
# 2. POWER: planted 0.5-logit drift on one item, complete data
# ---------------------------------------------------------------------
deltaB_drift <- delta0
deltaB_drift[drift_item_idx] <- deltaB_drift[drift_item_idx] + drift_size

set.seed(2)
n_pow <- 150
pow_drift_item <- logical(n_pow)
pow_other_flags <- vector("list", n_pow)
t0 <- Sys.time()
for (r in seq_len(n_pow)) {
  eq <- run_pair(300000 + r, 400000 + r, delta0, deltaB_drift)
  row <- eq$table[eq$table$item == drift_item, ]
  pow_drift_item[r] <- isTRUE(row$drift)
  pow_other_flags[[r]] <- eq$table$drift[eq$table$item != drift_item]
}
t1 <- Sys.time()
log("equate_tests POWER complete: %d reps in %.1fs", n_pow, as.numeric(t1 - t0, units = "secs"))
pow_rate <- mean(pow_drift_item)
pow_mcse <- sqrt(pow_rate * (1 - pow_rate) / n_pow)
log("  planted-item (0.5 logit) detection rate: %.4f (MC se=%.4f)", pow_rate, pow_mcse)
of <- unlist(pow_other_flags); of <- of[!is.na(of)]
other_rate <- mean(of)
log("  other-item false-flag rate under planted drift: %.4f (n=%d)", other_rate, length(of))
results$equate_power <- list(rate = pow_rate, n = n_pow, mcse = pow_mcse, other_rate = other_rate)

# ---------------------------------------------------------------------
# 3. MCAR 25% missing: null + power
# ---------------------------------------------------------------------
set.seed(3)
n_mcar_null <- 40  # MCAR fits are ~7s/rep (much slower convergence than complete
                    # data); reps cut relative to the complete-data condition to
                    # fit the time budget -- MC error reported accordingly
mcar_null_flags <- vector("list", n_mcar_null)
t0 <- Sys.time()
for (r in seq_len(n_mcar_null)) {
  eq <- run_pair(500000 + r, 600000 + r, delta0, delta0, missing = 0.25)
  mcar_null_flags[[r]] <- eq$table$drift
}
t1 <- Sys.time()
log("equate_tests NULL (25%% MCAR) complete: %d reps in %.1fs", n_mcar_null, as.numeric(t1 - t0, units = "secs"))
mnf <- unlist(mcar_null_flags); mnf <- mnf[!is.na(mnf)]
mcar_null_rate <- mean(mnf)
mcar_null_mcse <- sqrt(mcar_null_rate * (1 - mcar_null_rate) / length(mnf))
log("  item-level drift flag rate under null (MCAR): %.4f (n_tests=%d, MC se=%.4f)", mcar_null_rate, length(mnf), mcar_null_mcse)
results$equate_null_mcar <- list(rate = mcar_null_rate, n = length(mnf), mcse = mcar_null_mcse)

set.seed(4)
n_mcar_pow <- 25
mcar_pow_drift_item <- logical(n_mcar_pow)
t0 <- Sys.time()
for (r in seq_len(n_mcar_pow)) {
  eq <- run_pair(700000 + r, 800000 + r, delta0, deltaB_drift, missing = 0.25)
  row <- eq$table[eq$table$item == drift_item, ]
  mcar_pow_drift_item[r] <- isTRUE(row$drift)
}
t1 <- Sys.time()
log("equate_tests POWER (25%% MCAR) complete: %d reps in %.1fs", n_mcar_pow, as.numeric(t1 - t0, units = "secs"))
mcar_pow_rate <- mean(mcar_pow_drift_item)
mcar_pow_mcse <- sqrt(mcar_pow_rate * (1 - mcar_pow_rate) / n_mcar_pow)
log("  planted-item detection rate (MCAR): %.4f (MC se=%.4f)", mcar_pow_rate, mcar_pow_mcse)
results$equate_power_mcar <- list(rate = mcar_pow_rate, n = n_mcar_pow, mcse = mcar_pow_mcse)

# ---------------------------------------------------------------------
# 4. Bank route: refusal when cov_location is not attached
# ---------------------------------------------------------------------
dRef <- simulate_rasch(2000, I, difficulty = delta0, seed = 999001)
fitRef <- rasch(dRef, id = "id")
bank_plain <- data.frame(item = fitRef$items$item, location = fitRef$items$location,
                          se = fitRef$items$se, max = fitRef$items$max)
dCur <- simulate_rasch(N, I, difficulty = delta0, seed = 999002)
fitCur <- rasch(dCur, id = "id")
eq_refuse <- equate_tests(fitCur, bank_plain, independent = TRUE)
refusal_ok <- !isTRUE(eq_refuse$inferential) &&
  grepl("cov_location", eq_refuse$note, fixed = TRUE)
log("Bank route WITHOUT cov_location: inferential=%s (expect FALSE); note mentions cov_location: %s",
    eq_refuse$inferential, grepl("cov_location", eq_refuse$note %||% "", fixed = TRUE))
results$bank_refusal <- list(inferential = eq_refuse$inferential, note = eq_refuse$note)

# zero-SE fixed bank: should be treated as identified without cov_location
bank_fixed <- bank_plain; bank_fixed$se <- 0
eq_fixed <- equate_tests(fitCur, bank_fixed, independent = TRUE)
log("Bank route with SE=0 (fixed bank): inferential=%s (expect TRUE)", eq_fixed$inferential)
results$bank_fixed <- list(inferential = eq_fixed$inferential)

# ---------------------------------------------------------------------
# 5. Bank route WITH cov_location attached (diagonal covariance): null + power
# ---------------------------------------------------------------------
bank_cov <- diag(fitRef$items$se^2)
rownames(bank_cov) <- colnames(bank_cov) <- fitRef$items$item
bank <- bank_plain
attr(bank, "cov_location") <- bank_cov

set.seed(5)
n_bank_null <- 150
bank_null_flags <- vector("list", n_bank_null)
t0 <- Sys.time()
for (r in seq_len(n_bank_null)) {
  d <- simulate_rasch(N, I, difficulty = delta0, seed = 1000000 + r)
  fit <- rasch(d, id = "id")
  eq <- equate_tests(fit, bank, independent = TRUE)
  bank_null_flags[[r]] <- eq$table$drift
}
t1 <- Sys.time()
log("Bank route NULL (with cov_location): %d reps in %.1fs", n_bank_null, as.numeric(t1 - t0, units = "secs"))
bnf <- unlist(bank_null_flags); bnf <- bnf[!is.na(bnf)]
bank_null_rate <- mean(bnf)
bank_null_mcse <- sqrt(bank_null_rate * (1 - bank_null_rate) / length(bnf))
log("  item-level drift flag rate under null (bank route): %.4f (n_tests=%d, MC se=%.4f)",
    bank_null_rate, length(bnf), bank_null_mcse)
results$bank_null <- list(rate = bank_null_rate, n = length(bnf), mcse = bank_null_mcse)

set.seed(6)
n_bank_pow <- 100
bank_pow_drift_item <- logical(n_bank_pow)
t0 <- Sys.time()
for (r in seq_len(n_bank_pow)) {
  d <- simulate_rasch(N, I, difficulty = deltaB_drift, seed = 1100000 + r)
  fit <- rasch(d, id = "id")
  eq <- equate_tests(fit, bank, independent = TRUE)
  row <- eq$table[eq$table$item == drift_item, ]
  bank_pow_drift_item[r] <- isTRUE(row$drift)
}
t1 <- Sys.time()
log("Bank route POWER (with cov_location): %d reps in %.1fs", n_bank_pow, as.numeric(t1 - t0, units = "secs"))
bank_pow_rate <- mean(bank_pow_drift_item)
bank_pow_mcse <- sqrt(bank_pow_rate * (1 - bank_pow_rate) / n_bank_pow)
log("  planted-item detection rate (bank route): %.4f (MC se=%.4f)", bank_pow_rate, bank_pow_mcse)
results$bank_power <- list(rate = bank_pow_rate, n = n_bank_pow, mcse = bank_pow_mcse)

saveRDS(results, "tools/simval/round1/equating/equate_results.rds")
log("DONE")
