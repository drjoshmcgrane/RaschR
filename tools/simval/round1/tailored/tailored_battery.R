suppressWarnings(pkgload::load_all(".", quiet = TRUE))
`%||%` <- function(a, b) if (is.null(a)) b else a
log <- function(...) cat(sprintf(...), "\n")

results <- list()
I <- 8; N <- 250
d0 <- seq(-2, 2.5, length.out = I)

# ---------------------------------------------------------------------
# 1. NULL (no guessing): one run, boot_reps = 999
# ---------------------------------------------------------------------
set.seed(21)
dat_null <- simulate_rasch(N, I, difficulty = d0, seed = 21)
fit_null <- rasch(dat_null, id = "id")
t0 <- Sys.time()
ta_null <- tailored_analysis(fit_null, chance = 0.25, se_method = "bootstrap", boot_reps = 999)
t1 <- Sys.time()
log("tailored NULL bootstrap (999 reps): %.1fs, %d responses removed, %d anchor items",
    as.numeric(t1 - t0, units = "secs"), ta_null$n_removed, length(ta_null$anchor_items))
n_sig_null <- sum(ta_null$table$significant %in% TRUE)
log("  items flagged significant under null: %d / %d", n_sig_null, nrow(ta_null$table))
print(ta_null$table[, c("item", "initial", "tailored", "origin_equated", "shift", "p_adj", "significant")])
results$tailored_null <- list(n_sig = n_sig_null, n_items = nrow(ta_null$table),
                               n_removed = ta_null$n_removed, boot_used = ta_null$boot_reps_used)

# ---------------------------------------------------------------------
# 2. p-floor warning: demonstrate it fires with too few boot reps (50),
#    and does NOT fire once boot_reps is large enough (399, floor = 16/400=0.04)
# ---------------------------------------------------------------------
floor_warned_low <- FALSE
withCallingHandlers(
  tailored_analysis(fit_null, chance = 0.25, se_method = "bootstrap", boot_reps = 50),
  warning = function(w) {
    if (grepl("smallest achievable", conditionMessage(w)))
      floor_warned_low <<- TRUE
    invokeRestart("muffleWarning")
  })
log("p-floor warning fires at boot_reps=50 (floor 2*8/51=0.314>0.05): %s", floor_warned_low)

floor_warned_ok <- FALSE
withCallingHandlers(
  tailored_analysis(fit_null, chance = 0.25, se_method = "bootstrap", boot_reps = 399),
  warning = function(w) {
    if (grepl("smallest achievable", conditionMessage(w)))
      floor_warned_ok <<- TRUE
    invokeRestart("muffleWarning")
  })
log("p-floor warning at boot_reps=399 (floor 16/400=0.040<0.05): %s (expect FALSE)", floor_warned_ok)
results$tailored_floor <- list(warned_low = floor_warned_low, warned_ok = floor_warned_ok)

# ---------------------------------------------------------------------
# 3. POWER: planted guessing on the hardest items (chance responding when
#    the item is hard), using simulate_rasch's built-in `guessing` floor,
#    which concentrates its effect on low-P (hard-item) cells -- i.e.
#    difficulty-targeted chance responding.
# ---------------------------------------------------------------------
guess_vec <- rep(0, I); guess_vec[c(I - 1, I)] <- 0.25   # 2 hardest items guessed
set.seed(22)
dat_pow <- simulate_rasch(N, I, difficulty = d0, guessing = guess_vec, seed = 22)
fit_pow <- rasch(dat_pow, id = "id")
ta_pow_desc <- tailored_analysis(fit_pow, chance = 0.25, se_method = "none")
log("tailored POWER descriptive: %d responses removed", ta_pow_desc$n_removed)
print(ta_pow_desc$table[, c("item", "initial", "tailored", "origin_equated", "shift", "removed")])
guessed_items <- sprintf("I%02d", c(I - 1, I))
shift_guessed <- ta_pow_desc$table$shift[ta_pow_desc$table$item %in% guessed_items]
shift_other <- ta_pow_desc$table$shift[!ta_pow_desc$table$item %in% guessed_items]
log("  mean shift on guessed items: %.3f | mean shift on other items: %.3f",
    mean(shift_guessed), mean(shift_other))
results$tailored_power_desc <- list(shift_guessed = shift_guessed, shift_other = shift_other,
                                     n_removed = ta_pow_desc$n_removed)

t0 <- Sys.time()
ta_pow_boot <- tailored_analysis(fit_pow, chance = 0.25, se_method = "bootstrap", boot_reps = 399)
t1 <- Sys.time()
log("tailored POWER bootstrap (399 reps): %.1fs", as.numeric(t1 - t0, units = "secs"))
print(ta_pow_boot$table[, c("item", "shift", "se", "p_adj", "significant")])
sig_guessed <- ta_pow_boot$table$significant[ta_pow_boot$table$item %in% guessed_items]
sig_other <- ta_pow_boot$table$significant[!ta_pow_boot$table$item %in% guessed_items]
log("  guessed items flagged significant: %d/%d | other items flagged: %d/%d",
    sum(sig_guessed %in% TRUE), length(sig_guessed),
    sum(sig_other %in% TRUE), length(sig_other))
results$tailored_power_boot <- list(sig_guessed = sig_guessed, sig_other = sig_other,
                                     boot_used = ta_pow_boot$boot_reps_used)

# ---------------------------------------------------------------------
# 4. MCAR 25% missing on top of the null (no-guessing) scenario: confirm
#    tailored_analysis still runs and gives no significant shifts
# ---------------------------------------------------------------------
set.seed(23)
dat_mcar <- simulate_rasch(N, I, difficulty = d0, missing = 0.25, seed = 23)
fit_mcar <- rasch(dat_mcar, id = "id")
mcar_ok <- TRUE; mcar_err <- NULL
ta_mcar <- tryCatch(
  tailored_analysis(fit_mcar, chance = 0.25, se_method = "bootstrap", boot_reps = 399),
  error = function(e) { mcar_err <<- conditionMessage(e); mcar_ok <<- FALSE; NULL })
if (mcar_ok) {
  n_sig_mcar <- sum(ta_mcar$table$significant %in% TRUE)
  log("tailored NULL under 25%% MCAR: ran OK, %d/%d items significant, %d responses removed",
      n_sig_mcar, nrow(ta_mcar$table), ta_mcar$n_removed)
  results$tailored_mcar <- list(ran = TRUE, n_sig = n_sig_mcar, n_items = nrow(ta_mcar$table))
} else {
  log("tailored NULL under 25%% MCAR: refused/error -- %s", mcar_err)
  results$tailored_mcar <- list(ran = FALSE, err = mcar_err)
}

saveRDS(results, "tools/simval/round1/tailored/tailored_results.rds")
log("DONE")
