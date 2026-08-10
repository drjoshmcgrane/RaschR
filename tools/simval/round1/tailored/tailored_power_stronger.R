suppressWarnings(pkgload::load_all("."), quiet = TRUE))
log <- function(...) cat(sprintf(...), "\n")
I <- 8; N <- 400
d0 <- seq(-2, 2.5, length.out = I)
guess_vec <- rep(0, I); guess_vec[c(I - 1, I)] <- 0.35   # stronger planted guessing
set.seed(42)
dat_pow <- simulate_rasch(N, I, difficulty = d0, guessing = guess_vec, seed = 42)
fit_pow <- rasch(dat_pow, id = "id")
ta_pow_desc <- tailored_analysis(fit_pow, chance = 0.25, se_method = "none")
guessed_items <- sprintf("I%02d", c(I - 1, I))
shift_guessed <- ta_pow_desc$table$shift[ta_pow_desc$table$item %in% guessed_items]
shift_other <- ta_pow_desc$table$shift[!ta_pow_desc$table$item %in% guessed_items]
log("descriptive: n_removed=%d, mean shift guessed=%.3f, mean shift other=%.3f",
    ta_pow_desc$n_removed, mean(shift_guessed), mean(shift_other))
t0 <- Sys.time()
ta_pow_boot <- tailored_analysis(fit_pow, chance = 0.25, se_method = "bootstrap", boot_reps = 399)
t1 <- Sys.time()
log("bootstrap (399 reps): %.1fs", as.numeric(t1 - t0, units = "secs"))
print(ta_pow_boot$table[, c("item", "shift", "se", "p_adj", "significant")])
sig_guessed <- ta_pow_boot$table$significant[ta_pow_boot$table$item %in% guessed_items]
sig_other <- ta_pow_boot$table$significant[!ta_pow_boot$table$item %in% guessed_items]
log("guessed items flagged: %d/%d | other items flagged: %d/%d",
    sum(sig_guessed %in% TRUE), length(sig_guessed), sum(sig_other %in% TRUE), length(sig_other))
saveRDS(list(desc = ta_pow_desc$table, boot = ta_pow_boot$table),
        "tools/simval/round1/tailored/power_stronger.rds")
log("DONE")
