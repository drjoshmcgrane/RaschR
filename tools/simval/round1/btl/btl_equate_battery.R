suppressWarnings(pkgload::load_all("."), quiet = TRUE))
`%||%` <- function(a, b) if (is.null(a)) b else a

log <- function(...) cat(sprintf(...), "\n")

# custom BTL generator with an EXPLICIT object-location vector (simulate_btl
# only lets us control the SPREAD of an evenly spaced sum-zero scale, not
# an individual object's location, which we need to plant object drift
# between two independently generated panels)
sim_btl_beta <- function(beta, n_judges = 12, reps_per_pair = 25, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  objs <- names(beta); K <- length(objs)
  jids <- sprintf("J%d", seq_len(n_judges))
  pr <- t(utils::combn(objs, 2))
  d <- data.frame(object_a = rep(pr[, 1], each = reps_per_pair),
                   object_b = rep(pr[, 2], each = reps_per_pair),
                   stringsAsFactors = FALSE)
  d$judge <- sample(jids, nrow(d), TRUE)
  lp <- beta[d$object_a] - beta[d$object_b]
  resp <- as.integer(stats::runif(nrow(d)) < stats::plogis(lp))
  d$winner <- ifelse(resp == 1L, d$object_a, d$object_b)
  rownames(d) <- NULL
  d
}

K <- 8
beta0 <- setNames(as.numeric(scale(seq_len(K))) * 1, sprintf("O%d", seq_len(K)))
drift_obj <- "O4"
drift_size <- 0.5
beta_drift <- beta0; beta_drift[drift_obj] <- beta_drift[drift_obj] + drift_size

results <- list()

fit_pair <- function(betaA, betaB, seedA, seedB) {
  dA <- sim_btl_beta(betaA, seed = seedA)
  dB <- sim_btl_beta(betaB, seed = seedB)
  f1 <- btl(dA, "object_a", "object_b", winner = "winner", judge = "judge")
  f2 <- btl(dB, "object_a", "object_b", winner = "winner", judge = "judge")
  list(f1 = f1, f2 = f2)
}

# ---------------------------------------------------------------------
# 1. NULL calibration: two independent panels, same true object locations
# ---------------------------------------------------------------------
set.seed(11)
n_null <- 300
null_flags <- vector("list", n_null)
t0 <- Sys.time()
for (r in seq_len(n_null)) {
  p <- fit_pair(beta0, beta0, 10000 + r, 20000 + r)
  eq <- btl_equate(p$f1, p$f2, independent = TRUE)
  null_flags[[r]] <- eq$table$drifting
}
t1 <- Sys.time()
log("btl_equate NULL complete: %d reps in %.1fs", n_null, as.numeric(t1 - t0, units = "secs"))
nf <- unlist(null_flags); nf <- nf[!is.na(nf)]
null_rate <- mean(nf); null_mcse <- sqrt(null_rate * (1 - null_rate) / length(nf))
log("  object-level drift flag rate under null: %.4f (n_tests=%d, MC se=%.4f)", null_rate, length(nf), null_mcse)
results$btl_null <- list(rate = null_rate, n = length(nf), mcse = null_mcse)

# ---------------------------------------------------------------------
# 2. POWER: planted 0.5-logit drift on one object
# ---------------------------------------------------------------------
set.seed(12)
n_pow <- 150
pow_drift <- logical(n_pow)
pow_other <- vector("list", n_pow)
t0 <- Sys.time()
for (r in seq_len(n_pow)) {
  p <- fit_pair(beta0, beta_drift, 30000 + r, 40000 + r)
  eq <- btl_equate(p$f1, p$f2, independent = TRUE)
  row <- eq$table[eq$table$object == drift_obj, ]
  pow_drift[r] <- isTRUE(row$drifting)
  pow_other[[r]] <- eq$table$drifting[eq$table$object != drift_obj]
}
t1 <- Sys.time()
log("btl_equate POWER complete: %d reps in %.1fs", n_pow, as.numeric(t1 - t0, units = "secs"))
pow_rate <- mean(pow_drift); pow_mcse <- sqrt(pow_rate * (1 - pow_rate) / n_pow)
log("  planted-object (0.5 logit) detection rate: %.4f (MC se=%.4f)", pow_rate, pow_mcse)
po <- unlist(pow_other); po <- po[!is.na(po)]
other_rate <- mean(po)
log("  other-object false-flag rate under planted drift: %.4f (n=%d)", other_rate, length(po))
results$btl_power <- list(rate = pow_rate, n = n_pow, mcse = pow_mcse, other_rate = other_rate)

# ---------------------------------------------------------------------
# 3. Refusal: independence not stated between two fitted btl calibrations
# ---------------------------------------------------------------------
p <- fit_pair(beta0, beta0, 90001, 90002)
eq_refuse <- btl_equate(p$f1, p$f2)  # independent = NULL (default)
refusal_ok <- !isTRUE(eq_refuse$inferential) &&
  any(grepl("independence", eq_refuse$notes, fixed = TRUE))
log("btl_equate refusal (independent unstated): inferential=%s, notes mention independence: %s",
    eq_refuse$inferential, any(grepl("independence", eq_refuse$notes)))
results$btl_refusal <- list(inferential = eq_refuse$inferential, notes = paste(eq_refuse$notes, collapse = " | "))

# dependent explicitly (independent = FALSE): also withheld
eq_dep <- btl_equate(p$f1, p$f2, independent = FALSE)
log("btl_equate with independent=FALSE: inferential=%s (expect FALSE)", eq_dep$inferential)
results$btl_dependent <- list(inferential = eq_dep$inferential)

# ---------------------------------------------------------------------
# 4. Non-convergent fit is refused outright
# ---------------------------------------------------------------------
tiny <- sim_btl_beta(setNames(c(-1, 1), c("Z1", "Z2")), n_judges = 2, reps_per_pair = 1, seed = 5)
bad_ok <- FALSE
err_msg <- NULL
tryCatch({
  f_bad <- btl(tiny, "object_a", "object_b", winner = "winner", judge = "judge")
  f_bad$converged <- FALSE  # force the guard path even if it happened to converge
  btl_equate(p$f1, f_bad, independent = TRUE)
}, error = function(e) { err_msg <<- conditionMessage(e); bad_ok <<- TRUE })
log("btl_equate refuses a non-converged fit: %s (%s)", bad_ok, err_msg %||% "")
results$btl_nonconverged_refusal <- list(refused = bad_ok, msg = err_msg)

saveRDS(results, "tools/simval/round1/btl/btl_results.rds")
log("DONE")
