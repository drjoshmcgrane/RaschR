suppressWarnings(pkgload::load_all("."), quiet=TRUE))

set.seed(1)
REPS <- 300
N <- 500; I <- 10

# ---- NULL: no planted dependence ----
null_res <- lapply(seq_len(REPS), function(r) {
  s <- simulate_rasch(n_persons = N, n_items = I, seed = 1000 + r)
  f <- rasch(s)
  rc <- residual_correlations(f)                 # default flag = NULL
  rc_flagged_default <- nrow(rc$flagged)          # should always be 0 (no threshold)
  rc_heur <- residual_correlations(f, flag = 0.2) # heuristic screening flag
  list(n_flag_default = rc_flagged_default,
       n_flag_heur = nrow(rc_heur$flagged),
       max_abs_q3star = max(abs(rc$pairs$q3_star)),
       sd_q3star = sd(rc$pairs$q3_star))
})

n_flag_default <- sapply(null_res, `[[`, "n_flag_default")
n_flag_heur    <- sapply(null_res, `[[`, "n_flag_heur")
max_abs        <- sapply(null_res, `[[`, "max_abs_q3star")
sd_q3star      <- sapply(null_res, `[[`, "sd_q3star")

cat("=== NULL condition (no planted dependence) ===\n")
cat(sprintf("Default flag=NULL: flags returned in %d/%d reps (should be 0 always)\n",
            sum(n_flag_default > 0), REPS))
cat(sprintf("Heuristic flag=0.2: mean flagged pairs per rep = %.3f (any-pair false-alarm rate = %.3f)\n",
            mean(n_flag_heur), mean(n_flag_heur > 0)))
cat(sprintf("q3_star spread: mean(max|q3*|) = %.3f, mean(sd(q3*)) = %.3f\n",
            mean(max_abs), mean(sd_q3star)))

# ---- POWER: planted dependent pair (items 5,6), strength = 1 ----
POWREPS <- 200
pow_res <- lapply(seq_len(POWREPS), function(r) {
  s <- simulate_rasch(n_persons = N, n_items = I,
                       dependence = list(pairs = list(c(5, 6)), strength = 3),
                       seed = 2000 + r)
  f <- rasch(s)
  rc <- residual_correlations(f, flag = 0.2)
  pr <- rc$pairs
  target <- pr[(pr$item_a == "I05" & pr$item_b == "I06") |
               (pr$item_a == "I06" & pr$item_b == "I05"), ]
  list(target_q3 = target$q3, target_q3star = target$q3_star,
       target_flagged = isTRUE(target$flagged),
       target_is_top = which(pr$item_b == pr$item_b)[1] # placeholder unused
  )
})
target_q3     <- sapply(pow_res, `[[`, "target_q3")
target_q3star <- sapply(pow_res, `[[`, "target_q3star")
target_flag_rate <- mean(sapply(pow_res, `[[`, "target_flagged"))

cat("\n=== POWER condition (planted dependence I05-I06, strength=3) ===\n")
cat(sprintf("Target pair mean Q3 = %.3f, mean Q3* = %.3f\n", mean(target_q3), mean(target_q3star)))
cat(sprintf("Target pair flagged (heuristic 0.2) in %.3f of reps (n=%d), MC error ~ %.3f\n",
            target_flag_rate, POWREPS, sqrt(target_flag_rate*(1-target_flag_rate)/POWREPS)))

mc_err_null <- sqrt(mean(n_flag_heur>0)*(1-mean(n_flag_heur>0))/REPS)
cat(sprintf("\nNull any-pair false-alarm MC error ~ %.3f\n", mc_err_null))
