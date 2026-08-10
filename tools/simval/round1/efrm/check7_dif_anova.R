suppressWarnings(pkgload::load_all("."), quiet=TRUE))
t0 <- Sys.time()

n_per_group <- 250
items_per_set <- 8
n_reps <- 80
alpha_level <- 0.05

make_data <- function(shift, seed) {
  d <- simulate_efrm(n_per_group, items_per_set, n_sets = 2, n_groups = 2,
                      set_unit_ratio = 1.3, group_unit_ratio = 1, seed = seed)
  tr <- attr(d, "truth")
  N <- nrow(d)
  set.seed(seed + 500000L)
  cohort <- factor(sample(c("C1", "C2"), N, replace = TRUE))
  d$cohort <- cohort
  alpha_true <- setNames(tr$alpha, names(tr$item_sets)); phi_true <- tr$phi
  set_of_item <- setNames(rep(names(tr$item_sets), each = items_per_set),
                          unlist(tr$item_sets))
  theta <- tr$theta
  grp_idx <- as.integer(d$group)
  rho <- unname(alpha_true[set_of_item["S1I01"]]) * phi_true[grp_idx]
  delta_i <- tr$difficulty["S1I01"]
  shift_vec <- ifelse(cohort == "C2", shift, 0)
  d$S1I01 <- rbinom(N, 1, plogis(rho * (theta - delta_i - shift_vec)))
  list(d = d, tr = tr)
}

run_condition <- function(shift, n_reps, seed0) {
  p_g1 <- rep(NA_real_, n_reps); p_g2 <- rep(NA_real_, n_reps)
  frame_excluded <- rep(NA, n_reps)
  n_ok <- 0L
  for (r in seq_len(n_reps)) {
    dd <- make_data(shift, seed0 + r)
    fit <- tryCatch(rasch_efrm(dd$d, item_sets = dd$tr$item_sets, groups = "group",
                                factors = "cohort", se_method = "hybrid"),
                     error = function(e) NULL)
    if (is.null(fit)) next
    da <- tryCatch(dif_anova(fit, factors = c("group", "cohort")),
                    error = function(e) NULL)
    if (is.null(da)) next
    n_ok <- n_ok + 1L
    frame_excluded[r] <- any(grepl("frame factor", da$notes)) && !("group" %in% da$summary$term)
    s <- da$summary
    row1 <- s[s$item == "S1I01:g1" & s$term == "cohort", ]
    row2 <- s[s$item == "S1I01:g2" & s$term == "cohort", ]
    if (nrow(row1)) p_g1[r] <- row1$p_uniform_adj[1]
    if (nrow(row2)) p_g2[r] <- row2$p_uniform_adj[1]
  }
  list(p_g1 = p_g1, p_g2 = p_g2, frame_excluded = frame_excluded, n_ok = n_ok)
}

mc_err <- function(p_hat, n) sqrt(p_hat * (1 - p_hat) / n)
report <- function(label, res) {
  f1 <- mean(res$p_g1 < alpha_level, na.rm = TRUE); n1 <- sum(is.finite(res$p_g1))
  f2 <- mean(res$p_g2 < alpha_level, na.rm = TRUE); n2 <- sum(is.finite(res$p_g2))
  cat(sprintf("[%s] n_ok=%d  frame-excluded-in-all-fits=%s\n  S1I01:g1 uniformDIF flag=%.3f (mcse %.3f, n=%d)\n  S1I01:g2 uniformDIF flag=%.3f (mcse %.3f, n=%d)\n",
              label, res$n_ok, all(res$frame_excluded, na.rm = TRUE),
              f1, mc_err(f1, n1), n1, f2, mc_err(f2, n2), n2))
}

r_null <- run_condition(0, n_reps, seed0 = 70000)
report("null (no cohort DIF)", r_null)

r_power <- run_condition(1.1, n_reps, seed0 = 80000)
report("power (cohort DIF shift=1.1 on S1I01)", r_power)

saveRDS(list(null = r_null, power = r_power),
        "tools/simval/round1/efrm/check7_results.rds")
cat("total elapsed:", as.numeric(Sys.time() - t0, units = "secs"), "s\n")
