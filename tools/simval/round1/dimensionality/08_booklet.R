suppressWarnings(pkgload::load_all(".", quiet=TRUE))

# Booklet (structural missing) design: 12-item bank, 3 booklets of 6 items,
# chain-linked with a 2-item anchor between adjacent booklets:
#   Booklet A: I01-I06
#   Booklet B: I05-I10
#   Booklet C: I09-I12, I01-I02
# Items I03/I04 (A-only) and I07/I08 (B-only) never co-occur with I11/I12
# (C-only) -- a structurally disjoint pair for the residual-PCA refusal test.
booklets <- list(A = paste0("I", sprintf("%02d", 1:6)),
                 B = paste0("I", sprintf("%02d", 5:10)),
                 C = paste0("I", sprintf("%02d", c(9:12, 1:2))))

to_booklet <- function(full, assign) {
  X <- full
  cn <- colnames(X)
  for (p in seq_len(nrow(X))) {
    given <- booklets[[assign[p]]]
    X[p, setdiff(cn, given)] <- NA
  }
  X
}

gen_booklet <- function(N, seed, dependence = NULL, second_dim = NULL) {
  s <- simulate_rasch(n_persons = N, n_items = 12, dependence = dependence,
                       second_dim = second_dim, seed = seed)
  set.seed(seed + 500000)
  assign <- sample(names(booklets), N, replace = TRUE)
  Xb <- to_booklet(as.matrix(s[paste0("I", sprintf("%02d", 1:12))]), assign)
  out <- data.frame(id = s$id, Xb, check.names = FALSE)
  out
}

cat("############ BOOKLET (structural missing) condition ############\n\n")

# ---- residual_pca refusal on the disjoint pair ----
d1 <- gen_booklet(900, seed = 31001)
f1 <- rasch(d1)
pc_res <- tryCatch(residual_pca(f1), error = function(e) e)
cat("[residual_pca on booklet design]\n")
if (inherits(pc_res, "error")) {
  cat("  REFUSED (expected -- disjoint item columns):", conditionMessage(pc_res), "\n")
} else {
  cat("  did NOT refuse -- first_eigen =", pc_res$first_eigen, "\n")
}

# ---- dimensionality_test on booklet design (should reuse the pca refusal path) ----
dt_res <- dimensionality_test(f1)
cat("\n[dimensionality_test on booklet design]\n")
cat("  multidimensional =", dt_res$multidimensional,
    "; note =", if (!is.null(dt_res$note)) dt_res$note else "(none)", "\n")

# ---- Q3 / residual_correlations: works on pairwise-complete residuals (no refusal needed) ----
rc_res <- tryCatch(residual_correlations(f1, flag = 0.2), error = function(e) e)
cat("\n[residual_correlations on booklet design]\n")
if (inherits(rc_res, "error")) {
  cat("  REFUSED:", conditionMessage(rc_res), "\n")
} else {
  cat("  ran fine; n flagged pairs =", nrow(rc_res$flagged),
      "; n pairs with any overlap =", nrow(rc_res$pairs), "(of", 12*11/2, "possible)\n")
}

# ---- dependence_magnitude for a pair that DOES co-occur (I05,I06, within booklet A) ----
REPS <- 60
null_d <- numeric(REPS); null_p <- numeric(REPS)
for (r in seq_len(REPS)) {
  d <- gen_booklet(900, seed = 32000 + r)
  f <- rasch(d)
  dm <- tryCatch(dependence_magnitude(f, dependent = "I06", independent = "I05"),
                 error = function(e) e)
  if (inherits(dm, "error")) { null_d[r] <- NA; null_p[r] <- NA; next }
  null_d[r] <- dm$d; null_p[r] <- dm$p
}
ok <- !is.na(null_p)
cat(sprintf("\n[dependence_magnitude null on booklet, co-occurring pair I05-I06] mean d=%.4f, reject@0.05=%.3f (n=%d usable/%d, MCerr=%.3f)\n",
            mean(null_d, na.rm = TRUE), mean(null_p[ok] < 0.05), sum(ok), REPS,
            sqrt(0.05*0.95/sum(ok))))

POWREPS <- 40
pow_d <- numeric(POWREPS)
for (r in seq_len(POWREPS)) {
  d <- gen_booklet(900, seed = 33000 + r, dependence = list(pairs = list(c(5, 6)), strength = 3))
  f <- rasch(d)
  dm <- tryCatch(dependence_magnitude(f, dependent = "I06", independent = "I05"),
                 error = function(e) e)
  pow_d[r] <- if (inherits(dm, "error")) NA else dm$d
}
cat(sprintf("[dependence_magnitude power on booklet, planted I05-I06 strength=3] mean d=%.3f (n usable=%d/%d)\n",
            mean(pow_d, na.rm = TRUE), sum(!is.na(pow_d)), POWREPS))

# ---- combine_items + spread_test on a co-occurring dependent pair ----
d2 <- gen_booklet(900, seed = 34001, dependence = list(pairs = list(c(5, 6)), strength = 3))
f2 <- rasch(d2)
fc2 <- tryCatch(combine_items(f2, list(c("I05", "I06"))), error = function(e) e)
cat("\n[combine_items on booklet, co-occurring dependent pair I05-I06]\n")
if (inherits(fc2, "error")) {
  cat("  REFUSED:", conditionMessage(fc2), "\n")
} else {
  st2 <- spread_test(fc2)
  print(st2)
}

# ---- combine_items on a NON-co-occurring pair (I03 [A-only], I11 [C-only]) -- should refuse ----
fc3 <- tryCatch(combine_items(f1, list(c("I03", "I11"))), error = function(e) e)
cat("\n[combine_items on booklet, NON-co-occurring pair I03-I11 (structural)]\n")
if (inherits(fc3, "error")) {
  cat("  REFUSED (expected):", conditionMessage(fc3), "\n")
} else {
  cat("  did NOT refuse -- n =", nrow(fc3$X), "\n")
}
