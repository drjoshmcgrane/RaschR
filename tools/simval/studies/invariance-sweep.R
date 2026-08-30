# Structural invariance sweep: identities every fit must satisfy exactly,
# checked across randomised designs. These are not calibration questions --
# a violation here is a bug, not a property of an estimator. Covers:
# probability identities of the response model; invariance of every reported
# number under item relabelling, item-column permutation, and person-row
# permutation; the sum-zero identification; score-table monotonicity;
# internal consistency of the assembled tables; and for paired comparisons,
# invariance under reversing the presentation of every comparison.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
NREP <- as.integer(Sys.getenv("SV_REPS", "12"))
fails <- character(0)
ok <- function(cond, what, seed) {
  if (!isTRUE(cond)) fails <<- c(fails, sprintf("[seed %d] %s", seed, what))
}
num_eq <- function(a, b, tol = 1e-8) {
  isTRUE(all.equal(a, b, tolerance = tol, check.attributes = FALSE))
}

skipped <- 0L
for (r in seq_len(NREP)) {
  seed <- 40000 + r
  set.seed(seed)
  # draw every design parameter BEFORE the simulate call: simulate_rasch
  # reseeds from `seed` before its arguments are (lazily) evaluated, so an
  # inline sample() would draw from the reseeded stream and the design could
  # not be reproduced outside the loop
  model <- sample(c("dichotomous", "PCM", "RSM"), 1)
  N <- sample(150:400, 1); L <- sample(6:12, 1)
  nc <- sample(3:5, 1); miss <- sample(c(0, 0.05), 1)
  d <- simulate_rasch(N, L, model = model, n_categories = nc,
                      missing = miss, seed = seed)
  X <- as.matrix(d[, grep("^I", names(d))])
  f <- tryCatch(rasch(d, id = "id",
                      model = if (model == "RSM") "RSM" else "PCM"),
                error = function(e) e)
  if (inherits(f, "error")) {
    # a documented refusal (an RSM whose data lost a top category has
    # unequal maxima) is the package answering correctly, not a failure
    if (grepl("RSM requires equal max score", conditionMessage(f))) {
      skipped <- skipped + 1L; next
    }
    fails <- c(fails, sprintf("[seed %d] fit refused: %s", seed,
                              conditionMessage(f)))
    next
  }

  # --- response-model identities ----------------------------------------
  th <- seq(-3, 3, by = 1.5)
  for (i in seq_len(ncol(X))) {
    tau <- f$tau_list[[i]]
    for (t in th) {
      mo <- item_moments(t, tau)
      ok(num_eq(sum(mo$P), 1, 1e-12), sprintf("P sums to 1 (item %d)", i), seed)
      ok(all(mo$P >= 0), sprintf("P nonnegative (item %d)", i), seed)
      ok(num_eq(mo$E, sum(seq_along(mo$P) * mo$P - mo$P), 1e-10) ||
         num_eq(mo$E, sum((seq_along(mo$P) - 1) * mo$P), 1e-10),
         sprintf("E = sum k P (item %d)", i), seed)
      ok(mo$V >= -1e-12, sprintf("V nonnegative (item %d)", i), seed)
    }
  }

  # --- identification and table consistency -----------------------------
  # identification: unanchored fits centre the ITEM LOCATIONS at zero
  if (is.null(f$refit_spec$anchors))
    ok(num_eq(mean(f$items$location), 0, 1e-6),
       "item locations centre at zero", seed)
  ok(num_eq(f$total_chisq, sum(f$items$chisq, na.rm = TRUE), 1e-8),
     "total chisq is the sum over items", seed)
  usable <- is.finite(f$items$p)
  ok(num_eq(f$items$p_adj[usable],
            p.adjust(f$items$p[usable], method = "holm")),
     "p_adj is Holm over the items", seed)
  ok(f$psi$PSI >= 0 && f$psi$PSI <= 1, "PSI in [0,1]", seed)
  sc <- f$score_table
  if (!is.null(sc) && all(c("score", "theta") %in% names(sc))) {
    fin <- is.finite(sc$theta)
    ok(!is.unsorted(sc$theta[fin]), "score table monotone in raw score", seed)
  }

  # --- item relabelling leaves every number unchanged -------------------
  d2 <- d
  new <- sprintf("Q%02d", rev(seq_len(L)))          # reversed, different names
  names(d2)[match(colnames(X), names(d2))] <- new
  f2 <- rasch(d2, id = "id", model = if (model == "RSM") "RSM" else "PCM",
              items = new)
  ok(num_eq(f2$items$location, f$items$location) &&
     num_eq(f2$items$chisq, f$items$chisq) &&
     num_eq(f2$person$theta, f$person$theta),
     "item relabelling changes nothing but names", seed)

  # --- item-column permutation permutes the item table ------------------
  pi <- sample(L)
  keep <- intersect(c("id", "group"), names(d))
  d3 <- d[, c(setdiff(keep, "group"), colnames(X)[pi],
              intersect("group", keep))]
  f3 <- rasch(d3, id = "id", model = if (model == "RSM") "RSM" else "PCM")
  m3 <- match(f$items$item, f3$items$item)
  ok(num_eq(f3$items$location[m3], f$items$location) &&
     num_eq(f3$items$infit_ms[m3], f$items$infit_ms) &&
     num_eq(f3$person$theta, f$person$theta),
     "column permutation permutes items, persons untouched", seed)

  # --- person-row permutation permutes persons, items untouched ---------
  pr <- sample(N)
  f4 <- rasch(d[pr, ], id = "id", model = if (model == "RSM") "RSM" else "PCM")
  m4 <- match(f$person$id, f4$person$id)
  ok(num_eq(f4$items$location, f$items$location) &&
     num_eq(f4$items$chisq, f$items$chisq) &&
     num_eq(f4$person$theta[m4], f$person$theta),
     "row permutation permutes persons, items untouched", seed)
}

# --- paired comparisons: presentation reversal -------------------------
for (r in seq_len(max(3L, NREP %/% 4L))) {
  seed <- 41000 + r
  d <- simulate_btl(7, 10, reps_per_pair = 10, seed = seed)
  f <- tryCatch(btl(d, "object_a", "object_b", winner = "winner",
                    judge = "judge"), error = function(e) NULL)
  if (is.null(f)) { fails <- c(fails, sprintf("[seed %d] btl refused", seed)); next }
  # swap the presented order of every comparison and flip the winner: the
  # design is the same set of judgements, so every location must survive
  d5 <- d
  d5$object_a <- d$object_b; d5$object_b <- d$object_a
  f5 <- btl(d5, "object_a", "object_b", winner = "winner", judge = "judge")
  m5 <- match(f$objects$object, f5$objects$object)
  ok(num_eq(f5$objects$location[m5], f$objects$location, 1e-8),
     "reversing every presentation leaves object locations", seed)
  ok(num_eq(sum(f$objects$location), 0, 1e-8), "object locations sum to zero", seed)
}

cat(sprintf("\n%d design(s) swept, %d skipped on documented refusals\n",
            NREP, skipped))
if (length(fails)) {
  cat("FAILURES:\n"); cat(paste(" ", fails), sep = "\n"); quit(status = 1)
} else cat("every invariance held\n")
