# STUDY: resolve-versus-drop
#
# When an item breaks frame invariance, is it better removed or resolved?
#
# misfit-repair.csv established that dropping a flagged item restores a
# planted unit ratio wherever the item is found, so detection rather than
# repair is the binding constraint. It did not ask what the repair costs.
# Dropping takes the item out of every frame, so it stops contributing to
# any person's measure: eleven items where there were twelve, for everyone,
# including the frames the item behaved perfectly well in.
#
# Resolving is the alternative the DIF literature uses. Give the item a
# separate location in each frame and it stops linking the frames -- which
# is what the diagnosis says is wrong with it -- while continuing to measure
# the person within their own frame. Strictly more information is retained.
#
# The question is whether that extra information is worth having, so the
# study measures both the unit ratio and the person side:
#
#   clean      no departure planted -- the reference
#   damaged    one item shifted in frame 2, no repair
#   dropped    drop_items()
#   resolved   resolve_frames()
#
# Reported per cell: the recovered group-unit ratio against the planted
# 1.40, and the person-side cost -- the mean standard error of the person
# measures and their correlation with the person locations that generated
# them, which is the quantity dropping an item can only reduce.
#
# 200 replicates, 500 persons per frame, 12 items. Serial, ~25 min.
#   Rscript tools/simval/studies/resolve-versus-drop.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "resolve-versus-drop"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()

ratio <- 1.40
lt <- log(ratio)
u <- c(ratio^-0.5, ratio^0.5)
N <- 500L
K <- 12L
delta <- seq(-1.6, 1.6, length.out = K)
items <- sprintf("I%02d", seq_len(K))
BAD <- "I04"
SHIFT <- 1.2                      # logits, in frame 2 only
R <- 200L

gen <- function(seed, shift) {
  set.seed(seed)
  th <- vector("list", 2L)
  X <- vector("list", 2L)
  for (g in 1:2) {
    th[[g]] <- stats::rnorm(N, 0, 1.3)
    sh <- rep(0, K)
    if (g == 2L) sh[match(BAD, items)] <- shift
    X[[g]] <- vapply(seq_len(K), function(i)
      stats::rbinom(N, 1, stats::plogis(u[g] * (th[[g]] - delta[i] - sh[i]))),
      numeric(N))
  }
  M <- rbind(X[[1]], X[[2]])
  colnames(M) <- items
  list(d = data.frame(id = sprintf("P%05d", seq_len(2L * N)), M,
                      group = rep(c("g1", "g2"), each = N),
                      check.names = FALSE),
       theta = c(th[[1]], th[[2]]))
}

fit_it <- function(d)
  tryCatch(rasch_efrm(d, items = items, item_sets = list(set1 = items),
                      groups = "group", id = "id", boot_reps = 0),
           error = function(e) NULL)

phi_ratio <- function(f) {
  p <- f$phi_table$phi[match(c("g1", "g2"), f$phi_table$group)]
  log(p[2] / p[1])
}

# person-side cost: what the repair leaves behind for measuring people
person_cost <- function(f, truth) {
  ok <- !f$person$extreme & is.finite(f$person$theta)
  th <- f$person$theta[ok]
  tr <- truth[match(f$person$id[ok], sprintf("P%05d", seq_along(truth)))]
  c(se = mean(f$person$se[ok], na.rm = TRUE),
    r = suppressWarnings(stats::cor(th, tr)))
}

run_cell <- function(label, shift, repair) {
  lr <- rep(NA_real_, R); se <- rep(NA_real_, R); rr <- rep(NA_real_, R)
  n_ref <- 0L
  for (r in seq_len(R)) {
    z <- gen(600000 + r, shift)
    f <- fit_it(z$d)
    if (is.null(f)) { n_ref <- n_ref + 1L; next }
    g <- if (is.null(repair)) f else
      tryCatch(repair(f), error = function(e) NULL)
    if (is.null(g)) { n_ref <- n_ref + 1L; next }
    lr[r] <- phi_ratio(g)
    pc <- person_cost(g, z$theta)
    se[r] <- pc[["se"]]; rr[r] <- pc[["r"]]
  }
  ok <- is.finite(lr)
  add(label, "log phi ratio", sum(ok),
      bias = mean(lr[ok]) - lt, emp_sd = stats::sd(lr[ok]),
      n_attempted = R, n_refused = n_ref,
      notes = c(sprintf("recovered %.3f vs planted %.3f",
                        exp(mean(lr[ok])), ratio),
                sprintf("mean person se %.4f", mean(se[ok], na.rm = TRUE)),
                sprintf("correlation with the generating locations %.4f",
                        mean(rr[ok], na.rm = TRUE))))
  cat(sprintf("[%s] %-34s %.3f  (bias %+.4f)  se %.4f  r %.4f  n=%d\n",
              format(Sys.time(), "%H:%M"), label, exp(mean(lr[ok])),
              mean(lr[ok]) - lt, mean(se[ok], na.rm = TRUE),
              mean(rr[ok], na.rm = TRUE), sum(ok)))
}

run_cell("clean", 0, NULL)
run_cell("damaged, no repair", SHIFT, NULL)
run_cell("dropped", SHIFT, function(f) drop_items(f, BAD, boot_reps = 0))
run_cell("resolved", SHIFT, function(f) resolve_frames(f, BAD, boot_reps = 0))

sv_write(do.call(rbind, rows), "resolve-versus-drop")
cat(sprintf("TOTAL elapsed: %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))
