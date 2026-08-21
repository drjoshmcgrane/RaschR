# STUDY: mfrm-pooled-dif
#
# Null calibration of item-level DIF after pooling MFRM item-by-rater residual
# cells. The imbalanced design gives one person group only two of six raters,
# so the number and composition of observed cells differ sharply by group.
# Run from the package root. Set SV_CORES on Unix-like systems to parallelise.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "mfrm-pooled-dif"
R <- 1000L
N <- 200L; I <- 5L; J <- 6L; m <- 3L
items <- paste0("I", seq_len(I)); raters <- paste0("R", seq_len(J))
delta <- setNames(seq(-1, 1, length.out = I), items)
severity <- setNames(as.numeric(scale(c(-1.1, -0.5, -0.1, 0.2, 0.6, 0.9))) * 0.6,
                     raters)
base_tau <- c(-1.2, 0, 1.2)

generate_null <- function(seed) {
  set.seed(seed)
  theta <- rnorm(N, 0, 1.2)
  grid <- expand.grid(p = seq_len(N), i = seq_len(I), r = seq_len(J))
  score <- integer(nrow(grid))
  for (i in seq_len(I)) for (j in seq_len(J)) {
    rows <- grid$i == i & grid$r == j
    score[rows] <- rasch:::.sim_item(
      theta[grid$p[rows]], base_tau + delta[i] + severity[j]
    )
  }
  data.frame(person = sprintf("P%03d", grid$p), item = items[grid$i],
             rater = raters[grid$r], score = score,
             stringsAsFactors = FALSE)
}

one <- function(r, imbalanced, seed0) {
  # Item and rater parameters are fixed over replicates. Only the sampled
  # persons and their responses change, so the two allocation cells assess
  # the same measurement design.
  d <- generate_null(seed0 + r)
  pn <- match(d$person, unique(d$person))
  d$group <- ifelse(pn <= 100L, "A", "B")
  if (imbalanced)
    d <- d[!(d$group == "B" & d$rater %in% paste0("R", 3:6)), , drop = FALSE]
  f <- tryCatch(rasch_mfrm(d, person = "person", item = "item",
                           score = "score", facets = "rater",
                           factors = "group"), error = function(e) NULL)
  if (is.null(f)) return(c(any = NA, prop = NA, refused = 1, nonconv = 0))
  if (!isTRUE(f$est$converged))
    return(c(any = NA, prop = NA, refused = 0, nonconv = 1))
  z <- tryCatch(dif_anova(f, factors = "group", pool_facets = TRUE),
                error = function(e) NULL)
  if (is.null(z)) return(c(any = NA, prop = NA, refused = 1, nonconv = 0))
  flag <- z$summary$uniform_DIF | z$summary$nonuniform_DIF
  c(any = any(flag), prop = mean(flag), refused = 0, nonconv = 0)
}

run_cell <- function(imbalanced, seed0) {
  cores <- suppressWarnings(as.integer(Sys.getenv("SV_CORES", "1")))
  if (!is.finite(cores) || cores < 1L) cores <- 1L
  z <- if (cores > 1L && .Platform$OS.type != "windows")
    parallel::mclapply(seq_len(R), one, imbalanced = imbalanced, seed0 = seed0,
                       mc.cores = cores, mc.set.seed = FALSE)
  else lapply(seq_len(R), one, imbalanced = imbalanced, seed0 = seed0)
  z <- do.call(rbind, z); ok <- is.finite(z[, "any"])
  scenario <- if (imbalanced)
    "group A: 6 raters; group B: 2 raters" else "both groups: 6 raters"
  rbind(
    sv_row(STUDY, scenario, "familywise null DIF rejection", sum(ok),
           familywise = mean(z[ok, "any"]), n_attempted = R,
           n_refused = sum(z[, "refused"]),
           n_nonconv = sum(z[, "nonconv"])),
    sv_row(STUDY, scenario, "mean per-item null DIF flag proportion", sum(ok),
           type1 = mean(z[ok, "prop"]), n_attempted = R,
           n_refused = sum(z[, "refused"]),
           n_nonconv = sum(z[, "nonconv"]),
           mc_override = list(type1 = sd(z[ok, "prop"]) / sqrt(sum(ok))),
           notes = "Monte Carlo SE uses per-replicate item proportions")
  )
}

sv_write(rbind(run_cell(FALSE, 310000L), run_cell(TRUE, 320000L)), STUDY)
