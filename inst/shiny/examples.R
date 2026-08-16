# Bundled example data used by the graphical interface. Kept in a separate
# source file so the exact datasets can also be reconstructed by the R-code
# disclosures without starting the app.

# --- demo data: 10 polytomous items, one disordered, DIF on Q05 -------------
.demo_data <- function(seed = 11, Np = 1200) {
  set.seed(seed)
  simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }
  mvec <- rep(c(2, 3), length.out = 10)
  tau_true <- lapply(mvec, function(m) sort(rnorm(m, 0, 0.9)))
  tau_true[[2]] <- c(1.2, -1.3, 0.6)                       # disordered item
  th <- rnorm(Np, 0, 1.4)
  grp <- rep(c("reference", "focal"), each = Np / 2)
  sex <- sample(c("female", "male"), Np, replace = TRUE)
  X <- sapply(seq_along(mvec), function(i) {
    sft <- if (i == 5) ifelse(grp == "focal", 0.9, 0) else numeric(Np)  # uniform DIF
    sapply(seq_len(Np), function(n) sample(0:mvec[i], 1, prob = simP(th[n] - sft[n], tau_true[[i]])))
  })
  colnames(X) <- sprintf("Q%02d", seq_along(mvec))
  data.frame(person_id = sprintf("P%04d", seq_len(Np)), X,
             group = grp, sex = sex, check.names = FALSE)
}

# dichotomous demo: 15 multiple-choice items (raw A-D responses), DIF planted
# on I05 by group, and I07 deliberately miskeyed (true correct C, key says A)
.demo_dich <- function(seed = 41, Np = 1000) {
  set.seed(seed)
  d <- seq(-2, 2, length.out = 15)
  grp <- rep(c("reference", "focal"), each = Np / 2)
  sex <- sample(c("female", "male"), Np, replace = TRUE)
  th <- rnorm(Np, 0, 1.3)
  X <- sapply(seq_along(d), function(i) {
    sft <- if (i == 5) ifelse(grp == "focal", 0.8, 0) else 0
    correct <- if (i == 7) "C" else "A"
    ok <- rbinom(Np, 1, plogis(th - d[i] - sft))
    ifelse(ok == 1, correct,
           sample(setdiff(c("A", "B", "C", "D"), correct), Np, replace = TRUE))
  })
  colnames(X) <- sprintf("I%02d", seq_along(d))
  data.frame(person_id = sprintf("P%04d", seq_len(Np)), X,
             group = grp, sex = sex, check.names = FALSE)
}

# the demo key: all "A" (so I07 is the discoverable miskey)
.demo_dich_key <- function()
  setNames(rep("A", 15), sprintf("I%02d", 1:15))

# paired-comparison demo: 8 essays compared pairwise by 10 judges, with
# judge J09 answering at random (discoverable in the judge fit table).
# Besides the winner column it carries a polytomous `preference` column (four
# ordered categories) simulated from the same object locations, so the
# polytomous-response role can be pointed at it, and a `margin` column ("a
# little" < "much") derived from the preference for the winner + margin
# entry path.
.demo_btl <- function(seed = 47, reps = 26) {
  set.seed(seed)
  # a moderate object spread: extreme objects have near-zero residual variance
  # in paired comparisons, which would manufacture spurious DIF, so the range
  # is kept modest and the planted DIF is put on the central objects
  beta <- setNames(seq(-1.0, 1.0, length.out = 8), sprintf("E%02d", 1:8))
  pr <- t(utils::combn(names(beta), 2))
  d <- data.frame(object_a = rep(pr[, 1], each = reps),
                  object_b = rep(pr[, 2], each = reps),
                  stringsAsFactors = FALSE)
  d$judge <- sprintf("J%02d", sample(1:10, nrow(d), replace = TRUE))
  # two judge factors (each constant within judge) so DIF can be modelled by
  # one factor or several jointly: panel splits judges 1-5 / 6-10, experience
  # splits the odd / even judges independently
  d$panel <- ifelse(d$judge %in% sprintf("J%02d", 1:5), "panel A", "panel B")
  d$experience <- ifelse(d$judge %in% sprintf("J%02d", c(1, 3, 5, 7, 9)),
                         "expert", "novice")
  # judgment order: process each judge's comparisons in sequence so the
  # within-judge history (exposure) is well defined
  d <- d[sample(nrow(d)), ]
  d$t <- ave(seq_len(nrow(d)), d$judge, FUN = seq_along)
  d <- d[order(d$judge, d$t), ]
  rownames(d) <- NULL

  # Several signals are built in so the diagnostics have something to find.
  # (1) Exposure: an object already met by the judge gains `expo` logits (a
  # seen-before advantage). (2) Panel DIF: panel A over-rewards E04, so its
  # location differs by panel. (3) Experience DIF: experts over-reward E05, a
  # second, independent judge factor. The two factors point at different
  # objects so each is cleanly attributable. Judge J09 answers at random (a
  # misfitting judge), as before.
  expo <- 0.7
  dif_panel <- "E04"; dif_exp_obj <- "E05"
  dif <- 1.2         # panel effect
  # larger than the panel effect: expert J09 answers at random (diluting it),
  # and the dependence-adjusted DIF screen absorbs the share of a sequential
  # effect that the carry-over covariate can carry
  dif_exp <- 1.8
  tau <- c(-1.1, 0, 1.1)
  lev <- c("much worse", "a little worse", "a little better", "much better")
  seen <- new.env(parent = emptyenv())
  winner <- character(nrow(d)); pref <- integer(nrow(d))
  for (r in seq_len(nrow(d))) {
    j <- d$judge[r]; a <- d$object_a[r]; b <- d$object_b[r]
    ba <- beta[[a]]; bb <- beta[[b]]
    if (d$panel[r] == "panel A") {
      ba <- ba + dif * (a == dif_panel); bb <- bb + dif * (b == dif_panel)
    }
    if (d$experience[r] == "expert") {
      ba <- ba + dif_exp * (a == dif_exp_obj); bb <- bb + dif_exp * (b == dif_exp_obj)
    }
    ba <- ba + expo * isTRUE(get0(paste(j, a), seen, ifnotfound = FALSE))
    bb <- bb + expo * isTRUE(get0(paste(j, b), seen, ifnotfound = FALSE))
    if (j == "J09") { p <- 0.5; Pp <- rep(0.25, 4) }
    else { p <- plogis(ba - bb); Pp <- item_moments(ba - bb, tau)$P }
    winner[r] <- if (runif(1) < p) a else b
    pref[r] <- sample.int(4, 1, prob = Pp)
    assign(paste(j, a), TRUE, seen); assign(paste(j, b), TRUE, seen)
  }
  d$winner <- winner
  d$preference <- factor(lev[pref], levels = lev, ordered = TRUE)
  # margin of win as an ordered factor (extreme categories are "much" wins)
  d$margin <- factor(ifelse(d$preference %in% c("much worse", "much better"),
                            "much", "a little"),
                     levels = c("a little", "much"), ordered = TRUE)
  d <- d[sample(nrow(d)), ]   # present in random row order (t keeps the order)
  rownames(d) <- NULL
  d
}

# rating scale demo: common step structure, item locations vary
.demo_rsm <- function(seed = 51, Np = 1000) {
  set.seed(seed)
  simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }
  loc <- seq(-1.2, 1.2, length.out = 8)
  step <- c(-0.9, 0.0, 0.9)
  grp <- rep(c("reference", "focal"), each = Np / 2)
  th <- rnorm(Np, 0, 1.3)
  X <- sapply(loc, function(b) sapply(th, function(t)
    sample(0:3, 1, prob = simP(t, b + step))))
  colnames(X) <- sprintf("R%02d", seq_along(loc))
  data.frame(person_id = sprintf("P%04d", seq_len(Np)), X,
             group = grp, check.names = FALSE)
}

# rated (MFRM) demo, wide layout: 5 item columns, 6 raters (one erratic),
# incomplete design — one row per person-by-rater combination. The responses
# are simulated in long form (same structure and seed as always) and
# reshaped, so results are unchanged.
.demo_mfrm <- function(seed = 21, Np = 250) {
  set.seed(seed)
  simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }
  persons <- sprintf("P%04d", seq_len(Np)); raters <- paste0("Rater_", 1:6)
  th <- setNames(rnorm(Np, 0, 1.3), persons)
  rho <- setNames(c(-0.9, -0.4, -0.1, 0.1, 0.4, 0.9), raters)
  tau <- list(Essay = c(-1.2, 0.2, 1.1), Argument = c(-0.8, 0.5, 1.3),
              Evidence = c(-1.5, -0.2, 0.9), Style = c(-0.6, 0.4, 1.2),
              Mechanics = c(-1.0, 0.0, 1.0))
  d <- expand.grid(person = persons, item = names(tau), rater = raters,
                   stringsAsFactors = FALSE)
  seen <- unlist(lapply(persons, function(p) paste(p, sample(raters, 3))))
  d <- d[paste(d$person, d$rater) %in% seen, ]
  d$score <- mapply(function(p, i, r) {
    if (r == "Rater_6" && runif(1) < 0.2) return(sample(0:3, 1))  # erratic rater
    sample(0:3, 1, prob = simP(th[p], tau[[i]] + rho[r]))
  }, d$person, d$item, d$rater)
  # wide: one row per person-by-rater with one column per item
  w <- reshape(d, idvar = c("person", "rater"), timevar = "item",
               v.names = "score", direction = "wide")
  names(w) <- sub("^score\\.", "", names(w))
  w <- w[order(w$person, w$rater), c("person", "rater", names(tau))]
  rownames(w) <- NULL
  w
}

# frames demo: 2 person groups x 3 item sets with distinct units
.demo_efrm <- function(seed = 31, per_g = 350) {
  set.seed(seed)
  simP <- function(th, tau, r) { x <- 0:length(tau); p <- exp(r * (x * th - c(0, cumsum(tau)))); p / sum(p) }
  glev <- c("year5", "year7"); grp <- rep(glev, each = per_g); Np <- length(grp)
  phi <- c(year5 = 0.8, year7 = 1.25)
  sets <- rep(c("Number", "Algebra", "Space"), each = 6)
  alpha <- c(Number = 0.75, Algebra = 1.0, Space = 4 / 3)
  th <- rnorm(Np, 0, 1.3) + ifelse(grp == "year7", 0.5, 0)
  d <- as.numeric(sapply(c(-0.3, 0.1, 0.2), function(m) m + seq(-1.2, 1.2, length.out = 6)))
  X <- sapply(seq_along(sets), function(i) sapply(seq_len(Np), function(n)
    sample(0:2, 1, prob = simP(th[n], d[i] + c(-0.5, 0.5),
                               alpha[sets[i]] * phi[grp[n]]))))
  colnames(X) <- sprintf("%s_%02d", sets, seq_along(sets))
  data.frame(person_id = sprintf("P%04d", seq_len(Np)), X, year_group = grp,
             check.names = FALSE)
}
