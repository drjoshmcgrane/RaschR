# Extended conditional-null bootstrap validation for ordinary DIF. This study
# covers PCM and RSM responses, a three-level person factor, simultaneous
# correlated person factors, and partial-null families. The bootstrap samples
# response patterns conditional on each person's sufficient raw score, refits
# the model, and reconstructs the complete DIF analysis.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "100"))
NBOOT <- as.integer(Sys.getenv("SV_BOOT", "99"))
NCORE <- max(1L, as.integer(Sys.getenv("SV_CORES", "1")))
ONLY <- Sys.getenv("SV_ONLY", "")
L <- 6L

.pattern_cache <- new.env(parent = emptyenv())

pattern_structure <- function(m) {
  key <- paste(m, collapse = "-")
  if (exists(key, envir = .pattern_cache, inherits = FALSE))
    return(get(key, envir = .pattern_cache, inherits = FALSE))
  p <- as.matrix(expand.grid(lapply(m, function(x) 0:x)))
  storage.mode(p) <- "integer"
  ans <- list(pattern = p, by_score = split(seq_len(nrow(p)), rowSums(p)))
  assign(key, ans, envir = .pattern_cache)
  ans
}

# Under a PCM/RSM calibration, conditional on total score r,
#   Pr(X = x | R = r) propto exp{-sum_i sum_{k <= x_i} tau_ik}.
# The person location cancels. This reduces to the dichotomous conditional
# sampler when every item has maximum score one.
conditional_sample <- function(scores, tau_list) {
  m <- lengths(tau_list)
  st <- pattern_structure(m)
  p <- st$pattern
  cost <- numeric(nrow(p))
  for (i in seq_along(tau_list)) {
    cum <- c(0, cumsum(tau_list[[i]]))
    cost <- cost + cum[p[, i] + 1L]
  }
  X <- matrix(0L, length(scores), length(m))
  for (s in sort(unique(scores))) {
    take <- which(scores == s)
    j0 <- st$by_score[[as.character(s)]]
    if (!length(j0)) stop("no conditional response pattern for score ", s)
    if (length(j0) == 1L) {
      X[take, ] <- matrix(p[j0, ], length(take), length(m), byrow = TRUE)
    } else {
      lw <- -cost[j0]
      w <- exp(lw - max(lw))
      j <- sample(j0, length(take), replace = TRUE, prob = w)
      X[take, ] <- p[j, , drop = FALSE]
    }
  }
  colnames(X) <- names(tau_list)
  X
}

draw_item <- function(theta, tau, uniform = 0, nonuniform = 0,
                      centre = 0) {
  x <- 0:length(tau)
  cum <- c(0, cumsum(tau))
  eta <- outer(theta, x) - matrix(cum, length(theta), length(x), byrow = TRUE)
  # A uniform location shift adds the same shift to every threshold. The
  # non-uniform departure changes the category slope and crosses at `centre`,
  # avoiding an arbitrary average location effect in the focal group.
  eta <- eta + outer(nonuniform * (theta - centre) - uniform, x)
  eta <- eta - apply(eta, 1L, max)
  pr <- exp(eta)
  pr <- pr / rowSums(pr)
  cs <- t(apply(pr, 1L, cumsum))
  as.integer(rowSums(runif(length(theta)) > cs))
}

truth_thresholds <- function(model, n_categories) {
  m <- n_categories - 1L
  locations <- seq(-1.25, 1.25, length.out = L)
  shared <- if (m == 1L) 0 else seq(-1.05, 1.05, length.out = m)
  lapply(seq_len(L), function(i) {
    step <- if (model == "PCM") {
      z <- shared * (0.82 + 0.07 * i) +
        0.12 * sin(seq_len(m) * (i + 1L))
      z - mean(z)
    } else shared
    locations[i] + step
  })
}

make_design <- function(kind) {
  if (kind == "balanced") {
    factors <- data.frame(group = factor(rep(c("A", "B"), each = 300L)))
    theta_mean <- ifelse(factors$group == "B", 0, 0)
    target <- factors$group == "B"; centre <- 0
  } else if (kind == "imbalanced") {
    factors <- data.frame(group = factor(c(rep("A", 120L), rep("B", 480L))))
    theta_mean <- ifelse(factors$group == "B", 0.8, 0)
    target <- factors$group == "B"; centre <- 0.8
  } else if (kind == "three-level") {
    factors <- data.frame(group = factor(c(rep("A", 100L), rep("B", 200L),
                                           rep("C", 300L))))
    theta_mean <- c(A = -0.4, B = 0, C = 0.6)[as.character(factors$group)]
    target <- factors$group == "C"; centre <- 0.6
  } else if (kind == "multifactor") {
    cell <- rep(c("A.X", "A.Y", "B.X", "B.Y"), c(80L, 40L, 120L, 360L))
    sp <- strsplit(cell, ".", fixed = TRUE)
    factors <- data.frame(
      group = factor(vapply(sp, `[`, "", 1L)),
      region = factor(vapply(sp, `[`, "", 2L)))
    theta_mean <- 0.8 * (factors$group == "B") +
      0.4 * (factors$region == "Y")
    target <- factors$group == "B"; centre <- 1.1
  } else stop("unknown design")
  list(factors = factors, theta = rnorm(nrow(factors)) + theta_mean,
       target = target, centre = centre)
}

simulate_scenario <- function(sc) {
  des <- make_design(sc$design)
  tau <- truth_thresholds(sc$model, sc$n_categories)
  X <- matrix(0L, nrow(des$factors), L,
              dimnames = list(NULL, paste0("I", seq_len(L))))
  for (i in seq_len(L)) {
    u <- nu <- numeric(nrow(X))
    if (i == 3L) {
      u[des$target] <- sc$uniform
      nu[des$target] <- sc$nonuniform
    }
    X[, i] <- draw_item(des$theta, tau[[i]], u, nu, des$centre)
  }
  list(X = X, factors = des$factors)
}

fit_scenario <- function(X, factors, model) {
  dat <- data.frame(X, factors, check.names = FALSE)
  fit_model <- if (model == "dichotomous") "PCM" else model
  tryCatch(suppressWarnings(rasch(
    dat, model = fit_model, items = colnames(X), factors = names(factors))),
    error = function(e) NULL)
}

dif_statistics <- function(fit, effects) {
  da <- tryCatch(suppressWarnings(dif_anova(fit, effects = effects)),
                 error = function(e) NULL)
  if (is.null(da)) return(NULL)
  s <- da$summary
  u <- data.frame(key = paste(s$item, s$term, "uniform", sep = "|"),
                  item = s$item, term = s$term, kind = "uniform",
                  F = s$F_uniform, p = s$p_uniform,
                  p_adj = s$p_uniform_adj)
  n <- data.frame(key = paste(s$item, s$term, "nonuniform", sep = "|"),
                  item = s$item, term = s$term, kind = "nonuniform",
                  F = s$F_nonuniform, p = s$p_nonuniform,
                  p_adj = s$p_nonuniform_adj)
  out <- rbind(u, n)
  out[is.finite(out$F) & is.finite(out$p), , drop = FALSE]
}

one_replicate <- function(seed, sc) {
  set.seed(seed)
  dat <- simulate_scenario(sc)
  fit <- fit_scenario(dat$X, dat$factors, sc$model)
  if (is.null(fit) || !isTRUE(fit$est$converged)) return(NULL)
  obs <- dif_statistics(fit, sc$effects)
  if (is.null(obs) || !nrow(obs)) return(NULL)
  scores <- rowSums(dat$X)
  tau <- fit$tau_list[colnames(dat$X)]
  fb <- pb <- matrix(NA_real_, NBOOT, nrow(obs))
  got <- attempted <- 0L
  while (got < NBOOT && attempted < ceiling(1.35 * NBOOT)) {
    attempted <- attempted + 1L
    xb <- conditional_sample(scores, tau)
    if (!identical(as.integer(rowSums(xb)), as.integer(scores)))
      stop("conditional sampler did not preserve a raw score")
    bf <- fit_scenario(xb, dat$factors, sc$model)
    if (is.null(bf) || !isTRUE(bf$est$converged)) next
    bs <- dif_statistics(bf, sc$effects)
    if (is.null(bs)) next
    j <- match(obs$key, bs$key)
    if (anyNA(j) || any(!is.finite(bs$F[j])) || any(!is.finite(bs$p[j]))) next
    got <- got + 1L
    fb[got, ] <- bs$F[j]
    pb[got, ] <- bs$p[j]
  }
  if (got < NBOOT) return(NULL)

  # Use each refit's upper-tail reference probability because sparse class
  # intervals can change the term degrees of freedom between replicates.
  p_boot_raw <- (1 + colSums(sweep(pb, 2L, obs$p, "<="))) / (NBOOT + 1)
  min_p <- apply(pb, 1L, min)
  p_boot_family <- (1 + vapply(obs$p, function(p) sum(min_p <= p), 0L)) /
    (NBOOT + 1)
  current <- obs$p_adj < 0.05
  boot <- p_boot_family < 0.05
  target_key <- paste("I3", "group", sc$target_kind, sep = "|")
  target <- obs$key == target_key
  unaffected <- obs$item != "I3"
  other_factor <- if (sc$design == "multifactor") obs$term == "region" else
    rep(FALSE, nrow(obs))

  c(current_raw_uniform = mean(obs$p[obs$kind == "uniform"] < 0.05),
    bootstrap_raw_uniform = mean(p_boot_raw[obs$kind == "uniform"] < 0.05),
    current_raw_nonuniform = mean(obs$p[obs$kind == "nonuniform"] < 0.05),
    bootstrap_raw_nonuniform = mean(p_boot_raw[obs$kind == "nonuniform"] < 0.05),
    current_fwer = any(current), bootstrap_fwer = any(boot),
    current_target = if (any(target)) any(current[target]) else NA,
    bootstrap_target = if (any(target)) any(boot[target]) else NA,
    current_unaffected_fwer = any(current[unaffected]),
    bootstrap_unaffected_fwer = any(boot[unaffected]),
    current_other_factor_fwer = if (any(other_factor)) any(current[other_factor]) else NA,
    bootstrap_other_factor_fwer = if (any(other_factor)) any(boot[other_factor]) else NA)
}

run_condition <- function(sc) {
  seeds <- sample.int(.Machine$integer.max, NREP)
  z <- parallel::mclapply(seeds, one_replicate, sc = sc, mc.cores = NCORE,
                          mc.preschedule = TRUE)
  refused <- sum(vapply(z, is.null, TRUE))
  z <- z[!vapply(z, is.null, TRUE)]
  if (!length(z)) stop("no replicate completed in ", sc$name)
  z <- do.call(rbind, z)
  n <- nrow(z)
  null <- sc$uniform == 0 && sc$nonuniform == 0
  effect <- max(abs(sc$uniform), abs(sc$nonuniform))
  row <- function(quantity, col, field, note) {
    args <- list(study = "dif-conditional-bootstrap-extended",
      scenario = sc$name, quantity = quantity, n_reps = n,
      n_attempted = NREP, n_refused = refused, n_nonconv = 0L,
      effect = effect, notes = paste0(sc$model, "; ", sc$design, "; ",
        NBOOT, " conditional refits; ", note))
    args[[field]] <- mean(z[, col], na.rm = TRUE)
    if (field == "type1") args$mc_override <- list(
      type1 = stats::sd(z[, col], na.rm = TRUE) / sqrt(sum(is.finite(z[, col]))))
    do.call(sv_row, args)
  }
  if (null) rbind(
    row("current uniform item-wise Type I", "current_raw_uniform", "type1",
        "hybrid reference"),
    row("bootstrap uniform item-wise Type I", "bootstrap_raw_uniform", "type1",
        "term-specific empirical tail"),
    row("current non-uniform item-wise Type I", "current_raw_nonuniform", "type1",
        "residual-ANOVA reference"),
    row("bootstrap non-uniform item-wise Type I", "bootstrap_raw_nonuniform", "type1",
        "term-specific empirical tail"),
    row("current Holm FWER", "current_fwer", "familywise",
        "complete item-by-term family"),
    row("bootstrap minimum-p FWER", "bootstrap_fwer", "familywise",
        "complete item-by-term family"))
  else rbind(
    row("current target power", "current_target", "power",
        paste("I3", sc$target_kind, "group term")),
    row("bootstrap target power", "bootstrap_target", "power",
        paste("I3", sc$target_kind, "group term")),
    row("current unaffected-item FWER", "current_unaffected_fwer", "familywise",
        "all terms on I1-I2 and I4-I6"),
    row("bootstrap unaffected-item FWER", "bootstrap_unaffected_fwer", "familywise",
        "all terms on I1-I2 and I4-I6"),
    if (sc$design == "multifactor")
      row("current non-target-factor FWER", "current_other_factor_fwer", "familywise",
          "region terms across every item") else NULL,
    if (sc$design == "multifactor")
      row("bootstrap non-target-factor FWER", "bootstrap_other_factor_fwer", "familywise",
          "region terms across every item") else NULL)
}

scenario <- function(model, design, departure = "null", effect = 0) {
  list(name = paste(model, design, departure, format(effect, trim = TRUE)),
       model = model, design = design,
       n_categories = if (model == "dichotomous") 2L else 4L,
       uniform = if (departure == "uniform") effect else 0,
       nonuniform = if (departure == "nonuniform") effect else 0,
       target_kind = if (departure == "nonuniform") "nonuniform" else "uniform",
       effects = "main")
}

scenarios <- list(
  scenario("PCM", "balanced"), scenario("PCM", "imbalanced"),
  scenario("PCM", "imbalanced", "uniform", 0.6),
  scenario("PCM", "imbalanced", "nonuniform", 0.7),
  scenario("PCM", "imbalanced", "nonuniform", 1.4),
  scenario("RSM", "balanced"), scenario("RSM", "imbalanced"),
  scenario("RSM", "imbalanced", "uniform", 0.6),
  scenario("RSM", "imbalanced", "nonuniform", 0.7),
  scenario("RSM", "imbalanced", "nonuniform", 1.4),
  scenario("dichotomous", "three-level"),
  scenario("dichotomous", "three-level", "uniform", 0.6),
  scenario("dichotomous", "three-level", "nonuniform", 0.7),
  scenario("dichotomous", "three-level", "nonuniform", 1.4),
  scenario("dichotomous", "multifactor"),
  scenario("dichotomous", "multifactor", "uniform", 0.6),
  scenario("dichotomous", "multifactor", "nonuniform", 0.7),
  scenario("dichotomous", "multifactor", "nonuniform", 1.4))
if (nzchar(ONLY)) scenarios <- scenarios[vapply(scenarios, function(x)
  grepl(ONLY, x$name, fixed = TRUE), TRUE)]
if (!length(scenarios)) stop("SV_ONLY matched no scenario")

set.seed(8.31e7)
rows <- list()
for (i in seq_along(scenarios)) {
  message(sprintf("[%d/%d] %s", i, length(scenarios), scenarios[[i]]$name))
  rows[[i]] <- run_condition(scenarios[[i]])
  sv_write(do.call(rbind, rows), "dif-conditional-bootstrap-extended")
}
