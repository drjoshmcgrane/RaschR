# Partial-alternative operating characteristics of the public DIF bootstrap in
# the structural model families. One item or object has uniform DIF. Power for
# that member and familywise error among the remaining invariant members are
# calculated separately, because a fitted-global-null adjustment need not give
# strong control once the calibration or matching scores are contaminated.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "100"))
NBOOT <- as.integer(Sys.getenv("SV_BOOT", "99"))
NCORE <- max(1L, as.integer(Sys.getenv("SV_CORES", "4")))
ONLY <- Sys.getenv("SV_ONLY", "")
EFFECTS <- as.numeric(strsplit(Sys.getenv("SV_EFFECTS", ".5,.9"), ",",
                               fixed = TRUE)[[1]])
if (!length(EFFECTS) || any(!is.finite(EFFECTS)) || any(EFFECTS <= 0))
  stop("SV_EFFECTS must contain positive finite comma-separated shifts")

balanced <- function(ids, levels = c("A", "B"))
  stats::setNames(rep(levels, length.out = length(ids)), ids)

make_explanatory <- function(seed, effect) {
  d <- simulate_rasch(
    n_persons = 240, n_items = 6, difficulty = c(-1.2, 1.2),
    dif = list(items = "I03", uniform = effect), n_groups = 2, seed = seed)
  items <- sprintf("I%02d", seq_len(6L))
  q <- data.frame(item = items, feature = seq(-1.2, 1.2, length.out = 6L))
  fit <- rasch_explanatory(
    d[items], q, ~ feature, id = d$id,
    factors = data.frame(group = d$group), n_groups = 2)
  list(fit = fit, dif = dif_anova(fit, n_groups = 2, sizes = FALSE),
       target = "I03", unit = "item")
}

make_mfrm <- function(seed, effect) {
  d <- simulate_mfrm(100, 4, 3, n_categories = 3, seed = seed)
  tr <- attr(d, "truth")
  ids <- unique(as.character(d$person))
  group <- balanced(ids)
  d$group <- unname(group[as.character(d$person)])

  # Regenerate the affected item from the same fixed item, rater and person
  # truths, shifting only group B. This leaves the many-facet structure intact.
  set.seed(seed + 100000L)
  affected <- d$item == "I2" & d$group == "B"
  base_tau <- .sim_thresholds(0, 2L, 1.2)
  for (rater in names(tr$severity)) {
    ii <- which(affected & d$rater == rater)
    if (!length(ii)) next
    th <- tr$theta[match(d$person[ii], tr$person_id)] + effect
    tau <- base_tau + tr$difficulty["I2"] + tr$severity[rater]
    d$score[ii] <- .sim_item(th, tau)
  }
  fit <- rasch_mfrm(
    d, "person", "item", "score", facets = "rater",
    factors = "group", n_groups = 2)
  list(fit = fit, dif = dif_anova(fit, n_groups = 2, sizes = FALSE),
       target = "I2", unit = "item")
}

make_efrm <- function(seed, effect) {
  d <- simulate_efrm(
    n_per_group = 120, items_per_set = 5, n_sets = 2, n_groups = 2,
    set_unit_ratio = 1.25, group_unit_ratio = 1.2, seed = seed)
  tr <- attr(d, "truth")
  site <- character(nrow(d))
  for (g in levels(d$group)) {
    ii <- which(d$group == g)
    site[ii] <- unname(balanced(as.character(d$id[ii])))
  }
  d$site <- factor(site, levels = c("A", "B"))

  # Plant DIF in an external person factor while retaining the fitted frame
  # units. A positive effect makes the target easier for site B.
  target <- "S1I03"
  set.seed(seed + 200000L)
  s <- match(target, unlist(tr$item_sets))
  set_index <- which(vapply(tr$item_sets, function(z) target %in% z, TRUE))
  gidx <- match(as.character(d$group), names(tr$phi))
  rho <- tr$alpha[set_index] * tr$phi[gidx]
  shift <- ifelse(d$site == "B", effect, 0)
  d[[target]] <- rbinom(nrow(d), 1L, stats::plogis(
    rho * (tr$theta + shift - tr$difficulty[target])))

  fit <- rasch_efrm(
    d, item_sets = tr$item_sets, groups = "group", id = "id",
    factors = "site", n_groups = 2, boot_reps = 0)
  list(fit = fit, dif = dif_anova(fit, n_groups = 2, sizes = FALSE),
       target = target, unit = "item")
}

make_btl <- function(seed, effect) {
  d <- simulate_btl(5, 24, reps_per_pair = 40, seed = seed)
  tr <- attr(d, "truth")
  judges <- unique(as.character(d$judge))
  group <- balanced(judges)
  target <- "O3"
  set.seed(seed + 300000L)
  ga <- unname(group[as.character(d$judge)]) == "B"
  lp <- tr$location[d$object_a] - tr$location[d$object_b] +
    effect * ga * ((d$object_a == target) - (d$object_b == target))
  response <- rbinom(nrow(d), 1L, stats::plogis(lp))
  d$winner <- ifelse(response == 1L, d$object_a, d$object_b)
  fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  dif <- btl_dif(fit, group, min_n = 10)
  list(fit = fit, dif = dif, target = target, unit = "object")
}

builders <- list(
  explanatory = make_explanatory,
  mfrm = make_mfrm,
  efrm = make_efrm,
  btl = make_btl)
if (nzchar(ONLY)) {
  wanted <- trimws(strsplit(ONLY, ",", fixed = TRUE)[[1]])
  builders <- builders[intersect(wanted, names(builders))]
}
if (!length(builders)) stop("SV_ONLY selected no structural DIF scenario")

failure <- function(status, stage, condition = NULL) list(
  status = status, stage = stage,
  reason = if (is.null(condition)) "" else conditionMessage(condition),
  B_attempted = condition$B %||% 0L,
  B_used = condition$B_used %||% 0L,
  B_nonconverged = condition$B_nonconverged %||% 0L,
  B_errors = condition$B_errors %||% 0L)

one <- function(seed, builder, effect) {
  made <- tryCatch(builder(seed, effect), error = identity)
  if (inherits(made, "error"))
    return(failure(if (inherits(made, "rasch_refusal")) "refused" else
      "error", "observed fit or DIF", made))
  fit <- made$fit
  converged <- if (inherits(fit, "rasch_btl")) fit$converged else
    fit$est$converged
  if (!isTRUE(converged)) return(failure("nonconverged", "observed fit"))

  db <- tryCatch(suppressWarnings(dif_bootstrap(
    fit, made$dif, B = NBOOT, workers = 1L, seed = seed + 400000L)),
    error = identity)
  if (inherits(db, "error"))
    return(failure(if (inherits(db, "rasch_refusal")) "refused" else
      "error", "DIF bootstrap", db))

  ids <- as.character(made$dif$terms[[made$unit]])
  observed <- is.finite(made$dif$terms$p_adj) & made$dif$terms$p_adj < .05
  booted <- is.finite(db$terms$p_boot_adj) & db$terms$p_boot_adj < .05
  target <- ids == made$target
  list(
    status = "ok", stage = "analysed", reason = "",
    primary_power = any(observed[target]),
    bootstrap_power = any(booted[target]),
    primary_clean_fwer = any(observed[!target]),
    bootstrap_clean_fwer = any(booted[!target]),
    B_attempted = db$B, B_used = db$B_used,
    B_nonconverged = db$B_nonconverged, B_errors = db$B_errors)
}

set.seed(9.03e7)
rows <- list(); row_id <- 0L
for (kind in names(builders)) for (effect in EFFECTS) {
  scenario <- sprintf("%s_uniform_%0.2f", kind, effect)
  message(scenario)
  seeds <- sample.int(1e8, NREP)
  z <- parallel::mclapply(
    seeds, one, builder = builders[[kind]], effect = effect,
    mc.cores = NCORE, mc.preschedule = TRUE)
  status <- vapply(z, `[[`, "", "status")
  ok <- status == "ok"; zz <- z[ok]; n <- length(zz)
  if (!n) stop("no analysed replicate in ", scenario)
  val <- function(name) vapply(zz, `[[`, FALSE, name)
  pp <- val("primary_power"); bp <- val("bootstrap_power")
  pf <- val("primary_clean_fwer"); bf <- val("bootstrap_clean_fwer")
  b_attempted <- vapply(z, `[[`, 0, "B_attempted")
  b_used <- vapply(z, `[[`, 0, "B_used")
  b_nonconv <- vapply(z, `[[`, 0, "B_nonconverged")
  b_errors <- vapply(z, `[[`, 0, "B_errors")
  common <- list(
    study = "dif-bootstrap-structural-alternative", scenario = scenario,
    n_reps = n, n_attempted = NREP,
    n_refused = sum(status == "refused"),
    n_nonconv = sum(status == "nonconverged"),
    n_error = sum(status == "error"),
    n_boot_attempted = sum(b_attempted), n_boot_used = sum(b_used),
    n_boot_nonconv = sum(b_nonconv), n_boot_errors = sum(b_errors),
    effect = effect,
    notes = paste0("one uniform-DIF member; ", NBOOT,
      " fitted-global-null refits per analysed dataset; power for the affected ",
      "member and FWER among invariant members are separate"))
  mk <- function(quantity, ...) do.call(
    sv_row, c(common, list(quantity = quantity), list(...)))
  row_id <- row_id + 1L
  rows[[row_id]] <- rbind(
    mk("primary affected-unit power", power = mean(pp)),
    mk("bootstrap affected-unit power", power = mean(bp)),
    mk("primary invariant-unit FWER", familywise = mean(pf)),
    mk("bootstrap invariant-unit FWER", familywise = mean(bf)),
    mk("paired bootstrap-minus-primary power", bias = mean(bp - pp),
       mean_se = stats::sd(bp - pp) / sqrt(n), se_ratio = NA_real_),
    mk("paired bootstrap-minus-primary invariant FWER", bias = mean(bf - pf),
       mean_se = stats::sd(bf - pf) / sqrt(n), se_ratio = NA_real_),
    mk("mean usable bootstrap proportion", bias = mean(b_used[ok] / NBOOT),
       emp_sd = stats::sd(b_used[ok] / NBOOT)))
  sv_write(do.call(rbind, rows), "dif-bootstrap-structural-alternative")
}
