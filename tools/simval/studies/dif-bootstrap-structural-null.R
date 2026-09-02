# Null calibration of the DIF bootstrap in the structural model families.
# Every outer dataset is fitted through the public model and DIF functions;
# dif_bootstrap() then repeats the complete declared family under its fitted
# null. Performance is conditional on datasets that pass both the observed
# analysis and the bootstrap's predeclared usable-refit rule.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "50"))
NBOOT <- as.integer(Sys.getenv("SV_BOOT", "99"))
NCORE <- max(1L, as.integer(Sys.getenv("SV_CORES", "4")))
ONLY <- Sys.getenv("SV_ONLY", "")

balanced <- function(ids, levels = c("A", "B")) {
  z <- rep(levels, length.out = length(ids))
  stats::setNames(sample(z, length(z), replace = FALSE), ids)
}

make_explanatory <- function(seed) {
  set.seed(seed)
  N <- 240L
  group <- factor(rep(c("A", "B"), each = N / 2L))
  theta <- rnorm(N) + ifelse(group == "B", .45, 0)
  feature <- seq(-1.2, 1.2, length.out = 6L)
  X <- vapply(feature, function(d)
    rbinom(N, 1L, plogis(theta - d)), integer(N))
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  q <- data.frame(item = colnames(X), feature = feature)
  fit <- rasch_explanatory(
    X, q, ~ feature, factors = data.frame(group = group), n_groups = 2)
  list(fit = fit, dif = dif_anova(fit, n_groups = 2))
}

make_mfrm <- function(seed) {
  d <- simulate_mfrm(100, 4, 3, n_categories = 3, seed = seed)
  ids <- unique(as.character(d$person))
  set.seed(seed + 100000L)
  g <- balanced(ids)
  d$group <- unname(g[as.character(d$person)])
  fit <- rasch_mfrm(
    d, "person", "item", "score", facets = "rater",
    factors = "group", n_groups = 2)
  list(fit = fit, dif = dif_anova(fit, n_groups = 2))
}

make_efrm <- function(seed) {
  d <- simulate_efrm(
    n_per_group = 120, items_per_set = 5, n_sets = 2, n_groups = 2,
    set_unit_ratio = 1.25, group_unit_ratio = 1.2, seed = seed)
  tr <- attr(d, "truth")
  set.seed(seed + 200000L)
  site <- character(nrow(d))
  for (g in unique(as.character(d$group))) {
    ii <- which(as.character(d$group) == g)
    site[ii] <- unname(balanced(as.character(d$id[ii])))
  }
  d$site <- site
  fit <- rasch_efrm(
    d, item_sets = tr$item_sets, groups = "group", id = "id",
    factors = "site", n_groups = 2, boot_reps = 0)
  list(fit = fit, dif = dif_anova(fit, n_groups = 2))
}

make_btl <- function(seed, clustered = FALSE) {
  d <- simulate_btl(
    5, 24, reps_per_pair = 40,
    second_attribute = if (clustered) list(rho = .5) else NULL,
    seed = seed)
  fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  ids <- unique(as.character(fit$comparisons$judge))
  # In the clustered stress condition the simulator alternates its two
  # preference camps. This four-row pattern crosses the nominated group with
  # those camps exactly, retaining judge heterogeneity without true group DIF.
  group <- if (clustered) {
    number <- suppressWarnings(as.integer(sub("^J", "", ids)))
    stopifnot(all(is.finite(number)))
    stats::setNames(ifelse(number %% 4L %in% c(1L, 2L), "A", "B"), ids)
  }
  else {
    set.seed(seed + 300000L)
    balanced(ids)
  }
  dif <- btl_dif(fit, group, min_n = 10)
  list(fit = fit, dif = dif)
}

builders <- list(
  explanatory = make_explanatory,
  mfrm = make_mfrm,
  efrm = make_efrm,
  btl = function(seed) make_btl(seed, FALSE),
  btl_clustered = function(seed) make_btl(seed, TRUE))
if (nzchar(ONLY)) {
  wanted <- trimws(strsplit(ONLY, ",", fixed = TRUE)[[1]])
  builders <- builders[intersect(wanted, names(builders))]
}
if (!length(builders)) stop("SV_ONLY selected no structural DIF scenario")

one <- function(seed, builder) {
  made <- tryCatch(builder(seed), error = function(e) e)
  if (inherits(made, "error"))
    return(list(status = if (inherits(made, "rasch_refusal"))
      "refused" else "error"))
  fit <- made$fit
  converged <- if (inherits(fit, "rasch_btl")) fit$converged else
    fit$est$converged
  if (!isTRUE(converged)) return(list(status = "nonconverged"))
  db <- tryCatch(suppressWarnings(dif_bootstrap(
    fit, made$dif, B = NBOOT, workers = 1L, seed = seed + 400000L)),
    error = function(e) e)
  if (inherits(db, "error"))
    return(list(status = if (inherits(db, "rasch_refusal"))
      "refused" else "error"))
  fam <- is.finite(made$dif$terms$p_adj)
  list(
    status = "ok",
    primary_fwer = any(made$dif$terms$p_adj[fam] < .05),
    bootstrap_fwer = any(db$terms$p_boot_adj < .05, na.rm = TRUE),
    bootstrap_marginal = mean(db$terms$p_boot < .05, na.rm = TRUE),
    B_used = db$B_used, B_nonconverged = db$B_nonconverged,
    B_errors = db$B_errors)
}

set.seed(9.01e7)
rows <- list()
for (j in seq_along(builders)) {
  scenario <- names(builders)[j]
  message(sprintf("[%d/%d] %s", j, length(builders), scenario))
  seeds <- sample.int(1e8, NREP)
  z <- parallel::mclapply(
    seeds, one, builder = builders[[j]], mc.cores = NCORE,
    mc.preschedule = TRUE)
  status <- vapply(z, `[[`, "", "status")
  ok <- status == "ok"
  zz <- z[ok]
  n <- length(zz)
  if (!n) stop("no analysed replicate in ", scenario)
  val <- function(name) vapply(zz, `[[`, 0, name)
  primary <- val("primary_fwer")
  boot <- val("bootstrap_fwer")
  marginal <- val("bootstrap_marginal")
  used <- val("B_used")
  bnc <- val("B_nonconverged")
  berr <- val("B_errors")
  common <- list(
    study = "dif-bootstrap-structural-null", scenario = scenario,
    n_reps = n, n_attempted = NREP,
    n_refused = sum(status == "refused"),
    n_nonconv = sum(status == "nonconverged"),
    n_error = sum(status == "error"),
    n_boot_attempted = n * NBOOT, n_boot_used = sum(used),
    n_boot_nonconv = sum(bnc), n_boot_errors = sum(berr),
    notes = paste0(NBOOT, " fitted-null refits per analysed dataset; ",
      "performance conditional on the observed analysis and usable-refit rule"))
  mk <- function(quantity, ...) do.call(
    sv_row, c(common, list(quantity = quantity), list(...)))
  rows[[j]] <- rbind(
    mk("primary adjusted FWER", familywise = mean(primary)),
    mk("bootstrap minimum-p FWER", familywise = mean(boot)),
    mk("paired bootstrap-minus-primary FWER", bias = mean(boot - primary),
       mean_se = stats::sd(boot - primary) / sqrt(n), se_ratio = NA_real_),
    mk("bootstrap marginal term Type I", type1 = mean(marginal),
       mc_override = list(type1 = stats::sd(marginal) / sqrt(n))),
    mk("mean usable bootstrap proportion", bias = mean(used / NBOOT),
       emp_sd = stats::sd(used / NBOOT)))
  sv_write(do.call(rbind, rows), "dif-bootstrap-structural-null")
}
