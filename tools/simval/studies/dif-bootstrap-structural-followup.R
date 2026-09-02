# Follow-up null calibration for the two elevated structural cells. This run
# records where an outer dataset was lost and includes attempted bootstrap
# refits from a bootstrap refusal in the accounting.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "150"))
NBOOT <- as.integer(Sys.getenv("SV_BOOT", "99"))
NCORE <- max(1L, as.integer(Sys.getenv("SV_CORES", "4")))
ONLY <- Sys.getenv("SV_ONLY", "")

balanced <- function(ids, levels = c("A", "B")) {
  z <- rep(levels, length.out = length(ids))
  stats::setNames(sample(z, length(z), replace = FALSE), ids)
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
  rasch_efrm(
    d, item_sets = tr$item_sets, groups = "group", id = "id",
    factors = "site", n_groups = 2, boot_reps = 0)
}

make_btl <- function(seed) {
  d <- simulate_btl(5, 24, reps_per_pair = 40, seed = seed)
  fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  ids <- unique(as.character(fit$comparisons$judge))
  set.seed(seed + 300000L)
  attr(fit, "followup_group") <- balanced(ids)
  fit
}

builders <- list(efrm = make_efrm, btl = make_btl)
if (nzchar(ONLY)) {
  wanted <- trimws(strsplit(ONLY, ",", fixed = TRUE)[[1]])
  builders <- builders[intersect(wanted, names(builders))]
}
if (!length(builders)) stop("SV_ONLY selected no follow-up scenario")

failure <- function(status, stage, condition = NULL) {
  list(
    status = status, stage = stage,
    reason = if (is.null(condition)) "" else conditionMessage(condition),
    B_attempted = condition$B %||% 0L,
    B_used = condition$B_used %||% 0L,
    B_nonconverged = condition$B_nonconverged %||% 0L,
    B_errors = condition$B_errors %||% 0L)
}

one <- function(seed, builder, kind) {
  fit <- tryCatch(builder(seed), error = identity)
  if (inherits(fit, "error"))
    return(failure(if (inherits(fit, "rasch_refusal")) "refused" else
      "error", "observed fit", fit))
  converged <- if (inherits(fit, "rasch_btl")) fit$converged else
    fit$est$converged
  if (!isTRUE(converged))
    return(failure("nonconverged", "observed fit"))

  dif <- tryCatch(if (identical(kind, "btl"))
    btl_dif(fit, attr(fit, "followup_group"), min_n = 10) else
      dif_anova(fit, n_groups = 2), error = identity)
  if (inherits(dif, "error"))
    return(failure(if (inherits(dif, "rasch_refusal")) "refused" else
      "error", "observed DIF", dif))

  db <- tryCatch(suppressWarnings(dif_bootstrap(
    fit, dif, B = NBOOT, workers = 1L, seed = seed + 400000L)),
    error = identity)
  if (inherits(db, "error"))
    return(failure(if (inherits(db, "rasch_refusal")) "refused" else
      "error", "DIF bootstrap", db))

  fam <- is.finite(dif$terms$p_adj)
  list(
    status = "ok", stage = "analysed", reason = "",
    primary_fwer = any(dif$terms$p_adj[fam] < .05),
    bootstrap_fwer = any(db$terms$p_boot_adj < .05, na.rm = TRUE),
    bootstrap_marginal = mean(db$terms$p_boot < .05, na.rm = TRUE),
    B_attempted = db$B, B_used = db$B_used,
    B_nonconverged = db$B_nonconverged, B_errors = db$B_errors)
}

set.seed(9.02e7)
rows <- list()
for (j in seq_along(builders)) {
  scenario <- names(builders)[j]
  message(sprintf("[%d/%d] %s", j, length(builders), scenario))
  seeds <- sample.int(1e8, NREP)
  z <- parallel::mclapply(
    seeds, one, builder = builders[[j]], kind = scenario,
    mc.cores = NCORE, mc.preschedule = TRUE)
  status <- vapply(z, `[[`, "", "status")
  ok <- status == "ok"
  zz <- z[ok]
  n <- length(zz)
  if (!n) stop("no analysed replicate in ", scenario)
  val <- function(name) vapply(zz, `[[`, 0, name)
  primary <- val("primary_fwer")
  boot <- val("bootstrap_fwer")
  marginal <- val("bootstrap_marginal")
  b_attempted <- vapply(z, `[[`, 0, "B_attempted")
  b_used <- vapply(z, `[[`, 0, "B_used")
  b_nonconv <- vapply(z, `[[`, 0, "B_nonconverged")
  b_errors <- vapply(z, `[[`, 0, "B_errors")
  common <- list(
    study = "dif-bootstrap-structural-followup", scenario = scenario,
    n_reps = n, n_attempted = NREP,
    n_refused = sum(status == "refused"),
    n_nonconv = sum(status == "nonconverged"),
    n_error = sum(status == "error"),
    n_boot_attempted = sum(b_attempted), n_boot_used = sum(b_used),
    n_boot_nonconv = sum(b_nonconv), n_boot_errors = sum(b_errors),
    notes = paste0(NBOOT, " fitted-null refits per analysable dataset; ",
      "performance conditional on the observed analysis and usable-refit rule"))
  mk <- function(quantity, ...) {
    args <- c(common, list(quantity = quantity))
    extra <- list(...)
    args[names(extra)] <- extra
    do.call(sv_row, args)
  }
  perf <- rbind(
    mk("primary adjusted FWER", familywise = mean(primary)),
    mk("bootstrap minimum-p FWER", familywise = mean(boot)),
    mk("paired bootstrap-minus-primary FWER", bias = mean(boot - primary),
       mean_se = stats::sd(boot - primary) / sqrt(n), se_ratio = NA_real_),
    mk("bootstrap marginal term Type I", type1 = mean(marginal),
       mc_override = list(type1 = stats::sd(marginal) / sqrt(n))),
    mk("mean usable bootstrap proportion", bias = mean(val("B_used") / NBOOT),
       emp_sd = stats::sd(val("B_used") / NBOOT)))
  stage <- vapply(z, `[[`, "", "stage")
  lost <- which(!ok)
  detail <- if (!length(lost)) NULL else do.call(rbind, lapply(
    split(lost, paste(stage[lost], status[lost], sep = ": ")), function(ii) {
      reasons <- unique(vapply(z[ii], `[[`, "", "reason"))
      mk(paste0("outer loss: ", stage[ii[1]], " (", status[ii[1]], ")"),
         effect = length(ii), bias = length(ii) / NREP,
         notes = paste(reasons[nzchar(reasons)], collapse = " | "))
    }))
  rows[[j]] <- rbind(perf, detail)
  sv_write(do.call(rbind, rows), "dif-bootstrap-structural-followup")
}
