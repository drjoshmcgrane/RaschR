# Cross-family conformance for dif_bootstrap(). This checks implementation,
# not operating characteristics: each supported structural family is generated
# and refitted through an independently orchestrated public fitting call, then
# compared with the exported bootstrap's retained F and p matrices.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

dif_source_md5 <- unname(tools::md5sum("R/dif-bootstrap.R"))

B <- as.integer(Sys.getenv("SV_BOOT", "19"))
workers <- max(1L, as.integer(Sys.getenv("SV_CORES", "1")))

align_dif <- function(dif, observed, btl = FALSE) {
  family <- if (btl) dif$term_ids != "band" else
    !dif$term_ids %in% c("Residuals", "ci")
  id <- if (btl) "object" else "item"
  key <- .factor_keys(data.frame(
    unit = as.character(dif$terms[[id]][family]),
    term = dif$term_ids[family], stringsAsFactors = FALSE))
  at <- match(observed, key)
  if (length(key) != length(observed) || anyNA(at) || anyDuplicated(key) ||
      any(!is.finite(dif$terms$F_value[family][at])) ||
      any(!is.finite(dif$terms$p[family][at]))) return(NULL)
  list(F = dif$terms$F_value[family][at], p = dif$terms$p[family][at])
}

rasch_design_dif <- function(fit, design) dif_anova(
  fit, factors = design$factors, n_groups = design$n_groups,
  p_adjust = design$p_adjust, alpha = design$alpha,
  effects = design$effects, sizes = FALSE, id = design$id,
  within = design$within, pool_facets = design$pool_facets)

make_cases <- function() {
  set.seed(9801)
  N <- 160L; group <- factor(rep(c("A", "B"), each = N / 2L))
  theta <- rnorm(N)
  X <- sapply(seq(-1, 1, length.out = 6L), function(d)
    rbinom(N, 1L, plogis(theta - d)))
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  q <- data.frame(item = colnames(X), feature = seq_len(ncol(X)))
  ex <- rasch_explanatory(X, q, ~ feature,
    factors = data.frame(group = group), n_groups = 2)

  dm <- simulate_mfrm(n_persons = 110, n_items = 4, n_raters = 3,
    n_categories = 3, seed = 9802)
  ids <- unique(dm$person)
  dm$group <- rep(c("A", "B"), length.out = length(ids))[
    match(dm$person, ids)]
  mf <- rasch_mfrm(dm, "person", "item", "score", facets = "rater",
    factors = "group", n_groups = 2)

  de <- simulate_efrm(n_per_group = 110, items_per_set = 5, n_sets = 2,
    n_groups = 2, seed = 9803)
  et <- attr(de, "truth")
  de$site <- rep(c("A", "B"), length.out = nrow(de))
  ef <- rasch_efrm(de, item_sets = et$item_sets, groups = "group", id = "id",
    factors = "site", boot_reps = 0, n_groups = 2)

  dc <- simulate_btl(n_objects = 5, n_judges = 20, reps_per_pair = 40,
    seed = 9804)
  bt <- btl(dc, "object_a", "object_b", winner = "winner", judge = "judge")
  judges <- unique(bt$comparisons$judge)
  jg <- setNames(rep(c("A", "B"), length.out = length(judges)), judges)

  list(
    explanatory = list(fit = ex, dif = dif_anova(ex, n_groups = 2)),
    mfrm = list(fit = mf, dif = dif_anova(mf, n_groups = 2)),
    efrm = list(fit = ef, dif = dif_anova(ef, n_groups = 2)),
    btl = list(fit = bt,
      dif = btl_dif(bt, jg, objects = "O3", min_n = 10)))
}

manual_explanatory <- function(fit, design) {
  xb <- .fit_gen_conditional(fit$X, fit$tau_list,
    if (anyNA(fit$X)) is.na(fit$X) else NULL)
  spec <- fit$refit_spec
  bf <- rasch_explanatory(xb,
    predictors = fit$explanatory$source_predictors,
    formula = fit$explanatory$formula, level = fit$explanatory$level,
    id = fit$person$id, factors = fit$factors,
    n_groups = spec$n_groups %||% fit$n_groups,
    maxit = spec$maxit, tol = spec$tol)
  list(dif = rasch_design_dif(bf, design),
       preservation = identical(rowSums(xb, na.rm = TRUE),
                                rowSums(fit$X, na.rm = TRUE)))
}

manual_mfrm <- function(fit, design) {
  xb <- .fit_gen_conditional(fit$X, fit$tau_list,
    if (anyNA(fit$X)) is.na(fit$X) else NULL)
  vm <- fit$virtual_map; N <- nrow(xb); L <- ncol(xb)
  d <- data.frame(.person = rep(fit$person$id, times = L),
    .item = rep(vm$item, each = N), .score = as.vector(xb),
    stringsAsFactors = FALSE, check.names = FALSE)
  for (f in fit$facet_spec) d[[f]] <- rep(vm[[f]], each = N)
  spec <- fit$refit_spec
  bf <- rasch_mfrm(d, ".person", ".item", ".score",
    facets = fit$facet_spec, interaction = fit$interaction,
    factors = fit$factors, n_groups = spec$n_groups %||% fit$n_groups,
    maxit = spec$maxit, tol = spec$tol)
  list(dif = rasch_design_dif(bf, design),
       preservation = identical(rowSums(xb, na.rm = TRUE),
                                rowSums(fit$X, na.rm = TRUE)))
}

manual_efrm <- function(fit, design) {
  vm <- fit$virtual_map
  xb <- matrix(NA_integer_, nrow(fit$X), ncol(fit$X),
               dimnames = dimnames(fit$X))
  block <- .factor_keys(vm[, c("set", "group"), drop = FALSE])
  preserved <- TRUE
  for (bb in unique(block)) {
    cols <- which(block == bb); r <- unique(fit$disc[cols])
    z <- .fit_gen_conditional(fit$X[, cols, drop = FALSE],
      lapply(fit$tau_list[cols], `*`, r),
      is.na(fit$X[, cols, drop = FALSE]))
    xb[, cols] <- z
    preserved <- preserved && identical(rowSums(z, na.rm = TRUE),
      rowSums(fit$X[, cols, drop = FALSE], na.rm = TRUE))
  }
  base <- matrix(NA_integer_, nrow(xb), length(fit$set_of),
                 dimnames = list(NULL, names(fit$set_of)))
  for (it in names(fit$set_of)) {
    cols <- which(vm$item == it); z <- xb[, cols, drop = FALSE]
    hit <- which(rowSums(!is.na(z)) == 1L)
    if (length(hit)) base[hit, it] <- z[cbind(
      hit, max.col(!is.na(z[hit, , drop = FALSE]), ties.method = "first"))]
  }
  spec <- fit$refit_spec
  bf <- rasch_efrm(base, item_sets = fit$set_of,
    groups = fit$factors[[fit$frame_group[1L]]], id = fit$person$id,
    n_groups = spec$n_groups %||% fit$n_groups,
    maxit = spec$maxit, tol = spec$tol,
    min_link_persons = spec$min_link_persons,
    se_method = "hybrid", boot_reps = 0, workers = 1)
  list(dif = rasch_design_dif(bf, design), preservation = preserved)
}

manual_btl <- function(fit, design) {
  d <- .btl_boot_data(fit); spec <- fit$refit_spec
  args <- list(data = d, object_a = "object_a", object_b = "object_b",
    response = "response", count = "count",
    thresholds = spec$thresholds, position = isTRUE(spec$position),
    anchors = spec$anchors, maxit = spec$maxit, tol = spec$tol)
  if (any(!is.na(d$judge))) args$judge <- "judge"
  if (isTRUE(spec$has_order)) args$order <- "order"
  bf <- do.call(btl, args)
  bd <- do.call(btl_dif, c(list(fit = bf), design))
  source_pairs <- .btl_pair_key(fit$comparisons$object_a,
                                fit$comparisons$object_b)
  draw_pairs <- .btl_pair_key(d$object_a, d$object_b)
  list(dif = bd, preservation =
    setequal(unique(source_pairs), unique(draw_pairs)) &&
      setequal(unique(fit$comparisons$judge), unique(d$judge)))
}

manual_run <- function(case, seed, kind) {
  fit <- case$fit; dif <- case$dif; design <- dif$bootstrap_design
  btl <- identical(kind, "btl")
  family <- if (btl) dif$term_ids != "band" else
    !dif$term_ids %in% c("Residuals", "ci")
  id <- if (btl) "object" else "item"
  observed <- .factor_keys(data.frame(
    unit = as.character(dif$terms[[id]][family]),
    term = dif$term_ids[family], stringsAsFactors = FALSE))
  set.seed(seed); seeds <- sample.int(.Machine$integer.max, B)
  ans <- vector("list", B); preserved <- logical(B)
  for (b in seq_len(B)) {
    set.seed(seeds[b])
    z <- tryCatch(switch(kind,
      explanatory = manual_explanatory(fit, design),
      mfrm = manual_mfrm(fit, design),
      efrm = manual_efrm(fit, design),
      btl = manual_btl(fit, design)), error = function(e) NULL)
    if (is.null(z)) next
    preserved[b] <- isTRUE(z$preservation)
    ans[[b]] <- align_dif(z$dif, observed, btl)
  }
  keep <- !vapply(ans, is.null, logical(1))
  list(F = do.call(rbind, lapply(ans[keep], `[[`, "F")),
       p = do.call(rbind, lapply(ans[keep], `[[`, "p")),
       used = sum(keep), preservation_failures = sum(!preserved[keep]))
}

cases <- make_cases(); rows <- list()
for (j in seq_along(cases)) {
  kind <- names(cases)[j]; case <- cases[[j]]; seed <- 9900L + j
  public <- suppressWarnings(dif_bootstrap(
    case$fit, case$dif, B = B, workers = workers, seed = seed))
  manual <- manual_run(case, seed, kind)
  family <- if (identical(kind, "btl")) case$dif$term_ids != "band" else
    !case$dif$term_ids %in% c("Residuals", "ci")
  observed <- case$dif$terms[family, , drop = FALSE]
  min_p <- apply(manual$p, 1L, min)
  current_adj <- vapply(observed$p, function(p)
    (1 + sum(min_p <= p)) / (nrow(manual$p) + 1), numeric(1))
  legacy_raw <- vapply(seq_len(ncol(manual$F)), function(k)
    (1 + sum(manual$F[, k] >= observed$F_value[k])) /
      (nrow(manual$F) + 1), numeric(1))
  legacy_adj <- pmax(legacy_raw, current_adj)
  vals <- c(
    max_abs_F_difference = if (identical(dim(public$replicates$F),
      dim(manual$F))) max(abs(public$replicates$F - manual$F)) else Inf,
    max_abs_reference_p_difference = if (identical(dim(public$replicates$p),
      dim(manual$p))) max(abs(public$replicates$p - manual$p)) else Inf,
    retained_replicate_difference = public$B_used - manual$used,
    sufficient_score_or_design_failures = manual$preservation_failures,
    legacy_to_reference_adjusted_max_change = max(abs(
      legacy_adj - current_adj)),
    legacy_to_reference_flag_changes = sum(
      (legacy_adj < case$dif$alpha) != (current_adj < case$dif$alpha)))
  rows[[j]] <- do.call(rbind, lapply(names(vals), function(q)
    sv_row(study = "dif-bootstrap-model-conformance", scenario = kind,
      quantity = q, n_reps = 1L, n_attempted = 1L,
      n_refused = 0L, n_nonconv = 0L, n_error = 0L,
      n_boot_attempted = B, n_boot_used = public$B_used,
      n_boot_nonconv = public$B_nonconverged,
      n_boot_errors = public$B_errors, bias = unname(vals[q]),
      notes = paste("public result versus independently orchestrated public",
                    "fit; exact family and generated-design check;",
                    "R/dif-bootstrap.R md5", dif_source_md5))))
}
sv_write(do.call(rbind, rows), "dif-bootstrap-model-conformance")
