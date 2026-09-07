# Rerun the historical booklet null and the unequal-exposure regression
# design under the requested-interval bootstrap rule. No historical outputs
# are overwritten. Run from the package root with Rscript.
suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "100"))
CORES <- min(as.integer(Sys.getenv("SV_CORES", "4")), .rasch_available_workers())
stopifnot(is.finite(NREP), NREP >= 1L, is.finite(CORES), CORES >= 1L)
STUDY <- "item-fit-bootstrap-intervals"
OUTPUT <- if (NREP < 100L) paste0(STUDY, "-smoke") else STUDY
STATS <- c("chisq", "fit_resid", "infit_z", "outfit_z")
SCENARIOS <- c("linked booklets", "unequal exposure")
HARNESS_HASH <- unname(tools::md5sum("tools/simval/harness.R"))
R_HASHES <- tools::md5sum(sort(list.files("R", "[.]R$", full.names = TRUE)))

generate <- function(r, scenario) {
  n <- 400L
  if (scenario == "linked booklets") {
    # Identical generator and seeds to booklet_one() in item-fit-bootstrap.R.
    set.seed(5e6 + r)
    grp <- rep(c("low", "high"), each = n / 2)
    th <- rnorm(n, ifelse(grp == "low", -.75, .75), 1)
    X <- sapply(seq(-2, 2, length.out = 15),
                 function(dd) rbinom(n, 1, plogis(th - dd)))
    colnames(X) <- sprintf("I%02d", 1:15)
    X[grp == "low", 11:15] <- NA
    X[grp == "high", 1:5] <- NA
  } else {
    set.seed(743L + r)
    X <- matrix(rbinom(n * 8L, 1,
      plogis(outer(rnorm(n), seq(-1, 1, length.out = 8), "-"))), n, 8,
      dimnames = list(NULL, paste0("I", 1:8)))
    X[121:400, 7:8] <- NA_integer_
  }
  X
}

one <- function(r, scenario) {
  B <- if (scenario == "linked booklets") 600L else 399L
  seed <- if (scenario == "linked booklets") 6e6 + r else 9e6 + r
  out <- data.frame(scenario = scenario, replicate = r, B = B,
    seed = seed, status = "error", message = "", warnings = "",
    n_boot_attempted = 0L, n_boot_used = 0L, n_boot_nonconv = 0L,
    n_boot_errors = 0L, same_data_max_error = NA_real_,
    imposed_count_max_error = NA_real_, min_intervals = NA_integer_,
    max_intervals = NA_integer_, total = NA_real_)
  for (st in STATS) {
    out[[paste0(st, "_marginal")]] <- NA_real_
    out[[paste0(st, "_familywise")]] <- NA_real_
  }
  item_rows <- NULL
  warnings <- character(0)
  tryCatch(withCallingHandlers({
    X <- generate(r, scenario)
    fit <- rasch(X)
    if (!isTRUE(fit$est$converged)) {
      out$status <- "nonconverged"
    } else {
      stopifnot(is.null(fit$refit_spec$n_groups))
      same <- .fit_refit(X, fit$model, .refit_n_groups(fit), NULL, fit$m, 60, 1e-8)
      old <- .fit_refit(X, fit$model, fit$n_groups, NULL, fit$m, 60, 1e-8)
      out$same_data_max_error <- max(abs(same$chisq - fit$items$chisq))
      out$imposed_count_max_error <- max(abs(old$chisq - fit$items$chisq))
      stopifnot(out$same_data_max_error < 1e-8)
      intervals <- vapply(fit$ci_item,
        function(x) length(unique(x[!is.na(x)])), 0L)
      out$min_intervals <- min(intervals)
      out$max_intervals <- max(intervals)
      out$n_boot_attempted <- B
      bs <- fit_bootstrap(fit, B = B, seed = seed, workers = 1L)
      stopifnot(identical(bs$algorithm, "loo-maxt-2"))
      out$n_boot_used <- bs$B_used
      out$n_boot_nonconv <- bs$B_nonconverged
      out$n_boot_errors <- bs$B_errors
      out$total <- as.numeric(bs$total$chisq_p_boot < .05)
      item_rows <- data.frame(scenario = scenario, replicate = r,
        item = fit$items$item, intervals = intervals, bs$items[, -1L],
        check.names = FALSE)
      for (st in STATS) {
        p <- bs$items[[paste0(st, "_p_boot")]]
        adj <- bs$items[[paste0(st, "_p_boot_adj")]]
        # Do not turn an untestable family member into a non-rejection.
        if (all(is.finite(p)))
          out[[paste0(st, "_marginal")]] <- mean(p < .05)
        if (all(is.finite(adj)))
          out[[paste0(st, "_familywise")]] <- as.numeric(any(adj < .05))
      }
      out$status <- "ok"
    }
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    out$status <<- if (inherits(e, "rasch_refusal")) "refused" else "error"
    out$message <<- conditionMessage(e)
    fields <- c(n_boot_attempted = "B", n_boot_used = "B_used",
                n_boot_nonconv = "B_nonconverged", n_boot_errors = "B_errors")
    for (target in names(fields)) {
      # Failure conditions from the bootstrap retain its complete accounting.
      nm <- fields[[target]]
      if (length(e[[nm]]) == 1L && is.finite(e[[nm]])) out[[target]] <<- e[[nm]]
    }
  })
  out$warnings <- paste(unique(warnings), collapse = " | ")
  list(attempt = out, items = item_rows)
}

stamp <- function(d) {
  d$script <- .sv_prov$script
  d$script_md5 <- .sv_prov$hash
  d$package_sha <- .sv_prov$sha
  d$r_tree_md5 <- .sv_prov$rtree
  d$executed <- .sv_prov$date
  d$harness_md5 <- HARNESS_HASH
  d$algorithm <- "loo-maxt-2"
  d
}

summarise <- function(attempts) {
  rows <- list()
  for (scenario in SCENARIOS) {
    d <- attempts[attempts$scenario == scenario, ]
    if (!nrow(d)) next
    for (metric in c("total", as.vector(outer(STATS,
                            c("marginal", "familywise"), paste, sep = "_")))) {
      ok <- d$status == "ok" & is.finite(d[[metric]])
      values <- d[[metric]][ok]
      estimate <- if (length(values)) mean(values) else NA_real_
      mc <- if (length(values) > 1L) sd(values) / sqrt(length(values)) else NA_real_
      family <- grepl("familywise$", metric)
      args <- list(study = STUDY, scenario = scenario, quantity = metric,
        n_reps = sum(ok), n_attempted = nrow(d),
        n_refused = sum(d$status == "refused"),
        n_nonconv = sum(d$status == "nonconverged"),
        n_error = sum(d$status == "error"), n_withheld = 0L,
        n_metric_unavailable = sum(d$status == "ok" & !is.finite(d[[metric]])),
        n_boot_attempted = sum(d$n_boot_attempted), n_boot_used = sum(d$n_boot_used),
        n_boot_nonconv = sum(d$n_boot_nonconv), n_boot_errors = sum(d$n_boot_errors),
        notes = paste0("B = ", d$B[1], "; automatic per-item intervals; ",
          "rates conditional on successful fits with the complete metric family; ",
          "MCSE from independent dataset-level values; original booklet seeds retained"))
      args[[if (family) "familywise" else "type1"]] <- estimate
      args$mc_override <- setNames(list(mc), if (family) "familywise" else "type1")
      row <- do.call(sv_row, args)
      row$harness_md5 <- HARNESS_HASH
      row$algorithm <- "loo-maxt-2"
      rows[[length(rows) + 1L]] <- row
    }
  }
  do.call(rbind, rows)
}

attempts <- items <- list()
for (scenario in SCENARIOS) {
  for (batch in split(seq_len(NREP), ceiling(seq_len(NREP) / 8L))) {
    ans <- parallel::mclapply(batch, one, scenario = scenario,
                              mc.cores = CORES, mc.set.seed = FALSE)
    stopifnot(all(vapply(ans, is.list, logical(1))))
    attempts <- c(attempts, lapply(ans, `[[`, "attempt"))
    items <- c(items, lapply(ans, `[[`, "items"))
    d <- do.call(rbind, attempts)
    sv_write(stamp(d), paste0(OUTPUT, "-attempts"))
    it <- do.call(rbind, items)
    if (!is.null(it)) sv_write(stamp(it), paste0(OUTPUT, "-items"))
    sv_write(summarise(d), OUTPUT)
    cat(sprintf("%s: %d/%d datasets; %d total analysed, %d refused, %d errors\n",
      scenario, max(batch), NREP, sum(d$status == "ok"),
      sum(d$status == "refused"), sum(d$status == "error")))
    flush.console()
  }
}
stopifnot(identical(unname(tools::md5sum(.sv_prov$script)), .sv_prov$hash),
          identical(unname(tools::md5sum("tools/simval/harness.R")), HARNESS_HASH),
          identical(tools::md5sum(names(R_HASHES)), R_HASHES))
cat("Completed with matching study, harness and R-source hashes.\n")
