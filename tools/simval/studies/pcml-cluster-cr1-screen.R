# Finite-person-cluster screen for ordinary PCML equating. The production
# covariance uses the unscaled person-cluster score meat and the drift tests
# use Welch--Satterthwaite t references. This seed-paired study asks whether
# multiplying each repeated-person covariance by the common CR1 factor
# G/(G-1) improves the release target: Holm familywise control. It does not
# alter the estimator.

pkgload::load_all(".", quiet = TRUE)
source("tools/simval/harness.R")

draw_pcm_item <- function(theta, tau) {
  score <- 0:length(tau)
  lp <- outer(theta, score, "*") -
    matrix(c(0, cumsum(tau)), nrow = length(theta),
           ncol = length(score), byrow = TRUE)
  lp <- lp - apply(lp, 1L, max)
  p <- exp(lp)
  p <- p / rowSums(p)
  u <- runif(length(theta))
  as.integer(rowSums(u > t(apply(p, 1L, cumsum))))
}

make_form <- function(G, repeats, L, kind) {
  theta <- rep(rnorm(G), each = repeats)
  d <- seq(-1, 1, length.out = L)
  if (identical(kind, "dichotomous")) {
    pr <- plogis(outer(theta, d, "-"))
    X <- matrix(rbinom(length(pr), 1L, pr), nrow = length(theta), ncol = L)
  } else {
    X <- vapply(d, function(di)
      draw_pcm_item(theta, di + c(-0.7, 0.7)), integer(length(theta)))
  }
  colnames(X) <- paste0("I", seq_len(L))
  id <- rep(sprintf("P%03d", seq_len(G)), each = repeats)
  rasch(X, id = id)
}

one_rep <- function(seed, G, repeats, L, kind) {
  set.seed(seed)
  fits <- tryCatch(list(make_form(G, repeats, L, kind),
                        make_form(G, repeats, L, kind)),
                   error = function(e) e)
  if (inherits(fits, "error")) return(list(status = "error"))
  if (any(!vapply(fits, function(f) isTRUE(f$est$converged), logical(1))))
    return(list(status = "nonconv"))
  if (any(!vapply(fits, function(f) isTRUE(f$est$cluster_inference),
                  logical(1))))
    return(list(status = "withheld"))
  eq <- tryCatch(equate_tests(fits[[1L]], fits[[2L]], independent = TRUE),
                 error = function(e) e)
  if (inherits(eq, "error")) return(list(status = "error"))
  if (!isTRUE(eq$inferential)) return(list(status = "withheld"))
  if (anyNA(eq$table$p)) return(list(status = "unavailable"))

  cr1 <- G / (G - 1)
  p_t <- eq$table$p
  p_cr1 <- 2 * pt(-abs(eq$table$t / sqrt(cr1)), df = eq$table$df)
  list(status = "analysed",
       value = c(t_size = mean(p_t < 0.05),
                 cr1_size = mean(p_cr1 < 0.05),
                 t_fwer = any(p.adjust(p_t, "holm", n = L) < 0.05),
                 cr1_fwer = any(p.adjust(p_cr1, "holm", n = L) < 0.05)))
}

run_cell <- function(G, kind, reps = 500L, repeats = 10L, L = 5L) {
  seed0 <- if (identical(kind, "dichotomous")) 910000L else 920000L
  ans <- lapply(seq_len(reps), function(r)
    one_rep(seed0 + 10000L * G + r, G, repeats, L, kind))
  status <- vapply(ans, `[[`, character(1L), "status")
  keep <- status == "analysed"
  z <- do.call(cbind, lapply(ans[keep], `[[`, "value"))
  n <- ncol(z)
  accounting <- list(n_attempted = reps, n_refused = 0L,
                     n_nonconv = sum(status == "nonconv"),
                     n_error = sum(status == "error"),
                     n_withheld = sum(status == "withheld"),
                     n_metric_unavailable = sum(status == "unavailable"))
  scenario <- sprintf("%s; G=%d; 10 rows/person; 5 items",
                      kind, G)
  note <- paste(
    "Two independent null calibrations; current Welch t compared on the same",
    "fits with covariance multiplied by G/(G-1); Holm family is all 5 items")
  make_row <- function(metric, method, familywise = FALSE) {
    rate <- mean(z[metric, ])
    override <- if (!familywise)
      list(type1 = stats::sd(z[metric, ]) / sqrt(n)) else list()
    do.call(sv_row, c(list(
      study = "pcml-cluster-cr1-screen", scenario = scenario,
      quantity = if (familywise)
        paste(method, "Holm familywise rejection") else
        paste(method, "pooled raw rejection"),
      n_reps = n,
      type1 = if (familywise) NA_real_ else rate,
      familywise = if (familywise) rate else NA_real_,
      mc_override = override, notes = note), accounting))
  }
  rbind(make_row("t_size", "Welch t"),
        make_row("cr1_size", "G/(G-1) plus Welch t"),
        make_row("t_fwer", "Welch t", TRUE),
        make_row("cr1_fwer", "G/(G-1) plus Welch t", TRUE))
}

cells <- c(lapply(c(10L, 12L, 20L, 50L), function(G)
  list(G = G, kind = "dichotomous")),
  list(list(G = 20L, kind = "PCM")))
cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (length(cores) != 1L || !is.finite(cores) || cores < 1L) cores <- 4L
rows <- parallel::mclapply(cells, function(cell)
  run_cell(cell$G, cell$kind), mc.cores = min(length(cells), cores))
rows <- do.call(rbind, rows)
sv_write(rows, "pcml-cluster-cr1-screen")
print(rows[, c("scenario", "quantity", "n_reps", "type1",
               "mc_se_type1", "familywise", "mc_se_familywise",
               "n_nonconv", "n_error", "n_withheld")], row.names = FALSE)
