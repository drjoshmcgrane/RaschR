# STUDY: btl-dif-multicell-df
# Null calibration of a four-cell BTL-DIF interaction contrast. The resolved
# object fit is the same joint refit used by btl_dif(); fitting it directly
# avoids selecting contrasts through the preceding omnibus test.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

STUDY <- "btl-dif-multicell-df"
N_REP <- 500L
objects <- paste0("O", 1:4)
pairs <- t(utils::combn(objects, 2))
judges <- sprintf("J%02d", 1:48)
cells <- stats::setNames(rep(c("AX", "AY", "BX", "BY"), each = 12), judges)
beta <- stats::setNames(seq(-0.9, 0.9, length.out = length(objects)), objects)
copies <- paste0("O3 (", c("AX", "AY", "BX", "BY"), ")")
weights <- stats::setNames(c(0.5, -0.5, -0.5, 0.5), copies)

one_rep <- function(r) {
  set.seed(270000L + r)
  judge_object <- matrix(
    stats::rnorm(length(judges) * length(objects), 0, 0.45),
    length(judges), length(objects),
    dimnames = list(judges, objects))
  rows <- lapply(judges, function(j) {
    z <- data.frame(a = rep(pairs[, 1], each = 2L),
                    b = rep(pairs[, 2], each = 2L),
                    judge = j, stringsAsFactors = FALSE)
    eta <- beta[z$a] - beta[z$b] +
      judge_object[j, z$a] - judge_object[j, z$b]
    z$x <- stats::rbinom(nrow(z), 1L, stats::plogis(eta))
    z
  })
  d <- do.call(rbind, rows)
  cell <- unname(cells[d$judge])
  a <- ifelse(d$a == "O3", paste0("O3 (", cell, ")"), d$a)
  b <- ifelse(d$b == "O3", paste0("O3 (", cell, ")"), d$b)
  fit <- tryCatch(rasch:::.btl_graded(
    a, b, d$x, d$judge, rep(1, nrow(d)), c("0", "1"),
    maxit = 60L, tol = 1e-8, notes = character(0)),
    error = function(e) NULL)
  if (is.null(fit))
    return(c(estimate = NA, se = NA, df = NA, p = NA, p_old = NA,
             refused = 1, nonconv = 0))
  if (!isTRUE(fit$converged))
    return(c(estimate = NA, se = NA, df = NA, p = NA, p_old = NA,
             refused = 0, nonconv = 1))
  idx <- match(copies, fit$objects$object)
  if (anyNA(idx))
    return(c(estimate = NA, se = NA, df = NA, p = NA, p_old = NA,
             refused = 1, nonconv = 0))
  location <- stats::setNames(fit$objects$location[idx], copies)
  covariance <- fit$cov_beta[copies, copies, drop = FALSE]
  estimate <- sum(weights * location)
  se <- sqrt(drop(t(weights) %*% covariance %*% weights))
  effective <- vapply(c("AX", "AY", "BX", "BY"), function(cell) {
    use <- d$judge[cells[d$judge] == cell &
                     (d$a == "O3" | d$b == "O3")]
    workload <- table(use)
    sum(workload)^2 / sum(workload^2)
  }, 0)
  df <- rasch:::.btl_dif_contrast_df(effective)
  old_df <- sum(effective) - 2
  statistic <- estimate / se
  c(estimate = estimate, se = se, df = df,
    p = 2 * stats::pt(-abs(statistic), df),
    p_old = 2 * stats::pt(-abs(statistic), old_df),
    refused = 0, nonconv = 0)
}

draws <- do.call(rbind, lapply(seq_len(N_REP), one_rep))
ok <- is.finite(draws[, "p"])
n_ok <- sum(ok)
n_refused <- sum(draws[, "refused"])
n_nonconv <- sum(draws[, "nonconv"])
current_size <- mean(draws[ok, "p"] < 0.05)
pooled_size <- mean(draws[ok, "p_old"] < 0.05)
estimate <- draws[ok, "estimate"]
se <- draws[ok, "se"]
df <- draws[ok, "df"]

rows <- rbind(
  sv_row(
    STUDY, "balanced 2 x 2 judge cells; 12 judges per cell",
    "weakest-cell reference", n_ok,
    bias = mean(estimate), emp_sd = stats::sd(estimate), mean_se = mean(se),
    coverage95 = mean(abs(estimate) <= stats::qt(0.975, df) * se),
    type1 = current_size, n_attempted = N_REP,
    n_refused = n_refused, n_nonconv = n_nonconv,
    notes = "four-cell interaction; df = minimum effective judges minus one"),
  sv_row(
    STUDY, "balanced 2 x 2 judge cells; 12 judges per cell",
    "superseded pooled-count reference", n_ok,
    bias = mean(estimate), emp_sd = stats::sd(estimate), mean_se = mean(se),
    coverage95 = mean(abs(estimate) <= stats::qt(0.975, 46) * se),
    type1 = pooled_size, n_attempted = N_REP,
    n_refused = n_refused, n_nonconv = n_nonconv,
    notes = "diagnostic only; df = sum effective judges minus two"))

print(rows[, c("quantity", "n_reps", "type1", "mc_se_type1",
               "coverage95", "se_ratio")], row.names = FALSE)
sv_write(rows, "btl-dif-multicell-df")
