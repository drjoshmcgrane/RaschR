# STUDY: btl-equating-clustered
#
# Null familywise calibration of common-object drift tests when each BTL
# calibration uses a judge-clustered sandwich covariance. The two panels are
# independent; judge-by-object deviations induce dependence among each judge's
# comparisons. The current contrast-specific Welch-Satterthwaite reference is
# compared with the superseded normal reference on the same fitted statistics.
#
# Serial. Rscript tools/simval/studies/btl-equating-clustered.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "btl-equating-clustered"
R <- 1000L

common <- sprintf("O%02d", 1:5)
truth <- setNames(seq(-1.8, 1.8, length.out = 15),
                  c(common, sprintf("A%02d", 1:5), sprintf("B%02d", 1:5)))

panel_fit <- function(objects, prefix, seed) {
  set.seed(seed)
  judges <- sprintf("%sJ%02d", prefix, 1:12)
  pair <- t(utils::combn(objects, 2))
  d <- data.frame(object_a = rep(pair[, 1], each = 25),
                  object_b = rep(pair[, 2], each = 25),
                  stringsAsFactors = FALSE)
  d$judge <- rep(judges, length.out = nrow(d))
  d <- d[sample.int(nrow(d)), , drop = FALSE]
  dev <- matrix(stats::rnorm(length(judges) * length(objects), 0, 0.6),
                nrow = length(judges), dimnames = list(judges, objects))
  lp <- truth[d$object_a] - truth[d$object_b] +
    dev[cbind(d$judge, d$object_a)] - dev[cbind(d$judge, d$object_b)]
  d$winner <- ifelse(stats::runif(nrow(d)) < stats::plogis(lp),
                     d$object_a, d$object_b)
  btl(d, "object_a", "object_b", "winner", judge = "judge")
}

fw_t <- fw_normal <- rep(NA, R); n_ref <- n_nc <- 0L
for (r in seq_len(R)) {
  f1 <- tryCatch(panel_fit(c(common, sprintf("A%02d", 1:5)), "A",
                           810000 + 2L * r), error = function(e) NULL)
  f2 <- tryCatch(panel_fit(c(common, sprintf("B%02d", 1:5)), "B",
                           810001 + 2L * r), error = function(e) NULL)
  if (is.null(f1) || is.null(f2)) { n_ref <- n_ref + 1L; next }
  if (!isTRUE(f1$converged) || !isTRUE(f2$converged)) { n_nc <- n_nc + 1L; next }
  eq <- btl_equate(f1, f2, independent = TRUE)
  fw_t[r] <- any(eq$table$p_adj < 0.05, na.rm = TRUE)
  pn <- 2 * stats::pnorm(-abs(eq$table$t))
  fw_normal[r] <- any(stats::p.adjust(pn, "holm") < 0.05, na.rm = TRUE)
}
ok <- !is.na(fw_t)
rows <- rbind(
  sv_row(STUDY, "12 judges per independent panel; judge-object SD 0.6",
         "familywise null rejection: Welch-Satterthwaite + Holm", sum(ok),
         familywise = mean(fw_t[ok]), n_attempted = R,
         n_refused = n_ref, n_nonconv = n_nc,
         notes = "five common objects; 25 comparisons per object pair"),
  sv_row(STUDY, "same fitted draws (counterfactual)",
         "familywise null rejection: superseded normal + Holm", sum(ok),
         familywise = mean(fw_normal[ok]), n_attempted = R,
         n_refused = n_ref, n_nonconv = n_nc,
         notes = "same t statistics and clustered covariances; only the reference distribution changes"))
sv_write(rows, STUDY)
