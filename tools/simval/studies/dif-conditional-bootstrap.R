# Conditional null bootstrap for ordinary DIF. Each bootstrap preserves every
# person's sufficient total score, samples a response pattern under the fitted
# invariant Rasch calibration, refits the model, reconstructs its residual
# diagnostics, and recalculates both uniform and non-uniform terms.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "300"))
NBOOT <- as.integer(Sys.getenv("SV_BOOT", "199"))
NCORE <- max(1L, as.integer(Sys.getenv("SV_CORES", "1")))
L <- 8L

pattern_bank <- local({
  p <- as.matrix(expand.grid(rep(list(0:1), L)))
  split.data.frame(p, rowSums(p))
})

conditional_sample <- function(scores, delta) {
  X <- matrix(0L, length(scores), length(delta))
  for (s in sort(unique(scores))) {
    take <- which(scores == s)
    p <- pattern_bank[[as.character(s)]]
    if (nrow(p) == 1L) {
      X[take, ] <- rep(as.integer(p[1L, ]), each = length(take))
    } else {
      lw <- -drop(as.matrix(p) %*% delta)
      w <- exp(lw - max(lw))
      j <- sample.int(nrow(p), length(take), replace = TRUE, prob = w)
      X[take, ] <- as.matrix(p[j, , drop = FALSE])
    }
  }
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  X
}

dif_statistics <- function(fit) {
  group <- factor(fit$factors$group)
  ci <- .dif_class_intervals(fit, .dif_n_groups(fit, group))
  Fv <- pv <- matrix(NA_real_, L, 2L,
    dimnames = list(colnames(fit$residuals), c("group", "group:ci")))
  for (i in seq_len(L)) {
    d <- data.frame(z = fit$residuals[, i], group = group, ci = ci)
    d <- d[stats::complete.cases(d), , drop = FALSE]
    ft <- .dif_type2(d, c("group", "ci", "group:ci"),
                     variance = "hc3", robust_terms = "group")
    j <- match(colnames(Fv), ft$term)
    Fv[i, ] <- ft$F_value[j]
    pv[i, ] <- ft$p[j]
  }
  list(F = Fv, p = pv)
}

simulate_data <- function(scenario) {
  imbalanced <- scenario != "balanced null"
  n_a <- if (imbalanced) 120L else 300L
  n_b <- if (imbalanced) 480L else 300L
  group <- factor(c(rep("A", n_a), rep("B", n_b)))
  theta <- rnorm(n_a + n_b) +
    if (imbalanced) ifelse(group == "B", 0.8, 0) else 0
  delta <- seq(-1.4, 1.4, length.out = L)
  eta <- outer(theta, delta, "-")
  if (scenario == "uniform DIF 0.6")
    eta[group == "B", 3L] <- eta[group == "B", 3L] - 0.6
  if (scenario == "non-uniform DIF 0.7")
    eta[group == "B", 3L] <- eta[group == "B", 3L] +
      0.7 * (theta[group == "B"] - 0.8)
  X <- matrix(rbinom(length(eta), 1, plogis(eta)), nrow(eta), ncol(eta))
  colnames(X) <- paste0("I", seq_len(L))
  list(X = X, group = group)
}

one_replicate <- function(seed, scenario) {
  set.seed(seed)
  dat <- simulate_data(scenario)
  fit <- tryCatch(rasch(data.frame(dat$X, group = dat$group),
                        factors = "group"), error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$est$converged)) return(NULL)
  obs <- dif_statistics(fit)
  if (any(!is.finite(obs$F)) || any(!is.finite(obs$p))) return(NULL)
  scores <- rowSums(dat$X)
  delta <- stats::setNames(fit$items$location, fit$items$item)[colnames(dat$X)]
  fb <- pb <- matrix(NA_real_, NBOOT, 2L * L)
  got <- 0L; attempted <- 0L
  while (got < NBOOT && attempted < ceiling(1.25 * NBOOT)) {
    attempted <- attempted + 1L
    xb <- conditional_sample(scores, delta)
    bf <- tryCatch(rasch(data.frame(xb, group = dat$group),
                         factors = "group"), error = function(e) NULL)
    if (is.null(bf) || !isTRUE(bf$est$converged)) next
    bs <- tryCatch(dif_statistics(bf), error = function(e) NULL)
    if (is.null(bs) || any(!is.finite(bs$F)) || any(!is.finite(bs$p))) next
    got <- got + 1L
    fb[got, ] <- as.vector(bs$F)
    pb[got, ] <- as.vector(bs$p)
  }
  if (got < NBOOT) return(NULL)
  fo <- as.vector(obs$F); po <- as.vector(obs$p)
  p_boot_raw <- (1 + colSums(sweep(fb, 2L, fo, ">="))) / (NBOOT + 1)
  min_p <- apply(pb, 1L, min)
  p_boot_family <- (1 + vapply(po, function(p) sum(min_p <= p), 0L)) /
    (NBOOT + 1)
  p_default <- stats::p.adjust(po, "holm")
  iu <- 3L; inn <- L + 3L
  c(default_raw_uniform = mean(obs$p[, "group"] < 0.05),
    bootstrap_raw_uniform = mean(p_boot_raw[seq_len(L)] < 0.05),
    default_raw_nonuniform = mean(obs$p[, "group:ci"] < 0.05),
    bootstrap_raw_nonuniform = mean(p_boot_raw[L + seq_len(L)] < 0.05),
    default_fwer = any(p_default < 0.05),
    bootstrap_fwer = any(p_boot_family < 0.05),
    default_uniform_i3 = p_default[iu] < 0.05,
    bootstrap_uniform_i3 = p_boot_family[iu] < 0.05,
    default_nonuniform_i3 = p_default[inn] < 0.05,
    bootstrap_nonuniform_i3 = p_boot_family[inn] < 0.05)
}

run_condition <- function(scenario) {
  seeds <- sample.int(.Machine$integer.max, NREP)
  z <- parallel::mclapply(seeds, one_replicate, scenario = scenario,
                          mc.cores = NCORE, mc.preschedule = TRUE)
  refused <- sum(vapply(z, is.null, TRUE))
  z <- z[!vapply(z, is.null, TRUE)]
  if (!length(z)) stop("no replicate completed in ", scenario)
  z <- do.call(rbind, z); n <- nrow(z)
  is_null <- grepl("null$", scenario)
  effect <- if (scenario == "uniform DIF 0.6") 0.6 else
    if (scenario == "non-uniform DIF 0.7") 0.7 else 0
  mk <- function(quantity, col, field, note) {
    args <- list(study = "dif-conditional-bootstrap", scenario = scenario,
      quantity = quantity, n_reps = n, n_attempted = NREP,
      n_refused = refused, n_nonconv = 0L, effect = effect,
      notes = paste0("conditional Rasch null; ", NBOOT,
                     " bootstrap refits; ", note))
    args[[field]] <- mean(z[, col])
    if (field == "type1") args$mc_override <- list(
      type1 = stats::sd(z[, col]) / sqrt(n))
    do.call(sv_row, args)
  }
  if (is_null) rbind(
    mk("default uniform item-wise Type I", "default_raw_uniform", "type1",
       "HC3 uniform reference"),
    mk("bootstrap uniform item-wise Type I", "bootstrap_raw_uniform", "type1",
       "term-specific empirical tail"),
    mk("default non-uniform item-wise Type I", "default_raw_nonuniform", "type1",
       "residual-ANOVA reference"),
    mk("bootstrap non-uniform item-wise Type I", "bootstrap_raw_nonuniform", "type1",
       "term-specific empirical tail"),
    mk("default Holm FWER", "default_fwer", "familywise",
       "one family over all item-by-term tests"),
    mk("bootstrap minimum-p FWER", "bootstrap_fwer", "familywise",
       "single-step family reference over all item-by-term tests"))
  else if (scenario == "uniform DIF 0.6") rbind(
    mk("default uniform power", "default_uniform_i3", "power",
       "Holm-adjusted I3 uniform term"),
    mk("bootstrap uniform power", "bootstrap_uniform_i3", "power",
       "minimum-p adjusted I3 uniform term"))
  else rbind(
    mk("default non-uniform power", "default_nonuniform_i3", "power",
       "Holm-adjusted I3 non-uniform term"),
    mk("bootstrap non-uniform power", "bootstrap_nonuniform_i3", "power",
       "minimum-p adjusted I3 non-uniform term"))
}

set.seed(8.29e7)
rows <- do.call(rbind, lapply(
  c("balanced null", "imbalanced ability null", "uniform DIF 0.6",
    "non-uniform DIF 0.7"), run_condition))
sv_write(rows, "dif-conditional-bootstrap")
