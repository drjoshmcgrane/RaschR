# STUDY: tailored-bootstrap-topup
#
# Power of the complete tailored-analysis bootstrap at two guessing effects.
# Automatic anchor selection is repeated inside every person-bootstrap draw.
# The eight-item design uses 399 draws, for which the smallest attainable
# Holm-adjusted probability is 2 * 8 / 400 = 0.04. Run from the package root;
# set SV_CORES on Unix-like systems to parallelise over generated samples.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "tailored-bootstrap-topup"
R <- 80L; N <- 300L; I <- 8L; B <- 399L
difficulty <- seq(-2, 2.5, length.out = I)
planted <- sprintf("I%02d", c(I - 1L, I))

one <- function(r, guessing, seed0) {
  gv <- rep(0, I); gv[(I - 1L):I] <- guessing
  d <- simulate_rasch(N, I, difficulty = difficulty, guessing = gv,
                      theta_mean = 0, theta_sd = 1, seed = seed0 + r)
  f <- tryCatch(rasch(d, id = "id"), error = function(e) NULL)
  if (is.null(f))
    return(c(detection = NA, any_detection = NA, clean_any = NA,
             refused = 1, nonconv = 0))
  if (!isTRUE(f$est$converged))
    return(c(detection = NA, any_detection = NA, clean_any = NA,
             refused = 0, nonconv = 1))
  z <- tryCatch(suppressWarnings(tailored_analysis(
    f, chance = 0.25, se_method = "bootstrap", boot_reps = B,
    seed = seed0 + 1000000L + r
  )), error = function(e) NULL)
  if (is.null(z))
    return(c(detection = NA, any_detection = NA, clean_any = NA,
             refused = 1, nonconv = 0))
  sig <- z$table$significant %in% TRUE
  is_planted <- z$table$item %in% planted
  c(detection = mean(sig[is_planted]), any_detection = any(sig[is_planted]),
    clean_any = any(sig[!is_planted]), refused = 0, nonconv = 0)
}

run_cell <- function(guessing, seed0) {
  cores <- suppressWarnings(as.integer(Sys.getenv("SV_CORES", "1")))
  if (!is.finite(cores) || cores < 1L) cores <- 1L
  z <- if (cores > 1L && .Platform$OS.type != "windows")
    parallel::mclapply(seq_len(R), one, guessing = guessing, seed0 = seed0,
                       mc.cores = cores, mc.set.seed = FALSE)
  else lapply(seq_len(R), one, guessing = guessing, seed0 = seed0)
  z <- do.call(rbind, z)
  ok <- is.finite(z[, "detection"])
  scenario <- sprintf(
    "N=%d; I=%d; two hardest items guess=%.2f; chance=.25; B=%d",
    N, I, guessing, B
  )
  rbind(
    sv_row(STUDY, scenario,
      "mean planted-item detection proportion (Holm adjusted)", sum(ok),
      power = mean(z[ok, "detection"]), effect = guessing,
      mc_override = list(
        power = sd(z[ok, "detection"]) / sqrt(sum(ok))
      ), n_attempted = R, n_refused = sum(z[, "refused"]),
      n_nonconv = sum(z[, "nonconv"]),
      notes = "Monte Carlo SE uses per-replicate planted-item proportions"),
    sv_row(STUDY, scenario,
      "probability of detecting at least one planted item", sum(ok),
      power = mean(z[ok, "any_detection"]), effect = guessing,
      n_attempted = R, n_refused = sum(z[, "refused"]),
      n_nonconv = sum(z[, "nonconv"])),
    sv_row(STUDY, scenario,
      "familywise false flag among the six unplanted items", sum(ok),
      familywise = mean(z[ok, "clean_any"]), effect = 0,
      n_attempted = R, n_refused = sum(z[, "refused"]),
      n_nonconv = sum(z[, "nonconv"]))
  )
}

sv_write(rbind(run_cell(0.15, 410000L), run_cell(0.30, 420000L)), STUDY)
