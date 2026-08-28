# rasch :: telling a refusal apart from a fault
# ===========================================================================
# Some analyses cannot be run on some designs, and saying so is a result
# rather than a failure. Residual principal components are undefined when
# item columns share no respondents, because a correlation between them does
# not exist to be decomposed; a factor that defines a frame cannot also be
# tested for differential functioning within that frame, because it takes one
# level among the persons who respond there. Both are conclusions the package
# draws deliberately, and both are worth reading.
#
# A caller cannot tell either apart from a genuine fault when both arrive as
# a bare error. The application showed these considered refusals in the same
# red as a crash, which reads as the package having broken. Signalling them
# with their own condition class lets a caller answer each in the right
# voice, and leaves anything unclassed to be treated as the fault it is.
# ===========================================================================

# Signal a deliberate refusal. Inherits from "error", so an unprepared caller
# still stops and every existing expectation on the message still holds.
.refuse <- function(...) {
  stop(structure(
    class = c("rasch_refusal", "error", "condition"),
    list(message = paste0(...), call = NULL)))
}

# Dichotomous ETS classification on the log-odds scale. `p` and `p_beyond`
# are the probabilities used for decisions; callers that report a family of
# comparisons pass their familywise-adjusted values here.
ETS_DELTA_PER_LOGIT <- 2.35

.ets_p_beyond <- function(difference, se) {
  a_cut <- 1.0 / ETS_DELTA_PER_LOGIT      # 0.43 as published
  d <- abs(difference)
  # C additionally requires the magnitude to be significantly beyond the A
  # ceiling. Shervish's interval-null p-value retains both normal tails.
  out <- stats::pnorm((-a_cut - d) / se) +
    stats::pnorm((a_cut - d) / se)
  out[is.na(difference) | is.na(se) | se <= 0] <- NA_real_
  out
}

.ets_category <- function(difference, se, p, alpha = 0.05,
                          p_beyond = NULL) {
  a_cut <- 1.0 / ETS_DELTA_PER_LOGIT      # 0.43 as published
  c_cut <- 1.5 / ETS_DELTA_PER_LOGIT      # 0.64 as published
  d <- abs(difference)
  sig <- !is.na(p) & p < alpha
  if (is.null(p_beyond)) p_beyond <- .ets_p_beyond(difference, se)
  beyond <- !is.na(se) & se > 0 & !is.na(p_beyond) & p_beyond < alpha
  out <- ifelse(!sig | d <= a_cut, "A",
                ifelse(d >= c_cut & beyond, "C", "B"))
  out[is.na(difference) | is.na(se)] <- NA_character_
  sign_c <- ifelse(is.na(difference) | out == "A", "",
                   ifelse(difference > 0, "+", "-"))
  ifelse(is.na(out), NA_character_, paste0(out, sign_c))
}
