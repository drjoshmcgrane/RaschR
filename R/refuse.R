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

# ---------------------------------------------------------------------------
# The ETS differential-functioning categories
# ===========================================================================
# ETS classifies an item's differential functioning as A (negligible), B
# (moderate) or C (large) from the Mantel-Haenszel common odds ratio on the
# delta scale, where MH D-DIF = -2.35 log(alpha_MH). Under the Rasch model
# that log odds ratio and the conditional difference in item location
# estimate the same quantity, both being conditional on the total score, so
# the delta thresholds convert exactly: 2.35 delta units to the logit.
#
#   A   not significant, or below 1.0 delta (0.426 logits)
#   C   at or above 1.5 delta (0.638 logits) AND significantly beyond 1.0
#   B   everything else
#
# The sign follows ETS practice: positive where the second level is
# favoured. The conversion is only available for a dichotomous item; the
# polytomous rule uses a standardised mean difference in the observed-score
# metric, which is a different statistic rather than a rescaling of this one.
# ---------------------------------------------------------------------------
ETS_DELTA_PER_LOGIT <- 2.35

.ets_category <- function(difference, se, p_adj, alpha = 0.05,
                          dichotomous = TRUE) {
  if (!isTRUE(dichotomous)) return(rep(NA_character_, length(difference)))
  a_cut <- 1.0 / ETS_DELTA_PER_LOGIT
  c_cut <- 1.5 / ETS_DELTA_PER_LOGIT
  d <- abs(difference)
  sig <- !is.na(p_adj) & p_adj < alpha
  # C additionally requires the magnitude to be significantly beyond the A
  # ceiling, which is a one-sided test against a non-zero null
  beyond <- !is.na(se) & se > 0 &
    (d - a_cut) / se > stats::qnorm(1 - alpha)
  out <- ifelse(!sig | d < a_cut, "A",
                ifelse(d >= c_cut & beyond, "C", "B"))
  out[is.na(difference) | is.na(se)] <- NA_character_
  sign_c <- ifelse(is.na(difference) | out == "A", "",
                   ifelse(difference > 0, "+", "-"))
  ifelse(is.na(out), NA_character_, paste0(out, sign_c))
}
