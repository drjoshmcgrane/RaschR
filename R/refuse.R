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
# (moderate) or C (large) from the Mantel-Haenszel common odds ratio, and
# the published scheme is stated directly on the log-odds scale: A below
# 0.43, C at or above 0.64 (Zieky 1993; Penfield 2007; Golia 2012). Those
# cut-values are the delta-metric thresholds of 1.0 and 1.5 divided by the
# 2.35 delta units that make up a logit.
#
# Linacre and Wright (1989) showed the Rasch and Mantel-Haenszel approaches
# rest on the same relative odds, so under the Rasch model log(alpha_MH)
# equals the difference in item location and the scheme applies to that
# difference directly. Golia (2012, section 2.2) carries this to the partial
# credit model: an item's difficulty decomposes as delta_G + tau_j, so where
# differential functioning shifts the location and leaves the thresholds
# alone -- uniform DIF, which is what a resolved location difference
# estimates -- the signed area is J times that shift, and the same 0.43 and
# 0.64 apply per threshold. So the categories are available for dichotomous
# and polytomous items alike, on one metric.
#
#   A   not significant, or below 0.43 logits
#   C   at or above 0.64 logits AND significantly beyond 0.43
#   B   everything else
#
# The sign follows ETS practice: positive where the second level is
# favoured. A separate ETS convention categorises polytomous items on a
# standardised mean difference at 0.17 and 0.25, but that is an
# observed-score statistic rather than this one, and Zwick, Thayer and
# Mazzeo (1997) record that ETS had no official polytomous policy.
# ---------------------------------------------------------------------------
ETS_DELTA_PER_LOGIT <- 2.35

.ets_category <- function(difference, se, p_adj, alpha = 0.05) {
  a_cut <- 1.0 / ETS_DELTA_PER_LOGIT      # 0.43 as published
  c_cut <- 1.5 / ETS_DELTA_PER_LOGIT      # 0.64 as published
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
