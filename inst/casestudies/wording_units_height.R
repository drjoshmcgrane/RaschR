# Case study: a real wording effect that leaves the person ordering alone
# ===========================================================================
# The companion case study, wording_units_selfesteem.R, finds a set-unit
# difference between positively and negatively worded items that is heavily
# attenuated when one anomalous item, Q8, is removed. This study asks whether
# the same sensitivity appears in a deliberately balanced wording design.
#
# The Height Inventory (Recka, 2018) asks 26 questions about self-perceived
# height, 13 worded in the tall direction ("I often have to ask others to
# reach something for me" is scored the other way) and 13 in the short
# direction and reverse-coded. The wording halves are balanced by design
# rather than by accident, the construct is unambiguous, and the sample is
# large. Under the extended frame of reference model the wording direction
# defines two item sets, each with its own unit: rho = alpha_set.
#
# Two questions follow, and they have different answers. Is the unit
# difference a property of the wording or of one badly behaved item? And
# does it change what we conclude about a person? The inventory carries a
# self-reported height in centimetres, collected on a different response
# format from the 26 rating items, which gives a criterion to check the
# second question against rather than arguing it.
#
# Run time: several minutes (4,499 respondents, 26 items, 27 model fits).
# ===========================================================================
library(rasch)
set.seed(25)

if (!requireNamespace("ShinyItemAnalysis", quietly = TRUE))
  stop("this case study needs the 'ShinyItemAnalysis' package for its data: ",
       "install.packages('ShinyItemAnalysis')")

H <- ShinyItemAnalysis::HeightInventory
tall <- names(H)[1:13]                  # worded "I am tall"
short <- names(H)[14:26]                # worded "I am short", reverse-coded
items <- c(tall, short)
ok <- complete.cases(H[, items]) & is.finite(H$HeightCM)
df <- data.frame(H[ok, items], gender = H$gender[ok], check.names = FALSE)
cm <- H$HeightCM[ok]
cat(sprintf("%d respondents with a complete inventory and a reported height\n",
            nrow(df)))

# equal-unit Rasch versus wording-set EFRM -----------------------------------
f0 <- rasch(df, items = items, factors = "gender")
set.seed(26)
f1 <- rasch_efrm(df, items = items, groups = rep("all", nrow(df)),
                 item_sets = list(tall = tall, short = short),
                 se_method = "hybrid", boot_reps = 300)
print(f1$alpha_table, digits = 3)
print(f1$efrm_vs_rasch$unit_tests, digits = 3)

# is the difference carried by one item? -------------------------------------
# In the self-esteem data the worst-fitting item departs three times as far
# as the next, and removing it removes the effect. Rank the items the same
# way here and refit without each in turn: if the direction belongs to the
# wording rather than to an item, no single removal should reverse it.
fr <- f1$items[order(abs(f1$items$fit_resid), decreasing = TRUE), ]
fr$item <- sub(":.*$", "", fr$item)
cat(sprintf("\nlargest fit residual %.1f (%s), next %.1f (%s), ratio %.2f\n",
            abs(fr$fit_resid[1]), fr$item[1], abs(fr$fit_resid[2]),
            fr$item[2], abs(fr$fit_resid[1] / fr$fit_resid[2])))

ratio <- function(f) {
  a <- f$alpha_table
  unname(a$alpha[a$set == "tall"] / a$alpha[a$set == "short"])
}
each <- vapply(items, function(it)
  tryCatch(ratio(drop_items(f1, it, boot_reps = 0)),
           error = function(e) NA_real_), numeric(1))
cat(sprintf("unit ratio with all items: %.3f\n", ratio(f1)))
cat(sprintf("across the %d single-item removals: %.3f to %.3f, %d reverse it\n",
            length(each), min(each, na.rm = TRUE), max(each, na.rm = TRUE),
            sum(each < 1, na.rm = TRUE)))

# cross-check against free slopes --------------------------------------------
# One unit per set is an average over that set. A generalized partial credit
# model frees a slope per item on the same respondents, so its slopes say
# whether the average conceals as much as it summarises.
if (requireNamespace("mirt", quietly = TRUE)) {
  m <- mirt::mirt(as.data.frame(df[, items]), 1, itemtype = "gpcm",
                  verbose = FALSE)
  a1 <- mirt::coef(m, simplify = TRUE)$items[, "a1"]
  gm <- function(z) exp(mean(log(z)))
  cat(sprintf(paste("\nfree slopes: tall %.3f (%.2f to %.2f),",
                    "short %.3f (%.2f to %.2f)\n"),
              gm(a1[tall]), min(a1[tall]), max(a1[tall]),
              gm(a1[short]), min(a1[short]), max(a1[short])))
  cat(sprintf("GPCM slope ratio %.3f against the frame model's %.3f\n",
              gm(a1[tall]) / gm(a1[short]), ratio(f1)))
  cat(sprintf("of the 13 steepest, %d are worded in the tall direction\n",
              sum(names(sort(a1, decreasing = TRUE))[1:13] %in% tall)))
} else {
  message("install the 'mirt' package for the free-slope cross-check: ",
          "install.packages('mirt')")
}

# does it change what we conclude about a person? ----------------------------
# The frame model does not leave the raw score sufficient -- persons with the
# same score receive different measures -- so it reorders people. Whether that
# reordering is an improvement is a question for a criterion, not for the fit.
keep <- !f0$person$extreme & !f1$person$extreme
t0 <- f0$person$theta[keep]; t1 <- f1$person$theta[keep]
y <- cm[keep]
cat(sprintf("\nmax measure spread within one raw score: %.3f logits\n",
            max(tapply(t1, f0$person$raw[keep],
                       function(z) diff(range(z))), na.rm = TRUE)))
cat(sprintf("persons moving more than 0.1 logits: %.1f%%\n",
            100 * mean(abs(t1 - t0) > 0.1)))

ci <- function(a, b) {
  r <- cor.test(a, b)
  sprintf("%.4f [%.4f, %.4f]", r$estimate, r$conf.int[1], r$conf.int[2])
}
cat(sprintf("correlation with reported height, equal unit  : %s\n", ci(t0, y)))
cat(sprintf("correlation with reported height, wording unit: %s\n", ci(t1, y)))

# the two correlations share a variable, so compare them as dependent
# correlations (Williams's t; matches psych::r.test)
r12 <- cor(t0, y); r13 <- cor(t1, y); r23 <- cor(t0, t1); n <- length(y)
Rd <- (1 - r12^2 - r13^2 - r23^2) + 2 * r12 * r13 * r23
tt <- (r12 - r13) * sqrt((n - 1) * (1 + r23) /
        (2 * Rd * (n - 1) / (n - 3) + (r12 + r13)^2 / 4 * (1 - r23)^3))
cat(sprintf("difference between them: t = %.3f on %d df, p = %.3f\n",
            tt, n - 3, 2 * pt(-abs(tt), n - 3)))
cat(sprintf("correlation between the two sets of measures: %.5f\n", r23))

# what the pair of case studies shows ----------------------------------------
# Here the unit difference is distributed across the wording sets rather than
# concentrated in one item. The two largest fit residuals are nearly equal,
# every single-item removal leaves the ratio above one (1.19 to 1.32), and a
# free-slope model fitted to the same people gives the same ratio, 1.258.
#
# And it does not matter for measuring a person. The frame model spreads
# people who share a raw score by up to half a logit and moves most of the
# sample, yet the two sets of measures correlate above 0.999 and predict
# reported height equally well -- if anything the equal-unit measure is
# marginally the better predictor. A unit difference decisive at p < 1e-30
# is not thereby consequential for person measurement.
#
# The comparison with the self-esteem study is therefore one of stability,
# not merely statistical significance. There, removing Q8 reduces the ratio
# from about 1.32 to 1.08; here no single removal removes the difference.
# What the frame model supplies in both cases is an account of how the items
# behave. Read it as a claim about persons only after checking against
# something outside the model.
