library(testthat)
library(rasch)

# CRAN's incoming pretest enforces a ten-minute overall budget, which the
# complete suite exceeds on the Windows builder. On CRAN a core subset runs,
# exercising every estimator and inference path once; the complete suite runs
# whenever NOT_CRAN is true, locally and in continuous integration on three
# operating systems.
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
  test_check("rasch")
} else {
  core <- c(
    "app-project", "btl-efrm", "btl-equating", "btl-targeting", "compare-ic",
    "dif-bootstrap", "dif-contrasts", "fit-residual", "format",
    "identification", "item-disc",
    "mc-scoring", "missing", "pcml-pc", "plots", "recovery", "resolve-frames",
    "statistical-validity", "wright-map"
  )
  test_check("rasch", filter = paste0("^(", paste(core, collapse = "|"), ")$"))
}
