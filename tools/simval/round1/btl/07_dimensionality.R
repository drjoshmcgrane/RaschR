suppressWarnings(pkgload::load_all(".", quiet=TRUE))

K <- 10; J <- 14; reps_pp <- 25

## --- (7c) shared-order verdict withheld: simulate_btl's own dependence-order
## layout puts every judge through (nearly) the same pair sequence (pair-major
## row layout, stable sort within judge) -- the package's own shared_order
## detector should catch this and withhold the verdict (a PASS, not an error).
n_shared_check <- 20
shared_withheld <- logical(n_shared_check)
for (i in seq_len(n_shared_check)) {
  d <- simulate_btl(n_objects = K, n_judges = J, reps_per_pair = reps_pp,
                     dependence = list(exposure = 0, carry_over = 0), seed = 30000 + i)
  f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge", order = "order")
  bd <- btl_dimensionality(f, reps = 60)
  shared_withheld[i] <- is.na(bd$leading_structured)
}
cat("== Shared-order verdict withheld (simulate_btl's default order layout) ==\n")
cat("verdict withheld (NA) in", sum(shared_withheld), "/", n_shared_check, "replicates\n")

## --- custom generator: genuinely per-judge-randomised order ---------------
## (simulate_btl's own order column is pair-major within judge -- shown above
## to trigger the shared-order guard -- so a TRUE randomised-order design
## needs a custom layer, per the brief)
randomise_order <- function(d) {
  d$order <- unlist(tapply(seq_len(nrow(d)), d$judge, function(ix) sample(seq_along(ix))),
                     use.names = FALSE)[order(order(d$judge))]
  d
}

## --- (7a) null flag rate with randomised orders ---------------------------
n_rep <- 150
flag_null <- rep(NA, n_rep)
for (i in seq_len(n_rep)) {
  ok <- tryCatch({
    d <- simulate_btl(n_objects = K, n_judges = J, reps_per_pair = reps_pp, seed = 31000 + i)
    d$judge_dummy <- d$judge  # already has judge; keep name clarity
    d <- randomise_order(d)
    f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge", order = "order")
    bd <- btl_dimensionality(f, reps = 60)
    flag_null[i] <- isTRUE(bd$leading_structured)
    TRUE
  }, error = function(e) FALSE)
}
cat("\n== btl_dimensionality null flag rate, randomised per-judge order (K=10,J=14) ==\n")
cat("flag rate (leading bimension above 95th pct reference):",
    round(mean(flag_null, na.rm = TRUE), 4), " n=", sum(!is.na(flag_null)),
    " MC err~", round(sqrt(0.05*0.95/sum(!is.na(flag_null))), 4), "\n")

## --- (7b) planted judge-camp second attribute: power ----------------------
n_rep2 <- 100
flag_camp <- rep(NA, n_rep2)
for (i in seq_len(n_rep2)) {
  ok <- tryCatch({
    d <- simulate_btl(n_objects = K, n_judges = J, reps_per_pair = reps_pp,
                       second_attribute = list(rho = 0.1), seed = 32000 + i)
    f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
    bd <- btl_dimensionality(f, reps = 60)
    flag_camp[i] <- isTRUE(bd$leading_structured)
    TRUE
  }, error = function(e) FALSE)
}
cat("\n== btl_dimensionality power, planted judge-camp second attribute (rho=0.1) ==\n")
cat("flag rate:", round(mean(flag_camp, na.rm = TRUE), 4), " n=", sum(!is.na(flag_camp)), "\n")
