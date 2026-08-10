suppressWarnings(pkgload::load_all("."), quiet = TRUE))
t_start <- Sys.time()

results <- list()
add_result <- function(name, expr_ok, msg) {
  results[[length(results) + 1]] <<- data.frame(guard = name, refused = expr_ok, message = msg)
}

try_expect_error <- function(name, expr, pattern) {
  out <- tryCatch({ force(expr); list(err = NULL) },
                   error = function(e) list(err = conditionMessage(e)))
  if (is.null(out$err)) {
    add_result(name, FALSE, "NO ERROR RAISED (expected one)")
  } else {
    hit <- grepl(pattern, out$err, ignore.case = TRUE)
    add_result(name, hit, out$err)
  }
}

## (a) within-set disconnection ---------------------------------------------
# set1 has 4 objects split into two disconnected pairs (A1-A2) and (A3-A4),
# never compared to each other; set2 is a normal fully-connected 4-object set
# providing the cross-set link.
set.seed(1)
mk_pairs <- function(objs, n = 15) {
  pr <- t(utils::combn(objs, 2))
  data.frame(object_a = rep(pr[, 1], n), object_b = rep(pr[, 2], n),
             stringsAsFactors = FALSE)
}
d_a1 <- mk_pairs(c("A1", "A2"))
d_a2 <- mk_pairs(c("A3", "A4"))
within1 <- rbind(d_a1, d_a2)
within2 <- mk_pairs(c("B1", "B2", "B3", "B4"))
cross <- expand.grid(object_a = c("A1", "A2", "A3", "A4"),
                      object_b = c("B1", "B2", "B3", "B4"),
                      stringsAsFactors = FALSE)
cross <- cross[rep(seq_len(nrow(cross)), 15), ]
allrows <- rbind(within1, within2, cross)
set.seed(2)
allrows$winner <- ifelse(runif(nrow(allrows)) < 0.5, allrows$object_a, allrows$object_b)
allrows$judge <- sample(sprintf("J%02d", 1:20), nrow(allrows), replace = TRUE)
try_expect_error("within-set disconnection",
  btl_efrm(allrows, "object_a", "object_b", "winner", "judge",
           panels = setNames(rep("panel1", 20), sprintf("J%02d", 1:20)),
           object_sets = list(set1 = c("A1", "A2", "A3", "A4"),
                               set2 = c("B1", "B2", "B3", "B4")),
           se_method = "conditional"),
  "not.*connected")

## (b) cross-set separation ---------------------------------------------------
# set1 always beats set2 in every cross-set comparison (deterministic), while
# each set is internally well connected and balanced.
set.seed(3)
within1b <- mk_pairs(c("A1", "A2", "A3", "A4"))
within2b <- mk_pairs(c("B1", "B2", "B3", "B4"))
within1b$winner <- ifelse(runif(nrow(within1b)) < 0.5, within1b$object_a, within1b$object_b)
within2b$winner <- ifelse(runif(nrow(within2b)) < 0.5, within2b$object_a, within2b$object_b)
crossb <- expand.grid(object_a = c("A1", "A2", "A3", "A4"),
                       object_b = c("B1", "B2", "B3", "B4"),
                       stringsAsFactors = FALSE)
crossb <- crossb[rep(seq_len(nrow(crossb)), 15), ]
crossb$winner <- crossb$object_a       # A always wins: deterministic separation
allb <- rbind(within1b, within2b, crossb)
allb$judge <- sample(sprintf("J%02d", 1:20), nrow(allb), replace = TRUE)
try_expect_error("cross-set separation",
  btl_efrm(allb, "object_a", "object_b", "winner", "judge",
           panels = setNames(rep("panel1", 20), sprintf("J%02d", 1:20)),
           object_sets = list(set1 = c("A1", "A2", "A3", "A4"),
                               set2 = c("B1", "B2", "B3", "B4")),
           se_method = "conditional"),
  "separat")

## (c) degenerate set unit -> alpha reported NA ------------------------------
# set2's four objects are exactly tied within-set (50/50 balanced outcomes,
# no spread at all) but the cross-set contests still separate set2 from set1
# with a clear (non-degenerate) margin, giving kappa[2] identification but
# nothing for alpha[2] to scale.
set.seed(4)
within1c <- mk_pairs(c("A1", "A2", "A3", "A4"), n = 30)
within1c$winner <- ifelse((seq_len(nrow(within1c)) %% 2) == 0,
                           within1c$object_a, within1c$object_b)  # exact 50/50 per pair... but need spread for set1
# give set1 real spread: bias by rank
pr1 <- t(utils::combn(c("A1", "A2", "A3", "A4"), 2))
rank1 <- setNames(4:1, c("A1", "A2", "A3", "A4"))
within1c <- do.call(rbind, lapply(seq_len(nrow(pr1)), function(i) {
  oa <- pr1[i, 1]; ob <- pr1[i, 2]
  p <- plogis(1.5 * (rank1[[oa]] - rank1[[ob]]))
  n <- 40
  data.frame(object_a = oa, object_b = ob,
             winner = ifelse(runif(n) < p, oa, ob), stringsAsFactors = FALSE)
}))
within2c <- t(utils::combn(c("B1", "B2", "B3", "B4"), 2))
within2c <- do.call(rbind, lapply(seq_len(nrow(within2c)), function(i) {
  oa <- within2c[i, 1]; ob <- within2c[i, 2]
  n <- 40
  data.frame(object_a = oa, object_b = ob,
             winner = ifelse((seq_len(n) %% 2) == 0, oa, ob), stringsAsFactors = FALSE)
}))
crossc <- expand.grid(object_a = c("A1", "A2", "A3", "A4"),
                       object_b = c("B1", "B2", "B3", "B4"),
                       stringsAsFactors = FALSE)
crossc <- crossc[rep(seq_len(nrow(crossc)), 20), ]
set.seed(5)
crossc$winner <- ifelse(runif(nrow(crossc)) < 0.7, crossc$object_a, crossc$object_b)
allc <- rbind(within1c, within2c, crossc)
allc$judge <- sample(sprintf("J%02d", 1:24), nrow(allc), replace = TRUE)
fit_c <- tryCatch(btl_efrm(allc, "object_a", "object_b", "winner", "judge",
                            panels = setNames(rep("panel1", 24), sprintf("J%02d", 1:24)),
                            object_sets = list(set1 = c("A1", "A2", "A3", "A4"),
                                                set2 = c("B1", "B2", "B3", "B4")),
                            se_method = "conditional"),
                   error = function(e) e)
if (inherits(fit_c, "error")) {
  add_result("degenerate set unit -> alpha NA", FALSE,
             paste("ERROR (expected a fit with alpha[set2]=NA):", conditionMessage(fit_c)))
} else {
  alpha2 <- fit_c$alpha_table$alpha[fit_c$alpha_table$set == "set2"]
  has_note <- any(grepl("unidentified", fit_c$notes))
  ok <- is.na(alpha2) && has_note
  add_result("degenerate set unit -> alpha NA", ok,
             sprintf("alpha[set2] = %s; notes: %s", format(alpha2),
                     paste(fit_c$notes, collapse = " | ")))
}

## (d) judges in two panels --------------------------------------------------
set.seed(6)
within1d <- mk_pairs(c("A1", "A2", "A3", "A4"))
within2d <- mk_pairs(c("B1", "B2", "B3", "B4"))
crossd <- expand.grid(object_a = c("A1", "A2", "A3", "A4"),
                       object_b = c("B1", "B2", "B3", "B4"),
                       stringsAsFactors = FALSE)
crossd <- crossd[rep(seq_len(nrow(crossd)), 15), ]
alld <- rbind(within1d, within2d, crossd)
alld$winner <- ifelse(runif(nrow(alld)) < 0.5, alld$object_a, alld$object_b)
alld$judge <- sample(sprintf("J%02d", 1:20), nrow(alld), replace = TRUE)
alld$panel <- ifelse(as.integer(sub("J", "", alld$judge)) <= 10, "panel1", "panel2")
# corrupt: J01's panel flips between panel1 and panel2 across rows
flip <- alld$judge == "J01"
alld$panel[flip] <- ifelse(seq_len(sum(flip)) %% 2 == 0, "panel1", "panel2")
try_expect_error("judge assigned to two panels",
  btl_efrm(alld, "object_a", "object_b", "winner", "judge", panels = "panel",
           object_sets = list(set1 = c("A1", "A2", "A3", "A4"),
                               set2 = c("B1", "B2", "B3", "B4")),
           se_method = "conditional"),
  "more than one panel")

## (e) bonus: cross-set comparisons touch only 1 object of a set (structural)
set.seed(7)
within1e <- mk_pairs(c("A1", "A2", "A3", "A4"))
within2e <- mk_pairs(c("B1", "B2", "B3", "B4"))
crosse <- data.frame(object_a = rep("A1", 60), object_b = "B1", stringsAsFactors = FALSE)
alle <- rbind(within1e, within2e, crosse)
alle$winner <- ifelse(runif(nrow(alle)) < 0.5, alle$object_a, alle$object_b)
alle$judge <- sample(sprintf("J%02d", 1:20), nrow(alle), replace = TRUE)
try_expect_error("cross-set touches <2 objects of a set (structural)",
  btl_efrm(alle, "object_a", "object_b", "winner", "judge",
           panels = setNames(rep("panel1", 20), sprintf("J%02d", 1:20)),
           object_sets = list(set1 = c("A1", "A2", "A3", "A4"),
                               set2 = c("B1", "B2", "B3", "B4")),
           se_method = "conditional", min_link = 5),
  "unidentified")

## (f) bonus: sets not reachable via cross-set pairs (min_link too high) ----
set.seed(8)
within1f <- mk_pairs(c("A1", "A2", "A3", "A4"))
within2f <- mk_pairs(c("B1", "B2", "B3", "B4"))
crossf <- expand.grid(object_a = c("A1", "A2", "A3", "A4"),
                       object_b = c("B1", "B2", "B3", "B4"),
                       stringsAsFactors = FALSE)
crossf <- crossf[rep(seq_len(nrow(crossf)), 2), ]     # only 32 cross comps total
allf <- rbind(within1f, within2f, crossf)
allf$winner <- ifelse(runif(nrow(allf)) < 0.5, allf$object_a, allf$object_b)
allf$judge <- sample(sprintf("J%02d", 1:20), nrow(allf), replace = TRUE)
try_expect_error("sets unreachable given min_link",
  btl_efrm(allf, "object_a", "object_b", "winner", "judge",
           panels = setNames(rep("panel1", 20), sprintf("J%02d", 1:20)),
           object_sets = list(set1 = c("A1", "A2", "A3", "A4"),
                               set2 = c("B1", "B2", "B3", "B4")),
           se_method = "conditional", min_link = 100),
  "not reachable")

res <- do.call(rbind, results)
print(res, row.names = FALSE)
saveRDS(res, "tools/simval/round1/btl-efrm/04_guards.rds")
cat("\ntotal time (s):", as.numeric(Sys.time() - t_start, units = "secs"), "\n")
