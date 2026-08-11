suppressWarnings(pkgload::load_all(".", quiet = TRUE))
t_start <- Sys.time()

if (!requireNamespace("psychotools", quietly = TRUE))
  stop("psychotools not available")
data("GermanParties2009", package = "psychotools")
gp <- GermanParties2009

m <- as.matrix(gp$preference)
pair_names <- strsplit(colnames(m), ":", fixed = TRUE)
d6 <- do.call(rbind, lapply(seq_along(pair_names), function(j) {
  a <- pair_names[[j]][1]; b <- pair_names[[j]][2]
  data.frame(object_a = a, object_b = b,
             winner   = ifelse(m[, j] == 1, a, b),
             judge    = sprintf("R%03d", seq_len(nrow(m))),
             crisis   = as.character(gp$crisis),
             stringsAsFactors = FALSE)
}))
d <- d6[d6$object_a != "none" & d6$object_b != "none", ]

blocs <- list(left = c("Linke", "Gruene", "SPD"), right = c("CDU/CSU", "FDP"))
crisis_of <- setNames(d$crisis[!duplicated(d$judge)], d$judge[!duplicated(d$judge)])

set.seed(2009)
f <- tryCatch(btl_efrm(d, "object_a", "object_b", "winner", "judge",
                        panels = crisis_of, object_sets = blocs, boot_reps = 80),
              error = function(e) e)

if (inherits(f, "error")) {
  cat("FIT FAILED:", conditionMessage(f), "\n")
} else {
  cat("Fit converged:", f$converged, "\n")
  cat("Notes:\n"); print(f$notes)
  cat("\nalpha_table:\n"); print(f$alpha_table, digits = 3)
  cat("\nkappa_table:\n"); print(f$kappa_table, digits = 3)
  cat("\nphi_table:\n"); print(f$phi_table, digits = 3)
  right_alpha <- f$alpha_table$alpha[f$alpha_table$set == "right"]
  has_screen_note <- any(grepl("no stable|unidentified", f$notes))
  cat("\nright-bloc alpha (expect NA, one internal pair split near-even):", right_alpha, "\n")
  cat("screen note present:", has_screen_note, "\n")
  cat("PASS (fits + screen note as documented):", isTRUE(f$converged) && is.na(right_alpha) && has_screen_note, "\n")
}
cat("\ntotal time (s):", as.numeric(Sys.time() - t_start, units = "secs"), "\n")
