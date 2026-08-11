suppressWarnings(pkgload::load_all(".", quiet=TRUE))
source("tools/simval/round1/dif/gen_repeated.R")
# (b) structural differential missingness, no "missing=" argument used
n_items_b <- 20; skip_idx <- 16:20; N_b <- 600
d <- simulate_rasch(N_b, n_items_b, n_groups = 2, seed = 1)
d[d$group == "g1", sprintf("I%02d", skip_idx)] <- NA
t0 <- Sys.time()
fit <- rasch(d, id = "id", factors = "group")
t1 <- Sys.time()
da <- dif_anova(fit)
t2 <- Sys.time()
cat("(b) fit time:", as.numeric(t1-t0,units="secs"), " dif_anova time:", as.numeric(t2-t1,units="secs"), "\n")

# (c) mixed design with dropout
N_c <- 320; n_items_c <- 8; common_shift_item <- 5
set.seed(1)
long <- gen_repeated(N = N_c, n_items = n_items_c, g_levels = 2, occ_levels = 2,
                      occ_dif_item = common_shift_item, occ_dif_shift = 1.0, seed = 1)
is_g2_t2 <- long$group == "g2" & long$occasion == "T2"
is_g1_t2 <- long$group == "g1" & long$occasion == "T2"
drop_g2 <- is_g2_t2 & (runif(nrow(long)) < 0.55)
drop_g1 <- is_g1_t2 & (runif(nrow(long)) < 0.10)
item_cols <- sprintf("I%02d", seq_len(n_items_c))
long[drop_g2 | drop_g1, item_cols] <- NA
keep <- rowSums(!is.na(long[, item_cols])) > 0
long <- long[keep, ]
t0 <- Sys.time()
fit <- rasch(long, id = "person", factors = c("group", "occasion"))
t1 <- Sys.time()
da <- dif_anova(fit, within = "occasion", effects = "factorial")
t2 <- Sys.time()
cat("(c) fit time:", as.numeric(t1-t0,units="secs"), " dif_anova time:", as.numeric(t2-t1,units="secs"), "\n")
