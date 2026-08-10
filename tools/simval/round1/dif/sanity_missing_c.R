suppressWarnings(pkgload::load_all("."), quiet=TRUE))
source("gen_repeated.R")
N_c <- 320; n_items_c <- 8; common_shift_item <- 5
seed <- 1
set.seed(seed)
long <- gen_repeated(N = N_c, n_items = n_items_c, g_levels = 2, occ_levels = 2,
                      occ_dif_item = common_shift_item, occ_dif_shift = 1.0, seed = seed)
is_g2_t2 <- long$group == "g2" & long$occasion == "T2"
is_g1_t2 <- long$group == "g1" & long$occasion == "T2"
drop_g2 <- is_g2_t2 & (runif(nrow(long)) < 0.55)
drop_g1 <- is_g1_t2 & (runif(nrow(long)) < 0.10)
item_cols <- sprintf("I%02d", seq_len(n_items_c))
long[drop_g2 | drop_g1, item_cols] <- NA
keep <- rowSums(!is.na(long[, item_cols])) > 0
long <- long[keep, ]
cat("rows kept:", nrow(long), "of", N_c*2, "\n")
fit <- rasch(long, id = "person", factors = c("group", "occasion"))
da <- dif_anova(fit, within = "occasion", effects = "factorial")
print(da$summary[, c("item","term","uniform_DIF","nonuniform_DIF")])
print(da$notes)
