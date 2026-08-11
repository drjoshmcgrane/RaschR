suppressWarnings(pkgload::load_all(".", quiet=TRUE))
source("tools/simval/round1/dif/gen_repeated.R")

# within-person only (occasion 2 levels), N=300 persons
t0 <- Sys.time()
d1 <- gen_repeated(N = 300, n_items = 8, g_levels = 1, occ_levels = 2, seed = 1)
fit1 <- rasch(d1, id = "person", factors = "occasion")
da1 <- dif_anova(fit1, within = "occasion")
t1 <- Sys.time()
cat("within-only (occ=2) fit+dif time:", as.numeric(t1-t0,units="secs"), "\n")

# within occ=3 (GG)
t0 <- Sys.time()
d2 <- gen_repeated(N = 300, n_items = 8, g_levels = 1, occ_levels = 3, seed = 1)
fit2 <- rasch(d2, id = "person", factors = "occasion")
da2 <- dif_anova(fit2, within = "occasion")
t1 <- Sys.time()
cat("within-only (occ=3) fit+dif time:", as.numeric(t1-t0,units="secs"), "\n")

# mixed: group(2) x occasion(2), factorial
t0 <- Sys.time()
d3 <- gen_repeated(N = 300, n_items = 8, g_levels = 2, occ_levels = 2, seed = 1)
fit3 <- rasch(d3, id = "person", factors = c("group","occasion"))
da3 <- dif_anova(fit3, within = "occasion", effects = "factorial")
t1 <- Sys.time()
cat("mixed (2x2) factorial fit+dif time:", as.numeric(t1-t0,units="secs"), "\n")
print(unique(da3$summary$term))
