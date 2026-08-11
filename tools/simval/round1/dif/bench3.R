suppressWarnings(pkgload::load_all(".", quiet=TRUE))

t0 <- Sys.time()
d <- simulate_rasch(500, 20, dif = list(items = "I06", uniform = 0.7),
                     n_groups = 2, seed = 1)
fit <- rasch(d, id = "id", factors = "group")
ds <- dif_size(fit, "I06", by = "group")
t1 <- Sys.time()
cat("dif_size time:", as.numeric(t1-t0,units="secs"), "\n")
print(ds$pairs)

t0 <- Sys.time()
dc <- dif_contrasts(fit, items = "I06")
t1 <- Sys.time()
cat("dif_contrasts time:", as.numeric(t1-t0,units="secs"), "\n")
print(dc$table[dc$table$item=="I06",])

t0 <- Sys.time()
rd <- resolve_dif(fit)
t1 <- Sys.time()
cat("resolve_dif time:", as.numeric(t1-t0,units="secs"), "\n")
print(rd$splits)
print(rd$stopped)
