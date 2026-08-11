suppressWarnings(pkgload::load_all(".", quiet=TRUE))
t0 <- Sys.time()
d <- simulate_btl(n_objects=10, n_judges=14, reps_per_pair=25,
                   dependence=list(exposure=0, carry_over=0), seed=1)
f <- btl(d, "object_a","object_b", winner="winner", judge="judge", order="order")
bd <- btl_dimensionality(f, reps=100)
t1 <- Sys.time()
cat("time:", as.numeric(t1-t0,units="secs"), "\n")
print(bd$bimensions)
cat("leading_structured:", bd$leading_structured, "\n")
cat(bd$notes, sep="\n")

# camp
t0 <- Sys.time()
d2 <- simulate_btl(n_objects=10, n_judges=14, reps_per_pair=25,
                    second_attribute=list(rho=0.2), seed=2)
f2 <- btl(d2, "object_a","object_b", winner="winner", judge="judge")
bd2 <- btl_dimensionality(f2, reps=100)
t1 <- Sys.time()
cat("\ncamp time:", as.numeric(t1-t0,units="secs"), "\n")
print(bd2$bimensions)
cat("leading_structured:", bd2$leading_structured, "\n")
