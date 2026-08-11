suppressWarnings(pkgload::load_all(".", quiet=TRUE))
t0 <- Sys.time()
d <- simulate_btl_efrm(n_objects_per_set=6, n_sets=2, n_panels=2, n_judges_per_panel=6,
                        reps_within=20, reps_cross=20,
                        panel_units=c(1,1.3), set_units=c(1,1.4), set_origins=c(0,0.5), seed=1)
t1 <- Sys.time()
fit <- btl_efrm(d, "object_a","object_b","winner","judge","panel",
                 object_sets=attr(d,"truth")$object_sets, se_method="conditional")
t2 <- Sys.time()
cat("sim time:", as.numeric(t1-t0,units="secs"), "\n")
cat("conditional fit time:", as.numeric(t2-t1,units="secs"), "\n")
fitb <- btl_efrm(d, "object_a","object_b","winner","judge","panel",
                 object_sets=attr(d,"truth")$object_sets, se_method="judge_bootstrap", boot_reps=60)
t3 <- Sys.time()
cat("judge_bootstrap(60) fit time:", as.numeric(t3-t2,units="secs"), "\n")
fitp <- btl_efrm(d, "object_a","object_b","winner","judge","panel",
                 object_sets=attr(d,"truth")$object_sets, se_method="bootstrap", boot_reps=60)
t4 <- Sys.time()
cat("param bootstrap(60) fit time:", as.numeric(t4-t3,units="secs"), "\n")
print(fit$alpha_table)
print(fit$phi_table)
