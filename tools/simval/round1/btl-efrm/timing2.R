suppressWarnings(pkgload::load_all(".", quiet=TRUE))
d <- simulate_btl_efrm(n_objects_per_set=6, n_sets=2, n_panels=2, n_judges_per_panel=6,
                        reps_within=20, reps_cross=20,
                        panel_units=c(0.8,1.25), set_units=c(1,1.4), set_origins=c(0,0.5), seed=1)
t0<-Sys.time()
f<-btl_efrm(d,"object_a","object_b","winner","judge","panel",object_sets=attr(d,"truth")$object_sets,
            se_method="judge_bootstrap", boot_reps=100)
cat("judge_boot 100:", as.numeric(Sys.time()-t0,units="secs"),"\n")
t0<-Sys.time()
f2<-btl_efrm(d,"object_a","object_b","winner","judge","panel",object_sets=attr(d,"truth")$object_sets,
            se_method="bootstrap", boot_reps=100)
cat("param_boot 100:", as.numeric(Sys.time()-t0,units="secs"),"\n")
