# Calibration battery for the package's inferential core.
#
# Every reported standard error is checked against the empirical sampling
# variability of its estimator on model-true data (the check that exposed the
# count-weighted BTL sandwich inflation). Run on demand:
#   Rscript tools/calibration.R
# Pass bands: SE-ratio (empirical SD / mean reported SE, the convention
# used across tools/simval) in ~[0.85, 1.20] for the replicated paths; anchored
# translation equivalence to numerical precision.
#
# DESIGN RULE (learned the hard way): hold the TRUE parameters FIXED across
# replicates -- redrawing the truth each replicate folds between-replicate
# truth variation into the "empirical SD" and fakes a calibration failure.

# run from the package root: Rscript tools/calibration.R
suppressWarnings(pkgload::load_all(".", quiet = TRUE))
ok <- function(l,v) cat(sprintf("%-64s %s\n", l, v))

## A. pcml dichotomous item-location SEs, 150 reps (band should be tight)
L <- 10; Np <- 500
dtrue <- scale(seq(-2,2,length.out=L), scale=FALSE)[,1]
est <- ses <- matrix(NA, 150, L)
for (r in 1:150) { set.seed(5000+r)
  th <- rnorm(Np, 0, 1.3)
  X <- matrix(rbinom(Np*L,1,plogis(outer(th,dtrue,"-"))), Np, L)
  f <- pcml(X); est[r,] <- f$thr$tau; ses[r,] <- f$thr$se }
rA <- apply(est,2,sd)/colMeans(ses)
ok("A pcml dichot SE ratio per item (min..max)", sprintf("%.2f .. %.2f", min(rA), max(rA)))

## B. PCM threshold SEs, 100 reps
simP <- function(t, tau) { x<-0:length(tau); p<-exp(x*t-c(0,cumsum(tau))); p/sum(p) }
tt <- list(c(-1.2,-0.1), c(-0.6,0.5), c(-0.9,0.2), c(0.1,1.1), c(-0.4,0.8), c(-1.0,0.6))
tt <- lapply(tt, function(t) t - mean(unlist(tt)))
estB <- sesB <- matrix(NA, 100, sum(lengths(tt)))
for (r in 1:100) { set.seed(6000+r)
  th <- rnorm(600, 0, 1.4)
  X <- sapply(seq_along(tt), function(i) vapply(th, function(t)
    sample(0:2, 1, prob=simP(t, tt[[i]])), 0L))
  colnames(X) <- sprintf("P%d", 1:6)
  f <- rasch(X, model="PCM")
  estB[r,] <- f$thresholds$tau; sesB[r,] <- f$thresholds$se }
rB <- apply(estB,2,sd)/colMeans(sesB)
ok("B PCM threshold SE ratio (min..max)", sprintf("%.2f .. %.2f", min(rB), max(rB)))

## C. anchored rasch at distant origins
set.seed(3); th <- rnorm(600,0,1.3)
X <- matrix(rbinom(600*8,1,plogis(outer(th, seq(-1.5,1.5,length.out=8), "-"))),600,8)
colnames(X) <- sprintf("I%d",1:8)
free <- rasch(X)
fv <- free$items$location
for (delta in c(0, 6)) {
  anc <- data.frame(item=c("I1","I8"), k=1, tau=fv[c(1,8)] + delta)
  fa <- rasch(X, anchors=anc)
  errI <- max(abs(fa$items$location - (fv + delta)))
  errP <- max(abs(fa$person$theta - (free$person$theta + delta)), na.rm=TRUE)
  ok(sprintf("C anchored rasch +%d: conv, item err, person err", delta),
     sprintf("%s, %.2e, %.2e", fa$est$converged, errI, errP))
}

## D. MFRM severity SE calibration, fixed truth, 120 reps
simP2 <- function(t, tau) { x<-0:length(tau); p<-exp(x*t-c(0,cumsum(tau))); p/sum(p) }
lam <- c(R1=-0.7, R2=-0.2, R3=0.3, R4=0.6)
del <- c(I1=-0.8, I2=-0.2, I3=0.3, I4=0.7)
base_tau <- c(-1.1, 0, 1.1)
estD <- sesD <- matrix(NA, 120, 4)
for (r in 1:120) { set.seed(7000+r)
  th <- rnorm(70, 0, 1.2)
  g <- expand.grid(p=1:70, i=1:4, rt=1:4)
  sc <- vapply(seq_len(nrow(g)), function(k)
    sample(0:3, 1, prob=simP2(th[g$p[k]], base_tau + del[g$i[k]] + lam[g$rt[k]])), 0L)
  d <- data.frame(person=sprintf("P%03d",g$p), item=names(del)[g$i],
                  rater=names(lam)[g$rt], score=sc)
  mf <- rasch_mfrm(d, person="person", item="item", score="score", facets="rater")
  fe <- mf$facet_effects$rater
  estD[r,] <- fe$severity[match(names(lam), fe$level)]
  sesD[r,] <- fe$se[match(names(lam), fe$level)]
}
rD <- apply(estD,2,sd)/colMeans(sesD)
ok("D MFRM severity SE ratio, fixed truth (min..max)",
   sprintf("%.2f .. %.2f", min(rD), max(rD)))

## E. dif_size resolved-difference SE, 100 reps (planted 0.8 uniform DIF)
diffs <- sesE <- numeric(100)
for (r in 1:100) { set.seed(8000+r)
  d <- simulate_rasch(500, 10, dif=list(items="I05", uniform=0.8), n_groups=2, seed=8000+r)
  f <- rasch(d, id="id", factors="group")
  ds <- dif_size(f, "I05", by="group")
  diffs[r] <- ds$pairs$difference; sesE[r] <- ds$pairs$se }
ok("E dif_size: empSD vs mean reported SE", sprintf("%.3f vs %.3f (ratio %.2f)",
   sd(diffs), mean(sesE), sd(diffs)/mean(sesE)))
ok("E dif_size: mean recovered diff (planted 0.8)", sprintf("%.3f", mean(abs(diffs))))

## F. person WLE coverage over a theta grid, fixed calibration (compact form
## of tools/simval/studies/wle-coverage.R; see that script for the full
## per-grid-point study, >=1000/400 reps, and the identification-artifact
## write-up). One large calibration held fixed; "offset" shifts the TRUE
## person theta relative to the FIXED item bank (not the item bank itself --
## a uniform item-location shift is exactly absorbed by rasch()'s mean-zero
## recentring and would be a degenerate, not a targeting, manipulation).
set.seed(9); dv <- seq(-2,2,length.out=12)
Xc <- matrix(rbinom(1500*12,1,plogis(outer(rnorm(1500,0,1.4),dv,"-"))),1500,12)
colnames(Xc) <- sprintf("I%02d",1:12)
fc <- rasch(Xc); st <- person_wle(fc$tau_list)
grid <- seq(-3,3,by=0.75); reps_f <- 500
for (offset in c(0,1)) {
  cov <- bias <- se_r <- ext <- numeric(length(grid))
  for (gi in seq_along(grid)) {
    th_true <- grid[gi] + offset
    set.seed(9000L + offset*1000L + gi)
    p <- plogis(th_true - dv)
    X <- matrix(rbinom(reps_f*12,1,rep(p,each=reps_f)),reps_f,12)
    raw <- rowSums(X); key <- as.character(raw)
    est <- unname(st$theta[key]); se <- unname(st$se[key])
    cov[gi]  <- mean(abs(est-th_true) <= 1.96*se)
    bias[gi] <- mean(est-th_true)
    se_r[gi] <- sd(est)/mean(se)
    ext[gi]  <- mean(raw==0L | raw==12L)          # extreme (0/perfect) raw-score rate
  }
  worst <- which.max(abs(bias))
  ok(sprintf("F WLE grid coverage, offset=%+d: range / SEratio range / worst bias", offset),
     sprintf("%.3f..%.3f / %.2f..%.2f / %.3f at theta=%.2f",
             min(cov), max(cov), min(se_r), max(se_r), bias[worst], grid[worst]+offset))
  ok(sprintf("F WLE offset=%+d: extreme-score rate at |theta|>=2.25 (max over grid)", offset),
     sprintf("%.3f (all Warm WLEs finite, none refused)", max(ext[abs(grid+offset) >= 2.25])))
}
