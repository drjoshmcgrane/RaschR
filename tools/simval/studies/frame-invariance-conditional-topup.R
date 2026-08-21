# STUDY: frame-invariance-conditional-topup
#
# Fresh-seed top-up of the largest-sample conditional null cell after the
# main power study returned a combined Holm familywise rate close to the
# upper edge of Monte Carlo uncertainty. This isolates that cell rather than
# repeating the departure scenarios.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "frame-invariance-conditional-topup"

R <- as.integer(Sys.getenv("SV_FRAME_TOPUP_REPS", "2000"))
N <- 2000L
K <- 8L
items <- sprintf("I%02d", seq_len(K))
delta <- seq(-1.5, 1.5, length.out = K)
ratio <- 1.40
unit <- c(ratio^-0.5, ratio^0.5)
N_CORES <- as.integer(Sys.getenv("SV_CORES", "4"))
if (!is.finite(N_CORES) || N_CORES < 1L) N_CORES <- 1L

gen <- function(seed) {
  set.seed(seed)
  one <- function(g) {
    theta <- stats::rnorm(N, 0, 1.3)
    X <- vapply(seq_len(K), function(i)
      stats::rbinom(N, 1, stats::plogis(unit[g] * (theta - delta[i]))),
      numeric(N))
    colnames(X) <- items
    X
  }
  X <- rbind(one(1L), one(2L))
  data.frame(id = sprintf("P%05d", seq_len(2L * N)), X,
             group = rep(c("g1", "g2"), each = N), check.names = FALSE)
}

one_rep <- function(r) {
  f <- tryCatch(rasch_efrm(
    gen(1800000 + r), items = items, item_sets = list(set1 = items),
    groups = "group", id = "id", boot_reps = 0), error = function(e) NULL)
  if (is.null(f)) return(list(status = "refused"))
  if (!isTRUE(f$est$converged)) return(list(status = "nonconv"))
  z <- tryCatch(frame_invariance(f, se_method = "conditional"),
                error = function(e) NULL)
  if (is.null(z) || any(!is.finite(c(z$locations$p, z$locations$p_adj,
                                     z$discrimination$p,
                                     z$discrimination$p_adj))))
    return(list(status = "refused"))
  list(status = "analysed", location = z$locations,
       discrimination = z$discrimination)
}

rr <- if (.Platform$OS.type == "windows" || N_CORES == 1L)
  lapply(seq_len(R), one_rep) else
    parallel::mclapply(seq_len(R), one_rep,
      mc.cores = min(N_CORES, R), mc.preschedule = FALSE, mc.set.seed = FALSE)
status <- vapply(rr, `[[`, "", "status")
ok <- which(status == "analysed")
lp <- la <- dp <- da <- matrix(NA_real_, length(ok), K,
                                dimnames = list(NULL, items))
for (j in seq_along(ok)) {
  z <- rr[[ok[j]]]
  lp[j, z$location$item] <- z$location$p
  la[j, z$location$item] <- z$location$p_adj
  dp[j, z$discrimination$item] <- z$discrimination$p
  da[j, z$discrimination$item] <- z$discrimination$p_adj
}
complete <- stats::complete.cases(lp, la, dp, da)
lp <- lp[complete, , drop = FALSE]; la <- la[complete, , drop = FALSE]
dp <- dp[complete, , drop = FALSE]; da <- da[complete, , drop = FALSE]
n <- nrow(lp)
n_refused <- sum(status == "refused")
n_nonconv <- sum(status == "nonconv")
rows <- list()
put <- function(quantity, ..., mc_override = list())
  rows[[length(rows) + 1L]] <<- sv_row(
    STUDY, sprintf("null: unit ratio only | N = %d", N), quantity, n,
    ..., mc_override = mc_override, n_attempted = R,
    n_refused = n_refused, n_nonconv = n_nonconv)

for (channel in c("location", "discrimination")) {
  raw <- if (channel == "location") lp else dp
  holm <- if (channel == "location") la else da
  for (rule in c("raw", "combined Holm")) {
    M <- (if (rule == "raw") raw else holm) < 0.05
    per <- rowMeans(M)
    put(paste(channel, rule, "per-item type I error"),
        type1 = mean(per),
        mc_override = list(type1 = stats::sd(per) / sqrt(n)))
    put(paste(channel, rule, "channel familywise type I error"),
        familywise = mean(rowSums(M) > 0L))
  }
}
fw <- rowSums(cbind(la < 0.05, da < 0.05)) > 0L
put("combined Holm familywise type I error", familywise = mean(fw))
for (item in items)
  put(paste("discrimination raw type I error", item),
      type1 = mean(dp[, item] < 0.05))

sv_write(do.call(rbind, rows), "frame-invariance-conditional-topup")
cat(sprintf("analysed=%d/%d, combined Holm FWER=%.4f\n",
            n, R, mean(fw)))
