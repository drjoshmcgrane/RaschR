# rasch :: extended frame of reference for paired comparisons
# ===========================================================================
# The extended frame of reference model of Humphry (2005) and Humphry and
# Andrich (2008) states that the unit of the latent scale is a property of
# the frame of measurement, not a universal constant. rasch_efrm() fits that
# model for the persons-by-items design; btl_efrm() is this package's
# extension of the same idea to the Bradley-Terry-Luce family of paired
# comparisons (Bradley and Terry 1952; Luce 1959), where the frame is a
# judge-panel by object-set cell.
#
# Objects k are partitioned into sets s(k) (S sets), and judges j into panels
# g(j) (G panels). Each object has a true common-scale value
#
#     v_k = alpha_{s(k)} * beta_k + kappa_{s(k)},
#
# with beta_k the within-set (frame-unit) calibration location, alpha_s > 0
# the set unit, and kappa_s the set origin. A comparison judged in panel g
# carries the panel unit phi_g, and:
#
#   * WITHIN a set s (both objects in s):
#         logit P(a beats b) = phi_g * (beta_a - beta_b).
#     The origin kappa_s cancels, and the set unit alpha_s is confounded with
#     the spread of beta -- exactly as the item-set unit is absorbed in the
#     within-frame stage of the Rasch EFRM -- so it is not a parameter here.
#
#   * ACROSS sets (a in A, b in B):
#         logit P(a beats b) = phi_g * (v_a - v_b)
#                            = phi_g * (alpha_A beta_a - alpha_B beta_b
#                                        + kappa_A - kappa_B).
#     The cross-set comparisons place the two sets on one common scale and so
#     identify alpha and kappa.
#
# The lineage of a frame-dependent unit for comparative judgement is
# Thurstone's (1927) varying discriminal dispersion and the varying-precision
# paired-comparison models catalogued by David (1988); the measurement-unit
# reading of it is Humphry's. The paired-comparison form fitted here is this
# package's extension, stated for dichotomous winner data.
#
# Estimation is in two likelihood stages, mirroring the staged structure of
# rasch_efrm() but not its person-distribution link:
#
# Stage 1 (within frames): for each set the within-set comparisons, pooled
# over panels, fit the bilinear model logit = rho_{gs} (b_a - b_b) with
# b sum-zero and one reference panel fixed at rho = 1. This is the same
# constrained bilinear maximisation the Rasch EFRM performs in its stage 1;
# the ratios rho_{gs} = phi_g / phi_{ref(s)} estimate the panel units up to
# the set's reference panel, and are reconciled across sets by a
# precision-weighted least squares over the panel-by-set linking graph, then
# normalised to geometric mean one over panels. The reconciled reference-panel
# units put every set's b on the common panel scale, giving beta.
#
# Stage 2 (linking sets): with beta-hat and phi-hat fixed, the cross-set
# comparisons are a low-dimensional maximum likelihood in (log alpha, kappa)
# for the non-reference sets, solved by Newton with the analytic gradient and
# Hessian. The key theoretical point, and the reason this design is worth
# stating, is that the linking uses only comparison OUTCOMES: no distributional
# assumption about the objects is made, so the set units are identified WITHIN
# the conditional (person-free) framework. This is unlike the persons-by-items
# EFRM, whose item-set units are identified only from the person side (their
# distribution), a genuinely distributional step. The paired-comparison design
# supplies its own conditional link and needs no such assumption.
# ===========================================================================

# connected components of an undirected graph on 1..n given a two-column
# integer edge matrix; union-find, as used across the package's linking code
.btlef_components <- function(n, edges) {
  parent <- seq_len(n)
  find <- function(a) { while (parent[a] != a) a <- parent[a]; a }
  if (length(edges)) for (r in seq_len(nrow(edges))) {
    ra <- find(edges[r, 1]); rb <- find(edges[r, 2])
    if (ra != rb) parent[ra] <- rb
  }
  vapply(seq_len(n), find, 1L)
}

# Stage 1 bilinear solve for ONE frame set (or for the pooled single-unit
# model, called with one panel). Objects are indexed 1..K; panel is a
# character vector per comparison; the most-used panel is the reference
# (rho = 1). Parameters are the sum-zero location contrasts and the free
# panels' log discrimination. Fisher scoring with a step-halving line search
# gives the point estimate; a judge-clustered Godambe sandwich gives the
# covariance, exactly as in btl().
.btlef_stage1 <- function(ia, ib, y, panel, jd, K, maxit, tol, rho_fixed = NULL) {
  R <- length(ia)
  pcount <- table(panel)
  present <- names(pcount)
  ref <- present[which.max(as.integer(pcount))]     # reference panel: rho = 1
  free <- if (is.null(rho_fixed)) setdiff(present, ref) else character(0)
  Gf <- length(free)
  pf <- match(panel, free)                           # free-panel index, NA on ref
  rho0 <- if (is.null(rho_fixed)) rep(1, R) else as.numeric(rho_fixed[panel])
  B <- rbind(diag(K - 1L), rep(-1, K - 1L))          # K x (K-1) sum-zero map
  Bd <- B[ia, , drop = FALSE] - B[ib, , drop = FALSE]
  np <- (K - 1L) + Gf

  eval_th <- function(th) {
    bfree <- th[seq_len(K - 1L)]
    beta <- as.numeric(B %*% bfree)
    rho <- rho0
    if (Gf) {
      rr <- exp(th[(K - 1L) + seq_len(Gf)])
      ok <- !is.na(pf); rho[ok] <- rr[pf[ok]]
    }
    d <- beta[ia] - beta[ib]
    p <- plogis(rho * d)
    ll <- sum(ifelse(y == 1, log(pmax(p, 1e-300)), log(pmax(1 - p, 1e-300))))
    list(beta = beta, rho = rho, d = d, p = p, ll = ll)
  }
  design <- function(cur) {
    J <- Bd * cur$rho
    if (Gf) {
      Jr <- matrix(0, R, Gf)
      for (h in seq_len(Gf)) {
        sel <- which(pf == h); Jr[sel, h] <- cur$rho[sel] * cur$d[sel]
      }
      J <- cbind(J, Jr)
    }
    J
  }

  theta <- numeric(np); cur <- eval_th(theta)
  for (it in seq_len(maxit)) {
    J <- design(cur); u <- y - cur$p; av <- cur$p * (1 - cur$p)
    g <- crossprod(J, u); Fi <- crossprod(J, J * av)
    step <- tryCatch(solve(Fi, g),
                     error = function(e) solve(Fi + diag(1e-8, np), g))
    ms <- max(abs(step))                     # trust-region cap: a near-flat
    if (is.finite(ms) && ms > 5) step <- step * (5 / ms)   # log-rho direction
    lam <- 1; moved <- FALSE                 # cannot run away
    for (half in 1:30) {
      cand <- theta + lam * as.numeric(step); c2 <- eval_th(cand)
      if (is.finite(c2$ll) && c2$ll >= cur$ll - 1e-12) {
        theta <- cand; cur <- c2; moved <- TRUE; break
      }
      lam <- lam / 2
    }
    if (!moved) break
    if (max(abs(lam * step)) < tol) break
  }

  # judge-clustered Godambe sandwich (unclustered when every judge appears once)
  J <- design(cur); u <- y - cur$p; av <- cur$p * (1 - cur$p)
  Fi <- crossprod(J, J * av)
  # rank of the UNRIDGED information: panels that observe disjoint object
  # pairs leave the (location, log rho) system underdetermined -- the
  # ridged solve then lands somewhere on the flat manifold with a small
  # gradient, which the convergence criterion cannot distinguish from a
  # genuine optimum
  rc1 <- tryCatch(rcond(Fi), error = function(e) 0)
  rank_ok <- is.finite(rc1) && rc1 > 1e-10
  bread <- tryCatch(solve(Fi), error = function(e) solve(Fi + diag(1e-8, np)))
  Sr <- J * u
  Sc <- rowsum(Sr, jd)
  nc <- nrow(Sc)
  cluster_ok <- nc >= 10L && nc > np
  cov_theta <- if (nc > 1L) {
    bread %*% (crossprod(Sc) * nc / (nc - 1)) %*% bread
  } else {
    matrix(NA_real_, np, np)
  }
  # scale-free convergence: the gradient per comparison, invariant to the
  # number of comparisons -- an absolute threshold flags converged fits as
  # unconverged on large R (and would then misroute them into the screen);
  # a per-observation criterion also stays permissive at a boundary, where
  # the gradient vanishes but a Newton-decrement quadratic would not
  conv <- isTRUE(max(abs(crossprod(J, u))) < 1e-6 * R)

  cov_bb <- B %*% cov_theta[seq_len(K - 1L), seq_len(K - 1L), drop = FALSE] %*% t(B)
  se_beta <- sqrt(pmax(diag(cov_bb), 0))
  if (!cluster_ok) se_beta[] <- NA_real_
  rho_p <- if (is.null(rho_fixed)) setNames(rep(1, length(present)), present)
           else rho_fixed[present]
  cov_lrho <- matrix(0, Gf, Gf, dimnames = list(free, free))
  if (Gf) {
    li <- (K - 1L) + seq_len(Gf)
    rho_p[free] <- exp(theta[li])
    cov_lrho <- cov_theta[li, li, drop = FALSE]
    if (!cluster_ok) cov_lrho[] <- NA_real_
    dimnames(cov_lrho) <- list(free, free)
  }
  list(beta = cur$beta, se_beta = se_beta, p = cur$p, ll = cur$ll,
       ref = ref, panels = present, rho = rho_p, free = free,
       cov_lrho = cov_lrho, converged = conv, rank_ok = rank_ok,
       n_clusters = nc, n_parameters = np, cluster_ok = cluster_ok)
}

# Reconcile the per-set panel ratios into one set of panel units phi with
# geometric mean one, by generalised least squares on the panel graph. Each
# set contributes the observations log rho_{gs} = log phi_g - log phi_{ref(s)}
# for its free panels g, with the within-set covariance of those log-ratios
# carried across from stage 1; the observations are independent between sets,
# so the full observation covariance is block-diagonal by set. GLS on this
# gives both the precision-weighted point estimate (the reconciliation over
# the sets where a panel appears) and correctly correlated standard errors.
# Errors informatively when the panels are not connected through shared sets.
.btlef_reconcile_phi <- function(panels_u, blocks) {
  G <- length(panels_u)
  if (G == 1L)
    return(list(phi = setNames(1, panels_u),
                se_log_phi = setNames(NA_real_, panels_u),
                lphi = setNames(0, panels_u),
                cov_log_phi = matrix(0, 1L, 1L,
                  dimnames = list(panels_u, panels_u))))
  # flatten the per-set blocks into one observation vector, design and
  # block-diagonal covariance
  y <- numeric(0); pan <- ref <- character(0)
  Cov <- matrix(0, 0, 0)
  for (bk in blocks) {
    if (!length(bk$free)) next
    idx <- length(y) + seq_along(bk$free)
    y <- c(y, bk$lrho[bk$free]); pan <- c(pan, bk$free)
    ref <- c(ref, rep(bk$ref, length(bk$free)))
    cb <- bk$cov[bk$free, bk$free, drop = FALSE]
    cb <- cb + diag(1e-10, nrow(cb))                       # numerical floor
    Z <- matrix(0, nrow(Cov) + nrow(cb), ncol(Cov) + ncol(cb))
    if (nrow(Cov)) Z[seq_len(nrow(Cov)), seq_len(ncol(Cov))] <- Cov
    Z[nrow(Cov) + seq_len(nrow(cb)), ncol(Cov) + seq_len(ncol(cb))] <- cb
    Cov <- Z
  }
  if (!length(y))
    stop("the panels cannot be linked: no set contains comparisons from more ",
         "than one panel, so the panel units phi are unidentified")
  ei <- cbind(match(pan, panels_u), match(ref, panels_u))
  comp <- .btlef_components(G, ei)
  if (length(unique(comp)) > 1L)
    stop("the panel-by-set graph is not connected; panel units (phi) are ",
         "unidentified between: ",
         paste(tapply(panels_u, comp, paste, collapse = "+"), collapse = " | "))

  g0 <- panels_u[1]                                        # arbitrary anchor
  cols <- setdiff(panels_u, g0)
  X <- matrix(0, length(y), length(cols))
  for (r in seq_along(y)) {
    cp <- match(pan[r], cols); if (!is.na(cp)) X[r, cp] <- X[r, cp] + 1
    cr <- match(ref[r], cols); if (!is.na(cr)) X[r, cr] <- X[r, cr] - 1
  }
  W <- solve(Cov)                                          # GLS weight
  XtW <- t(X) %*% W
  covred <- solve(XtW %*% X)
  bred <- covred %*% (XtW %*% y)
  lphi <- setNames(numeric(G), panels_u); lphi[cols] <- bred
  cov_full <- matrix(0, G, G, dimnames = list(panels_u, panels_u))
  cov_full[cols, cols] <- covred
  A <- diag(G) - matrix(1 / G, G, G)                       # centre to geo-mean 1
  lphi_c <- as.numeric(A %*% lphi)
  cov_c <- A %*% cov_full %*% t(A)
  list(phi = setNames(exp(lphi_c), panels_u),
       se_log_phi = setNames(sqrt(pmax(diag(cov_c), 0)), panels_u),
       lphi = setNames(lphi_c, panels_u), cov_log_phi = cov_c)
}

# Stage 2: cross-set linking. With the frame locations beta and panel units
# phi held fixed, estimate (log alpha, kappa) for the non-reference sets by
# Newton on the cross-set comparison likelihood. Standard errors are the
# inverse observed information, conditional on stage 1 (the stage-1
# uncertainty is not propagated -- see the roxygen note).
.btlef_stage2 <- function(a, b, y, phg, sa, sb, bhat, sets_u, maxit, tol) {
  S <- length(sets_u); free <- sets_u[-1L]; nf <- S - 1L; np <- 2L * nf
  ba <- bhat[a]; bb <- bhat[b]
  fa <- match(sa, free); fb <- match(sb, free)             # NA on the reference set
  R <- length(y)

  eval_th <- function(th) {
    la <- th[seq_len(nf)]; kap <- th[nf + seq_len(nf)]
    alpha <- setNames(rep(1, S), sets_u); kappa <- setNames(rep(0, S), sets_u)
    alpha[free] <- exp(la); kappa[free] <- kap
    va <- alpha[sa] * ba + kappa[sa]; vb <- alpha[sb] * bb + kappa[sb]
    p <- plogis(phg * (va - vb))
    ll <- sum(ifelse(y == 1, log(pmax(p, 1e-300)), log(pmax(1 - p, 1e-300))))
    list(alpha = alpha, kappa = kappa, p = p, ll = ll)
  }
  # derivative design D (R x np): columns log alpha (free sets), then kappa
  design <- function(cur) {
    D <- matrix(0, R, np)
    for (j in seq_len(nf)) {
      s <- free[j]; aj <- cur$alpha[[s]]
      onA <- !is.na(fa) & fa == j; onB <- !is.na(fb) & fb == j
      D[, j] <- phg * (onA * aj * ba - onB * aj * bb)       # d eta / d log alpha_s
      D[, nf + j] <- phg * (onA - onB)                      # d eta / d kappa_s
    }
    D
  }

  solve_masked <- function(mask) {
    theta <- numeric(np); cur <- eval_th(theta)
    for (it in seq_len(maxit)) {
      Dm <- design(cur)[, mask, drop = FALSE]
      u <- y - cur$p; av <- cur$p * (1 - cur$p)
      g <- crossprod(Dm, u); Fi <- crossprod(Dm, Dm * av)
      sm <- tryCatch(solve(Fi, g),
                     error = function(e) solve(Fi + diag(1e-8, sum(mask)), g))
      step <- numeric(np); step[mask] <- as.numeric(sm)
      lam <- 1; moved <- FALSE
      for (half in 1:30) {
        cand <- theta + lam * step; c2 <- eval_th(cand)
        if (is.finite(c2$ll) && c2$ll >= cur$ll - 1e-12) {
          theta <- cand; cur <- c2; moved <- TRUE; break
        }
        lam <- lam / 2
      }
      if (!moved) break
      if (max(abs(lam * step)) < tol) break
    }
    cur
  }
  obs_info <- function(cur) {
    # observed information; the la-diagonal carries the curvature
    # d2 eta / d la^2
    D <- design(cur); u <- y - cur$p; av <- cur$p * (1 - cur$p)
    H <- crossprod(D, D * av)
    for (j in seq_len(nf)) {
      sj <- free[j]; aj <- cur$alpha[[sj]]
      onA <- !is.na(fa) & fa == j; onB <- !is.na(fb) & fb == j
      curv <- phg * (onA * aj * ba - onB * aj * bb)
      H[j, j] <- H[j, j] - sum(u * curv)
    }
    list(H = H, D = D, u = u)
  }

  cur <- solve_masked(rep(TRUE, np))
  oi <- obs_info(cur)

  # identification, classified by WHERE the information fails. A flat
  # direction confined to log-alpha columns is the degenerate-unit case:
  # the set's within-set locations are (near) indistinguishable, so its
  # unit has nothing to scale, but its origin kappa -- and with it the
  # placement of its objects -- is still identified. Refit with those
  # units fixed at the conventional 1 and report their alpha as NA. A
  # flat direction that loads on a kappa column means the set cannot be
  # PLACED at all: that is a structural failure of the cross-set design.
  alpha_unident <- setNames(rep(FALSE, nf), free)
  rank_ok <- TRUE; separated <- FALSE
  # Complete / quasi-complete separation of the cross-set comparisons:
  # every outcome is fitted at the boundary (p -> 0 or 1), so the
  # likelihood is unbounded and kappa / log-alpha run to the trust region
  # while the observed information stays finite -- the eigenvalue and rcond
  # checks below cannot see it (nothing is singular at the stopping point).
  # Detect it directly from the fitted probabilities: near-deterministic
  # prediction of essentially every cross-set outcome means the sets are
  # ordered by an unbounded margin and cannot be placed on a common scale.
  if (length(y) && mean(pmin(cur$p, 1 - cur$p) < 1e-4) > 0.99) {
    separated <- TRUE; rank_ok <- FALSE
  }
  eh <- eigen(oi$H, symmetric = TRUE)
  flat <- abs(eh$values) < max(abs(eh$values)) * 1e-10
  if (!separated && any(flat)) {
    V <- eh$vectors[, flat, drop = FALSE]
    la_load <- sqrt(rowSums(V[seq_len(nf), , drop = FALSE]^2))
    ka_load <- sqrt(rowSums(V[nf + seq_len(nf), , drop = FALSE]^2))
    if (any(ka_load > 1e-2)) rank_ok <- FALSE
    else {
      alpha_unident[la_load > 1e-2] <- TRUE
      mask <- rep(TRUE, np); mask[which(alpha_unident)] <- FALSE
      cur <- solve_masked(mask)
      oi <- obs_info(cur)
      kept <- which(mask)
      rc <- tryCatch(rcond(oi$H[kept, kept, drop = FALSE]),
                     error = function(e) 0)
      if (!(is.finite(rc) && rc > 1e-10)) rank_ok <- FALSE
    }
  }
  if (rank_ok) {
    rc <- tryCatch(rcond(oi$H), error = function(e) 0)
    if (!any(alpha_unident) && !(is.finite(rc) && rc > 1e-10))
      rank_ok <- FALSE
  }

  # covariance over the estimated parameters; fixed-by-convention units
  # carry zero rows (their uncertainty is not defined, and the reported
  # se_log_alpha is NA)
  kept <- which(!c(alpha_unident, rep(FALSE, nf)))
  cov <- matrix(0, np, np)
  cov[kept, kept] <- tryCatch(solve(oi$H[kept, kept, drop = FALSE]),
    error = function(e) {
      Dk <- oi$D[, kept, drop = FALSE]
      av <- cur$p * (1 - cur$p)
      solve(crossprod(Dk, Dk * av) + diag(1e-8, length(kept)))
    })
  # scale-free per-comparison gradient criterion (see .btlef_stage1)
  conv <- isTRUE(max(abs(crossprod(oi$D[, kept, drop = FALSE], oi$u))) <
                   1e-6 * length(y))
  se <- sqrt(pmax(diag(cov), 0))
  if (!rank_ok) se[] <- NA_real_
  alpha_rep <- cur$alpha
  alpha_rep[names(which(alpha_unident))] <- NA_real_
  se_la <- setNames(se[seq_len(nf)], free)
  se_la[alpha_unident] <- NA_real_
  list(alpha = alpha_rep, alpha_use = cur$alpha, kappa = cur$kappa,
       p = cur$p, ll = cur$ll, cov = cov,
       se_log_alpha = se_la,
       se_kappa = setNames(se[nf + seq_len(nf)], free),
       free = free, converged = conv, rank_ok = rank_ok,
       separated = separated, alpha_unident = alpha_unident)
}

# pooled log-of-mean-square fit residual over a set of comparisons, using the
# frame model's fitted probabilities (the paired-comparison residual logic of
# btl(): z = (y - p) / sqrt(p(1-p)), Andrich and Marais 2019 ch. 23)
.btlef_frame_fit <- function(y, p, f_cell = 1) {
  n <- length(y)
  if (n < 3L) return(list(fit_resid = NA_real_, df = NA_real_))
  V <- pmax(p * (1 - p), 1e-12)
  z2 <- (y - p)^2 / V
  mu4 <- p * (1 - p)^4 + (1 - p) * p^4
  c4v <- mu4 / V^2 - 1
  y2 <- sum(z2); f <- f_cell * n; v <- sum(c4v)
  list(fit_resid = if (v > 1e-8 && y2 > 0 && f > 0)
         f * (log(y2) - log(f)) / sqrt(v) else NA_real_, df = f)
}

# Core diagnostics on the fitted frame probabilities. These are calculated
# here, rather than borrowed from a single-unit BTL fit, so object and judge
# fit update when the frame structure updates the common-scale locations.
.btlef_diagnostics <- function(a, b, y, jd, pan, sa, sb, p, phi, alpha,
                               objects, n_parameters) {
  objs <- objects$object
  K <- length(objs)
  ia <- match(a, objs); ib <- match(b, objs)
  p <- pmin(pmax(as.numeric(p), 1e-12), 1 - 1e-12)
  V <- p * (1 - p)
  z2 <- (y - p)^2 / V
  mu4 <- p * (1 - p)^4 + (1 - p) * p^4
  c4v <- mu4 / V^2 - 1
  f_cell <- max((length(y) - n_parameters) / length(y), 0)
  pool <- function(sel) {
    n <- sum(sel)
    if (n < 3L || f_cell <= 0)
      return(list(infit_ms = NA_real_, outfit_ms = NA_real_,
                  fit_resid = NA_real_, df = NA_real_, n = n))
    y2 <- sum(z2[sel]); f <- f_cell * n
    wv <- sum(V[sel])
    infit <- if (wv > 1e-12) sum(z2[sel] * V[sel]) / (f_cell * wv)
             else NA_real_
    vv <- sum(c4v[sel])
    fr <- if (vv > 1e-8 && y2 > 0 && f > 0)
      f * (log(y2) - log(f)) / sqrt(vv) else NA_real_
    list(infit_ms = infit, outfit_ms = y2 / f, fit_resid = fr,
         df = f, n = n)
  }

  ofit <- lapply(seq_len(K), function(k) pool(ia == k | ib == k))
  wins <- vapply(seq_len(K), function(k)
    sum(y[ia == k]) + sum(1 - y[ib == k]), 0)
  objects$comparisons <- vapply(ofit, `[[`, 0, "n")
  objects$wins <- wins
  objects$infit_ms <- vapply(ofit, `[[`, 0, "infit_ms")
  objects$outfit_ms <- vapply(ofit, `[[`, 0, "outfit_ms")
  objects$fit_resid <- vapply(ofit, `[[`, 0, "fit_resid")
  objects$df_fit <- vapply(ofit, `[[`, 0, "df")

  ju <- sort(unique(jd))
  jfit <- lapply(ju, function(j) pool(jd == j))
  judges <- data.frame(
    judge = ju, n = vapply(jfit, `[[`, 0, "n"),
    infit_ms = vapply(jfit, `[[`, 0, "infit_ms"),
    outfit_ms = vapply(jfit, `[[`, 0, "outfit_ms"),
    fit_resid = vapply(jfit, `[[`, 0, "fit_resid"),
    df_fit = vapply(jfit, `[[`, 0, "df"), stringsAsFactors = FALSE)

  lo_first <- ia < ib
  key <- paste(pmin(ia, ib), pmax(ia, ib))
  obs_lo <- ifelse(lo_first, y, 1 - y)
  exp_lo <- ifelse(lo_first, p, 1 - p)
  n_pair <- tapply(rep(1, length(y)), key, sum)
  obs_m <- tapply(obs_lo, key, mean)
  exp_m <- tapply(exp_lo, key, mean)
  v_pair <- tapply(V, key, sum)
  zp <- tapply(obs_lo - exp_lo, key, sum) / sqrt(pmax(v_pair, 1e-12))
  idx <- do.call(rbind, strsplit(names(n_pair), " "))
  pairs <- data.frame(
    object_a = objs[as.integer(idx[, 1])],
    object_b = objs[as.integer(idx[, 2])], n = as.numeric(n_pair),
    obs_prop = as.numeric(obs_m), exp_prop = as.numeric(exp_m),
    residual = as.numeric(zp), chisq = as.numeric(zp)^2,
    stringsAsFactors = FALSE)
  used <- pairs$n >= 2L
  total_df <- sum(used) - n_parameters
  total_chisq <- if (total_df >= 1L) sum(pairs$chisq[used]) else NA_real_
  if (total_df < 1L) total_df <- NA_integer_

  # d eta / d(v_a - v_b): within-set comparisons use phi/alpha because
  # v = alpha * beta + kappa; cross-set comparisons use phi directly.
  slope <- unname(phi[pan])
  within <- sa == sb
  slope[within] <- slope[within] / unname(alpha[sa[within]])
  comparisons <- data.frame(
    object_a = a, object_b = b, response = y, weight = 1,
    judge = jd, panel = pan, set_a = sa, set_b = sb,
    expected = p, frame_slope = slope,
    information = slope^2 * V, stringsAsFactors = FALSE)

  list(objects = objects, judges = judges, pairs = pairs,
       comparisons = comparisons, total_chisq = total_chisq,
       total_df = total_df,
       total_p = if (is.finite(total_chisq))
         stats::pchisq(total_chisq, total_df, lower.tail = FALSE) else NA_real_)
}

#' Fit the extended frame of reference model for paired comparisons
#'
#' Fits paired comparisons when judges belong to panels and objects belong to
#' linked sets whose units or origins can differ. It combines the
#' Bradley--Terry--Luce model with Humphry's extended frame of reference
#' structure.
#'
#' @details
#' For object \eqn{k} in set \eqn{s}, let
#' \deqn{v_k=\alpha_s\beta_k+\kappa_s,}
#' where \eqn{\beta_k} is its within-set location, \eqn{\alpha_s>0} is the set
#' unit, and \eqn{\kappa_s} is the set origin. A comparison in panel \eqn{g}
#' has logit
#' \deqn{\phi_g(\beta_a-\beta_b)}
#' for objects in the same set, and
#' \deqn{\phi_g(v_a-v_b)}
#' for objects in different sets. Cross-set comparisons identify the common
#' scale. The first set fixes \eqn{\alpha=1} and \eqn{\kappa=0}; panel units
#' have geometric mean one.
#'
#' Estimation has two stages. Within-set comparisons estimate object locations
#' and panel-unit ratios. Weighted least squares reconciles the ratios over the
#' panel-by-set linking graph. Cross-set comparisons then estimate the set
#' units and origins. Unlike the person-by-item EFRM, this linking step uses
#' only comparison outcomes and does not require a distribution of persons.
#' The paired-comparison form is an extension of Humphry's model implemented in
#' this package.
#'
#' The default judge bootstrap resamples judges within panels and refits both
#' stages. The parametric bootstrap draws independent outcomes from the fitted
#' probabilities and uses normal and chi-square reference distributions.
#' \code{se_method = "conditional"} uses analytic stage-one
#' errors and conditions the linking errors on stage one; it is intended for
#' preliminary inspection. Its unit probabilities and omnibus tests are
#' withheld because it does not propagate stage-one uncertainty. Bootstrap
#' failures and boundary estimates are reported in \code{notes}.
#'
#' With one set, the model contains panel units only. With one set and one
#' panel, it reduces to \code{\link{btl}}. Omnibus Wald tests provide inference
#' for the unit families; individual contrasts are Holm-adjusted follow-ups.
#' Judge-bootstrap probabilities require at least six judges and 5.5 effective
#' judges in every contributing panel, and eight of each on a set link. The
#' support is returned in \code{unit_support}; estimates remain descriptive
#' when a probability is withheld.
#'
#' @param data A data frame with one comparison per row.
#' @param object_a,object_b Names of the columns holding the two compared
#'   objects.
#' @param winner Name of the winner column. A value must match one of the two
#'   objects in that row. \code{"tie"} and \code{"draw"} mark ties; other
#'   values are treated as missing.
#' @param judge Name of the judge column (clusters the stage-one standard
#'   errors and defines the panels when \code{panels} is a judge attribute).
#' @param panels Either the name of a judge-attribute column in \code{data} or
#'   a named vector mapping judge to panel.
#' @param object_sets A named list mapping set names to character vectors of
#'   object names; every compared object must belong to exactly one set.
#' @param response Not supported: this first implementation fits dichotomous
#'   winner data only. Supplying it raises an informative error.
#' @param ties \code{"drop"} (default, removed with a note) or \code{"error"}.
#' @param min_link Minimum number of cross-set comparisons a set pair must
#'   supply to be used for linking; sets not reachable from the reference set
#'   through sufficient cross-set pairs raise an error.
#' @param se_method Method used for standard errors. The default,
#'   \code{"judge_bootstrap"}, resamples judges within panels and retains
#'   dependence among a judge's comparisons. \code{"bootstrap"} instead
#'   draws independent outcomes from fitted probabilities. Both stages are
#'   refitted. \code{"conditional"} uses
#'   analytic stage-one standard errors for \code{beta} and \code{phi}, and
#'   inverse observed information for \code{alpha} and \code{kappa}
#'   conditional on the stage-one estimates. It is faster, but does not
#'   propagate stage-one uncertainty into the linking parameters; unit
#'   probabilities and omnibus tests are therefore withheld.
#' @param boot_reps Number of replicates for \code{se_method = "bootstrap"}
#'   or \code{"judge_bootstrap"}; at least 30 are required.
#' @param workers Number of judge-bootstrap workers. The default is four,
#'   reduced when the system limit is lower. The parametric bootstrap remains
#'   serial because its refits are inexpensive.
#' @param seed Optional bootstrap seed. The caller's random-number state is
#'   restored when estimation finishes.
#' @param progress Optional function called as \code{progress(stage, current,
#'   total)} during estimation.
#' @param cancel Optional zero-argument function checked between bootstrap
#'   batches. Returning \code{TRUE} stops with a \code{rasch_cancelled}
#'   condition.
#' @param maxit,tol Newton iteration cap and convergence tolerance.
#' @return An object of class \code{"rasch_btl_efrm"}. It contains the object
#'   estimates, group- and set-unit tables, origin shifts, omnibus unit tests,
#'   unit-specific judge support, frame definitions, convergence information,
#'   and analysis notes.
#' @references Andrich, D. (1978). Relationships between the Thurstone and
#'   Rasch approaches to item scaling. Applied Psychological Measurement,
#'   2(3), 451--462.
#'
#'   Bradley, R. A. and Terry, M. E. (1952). Rank analysis of
#'   incomplete block designs: I. The method of paired comparisons.
#'   Biometrika, 39, 324--345.
#'
#'   David, H. A. (1988). The Method of Paired Comparisons (2nd ed.). Griffin.
#'
#'   Humphry, S. M. (2005). Maintaining a common arbitrary unit in social
#'   measurement. PhD thesis, Murdoch University.
#'
#'   Humphry, S. M. (2012). Item set discrimination and the unit in the
#'   Rasch model. Journal of Applied Measurement, 13(2), 165--180.
#'
#'   Humphry, S. M. and Andrich, D. (2008). Understanding the unit in the
#'   Rasch model. Journal of Applied Measurement, 9(3), 249--264.
#'
#'   Luce, R. D. (1959). Individual Choice Behavior: A Theoretical
#'   Analysis. Wiley.
#'
#'   Thurstone, L. L. (1927). A law of comparative judgment. Psychological
#'   Review, 34, 273--286.
#' @seealso \code{\link{btl}}, \code{\link{rasch_efrm}},
#'   \code{\link{plot_btl_units}}, and \code{\link{simulate_btl_efrm}}.
#' @examples
#' \donttest{
#' d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
#'                        set_units = c(1, 1.4), set_origins = c(0, 0.8),
#'                        seed = 1)
#' fit <- btl_efrm(d, "object_a", "object_b", winner = "winner",
#'                 judge = "judge", panels = "panel",
#'                 object_sets = attr(d, "truth")$object_sets,
#'                 se_method = "conditional")
#' fit$alpha_table
#' }
#' @export
btl_efrm <- function(data, object_a, object_b, winner, judge, panels,
                     object_sets, response = NULL,
                     ties = c("drop", "error"), min_link = 20,
                     se_method = c("judge_bootstrap", "bootstrap",
                                   "conditional"),
                     boot_reps = 200, workers = 4L, seed = NULL,
                     progress = NULL, cancel = NULL,
                     maxit = 60, tol = 1e-8) {
  .check_column_names(data)
  ties <- match.arg(ties)
  se_method <- match.arg(se_method)
  if (length(boot_reps) != 1L || !is.finite(boot_reps) || boot_reps < 0L ||
      boot_reps != floor(boot_reps))
    stop("boot_reps must be one non-negative whole number")
  boot_reps <- as.integer(boot_reps)
  if (se_method %in% c("bootstrap", "judge_bootstrap") && boot_reps < 30L)
    stop("BTL-EFRM bootstrap inference needs at least 30 replicates")
  if (length(workers) != 1L || !is.finite(workers) || workers < 1L ||
      workers != floor(workers))
    stop("workers must be one positive whole number")
  workers <- as.integer(workers)
  workers <- if (se_method == "judge_bootstrap")
    min(workers, .rasch_available_workers(), boot_reps) else 1L
  if (workers > 1L &&
      !file.exists(system.file("DESCRIPTION", package = "rasch")))
    stop("parallel BTL-EFRM workers require an installed package; install ",
         "rasch before using workers above one")
  if (!is.null(progress) && !is.function(progress))
    stop("progress must be NULL or a function")
  if (!is.null(cancel) && !is.function(cancel))
    stop("cancel must be NULL or a function")
  if (!is.null(seed)) {
    if (length(seed) != 1L || !is.finite(seed) || seed < 0 ||
        seed != floor(seed) || seed > .Machine$integer.max)
      stop("seed must be NULL or one non-negative whole number within the integer range")
    seed <- as.integer(seed)
    old_seed <- .sim_seed_capture()
    on.exit(.sim_seed_restore(old_seed), add = TRUE)
    set.seed(seed)
  }
  report <- function(stage, current, total) {
    if (is.function(cancel) && isTRUE(cancel()))
      .rasch_cancelled("BTL-EFRM estimation")
    if (is.function(progress)) progress(stage, current, total)
    invisible(NULL)
  }
  if (!is.null(response))
    stop("btl_efrm fits dichotomous winner data only in this first ",
         "implementation; a polytomous `response` is not supported. Reduce the ",
         "polytomous margins to a winner, or use btl() for a single-frame polytomous ",
         "analysis.")
  data <- as.data.frame(data)
  for (col in c(object_a, object_b, winner, judge))
    if (!col %in% names(data)) stop("column not found: ", col)
  a <- trimws(as.character(data[[object_a]]))
  b <- trimws(as.character(data[[object_b]]))
  wn <- trimws(as.character(data[[winner]]))
  jd <- as.character(data[[judge]])

  # panels: a judge-attribute column, or a named judge -> panel vector
  if (length(panels) == 1L && is.character(panels) && panels %in% names(data)) {
    pan <- as.character(data[[panels]])
    # a panel is a judge attribute: one judge in two panels is a data error
    # (and the judge bootstrap would silently reclassify their rows)
    npan <- tapply(pan, jd, function(x) length(unique(x[!is.na(x)])))
    if (any(npan > 1L, na.rm = TRUE))
      stop("judge(s) assigned to more than one panel: ",
           paste(names(npan)[npan > 1L], collapse = ", "),
           "; a panel is a judge attribute and must be constant per judge")
  } else if (!is.null(names(panels)) && all(nzchar(names(panels)))) {
    pan <- unname(as.character(panels)[match(jd, names(panels))])
  } else {
    stop("`panels` must name a column of `data` or be a named vector ",
         "mapping judge to panel")
  }
  notes <- character(0)

  keep <- !is.na(a) & !is.na(b) & !is.na(wn) & !is.na(jd) & !is.na(pan) & a != b
  if (any(!keep)) {
    notes <- c(notes, sprintf("%d row(s) dropped (missing or self-comparison)",
                              sum(!keep)))
    a <- a[keep]; b <- b[keep]; wn <- wn[keep]; jd <- jd[keep]; pan <- pan[keep]
  }
  if (!length(a)) stop("no usable comparisons")

  y <- ifelse(wn == a, 1L, ifelse(wn == b, 0L, NA_integer_))
  is_tie <- is.na(y) & tolower(wn) %in% c("tie", "draw")
  miss <- is.na(y) & !is_tie
  if (any(miss)) {
    notes <- c(notes, sprintf(
      "%d row(s) with winner matching neither object treated as missing", sum(miss)))
    sel <- !miss
    a <- a[sel]; b <- b[sel]; y <- y[sel]; jd <- jd[sel]; pan <- pan[sel]
  }
  if (anyNA(y)) {
    nt <- sum(is.na(y))
    if (ties == "error") stop(nt, " tie(s) present; set ties = 'drop'")
    notes <- c(notes, sprintf("%d tie(s) dropped", nt))
    sel <- !is.na(y)
    a <- a[sel]; b <- b[sel]; y <- y[sel]; jd <- jd[sel]; pan <- pan[sel]
  }
  if (!length(a)) stop("no usable comparisons after cleaning")

  # --- object sets ----------------------------------------------------------
  if (!is.list(object_sets) || is.null(names(object_sets)) ||
      any(!nzchar(names(object_sets))))
    stop("`object_sets` must be a named list: set name -> object names")
  objs_all <- sort(unique(c(a, b)))
  set_of <- setNames(rep(NA_character_, length(objs_all)), objs_all)
  multi <- character(0)
  for (s in names(object_sets)) {
    hit <- intersect(as.character(object_sets[[s]]), objs_all)
    for (o in hit) {
      if (!is.na(set_of[o])) multi <- c(multi, o)
      set_of[o] <- s
    }
  }
  if (length(multi))
    stop("object(s) assigned to more than one set: ",
         paste(unique(multi), collapse = ", "))
  if (anyNA(set_of))
    stop("object(s) in the data not found in `object_sets` (every compared ",
         "object must belong to exactly one set): ",
         paste(objs_all[is.na(set_of)], collapse = ", "))
  sets_u <- sort(unique(set_of)); S <- length(sets_u)
  panels_u <- sort(unique(pan)); G <- length(panels_u)
  sa <- set_of[a]; sb <- set_of[b]
  within <- sa == sb

  if (any(table(set_of) < 2L))
    stop("every set needs at least two objects; offending set(s): ",
         paste(names(which(table(set_of) < 2L)), collapse = ", "))

  # --- structural checks that do not depend on the outcomes -----------------
  # panel linkage first (its failure is the more fundamental design flaw):
  # panels are connected when they share within-set comparisons in some set
  if (G > 1L) {
    pe <- unique(do.call(rbind, lapply(sets_u, function(s) {
      pg <- unique(pan[within & sa == s])
      if (length(pg) < 2L) return(NULL)
      t(combn(match(pg, panels_u), 2))
    })))
    pcomp <- .btlef_components(G, if (is.null(pe)) matrix(numeric(0), 0, 2) else pe)
    if (any(pcomp != pcomp[1]))
      stop("the panels cannot be linked: no set contains comparisons from ",
           "more than one panel connecting panel group(s) {",
           paste(panels_u[pcomp != pcomp[1]], collapse = ", "),
           "} to the rest; panel units phi are unidentified across them")
  }
  cross <- which(!within)
  n_cross <- data.frame(set_a = character(0), set_b = character(0),
                        n = integer(0), stringsAsFactors = FALSE)
  if (S > 1L) {
    if (!length(cross))
      stop("no cross-set comparisons: the sets cannot be linked to a common ",
           "scale (set units alpha and origins kappa are unidentified)")
    isa <- match(sa[cross], sets_u); isb <- match(sb[cross], sets_u)
    key <- paste(pmin(isa, isb), pmax(isa, isb))
    tab <- table(key); parts <- do.call(rbind, strsplit(names(tab), " ", fixed = TRUE))
    n_cross <- data.frame(set_a = sets_u[as.integer(parts[, 1])],
                          set_b = sets_u[as.integer(parts[, 2])],
                          n = as.integer(tab), stringsAsFactors = FALSE)
    rownames(n_cross) <- NULL
    used <- n_cross$n >= min_link
    edges <- cbind(match(n_cross$set_a[used], sets_u),
                   match(n_cross$set_b[used], sets_u))
    comp <- .btlef_components(S, edges)
    ref_comp <- comp[1]                                 # reference set = sets_u[1]
    if (any(comp != ref_comp))
      stop("set(s) not reachable from the reference set '", sets_u[1],
           "' through cross-set pairs with at least min_link = ", min_link,
           " comparisons: ", paste(sets_u[comp != ref_comp], collapse = ", "),
           " (increase the cross-set data or lower min_link)")
    # a set's UNIT alpha is identified by how its internal spread shows in
    # cross-set outcomes: cross-set comparisons touching only one of its
    # objects identify the origin kappa but leave alpha riding on nothing
    for (s in sets_u[-1]) {
      touched <- unique(c(a[cross][sa[cross] == s], b[cross][sb[cross] == s]))
      if (length(touched) < 2L)
        stop("cross-set comparisons touch only ", length(touched),
             " object(s) of set '", s, "': its unit (alpha) is ",
             "unidentified -- add cross-set comparisons involving at ",
             "least two of its objects")
    }
  }
  # within each set, the object comparison graph must be connected, or the
  # relative locations inside the set are unidentified (the stage-1 Newton
  # would land wherever the ridge sends it, exactly as in pcml())
  for (s in sets_u) {
    rows_s <- which(within & sa == s)
    os_s <- sort(names(set_of)[set_of == s])
    ia_s <- match(a[rows_s], os_s); ib_s <- match(b[rows_s], os_s)
    comp_s <- .btlef_components(length(os_s), cbind(ia_s, ib_s))
    if (length(unique(comp_s)) > 1L)
      stop("the within-set comparison graph of set '", s, "' is not ",
           "connected: relative locations are unidentified between {",
           paste(tapply(os_s, comp_s, paste, collapse = ", "),
                 collapse = "} and {"), "}")
  }

  # --- the two-stage estimator, callable on any outcome vector ---------------
  # (one function for the observed data and for every bootstrap replicate, so
  # the resampled pipeline is identical to the reported one)
  fit_once <- function(yy) {
    bhat <- setNames(rep(NA_real_, length(objs_all)), objs_all)
    se_bhat <- bhat
    ref_of_set <- setNames(rep(NA_character_, S), sets_u)
    within_p <- rep(NA_real_, length(a))              # frame-model fitted p
    blocks <- list()                                  # per-set panel-ratio blocks
    ll_within <- 0; s1_conv <- TRUE
    s1 <- list(); dropped <- character(0)
    for (s in sets_u) {
      rows <- which(within & sa == s)
      os <- sort(names(set_of)[set_of == s]); Ks <- length(os)
      ia <- match(a[rows], os); ib <- match(b[rows], os)
      if (anyNA(ia) || anyNA(ib) || any(is.na(match(os, unique(c(a[rows], b[rows]))))))
        stop("set '", s, "' has object(s) with no within-set comparison; ",
             "each object needs at least one comparison inside its own set")
      fit1 <- .btlef_stage1(ia, ib, yy[rows], pan[rows], jd[rows], Ks, maxit, tol)
      s1[[s]] <- list(fit = fit1, rows = rows, os = os, ia = ia, ib = ib)
      if (length(fit1$free) && !isTRUE(fit1$cluster_ok))
        stop("set '", s, "' has ", fit1$n_clusters,
             " independent judge clusters for ", fit1$n_parameters,
             " stage-one parameters. Estimating and precision-weighting its ",
             "panel-unit ratio needs at least 10 judges and more judges than ",
             "parameters; add independent judges or simplify the frame design")
      # A set carries information about the panel-unit ratios only through its
      # internal separation: when its within-set contests are all near-even
      # (or one-sided), the ratio log rho_gs is the quotient of two near-zero
      # (or unbounded) logits and its estimate diverges with a spurious
      # covariance. Screen such blocks out of the phi reconciliation rather
      # than let them poison it; the set is refit below at the reconciled
      # units, which the frame model says apply to it regardless.
      lr <- log(fit1$rho[fit1$free])
      usable <- isTRUE(fit1$converged) && isTRUE(fit1$rank_ok) &&
        (!length(fit1$free) ||
           (all(is.finite(lr)) && max(abs(lr)) < 4 &&
            all(is.finite(fit1$cov_lrho)) && all(diag(fit1$cov_lrho) > 0)))
      if (!usable && length(fit1$free)) { dropped <- c(dropped, s); next }
      if (!isTRUE(fit1$rank_ok))
        stop("the within-set information of set '", s, "' is singular: ",
             "its object locations are unidentified even with the panel ",
             "units fixed -- the within-set comparisons do not span the ",
             "objects")
      blocks[[s]] <- list(ref = fit1$ref, free = fit1$free,
                          lrho = setNames(lr, fit1$free),
                          cov = fit1$cov_lrho)
    }
    rec <- tryCatch(.btlef_reconcile_phi(panels_u, blocks), error = function(e) {
      if (length(dropped))
        stop(conditionMessage(e), "\n  (set(s) ", paste(dropped, collapse = ", "),
             " were excluded from the panel-unit reconciliation because their ",
             "within-set comparisons carry no stable panel-ratio information ",
             "-- their contests are near-even or one-sided, or their panels ",
             "observe disjoint object pairs, leaving the ratio ",
             "rank-deficient)", call. = FALSE)
      stop(e)
    })
    phi <- rec$phi
    for (s in sets_u) {
      rows <- s1[[s]]$rows; os <- s1[[s]]$os
      if (s %in% dropped) {
        # refit the set's locations with the panel units held at the
        # reconciled phi: beta comes out directly on the common scale
        fit1 <- .btlef_stage1(s1[[s]]$ia, s1[[s]]$ib, yy[rows], pan[rows],
                              jd[rows], length(os), maxit, tol,
                              rho_fixed = phi)
        if (!isTRUE(fit1$rank_ok))
          stop("set '", s, "': object locations are unidentified even ",
               "with the panel units held fixed -- the within-set ",
               "comparisons do not span the objects")
        s1[[s]]$fit <- fit1
        bhat[os] <- fit1$beta; se_bhat[os] <- fit1$se_beta
      } else {
        fit1 <- s1[[s]]$fit
        pr <- phi[[fit1$ref]]
        bhat[os] <- fit1$beta / pr; se_bhat[os] <- fit1$se_beta / pr
      }
      within_p[rows] <- fit1$p
      ref_of_set[s] <- fit1$ref
      ll_within <- ll_within + fit1$ll
      s1_conv <- s1_conv && isTRUE(fit1$converged)
    }
    # within-set fitted p on the common scale: logit = phi_g (bhat_a - bhat_b)
    within_p[within] <- plogis(phi[pan[within]] * (bhat[a[within]] - bhat[b[within]]))

    alpha <- setNames(rep(1, S), sets_u); kappa <- setNames(rep(0, S), sets_u)
    alpha_use <- alpha
    s2_alpha_unident <- setNames(logical(0), character(0))
    se_log_alpha <- setNames(rep(NA_real_, S), sets_u)
    se_kappa <- setNames(rep(NA_real_, S), sets_u)
    cov2 <- NULL; s2_conv <- TRUE; ll_cross <- 0; s2_rank_ok <- TRUE
    s2_separated <- FALSE
    p_all <- within_p
    if (S > 1L) {
      st2 <- .btlef_stage2(a[cross], b[cross], yy[cross], phi[pan[cross]],
                           sa[cross], sb[cross], bhat, sets_u, maxit, tol)
      alpha <- st2$alpha; alpha_use <- st2$alpha_use
      kappa <- st2$kappa; cov2 <- st2$cov
      s2_alpha_unident <- st2$alpha_unident
      se_log_alpha[st2$free] <- st2$se_log_alpha
      se_kappa[st2$free] <- st2$se_kappa
      s2_conv <- st2$converged && st2$rank_ok; ll_cross <- st2$ll
      s2_rank_ok <- st2$rank_ok; s2_separated <- st2$separated
      p_all[cross] <- st2$p
    }
    v <- alpha_use[set_of[objs_all]] * bhat[objs_all] +
      kappa[set_of[objs_all]]
    list(bhat = bhat, se_bhat = se_bhat, phi = phi,
         se_log_phi = rec$se_log_phi, cov_log_phi = rec$cov_log_phi,
         ref_of_set = ref_of_set,
         alpha = alpha, alpha_use = alpha_use, kappa = kappa, cov2 = cov2,
         s2_alpha_unident = s2_alpha_unident,
         se_log_alpha = se_log_alpha, se_kappa = se_kappa, v = v,
         within_p = within_p, p_all = p_all,
         ll_within = ll_within, ll_cross = ll_cross,
         dropped = dropped, s2_rank_ok = s2_rank_ok,
         s2_separated = s2_separated,
         converged = isTRUE(s1_conv && s2_conv))
  }

  report("two-stage fit", 0L, 1L)
  fit0 <- fit_once(y)
  report("two-stage fit", 1L, 1L)
  if (S > 1L && isTRUE(fit0$s2_separated))
    stop("the cross-set comparisons are (quasi-)completely separated: one ",
         "set beats the other in essentially every cross-set comparison, so ",
         "the sets are ordered by an unbounded margin and cannot be placed ",
         "on one scale (the set units alpha and origins kappa have no finite ",
         "estimate). Collect cross-set comparisons that some objects of the ",
         "weaker set sometimes win", call. = FALSE)
  if (S > 1L && !isTRUE(fit0$s2_rank_ok))
    stop("the cross-set information matrix is singular or ill-conditioned: ",
         "the cross-set comparisons cannot place the sets on one scale ",
         "(a flat direction of the information loads on a set origin ",
         "kappa) -- add cross-set comparisons that link every set")
  if (!isTRUE(fit0$converged))
    warning("BTL-EFRM estimation did NOT converge; increase maxit or inspect ",
            "the within- and cross-set comparison design", call. = FALSE)
  if (any(fit0$s2_alpha_unident))
    notes <- c(notes, paste0(
      "set unit(s) unidentified for ",
      paste(names(which(fit0$s2_alpha_unident)), collapse = ", "),
      ": their within-set locations are indistinguishable, so the unit ",
      "has nothing to scale; alpha is reported NA (fixed at 1 by ",
      "convention in the linked values) and the objects are placed ",
      "through the origin kappa alone"))
  bhat <- fit0$bhat; se_bhat <- fit0$se_bhat
  phi <- fit0$phi; ref_of_set <- fit0$ref_of_set
  alpha <- fit0$alpha; alpha_use <- fit0$alpha_use
  kappa <- fit0$kappa; cov2 <- fit0$cov2
  se_log_phi <- fit0$se_log_phi
  cov_log_phi <- fit0$cov_log_phi
  se_log_alpha <- fit0$se_log_alpha; se_kappa <- fit0$se_kappa
  v <- fit0$v; within_p <- fit0$within_p
  ll_within <- fit0$ll_within; ll_cross <- fit0$ll_cross

  # conditional (analytic) delta-method errors for the common-scale values
  se_v <- se_bhat[objs_all]                             # reference set: v = beta
  free <- sets_u[-1L]
  if (S > 1L) for (o in objs_all) {
    s <- set_of[[o]]; if (s == sets_u[1]) next
    j <- match(s, free); idx <- c(j, (S - 1L) + j)
    C2 <- cov2[idx, idx, drop = FALSE]
    gvec <- c(alpha_use[[s]] * bhat[[o]], 1)            # d v / d(log alpha, kappa)
    var_link <- drop(t(gvec) %*% C2 %*% gvec)
    se_v[[o]] <- sqrt(pmax(alpha_use[[s]]^2 * se_bhat[[o]]^2 + var_link, 0))
  }

  # --- bootstrap inference ---------------------------------------------------
  # Both bootstrap modes refit the full pipeline, carrying stage-one
  # uncertainty into the linking. Judge resampling is the default; the
  # parametric outcome bootstrap remains a model-based sensitivity analysis.
  boot_fail <- 0L
  cov_v <- NULL
  cov_draws <- function(D, idx) {
    V <- matrix(NA_real_, length(idx), length(idx))
    ok <- vapply(idx, function(j) all(is.finite(D[, j])), TRUE)
    if (any(ok))
      V[ok, ok] <- stats::cov(D[, idx[ok], drop = FALSE])
    V
  }
  if (se_method == "judge_bootstrap") {
    # resample JUDGES with replacement within each panel (the panel design
    # is fixed; judges are the sampling units), relabel the copies so
    # clusters stay distinct, and rerun the whole pipeline on each
    # resample: unlike the parametric bootstrap, this carries any
    # extra-model dependence within a judge's comparisons
    jd_rows <- split(seq_along(a), jd)
    pan_of_judge <- vapply(jd_rows, function(r) pan[r[1]], "")
    judges_by_panel <- split(names(jd_rows), pan_of_judge)
    nj_min <- min(lengths(judges_by_panel))
    if (nj_min < 2L)
      stop("judge bootstrap needs at least 2 judges in every panel ",
           "(resampling a single judge returns the same data every time, ",
           "so its SEs would be a spurious zero); use se_method = ",
           "'bootstrap' or 'conditional'")
    if (nj_min < 5L)
      notes <- c(notes, sprintf(
        "judge bootstrap: smallest panel has only %d judges; resampling so few is unstable and the SEs are rough", nj_min))
    report("judge bootstrap", 0L, boot_reps)
    takes <- vector("list", boot_reps)
    for (bb in seq_len(boot_reps)) {
      if (is.function(cancel) && isTRUE(cancel()))
        .rasch_cancelled("BTL-EFRM estimation")
      takes[[bb]] <- unlist(lapply(judges_by_panel, function(js)
        sample(js, length(js), replace = TRUE)), use.names = FALSE)
    }
    exp_ok <- is.finite(c(log(phi), log(alpha), kappa))
    one_judge_boot <- function(bb) {
      take <- takes[[bb]]
      idx <- unlist(jd_rows[take], use.names = FALSE)
      jd_new <- rep(paste0(take, "#", seq_along(take)),
                    lengths(jd_rows[take]))
      df_b <- data.frame(oa = a[idx], ob = b[idx],
                         win = ifelse(y[idx] == 1L, a[idx], b[idx]),
                         judge = jd_new, stringsAsFactors = FALSE)
      pmap_b <- setNames(pan[idx][!duplicated(jd_new)], unique(jd_new))
      fb <- tryCatch(suppressWarnings(
        utils::getFromNamespace("btl_efrm", "rasch")(
          df_b, "oa", "ob", "win", "judge", panels = pmap_b,
          object_sets = object_sets, ties = "drop",
          min_link = min_link, se_method = "conditional", workers = 1L,
          maxit = maxit, tol = tol)), error = function(e) NULL)
      if (is.null(fb) || !isTRUE(fb$converged)) return(NULL)
      lphi_b <- log(fb$phi_table$phi)[match(panels_u, fb$phi_table$panel)]
      la_b <- log(fb$alpha_table$alpha)[match(sets_u, fb$alpha_table$set)]
      ka_b <- fb$kappa_table$kappa[match(sets_u, fb$kappa_table$set)]
      bh_b <- fb$objects$beta_set[match(objs_all, fb$objects$object)]
      v_b <- fb$objects$v[match(objs_all, fb$objects$object)]
      # a unit the OBSERVED fit already reports NA (unidentified) is NA in
      # replicates too; only unexpected NAs mark a failed replicate
      if (anyNA(c(lphi_b, la_b, ka_b)[exp_ok])) return(NULL)
      c(lphi_b, la_b, ka_b, bh_b, v_b)
    }
    ans <- .rasch_boot_apply(
      boot_reps, one_judge_boot, workers,
      progress = function(current, total)
        report("judge bootstrap", current, total),
      cancel = cancel, label = "BTL-EFRM judge-bootstrap")
    good_draw <- !vapply(ans, is.null, logical(1))
    boot_fail <- sum(!good_draw)
    draws <- ans[good_draw]
    if (length(draws) < max(20L, ceiling(boot_reps / 2)))
      stop("judge bootstrap failed on ", boot_fail, " of ", boot_reps,
           " replicates; too few judges per panel for stable resampling -- ",
           "use se_method = 'bootstrap' or 'conditional'")
    D <- do.call(rbind, draws)
    colnames(D) <- c(paste0("log phi[", panels_u, "]"),
                     paste0("log alpha[", sets_u, "]"),
                     paste0("kappa[", sets_u, "]"),
                     paste0("beta[", objs_all, "]"), paste0("v[", objs_all, "]"))
    known_na <- !is.finite(c(log(phi), log(alpha), kappa,
                             bhat[objs_all], v))
    D[, known_na] <- NA_real_
    cov_log_phi <- cov_draws(D, seq_len(G))
    if (S > 1L) {
      i2 <- c(G + match(free, sets_u),
              G + S + match(free, sets_u))
      cov2 <- cov_draws(D, i2)
    }
    n_inf <- colSums(!is.finite(D)) * !known_na
    sds <- rep(NA_real_, ncol(D))
    ok_col <- n_inf == 0L & !known_na
    sds[ok_col] <- apply(D[, ok_col, drop = FALSE], 2, sd)
    if (any(n_inf > 0L))
      notes <- c(notes, paste0(
        "bootstrap: ", paste(sprintf("%s reached the boundary in %d of %d replicates",
                                     colnames(D)[n_inf > 0L], n_inf[n_inf > 0L],
                                     length(draws)), collapse = "; "),
        "; the parameter is weakly identified and its SE is reported as NA"))
    nO <- length(objs_all)
    iv <- G + 2L * S + nO + seq_len(nO)
    cov_v <- cov_draws(D, iv)
    dimnames(cov_v) <- list(objs_all, objs_all)
    se_log_phi <- setNames(sds[seq_len(G)], panels_u)
    se_log_alpha <- setNames(sds[G + seq_len(S)], sets_u)
    se_kappa <- setNames(sds[G + S + seq_len(S)], sets_u)
    se_log_alpha[sets_u[1]] <- NA_real_
    se_kappa[sets_u[1]] <- NA_real_
    se_bhat <- setNames(sds[G + 2L * S + seq_len(nO)], objs_all)
    se_v <- setNames(sds[G + 2L * S + nO + seq_len(nO)], objs_all)
    if (boot_fail > 0)
      notes <- c(notes, sprintf(
        "judge bootstrap: %d of %d replicates failed and were skipped",
        boot_fail, boot_reps))
  }
  if (se_method == "bootstrap") {
    p_hat <- pmin(pmax(fit0$p_all, 1e-8), 1 - 1e-8)
    draws <- list()
    report("parametric bootstrap", 0L, boot_reps)
    for (bb in seq_len(boot_reps)) {
      report("parametric bootstrap", bb - 1L, boot_reps)
      yb <- rbinom(length(p_hat), 1L, p_hat)
      fb <- tryCatch(fit_once(yb), error = function(e) NULL)
      if (is.null(fb) || !isTRUE(fb$converged)) { boot_fail <- boot_fail + 1L; next }
      draws[[length(draws) + 1L]] <- c(log(fb$phi), log(fb$alpha), fb$kappa,
                                       fb$bhat, fb$v)
    }
    report("parametric bootstrap", boot_reps, boot_reps)
    if (length(draws) < max(20L, ceiling(boot_reps / 2)))
      stop("parametric bootstrap failed on ", boot_fail, " of ", boot_reps,
           " replicates; the design is too sparse for stable resampling -- ",
           "add comparisons or use se_method = 'conditional'")
    D <- do.call(rbind, draws)
    colnames(D) <- c(paste0("log phi[", panels_u, "]"),
                     paste0("log alpha[", sets_u, "]"),
                     paste0("kappa[", sets_u, "]"),
                     paste0("beta[", objs_all, "]"), paste0("v[", objs_all, "]"))
    # a unit already reported NA in the observed fit (unidentified) is NA
    # in every replicate by the same logic: exclude it from the boundary
    # accounting so its absence is not misreported as boundary behaviour
    known_na <- !is.finite(c(log(phi), log(alpha), kappa,
                             bhat[objs_all], v))
    D[, known_na] <- NA_real_
    cov_log_phi <- cov_draws(D, seq_len(G))
    if (S > 1L) {
      i2 <- c(G + match(free, sets_u),
              G + S + match(free, sets_u))
      cov2 <- cov_draws(D, i2)
    }
    # a parameter that reaches its boundary in some replicates (a set unit
    # driven to zero when a resampled within-set order flips against the
    # cross-set evidence) has no normal sampling distribution: report NA
    # rather than a standard deviation over infinite draws
    n_inf <- colSums(!is.finite(D)) * !known_na
    sds <- rep(NA_real_, ncol(D))
    ok_col <- n_inf == 0L & !known_na
    sds[ok_col] <- apply(D[, ok_col, drop = FALSE], 2, sd)
    if (any(n_inf > 0L))
      notes <- c(notes, paste0(
        "bootstrap: ", paste(sprintf("%s reached the boundary in %d of %d replicates",
                                     colnames(D)[n_inf > 0L], n_inf[n_inf > 0L],
                                     length(draws)), collapse = "; "),
        "; the parameter is weakly identified and its SE is reported as NA"))
    nO <- length(objs_all)
    iv <- G + 2L * S + nO + seq_len(nO)
    cov_v <- cov_draws(D, iv)
    dimnames(cov_v) <- list(objs_all, objs_all)
    se_log_phi <- setNames(sds[seq_len(G)], panels_u)
    se_log_alpha <- setNames(sds[G + seq_len(S)], sets_u)
    se_kappa <- setNames(sds[G + S + seq_len(S)], sets_u)
    se_log_alpha[sets_u[1]] <- NA_real_                 # reference: fixed at 1 / 0
    se_kappa[sets_u[1]] <- NA_real_
    se_bhat <- setNames(sds[G + 2L * S + seq_len(nO)], objs_all)
    se_v <- setNames(sds[G + 2L * S + nO + seq_len(nO)], objs_all)
    if (boot_fail > 0)
      notes <- c(notes, sprintf(
        "bootstrap: %d of %d replicates failed and were skipped", boot_fail,
        boot_reps))
  }

  # --- equal-unit (single-unit) comparison ----------------------------------
  npar_frame <- (length(objs_all) - S) + (G - 1L) + 2L * (S - 1L)
  ll_frames <- ll_within + ll_cross
  single <- tryCatch(
    .btlef_stage1(match(a, objs_all), match(b, objs_all), y,
                  rep("all", length(a)), jd, length(objs_all), maxit, tol),
    error = function(e) NULL)
  ll_single <- if (is.null(single)) NA_real_ else single$ll
  equal_unit <- list(
    loglik_frames = ll_frames, loglik_single = ll_single,
    difference = if (is.na(ll_single)) NA_real_ else ll_frames - ll_single,
    two_delta_ll = if (is.na(ll_single)) NA_real_ else
      2 * (ll_frames - ll_single),
    parameters_frames = npar_frame,
    parameters_single = length(objs_all) - 1L,
    note = paste("descriptive composite-likelihood difference;",
                 "the omnibus Wald tests on the unit families carry the inference"))

  # --- structural tables ----------------------------------------------------
  # Judge-bootstrap inference is limited by the judges who contribute to the
  # particular unit. Count both raw and Kish-effective judges from their
  # comparison workloads. Sparse panels or cross-set links cannot borrow
  # denominator degrees of freedom from judges who informed other units.
  judge_support <- function(rows) {
    ww <- table(jd[rows]); ww <- as.numeric(ww[ww > 0])
    c(n_judges = length(ww), effective_judges = if (length(ww))
      sum(ww)^2 / sum(ww^2) else 0)
  }
  panel_support <- do.call(rbind, lapply(panels_u, function(g) {
    z <- judge_support(within & pan == g)
    data.frame(panel = g, n_judges = z[1L], effective_judges = z[2L],
               stringsAsFactors = FALSE)
  }))
  pa <- pmin(sa, sb); pb <- pmax(sa, sb)
  edge_key <- paste(pa, pb, sep = "\r")
  edge_levels <- unique(edge_key[!within])
  edge_support <- do.call(rbind, lapply(edge_levels, function(e) {
    z <- judge_support(!within & edge_key == e)
    ab <- strsplit(e, "\r", fixed = TRUE)[[1L]]
    data.frame(set_a = ab[1L], set_b = ab[2L], n_judges = z[1L],
               effective_judges = z[2L], stringsAsFactors = FALSE)
  }))
  if (is.null(edge_support)) edge_support <- data.frame(
    set_a = character(), set_b = character(), n_judges = numeric(),
    effective_judges = numeric())
  set_support <- do.call(rbind, lapply(sets_u, function(s) {
    rr <- edge_support$set_a == s | edge_support$set_b == s
    data.frame(set = s,
      n_judges = if (any(rr)) min(edge_support$n_judges[rr]) else 0,
      effective_judges = if (any(rr))
        min(edge_support$effective_judges[rr]) else 0,
      stringsAsFactors = FALSE)
  }))
  panel_ok <- panel_support$n_judges >= 6L &
    panel_support$effective_judges >= 5.5 - sqrt(.Machine$double.eps)
  set_ok <- set_support$n_judges >= 8L &
    set_support$effective_judges >= 8 - sqrt(.Machine$double.eps)
  if (se_method == "judge_bootstrap" && any(!panel_ok))
    notes <- c(notes, paste0(
      "panel-unit inference is withheld because panel(s) ",
      paste(panel_support$panel[!panel_ok], collapse = ", "),
      " have fewer than six judges or 5.5 effective judges"))
  if (se_method == "judge_bootstrap" &&
      (any(!panel_ok) || any(!set_ok[set_support$set %in% free])))
    notes <- c(notes, paste0(
      "set-unit and set-origin inference is withheld because a contributing ",
      "panel has fewer than six judges or 5.5 effective judges, or a link ",
      "has fewer than eight"))
  if (se_method == "judge_bootstrap") {
    pc <- panel_ok & panel_support$effective_judges < 8
    sc <- set_ok & set_support$effective_judges < 9.5 &
      set_support$set %in% free
    if (any(pc)) notes <- c(notes, paste0(
      "panel(s) ", paste(panel_support$panel[pc], collapse = ", "),
      " have 5.5--7.9 effective judges; interpret unit inference cautiously"))
    if (any(sc)) notes <- c(notes, paste0(
      "set link(s) ", paste(set_support$set[sc], collapse = ", "),
      " have 8.0--9.4 effective judges; interpret unit inference cautiously"))
  }
  # Judge-resampling inference uses a finite-sample t reference because the
  # judges are the independent sampling units. The parametric bootstrap draws
  # comparison outcomes independently conditional on the fitted design, so a
  # judge-based denominator is not its reference distribution.
  df_phi <- if (se_method == "judge_bootstrap" && all(panel_ok))
    max(floor(sum(panel_support$effective_judges)) - 1L, 1L)
    else if (se_method == "judge_bootstrap") NA_real_ else Inf
  z_phi <- log(phi) / se_log_phi
  phi_table <- data.frame(panel = panels_u, phi = unname(phi),
                          se_log_phi = unname(se_log_phi),
                          t = unname(z_phi), df = df_phi,
                          p = unname(2 * pt(-abs(z_phi), df_phi)),
                          stringsAsFactors = FALSE)
  df_set <- if (se_method == "judge_bootstrap") {
    z <- pmax(floor(set_support$effective_judges) - 1L, 1L)
    z[!set_ok | !all(panel_ok)] <- NA_real_; setNames(z, set_support$set)
  } else setNames(rep(Inf, S), sets_u)
  z_al <- log(alpha) / se_log_alpha
  alpha_table <- data.frame(set = sets_u, alpha = unname(alpha),
                            se_log_alpha = unname(se_log_alpha),
                            t = unname(z_al), df = unname(df_set[sets_u]),
                            p = unname(2 * pt(-abs(z_al), df_set[sets_u])),
                            stringsAsFactors = FALSE)
  z_ka <- kappa / se_kappa
  kappa_table <- data.frame(set = sets_u, kappa = unname(kappa),
                            se_kappa = unname(se_kappa),
                            t = unname(z_ka), df = unname(df_set[sets_u]),
                            p = unname(2 * pt(-abs(z_ka), df_set[sets_u])),
                            stringsAsFactors = FALSE)
  adjust_unit_table <- function(tab) {
    tab$p_adj <- NA_real_
    usable <- is.finite(tab$p)
    tab$p_adj[usable] <- stats::p.adjust(tab$p[usable], method = "holm")
    tab$significant <- ifelse(is.na(tab$p_adj), NA, tab$p_adj < 0.05)
    tab
  }
  phi_table <- adjust_unit_table(phi_table)
  alpha_table <- adjust_unit_table(alpha_table)
  kappa_table <- adjust_unit_table(kappa_table)
  if (identical(se_method, "conditional")) {
    # Conditional linking errors omit stage-one uncertainty. Null simulation
    # rejected 17.5% (phi) and 35.5% (alpha) at nominal 5%, so ordinary-looking
    # probabilities are not defensible. Retain estimates and their explicitly
    # conditional SEs for preliminary inspection, but withhold inference.
    withhold <- function(tab) {
      tab$t <- NA_real_; tab$df <- NA_real_; tab$p <- NA_real_
      tab$p_adj <- NA_real_; tab$significant <- NA
      tab
    }
    phi_table <- withhold(phi_table)
    alpha_table <- withhold(alpha_table)
    kappa_table <- withhold(kappa_table)
    notes <- c(notes, paste0(
      "unit probabilities and omnibus tests withheld for conditional standard ",
      "errors, which do not propagate stage-one uncertainty; use the judge ",
      "bootstrap for inference"))
  }

  # Hotelling-style F reference with judges as the independent units: the
  # unit covariances are judge-limited (the bootstrap resamples ~10-20
  # judges), and a chi-square reference on a covariance estimated from n
  # units rejects a true null at ~8.7% for the set origins in simulation
  # (550 replicates, 12 judges); W * (n - q) / (q (n - 1)) ~ F(q, n - q)
  # restores the nominal rate, exactly as for the MFRM interaction omnibus.
  wald_unit <- function(est, V, term, n_units = Inf, available = TRUE) {
    if (!length(est) || is.null(V)) return(NULL)
    ok <- is.finite(est) & is.finite(diag(V))
    if (!any(ok)) return(NULL)
    est <- est[ok]; V <- V[ok, ok, drop = FALSE]
    if (any(!is.finite(V))) return(NULL)
    ee <- eigen((V + t(V)) / 2, symmetric = TRUE)
    cutoff <- max(abs(ee$values)) * 1e-8
    use <- ee$values > cutoff
    if (!any(use)) return(NULL)
    Vinv <- ee$vectors[, use, drop = FALSE] %*%
      (t(ee$vectors[, use, drop = FALSE]) / ee$values[use])
    W <- drop(t(est) %*% Vinv %*% est)
    q <- sum(use)
    if (!available) {
      data.frame(term = term, df = q, df2 = NA_real_, wald = W,
                 f = NA_real_, p = NA_real_)
    } else if (is.infinite(n_units)) {
      data.frame(term = term, df = q, df2 = Inf, wald = W,
                 f = W / q,
                 p = stats::pchisq(W, q, lower.tail = FALSE))
    } else if (n_units > q) {
      Fs <- W * (n_units - q) / (q * (n_units - 1))
      data.frame(term = term, df = q, df2 = n_units - q, wald = W,
                 f = Fs, p = stats::pf(Fs, q, n_units - q,
                                       lower.tail = FALSE))
    } else
      data.frame(term = term, df = q, df2 = NA_real_, wald = W,
                 f = NA_real_, p = NA_real_)
  }
  omni_parts <- list()
  if (G > 1L) omni_parts[[length(omni_parts) + 1L]] <- wald_unit(
    log(phi), cov_log_phi, "panel units (phi)",
    n_units = if (se_method == "judge_bootstrap")
      floor(sum(panel_support$effective_judges)) else Inf,
    available = se_method != "judge_bootstrap" || all(panel_ok))
  if (S > 1L) {
    set_n <- if (se_method == "judge_bootstrap")
      floor(min(set_support$effective_judges[
        match(free, set_support$set)])) else Inf
    set_available <- se_method != "judge_bootstrap" ||
      (all(panel_ok) && all(set_ok[match(free, set_support$set)]))
    omni_parts[[length(omni_parts) + 1L]] <- wald_unit(
      log(alpha[free]),
      cov2[seq_along(free), seq_along(free), drop = FALSE],
      "set units (alpha)", n_units = set_n, available = set_available)
    omni_parts[[length(omni_parts) + 1L]] <- wald_unit(
      kappa[free],
      cov2[length(free) + seq_along(free),
           length(free) + seq_along(free), drop = FALSE],
      "set origins (kappa)", n_units = set_n, available = set_available)
  }
  unit_omnibus <- do.call(rbind, Filter(Negate(is.null), omni_parts))
  if (identical(se_method, "conditional") && !is.null(unit_omnibus)) {
    unit_omnibus$df2 <- NA_real_
    unit_omnibus$f <- NA_real_
    unit_omnibus$p <- NA_real_
  }

  objects <- data.frame(object = objs_all, set = unname(set_of[objs_all]),
                        location = unname(v), se = unname(se_v),
                        beta_set = unname(bhat[objs_all]),
                        se_beta = unname(se_bhat[objs_all]),
                        v = unname(v), se_v = unname(se_v),
                        stringsAsFactors = FALSE)
  rownames(objects) <- NULL

  # frames: one row per panel-by-set cell holding within-set comparisons
  fr <- list()
  f_cell <- max((length(y) - npar_frame) / length(y), 0)
  for (s in sets_u) for (g in panels_u) {
    rows <- which(within & sa == s & pan == g)
    if (!length(rows)) next
    ff <- .btlef_frame_fit(y[rows], within_p[rows], f_cell = f_cell)
    fr[[length(fr) + 1L]] <- data.frame(
      # within-set logit = phi_g (beta_a - beta_b) and beta = (v - kappa)/
      # alpha, so the discrimination applied to COMMON-SCALE differences
      # within this frame is phi/alpha (phi * alpha composed the map the
      # wrong way round)
      panel = g, set = s, rho = unname(phi[[g]] / alpha[[s]]),
      n_comparisons = length(rows),
      fit_resid = ff$fit_resid, df_fit = ff$df,
      stringsAsFactors = FALSE)
  }
  frames <- if (length(fr)) do.call(rbind, fr) else NULL
  if (!is.null(frames)) rownames(frames) <- NULL

  diag_frame <- .btlef_diagnostics(
    a, b, y, jd, pan, sa, sb, fit0$p_all, phi, alpha_use,
    objects, n_parameters = npar_frame)
  objects <- diag_frame$objects
  osi <- .psi(objects$location, objects$se)

  if (S == 1L)
    notes <- c(notes, "single set: panel-units model (set units alpha not estimated)")
  if (G == 1L && S == 1L)
    notes <- c(notes, "single panel and single set: reduces to btl()")
  if (length(fit0$dropped))
    notes <- c(notes, paste0(
      "set(s) ", paste(fit0$dropped, collapse = ", "), " carry no stable ",
      "panel-ratio information (within-set contests too close to even or ",
      "too one-sided); they were excluded from the phi reconciliation and ",
      "refit with the panel units held at the reconciled phi"))

  report("finalising", 1L, 1L)
  out <- list(objects = objects, phi_table = phi_table,
              alpha_table = alpha_table, kappa_table = kappa_table,
              unit_omnibus = unit_omnibus,
              unit_support = list(panel = panel_support, edge = edge_support,
                                  set = set_support,
                                  minimum_panel_judges = 6L,
                                  minimum_panel_effective_judges = 5.5,
                                  minimum_link_judges = 8L),
              frames = frames, equal_unit = equal_unit, n_cross = n_cross,
              sets = sets_u, panels = panels_u, reference_set = sets_u[1],
              n_comparisons = length(a),
              converged = fit0$converged,
              # rasch_btl-compatible core diagnostics, all evaluated under
              # the frame probabilities rather than copied from an equal-unit
              # fit. This lets the GUI and generic summaries use the active
              # common-scale calibration after frames are added.
              pairs = diag_frame$pairs, judges = diag_frame$judges,
              comparisons = diag_frame$comparisons,
              total_chisq = diag_frame$total_chisq,
              total_df = diag_frame$total_df, total_p = diag_frame$total_p,
              osi = osi, loglik = ll_frames, iterations = NA_integer_,
              n_parameters = npar_frame,
              clustered = TRUE, cov_beta = cov_v,
              thresholds = NULL, components = NULL,
              thr_structure = "dichotomous", m = 1L,
              categories = c("0", "1"), dependence = NULL,
              dependence_data = NULL,
              se_method = se_method,
              workers = workers, seed = seed,
              boot_reps = if (se_method %in% c("bootstrap", "judge_bootstrap"))
                boot_reps else NA_integer_,
              boot_reps_used = if (se_method %in%
                c("bootstrap", "judge_bootstrap")) length(draws) else NA_integer_,
              se_note = if (se_method == "judge_bootstrap")
                paste("standard errors from a judge-resampling bootstrap",
                      "(judges redrawn with replacement within panels, the",
                      "whole pipeline refitted): carries stage-one",
                      "uncertainty and any extra-model dependence within",
                      "judges")
              else if (se_method == "bootstrap")
                paste("standard errors from a parametric bootstrap of the",
                      "whole two-stage pipeline: the two-stage",
                      "estimates are unchanged, and their errors carry the",
                      "stage-one uncertainty into the linking")
              else
                paste("analytic standard errors for beta and phi; alpha and",
                      "kappa errors condition on the stage-one estimates and",
                      "do not represent total pipeline uncertainty; use",
                      "se_method = 'judge_bootstrap' for inference"),
              notes = notes)
  out <- .tag_tables(out)
  class(out) <- c("rasch_btl_efrm", "rasch_btl")
  out
}

#' @export
print.rasch_btl_efrm <- function(x, ...) {
  cat(sprintf(paste0("Bradley-Terry-Luce extended frame of reference: ",
                     "%d objects in %d set(s) x %d panel(s), %d comparisons\n"),
              nrow(x$objects), nrow(x$alpha_table), nrow(x$phi_table),
              x$n_comparisons))
  cat(sprintf("Two-stage maximum likelihood: %s; SEs %s\n",
              if (x$converged) "converged" else "NOT converged",
              if (identical(x$se_method, "judge_bootstrap"))
                sprintf("by judge-resampling bootstrap (B = %d)", x$boot_reps)
              else if (identical(x$se_method, "bootstrap"))
                sprintf("by parametric bootstrap (B = %d)", x$boot_reps)
              else "conditional (understate the linking uncertainty)"))
  if (nrow(x$alpha_table) == 1L)
    cat("Model: panel units only (single set; set units not estimated)\n")
  if (!is.null(x$unit_omnibus)) {
    cat("\nOmnibus Wald tests of equal units and origins:\n")
    print(.fmt_df(x$unit_omnibus), row.names = FALSE)
  }
  cat("\nPanel units (phi; exploratory Holm-adjusted contrasts):\n")
  print(.fmt_df(x$phi_table), row.names = FALSE)
  if (nrow(x$alpha_table) > 1L) {
    cat("\nSet units (alpha) and origins (kappa; reference set = ",
        x$reference_set, "):\n", sep = "")
    at <- merge(x$alpha_table[, c("set", "alpha", "se_log_alpha",
                                  "p_adj", "significant")],
                x$kappa_table[, c("set", "kappa", "se_kappa", "p_adj",
                                  "significant")],
                by = "set", sort = FALSE)
    names(at)[names(at) == "p_adj.x"] <- "p_adj_alpha"
    names(at)[names(at) == "significant.x"] <- "significant_alpha"
    names(at)[names(at) == "p_adj.y"] <- "p_adj_kappa"
    names(at)[names(at) == "significant.y"] <- "significant_kappa"
    print(.fmt_df(at), row.names = FALSE)
  }
  eu <- x$equal_unit
  if (!is.na(eu$difference))
    cat(sprintf(paste0("\nEqual-unit comparison: ll_frames - ll_single = ",
                       "%.3f (%s)\n"), eu$difference, eu$note))
  if (length(x$notes)) cat(sprintf("Notes: %s\n", paste(x$notes, collapse = "; ")))
  invisible(x)
}

# Object characteristic display for a frame-adjusted comparison fit. A single
# equal-unit logistic curve is not available: the slope depends on panel and,
# for within-set comparisons, on the object's set unit. Curves and observed
# points therefore share a colour within each fitted frame.
.plot_btl_efrm_icc <- function(fit, object, grid = NULL, min_n = 10) {
  ob <- fit$objects
  if (!object %in% ob$object) stop("no such object: ", object)
  vo <- ob$location[match(object, ob$object)]
  set_o <- ob$set[match(object, ob$object)]
  if (is.null(grid)) {
    rng <- range(ob$location) + c(-1, 1)
    grid <- seq(rng[1], rng[2], length.out = 201)
  }
  cm <- fit$comparisons
  take_a <- cm$object_a == object
  take_b <- cm$object_b == object
  d <- rbind(
    data.frame(opponent = cm$object_b[take_a], response = cm$response[take_a],
               weight = cm$weight[take_a], panel = cm$panel[take_a],
               set_opponent = cm$set_b[take_a],
               slope = cm$frame_slope[take_a]),
    data.frame(opponent = cm$object_a[take_b], response = 1 - cm$response[take_b],
               weight = cm$weight[take_b], panel = cm$panel[take_b],
               set_opponent = cm$set_a[take_b],
               slope = cm$frame_slope[take_b]))
  d <- d[d$opponent %in% ob$object, , drop = FALSE]
  d$frame <- ifelse(d$set_opponent == set_o,
                    paste(d$panel, "within", set_o),
                    paste(d$panel, "vs", d$set_opponent))
  sp <- split(seq_len(nrow(d)), .factor_keys(
    data.frame(frame = d$frame, opponent = d$opponent,
               stringsAsFactors = FALSE)))
  obs <- do.call(rbind, lapply(sp, function(ii) data.frame(
    frame = d$frame[ii[1]], opponent = d$opponent[ii[1]],
    loc = ob$location[match(d$opponent[ii[1]], ob$object)],
    mean = sum(d$weight[ii] * d$response[ii]) / sum(d$weight[ii]),
    n = sum(d$weight[ii]), slope = mean(d$slope[ii]))))
  obs <- obs[obs$n >= min_n, , drop = FALSE]
  frames <- unique(d$frame)
  slope <- vapply(frames, function(fr) mean(d$slope[d$frame == fr]), 0)
  cols <- setNames(rep(.rr$pal, length.out = length(frames)), frames)
  op <- .rr_canvas(range(grid), c(0, 1),
                   "Opponent common-scale location (logits)",
                   "Probability preferred",
                   sprintf("%s  (common-scale location %.3f)", object, vo))
  on.exit(par(op))
  for (fr in frames)
    lines(grid, stats::plogis(slope[[fr]] * (vo - grid)),
          col = cols[[fr]], lwd = 2.2)
  if (nrow(obs))
    points(obs$loc, obs$mean, pch = 21, bg = cols[obs$frame],
           col = "white", cex = 1.45, lwd = 1.1)
  .rr_legend("topright", c("Model", "Observed"), lwd = c(2.2, NA),
             pch = c(NA, 21), pt.bg = c(NA, .rr$blue),
             col = c(.rr$ink, "white"), pt.cex = 1.25)
  if (length(frames) <= 8L)
    .rr_legend("bottomleft", frames, lwd = 2.2, col = unname(cols), cex = .72)
  invisible(if (nrow(obs)) obs$opponent else character(0))
}

#' Plot the frame units of a paired-comparison EFRM fit
#'
#' Caterpillar plot of panel units \code{phi_g} and set units \code{alpha_s}
#' on the log scale, with 95 per cent intervals and unit one marked.
#'
#' @param fit A fitted object from \code{\link{btl_efrm}}.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' \donttest{
#' # see ?btl_efrm for a complete simulated example
#' }
#' @export
plot_btl_units <- function(fit) {
  if (!inherits(fit, "rasch_btl_efrm"))
    stop("plot_btl_units needs a rasch_btl_efrm fit")
  ph <- fit$phi_table; al <- fit$alpha_table
  rows <- rbind(
    data.frame(label = paste0("panel: ", ph$panel), kind = "panel",
               est = log(ph$phi), se = ph$se_log_phi, stringsAsFactors = FALSE),
    if (nrow(al) > 1L)
      data.frame(label = paste0("set: ", al$set), kind = "set",
                 est = log(al$alpha), se = al$se_log_alpha,
                 stringsAsFactors = FALSE))
  rows <- rows[order(rows$kind, rows$est), ]
  n <- nrow(rows)
  lo <- rows$est - 1.96 * rows$se; hi <- rows$est + 1.96 * rows$se
  colr <- ifelse(rows$kind == "panel", .rr$blue, .rr$purple)
  op <- par(mar = c(4.2, 9, 3.2, 1.5), mgp = c(2.5, 0.7, 0), tcl = -0.25,
            las = 1, col.axis = .rr$ink, col.lab = .rr$ink, col.main = .rr$ink,
            font.main = 2, cex.main = 1.15)
  on.exit(par(op))
  plot(NA, xlim = range(c(rows$est, lo, hi, 0), na.rm = TRUE) + c(-0.1, 0.1),
       ylim = c(0.5, n + 0.5),
       xlab = "log unit", ylab = "", axes = FALSE, main = "")
  abline(h = seq_len(n), col = .rr$grid, lwd = 0.8)
  abline(v = 0, lty = 2, col = .rr$soft)
  axis(1, col = .rr$grid, col.ticks = .rr$soft)
  axis(2, at = seq_len(n), labels = rows$label, cex.axis = 0.75,
       col = .rr$grid, col.ticks = NA)
  hs <- is.finite(rows$se)
  segments(lo[hs], which(hs), hi[hs], which(hs), lwd = 2.2,
           col = .rr$soft)
  points(rows$est, seq_len(n), pch = 21, cex = 1.5, bg = colr,
         col = "white", lwd = 1.2)
  .rr_legend("bottomright", c("panel unit (phi)", "set unit (alpha)"),
             pch = 21, pt.bg = c(.rr$blue, .rr$purple), col = "white",
             pt.cex = 1.2)
  invisible(NULL)
}
