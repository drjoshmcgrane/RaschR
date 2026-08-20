# rasch :: extended frame of reference model
# ===========================================================================
# Humphry's extended frame of reference model (Humphry 2005; Humphry &
# Andrich 2008). A frame F_sg is one item-set by person-group cell, with
# unit rho_sg = alpha_s * phi_g:
#
#   P(X_ni = x) prop exp( rho_sg * ( x*theta_n - sum_{h<=x} delta_ih ) )
#
# Within a frame all curves are parallel, so the partial credit model holds
# in the frame's natural unit and the pairwise conditional logic of
# Andrich & Luo (2003) applies unchanged within frames.
#
# Humphry (2005) states the model for dichotomous responses; the polytomous
# form above is this package's extension of it. It is characterised by the
# properties the model's logic requires: the unit multiplies the whole
# exponent, so (i) within every frame the partial credit model holds in the
# frame's natural unit (natural thresholds rho*(delta - c), parallel
# curves), which is what makes the pairwise conditional cancellation valid;
# (ii) the weighted score sum_i rho_i x_i remains sufficient for the person
# parameter; and (iii) it reduces exactly to the dichotomous statement when
# every item has two categories and to the ordinary PCM when every unit is
# one. Per-threshold discriminations would destroy (i) and with it Rasch
# measurement; a unit on only part of the exponent would destroy the
# natural-unit reading. A consequence worth noting when interpreting
# results: category widths in natural units scale with the frame unit.
#
# Estimation is in two routes, each assigned to the data structure that
# identifies it:
#
# Route 1 (person-free, within-frame pairwise conditional ML): the centred,
# alpha-absorbed set thresholds dtilde (sum-zero per set) and the person
# group units phi_g, identified because a set taken by two groups shows the
# same threshold pattern at two scales. The structural map
# tau_v = phi_g(v) * dtilde is bilinear, so estimation alternates two exact
# linear Newton steps and finishes with a joint polish for exact sandwich
# standard errors. Frame origins are NOT pairwise-identifiable (an additive
# constant per frame cancels at every pair total) and item-set units
# alpha_s are exactly confounded with the spread of the set's thresholds,
# so neither is a parameter of this stage.
#
# Route 2 (person-side linking): alpha_s and set locations mu_s from
# persons common to two sets. Within set s a person's frame estimate is
# u = alpha_s * (theta - mu_s), so the ratio of true-score variances over
# common persons estimates squared alpha ratios, and the offsets give the
# relative set locations; the log-ratios are reconciled over the linking
# graph by weighted least squares. The true-score variance is recovered by
# the truncated-score-moment correction
#     var_true = (var(u_hat) - E[w]) / E[g']^2,
# where g(u) and w(u) are the exact mean and variance of the WLE score map
# over the non-extreme score distribution at u (from the fitted
# thresholds), and E[.] over the person distribution is estimated by
# score-weight vectors phi solving E[phi(R)|u] ~ target(u) on a grid --
# unbiased for any person distribution because the raw score is
# sufficient. The naive correction var(u_hat) - mean(se^2) is badly
# calibrated on short tests: the reported SE^2 overstates the actual error
# variance and WLE shrinkage induces cov(u, error) < 0; at 8 dichotomous
# items per set both distortions exceed 50% of the variance and their
# imperfect cancellation biased the log unit ratio by +0.05.
#
# Identification: sum dtilde = 0 per set; sum_g log phi_g = 0;
# sum_s log alpha_s = 0; mean_s mu_s = 0. Together these give the
# product-of-units constraint over the frame grid with the arbitrary-unit
# origin at the mean set location. Cross-set item pairs carry no usable
# conditional information when units differ (the person parameter does not
# cancel) and are excluded from the pairwise stage.
# ===========================================================================

# Same-set pair filter: .pair_counts only returns pairs with overlapping
# persons, so cross-group pairs are already absent; cross-set pairs within a
# group must be removed because the person parameter does not cancel there.
.efrm_filter_pairs <- function(pairs, vmap) {
  Filter(function(pc) vmap$set[pc$i] == vmap$set[pc$j], pairs)
}

# Connected components by union-find; x is a list of integer edge pairs.
.efrm_components <- function(n, edges) {
  parent <- seq_len(n)
  find <- function(a) { while (parent[a] != a) a <- parent[a]; a }
  for (e in edges) {
    ra <- find(e[1]); rb <- find(e[2])
    if (ra != rb) parent[ra] <- rb
  }
  vapply(seq_len(n), find, 1L)
}

# Stage 1: block-coordinate Newton on the bilinear map tau = phi_g * dtilde,
# with a joint polish step and Godambe sandwich covariance.
.efrm_solve <- function(Xv, thr_v, m_v, vmap, pairs, drow, A_D,
                        maxit = 50, tol = 1e-7) {
  Mv <- nrow(thr_v)
  glevs <- sort(unique(vmap$group))
  G <- length(glevs)
  gidx <- match(vmap$group[thr_v$item], glevs)   # group of each virtual threshold
  Pd <- ncol(A_D)

  # one inner Newton block on a linear design tau = off + B beta
  newton_block <- function(B, off, beta, positive = FALSE, inner = 4) {
    glh <- .pcml_glh(drop(off + B %*% beta), thr_v, pairs, m_v)
    for (it in seq_len(inner)) {
      gb <- drop(crossprod(B, glh$g))
      Hb <- crossprod(B, glh$H %*% B)
      step <- tryCatch(solve(Hb, gb), error = function(e)
        solve(Hb - diag(1e-8, nrow(Hb)), gb))
      lam <- 1; ok <- FALSE; g2 <- glh
      for (half in 1:30) {
        cand <- beta - lam * step
        if (positive && any(cand <= 0)) { lam <- lam / 2; next }
        g2 <- .pcml_glh(drop(off + B %*% cand), thr_v, pairs, m_v)
        if (is.finite(g2$ll) && g2$ll >= glh$ll - 1e-12) { ok <- TRUE; break }
        lam <- lam / 2
      }
      if (!ok) break
      beta <- cand; glh <- g2
      if (max(abs(lam * step)) < tol) break
    }
    list(beta = beta, ll = glh$ll, glh = glh)
  }

  # start values: log-ratio least squares pooled over groups, phi = 1
  st <- .start_tau(Xv, thr_v)
  dt0 <- vapply(seq_len(max(drow)), function(d) mean(st[drow == d]), 0)
  # recentre per set
  for (s in unique(vmap$set)) {
    rows <- unique(drow[vmap$set[thr_v$item] == s])
    dt0[rows] <- dt0[rows] - mean(dt0[rows])
  }
  beta_d <- qr.coef(qr(A_D), dt0); beta_d[is.na(beta_d)] <- 0
  phi <- rep(1, G)

  ll <- -Inf
  for (outer in seq_len(maxit)) {
    # delta step: tau = B_d beta_d with B_d = phi_g(v) * A_D[drow(v), ]
    B_d <- A_D[drow, , drop = FALSE] * phi[gidx]
    res_d <- newton_block(B_d, 0, beta_d)
    beta_d <- res_d$beta
    dtil <- drop(A_D %*% beta_d)
    # phi step (skip when G = 1): tau = off + B_r phi[-1]
    if (G > 1L) {
      off <- dtil[drow] * (gidx == 1L)
      B_r <- matrix(0, Mv, G - 1L)
      for (g in 2:G) B_r[gidx == g, g - 1L] <- dtil[drow][gidx == g]
      res_r <- newton_block(B_r, off, phi[-1], positive = TRUE)
      phi <- c(1, res_r$beta)
      ll_new <- res_r$ll
    } else ll_new <- res_d$ll
    if (is.finite(ll) && abs(ll_new - ll) < tol * (abs(ll_new) + 1)) {
      ll <- ll_new; break
    }
    ll <- ll_new
  }

  # joint polish with the exact bilinear Hessian, then sandwich covariance
  tau_hat <- phi[gidx] * dtil[drow]
  glh <- .pcml_glh(tau_hat, thr_v, pairs, m_v)
  build_J <- function(dtil, phi) {
    J <- matrix(0, Mv, Pd + G - 1L)
    J[, seq_len(Pd)] <- A_D[drow, , drop = FALSE] * phi[gidx]
    if (G > 1L) for (g in 2:G)
      J[gidx == g, Pd + g - 1L] <- dtil[drow][gidx == g]
    J
  }
  struct_H <- function(J, glh) {
    H <- crossprod(J, glh$H %*% J)
    if (G > 1L) for (g in 2:G) {
      rows <- which(gidx == g)
      cb <- drop(crossprod(A_D[drow[rows], , drop = FALSE], glh$g[rows]))
      H[seq_len(Pd), Pd + g - 1L] <- H[seq_len(Pd), Pd + g - 1L] + cb
      H[Pd + g - 1L, seq_len(Pd)] <- H[Pd + g - 1L, seq_len(Pd)] + cb
    }
    H
  }
  for (polish in 1:10) {
    J <- build_J(dtil, phi)
    H <- struct_H(J, glh)
    g_full <- drop(crossprod(J, glh$g))
    step <- tryCatch(solve(H, g_full), error = function(e)
      solve(H - diag(1e-8, nrow(H)), g_full))
    cand_d <- beta_d - step[seq_len(Pd)]
    cand_p <- if (G > 1L) phi[-1] - step[Pd + seq_len(G - 1L)] else numeric(0)
    if (G > 1L && any(cand_p <= 0)) break
    dt_c <- drop(A_D %*% cand_d); phi_c <- c(1, cand_p)
    glh_c <- .pcml_glh(phi_c[gidx] * dt_c[drow], thr_v, pairs, m_v)
    if (!is.finite(glh_c$ll) || glh_c$ll < glh$ll - 1e-10) break
    beta_d <- cand_d; dtil <- dt_c; phi <- phi_c; glh <- glh_c
    if (max(abs(step)) < tol) break
  }
  tau_hat <- phi[gidx] * dtil[drow]

  J <- build_J(dtil, phi)
  H <- struct_H(J, glh)
  # identification of the group units: a (near-)null direction of the
  # UNRIDGED joint information that loads on a group's unit parameter
  # means the data cannot identify that unit -- the ridged solve would
  # land anywhere on the flat manifold and still show a small gradient.
  # (H is the negative expected Hessian here, so its eigenvalues are
  # nonnegative up to the bilinear cross term.)
  phi_unident <- setNames(rep(FALSE, G), glevs)
  if (G > 1L) {
    eh <- eigen(H, symmetric = TRUE)
    flat <- abs(eh$values) < max(abs(eh$values)) * 1e-10
    if (any(flat)) {
      V <- eh$vectors[, flat, drop = FALSE]
      lp_load <- sqrt(rowSums(V[Pd + seq_len(G - 1L), , drop = FALSE]^2))
      phi_unident[-1L][lp_load > 1e-2] <- TRUE
    }
  }
  Hinv <- tryCatch(solve(H), error = function(e) solve(H - diag(1e-8, nrow(H))))
  Jt <- .pcml_sandwich(Xv, thr_v, m_v, tau_hat, pairs)
  covb <- Hinv %*% crossprod(J, Jt %*% J) %*% Hinv
  conv <- max(abs(drop(crossprod(J, glh$g)))) < 1e-3

  # recentre to sum_g log phi = 0; dtilde absorbs the constant
  cc <- mean(log(phi))
  phi_c <- phi * exp(-cc); dtil_c <- dtil * exp(cc)
  # covariance of log phi (ref group has zero variance), centred
  cov_lp <- matrix(0, G, G)
  if (G > 1L) {
    idx <- Pd + seq_len(G - 1L)
    Dl <- diag(1 / phi[-1], G - 1L)
    cov_lp[2:G, 2:G] <- Dl %*% covb[idx, idx, drop = FALSE] %*% Dl
  }
  Ac <- diag(G) - matrix(1 / G, G, G)
  cov_lp <- Ac %*% cov_lp %*% t(Ac)
  # dtilde_c = exp(cc) A_D beta_d, where
  # cc = mean(log(phi)).  Because cc is estimated jointly with beta_d, the
  # delta-method Jacobian must include its derivatives with respect to every
  # free phi as well as the beta--phi cross-covariance. Treating cc as fixed
  # can materially misstate the common-unit item standard errors.
  J_dt <- matrix(0, nrow(A_D), Pd + G - 1L)
  J_dt[, seq_len(Pd)] <- exp(cc) * A_D
  if (G > 1L)
    J_dt[, Pd + seq_len(G - 1L)] <-
      outer(dtil_c, 1 / (G * phi[-1L]))
  cov_dt <- J_dt %*% covb %*% t(J_dt)

  # joint covariance of (dtilde, log phi): both are functions of the same
  # jointly estimated (beta_d, phi), so downstream propagation (the
  # linking-stage calibration redraws) must draw them together -- the
  # cross-covariance is part of the estimate, not an approximation choice
  J_lp <- matrix(0, G, Pd + G - 1L)
  if (G > 1L) J_lp[2:G, Pd + seq_len(G - 1L)] <- diag(1 / phi[-1], G - 1L)
  J_lp <- Ac %*% J_lp
  J_all <- rbind(J_dt, J_lp)
  cov_joint <- J_all %*% covb %*% t(J_all)

  dimnames(cov_lp) <- list(glevs, glevs)
  list(dtilde = dtil_c, phi = setNames(phi_c, glevs),
       se_log_phi = setNames(sqrt(pmax(diag(cov_lp), 0)), glevs),
       cov_log_phi = cov_lp, phi_unident = phi_unident,
       cov_dtilde = cov_dt, cov_joint = cov_joint,
       loglik = glh$ll, iterations = outer,
       converged = conv, gidx = gidx)
}

# Per-person correction moments for the set-unit linking. For each
# missing-data pattern the WLE score map W and the model score
# distribution are exact, so the mean g(u) and variance w(u) of W(R) over
# the non-extreme scores are exact functions of u given the thresholds.
# The linking needs E[w(u)] and E[g'(u)] over the (unobserved) person
# distribution; because the raw score is sufficient, score-weight vectors
# phi with E[phi(R) | u] ~ target(u) on a grid make mean(phi(R_i)) an
# estimate of E[target(u)] that is unbiased for ANY person distribution
# (both targets are bounded and tail-flat, so the ridge least-squares
# inversion is stable, unlike a direct representation of u^2). Returns
# per-person phi_w(R_i) and phi_g(R_i), NA for extreme or empty rows.
.person_link_moments <- function(X, tau_list, disc = 1) {
  N <- nrow(X)
  obs <- !is.na(X)
  pat <- apply(obs, 1, function(z) paste(which(z), collapse = ","))
  w_i <- g_i <- rep(NA_real_, N)
  for (key in unique(pat)) {
    cols <- as.integer(strsplit(key, ",", fixed = TRUE)[[1]])
    if (length(cols) < 2L) next
    tl <- tau_list[cols]
    maxr <- sum(vapply(tl, length, 1L))
    # a pattern needs a score range of at least 4 (>= 3 interior score
    # categories) to support the moment correction: at 2 interior
    # categories the map slope is degenerate and the estimator collapses
    # to a confidently wrong value (-0.27 log-ratio bias, near-zero
    # variance, in the limits study) rather than getting noisier
    if (maxr < 4L) next
    W <- person_wle(tl, disc = disc)$theta[as.character(0:maxr)]
    Wi <- W[-c(1L, maxr + 1L)]
    if (any(!is.finite(Wi))) next
    grid <- seq(min(Wi) - 3, max(Wi) + 3, length.out = 121L)
    # truncated score distribution at each grid point, by direct
    # convolution of the item category probabilities (same model as
    # person_wle; direct is faster than FFT at these lengths)
    B <- vapply(grid, function(u) {
      p <- 1
      for (tau in tl) {
        pi <- item_moments(u, tau, disc = disc)$P
        np <- numeric(length(p) + length(pi) - 1L)
        for (x in seq_along(pi)) {
          ix <- x:(x + length(p) - 1L)
          np[ix] <- np[ix] + p * pi[x]
        }
        p <- np
      }
      p <- pmax(p, 0)[-c(1L, maxr + 1L)]
      p / sum(p)
    }, numeric(maxr - 1L))
    gm <- colSums(B * Wi)
    wm <- colSums(B * Wi^2) - gm^2
    h <- grid[2L] - grid[1L]
    gp <- c(gm[2L] - gm[1L], (gm[-c(1L, 2L)] - gm[-c(120L, 121L)]) / 2,
            gm[121L] - gm[120L]) / h
    A2 <- t(B)
    M <- crossprod(A2) + diag(1e-6, maxr - 1L)
    phi_w <- drop(solve(M, crossprod(A2, wm)))
    phi_g <- drop(solve(M, crossprod(A2, gp)))
    sel <- which(pat == key)
    r <- rowSums(X[sel, cols, drop = FALSE])
    inner <- r >= 1 & r <= maxr - 1L
    w_i[sel[inner]] <- phi_w[r[inner]]
    g_i[sel[inner]] <- phi_g[r[inner]]
  }
  list(w = w_i, g = g_i)
}

# Stage 2: alpha_s and set locations mu_s from persons common to set pairs,
# with the truncated-score-moment correction (see .person_link_moments).
# Standard errors and the joint covariance of (log alpha, mu) come from a
# person bootstrap of the whole linking stage: each replicate resamples
# persons, recomputes the corrected variance ratios and offsets per set
# pair, and re-solves the linking least squares. This captures the
# nonlinearity of the correction and the reconciliation over the linking
# graph, which a leave-one-out jackknife of the raw log slope understates.
.efrm_link_sets <- function(u_mat, w_mat, g_mat, sets_u, min_link_persons,
                            boot_reps = 300, regen = NULL) {
  S <- ncol(u_mat)
  A <- rbind(diag(S - 1L), rep(-1, S - 1L))

  # one pass over a person index vector: per-edge corrected slopes and
  # offsets, then the two weighted least-squares solves. um/wm/gm default
  # to the point-estimate person matrices; the bootstrap can pass
  # regenerated ones (see below).
  link_once <- function(idx, um = u_mat, wm = w_mat, gm = g_mat,
                        hard = FALSE) {
    edges <- list(); ls_est <- off_est <- off_n <- numeric(0)
    for (a in seq_len(S - 1)) for (b in (a + 1):S) {
      ok <- idx[is.finite(um[idx, a]) & is.finite(um[idx, b]) &
                is.finite(wm[idx, a]) & is.finite(wm[idx, b]) &
                is.finite(gm[idx, a]) & is.finite(gm[idx, b])]
      if (length(ok) < min_link_persons) next
      u1 <- um[ok, a]; u2 <- um[ok, b]
      d1 <- mean(gm[ok, a]); d2 <- mean(gm[ok, b])
      if (!is.finite(d1) || !is.finite(d2) || d1 < 0.1 || d2 < 0.1) {
        if (hard)
          stop("the score maps of sets '", sets_u[a], "' and '", sets_u[b],
               "' are too flat to link (degenerate map slope)")
        next
      }
      v1 <- (var(u1) - mean(wm[ok, a])) / d1^2
      v2 <- (var(u2) - mean(wm[ok, b])) / d2^2
      if (!is.finite(v1) || !is.finite(v2) || v1 <= 0 || v2 <= 0) {
        if (hard)
          stop("too little true person variance to link sets '", sets_u[a],
               "' and '", sets_u[b], "'")
        next
      }
      ls <- 0.5 * (log(v2) - log(v1))                # log(alpha_b / alpha_a)
      edges[[length(edges) + 1L]] <- c(a, b)
      ls_est <- c(ls_est, ls)
      off_est <- c(off_est, mean(u2) - exp(ls) * mean(u1))  # alpha_b (mu_a - mu_b)
      off_n <- c(off_n, length(ok))
    }
    if (!length(edges)) {
      if (hard) stop("no set pairs share enough persons with informative ",
                     "score patterns to link the units: each person's ",
                     "pattern needs a score range of at least 4 within a ",
                     "set (at least four dichotomous items, or fewer ",
                     "polytomous ones) for the set-unit correction")
      return(NULL)
    }
    comp <- .efrm_components(S, edges)
    if (length(unique(comp)) > 1L) {
      if (hard)
        stop("item sets are not linked by common persons; relative units (alpha) ",
             "are unidentified between: ",
             paste(tapply(sets_u, comp, paste, collapse = "+"), collapse = " | "))
      return(NULL)
    }
    C <- matrix(0, length(edges), S)
    for (e in seq_along(edges)) { C[e, edges[[e]][2]] <- 1; C[e, edges[[e]][1]] <- -1 }
    sw <- sqrt(off_n); M <- C %*% A
    la <- drop(A %*% qr.coef(qr(M * sw), ls_est * sw)); la[is.na(la)] <- 0
    alpha <- exp(la)
    dmu <- off_est / alpha[vapply(edges, `[`, 1L, 2)]
    Cm <- matrix(0, length(edges), S)
    for (e in seq_along(edges)) { Cm[e, edges[[e]][1]] <- 1; Cm[e, edges[[e]][2]] <- -1 }
    mu <- drop(A %*% qr.coef(qr((Cm %*% A) * sw), dmu * sw)); mu[is.na(mu)] <- 0
    list(la = la, mu = mu, edges = edges, ls_est = ls_est, off_n = off_n)
  }

  N <- nrow(u_mat)
  point <- link_once(seq_len(N), hard = TRUE)

  # person bootstrap of the linking stage (skipped inside an outer
  # bootstrap). When a regen closure is supplied, each replicate also
  # redraws the within-frame thresholds and group units from their
  # estimated covariances and rebuilds the person estimates, so the
  # set-unit uncertainty carries the calibration noise that person
  # resampling alone cannot see: the set unit is a scale, and error in the
  # estimated threshold spread moves every person estimate's variance
  # coherently. Without this the hybrid log-alpha standard error
  # understates by ~20% in simulation and the unit tests reject at ~9-10%
  # instead of 5%.
  cov_link <- NULL
  if (boot_reps > 0) {
    # a pool of parameter draws cycled across the person resamples: each
    # replicate pairs one draw with one resample, so the replicate spread
    # still carries both variance components, at a fraction of the regen
    # cost of one draw per replicate
    n_draws <- if (is.null(regen)) 0L else {
      nd <- getOption("rasch.efrm_link_draws", max(50L, boot_reps %/% 5L))
      if (length(nd) != 1L || !is.finite(nd) || nd != floor(nd) || nd < 1)
        stop("option rasch.efrm_link_draws must be one positive whole number")
      as.integer(nd)
    }
    draws <- if (n_draws > 0L) vector("list", n_draws)
    reps <- matrix(NA_real_, boot_reps, 2L * S)
    for (r in seq_len(boot_reps)) {
      mats <- if (n_draws == 0L) NULL else {
        di <- ((r - 1L) %% n_draws) + 1L
        if (is.null(draws[[di]])) draws[[di]] <- regen()
        draws[[di]]
      }
      b <- if (is.null(mats)) link_once(sample.int(N, N, replace = TRUE))
           else link_once(sample.int(N, N, replace = TRUE),
                          um = mats$u, wm = mats$w, gm = mats$g)
      if (!is.null(b)) reps[r, ] <- c(b$la, b$mu)
    }
    reps <- reps[stats::complete.cases(reps), , drop = FALSE]
    if (nrow(reps) < 30)
      stop("the unit-linking bootstrap failed in most replicates; the linking ",
           "design is too weak for stable alpha estimation")
    cov_link <- stats::cov(reps)
  }

  list(alpha = setNames(exp(point$la), sets_u),
       se_log_alpha = if (is.null(cov_link)) setNames(rep(NA_real_, S), sets_u)
         else setNames(sqrt(pmax(diag(cov_link)[seq_len(S)], 0)), sets_u),
       mu = setNames(point$mu, sets_u),
       cov_link = cov_link,
       edges = data.frame(set_a = sets_u[vapply(point$edges, `[`, 1L, 1)],
                          set_b = sets_u[vapply(point$edges, `[`, 1L, 2)],
                          n = point$off_n, log_slope = point$ls_est))
}

# Person estimation under unequal frame units: the weighted score
# W = sum_i rho_i x_i is sufficient, so roots are solved once per
# missing-data pattern and unique weighted score.
.efrm_person_estimates <- function(X, tau_list, disc) {
  N <- nrow(X)
  obs <- !is.na(X)
  m <- vapply(tau_list, length, 1L)
  pat <- apply(obs, 1, function(z) paste(which(z), collapse = ","))
  theta <- se <- rep(NA_real_, N)
  raw <- rowSums(X, na.rm = TRUE); raw[rowSums(obs) == 0L] <- NA
  max_raw <- as.numeric(obs %*% m)
  Xw <- sweep(X, 2, disc, "*")
  W <- rowSums(Xw, na.rm = TRUE); W[rowSums(obs) == 0L] <- NA
  Wmax <- as.numeric(obs %*% (disc * m))
  extreme <- !is.na(W) & (W <= 1e-12 | W >= Wmax - 1e-12)

  for (key in unique(pat)) {
    cols <- as.integer(strsplit(key, ",", fixed = TRUE)[[1]])
    if (!length(cols)) next
    sel <- which(pat == key)
    r <- disc[cols]; tl <- tau_list[cols]
    for (Wu in unique(signif(W[sel], 12))) {
      who <- sel[signif(W[sel], 12) == Wu]
      g <- function(th) {
        mo <- lapply(seq_along(cols), function(j)
          item_moments(th, tl[[j]], disc = r[j]))
        E  <- vapply(mo, `[[`, 0, "E");  V <- vapply(mo, `[[`, 0, "V")
        m3 <- vapply(mo, `[[`, 0, "mu3")
        (Wu - sum(r * E)) + sum(r^3 * m3) / (2 * sum(r^2 * V))
      }
      root <- tryCatch(uniroot(g, c(-30, 30), tol = 1e-9)$root,
                       error = function(e) NA_real_)
      theta[who] <- root
      if (!is.na(root)) {
        V <- vapply(seq_along(cols), function(j)
          item_moments(root, tl[[j]], disc = r[j])$V, 0)
        se[who] <- 1 / sqrt(sum(r^2 * V))
      }
    }
  }
  data.frame(n_items = rowSums(obs), raw = raw, max_raw = max_raw,
             weighted_score = W, theta = theta, se = se, extreme = extreme)
}

#' Fit the extended frame of reference model
#'
#' Fits Humphry's extended frame of reference model, in which the unit can
#' differ across item-set by person-group frames. For item \eqn{i} in set
#' \eqn{s} and person \eqn{n} in group \eqn{g},
#' \deqn{P(X_{ni}=x)=\frac{\exp\{\rho_{sg}[x\theta_n-
#'   \sum_{k=1}^{x}\delta_{ik}]\}}
#'   {\sum_{y=0}^{m_i}\exp\{\rho_{sg}[y\theta_n-
#'   \sum_{k=1}^{y}\delta_{ik}]\}},\qquad
#'   \rho_{sg}=\alpha_s\phi_g.}
#'
#' @details
#' The partial credit model holds within each frame in its natural unit.
#' Person-group units \eqn{\phi_g} and centred set thresholds are estimated by
#' within-frame pairwise conditional maximum likelihood.
#'
#' The two kinds of unit are estimated by different routes, because the
#' design offers different information about each. A set taken by two
#' person groups gives the same items calibrated at two scales, so the
#' threshold pattern identifies \eqn{\phi_g} conditionally, without
#' reference to person estimates. Item sets partition the items, so no item
#' is calibrated at two scales and \eqn{\alpha_s} has no such channel: it is
#' identified person-side, from the dispersion of the estimates of persons
#' common to two sets. That route carries estimation error into the
#' quantity being compared, which is why it needs the correction described
#' next, and why the set units are less precise than the group units at the
#' same sample size. Placing an item in two sets would create the missing
#' channel but requires administering it twice to the same person, and the
#' resulting carry-over biases the recovered ratio far more than the error
#' it would avoid; overlapping sets are therefore refused.
#'
#' Item-set units \eqn{\alpha_s} and set locations are estimated from persons
#' common to linked sets by true-score variance ratios, with the true-score variance
#' recovered by a truncated-score-moment correction: the mean and variance
#' of the weighted likelihood score map over the non-extreme scores are
#' exact functions of the person location given the fitted thresholds, and
#' their person-distribution expectations are estimated through score
#' weights that are unbiased for any person distribution because the raw
#' score is sufficient. The naive \eqn{var(\hat u) - mean(SE^2)} correction
#' is badly calibrated on short tests (the reported error variance
#' overstates the actual one and weighted likelihood shrinkage makes the
#' errors covary negatively with the locations); its residual distortion
#' biased the log unit ratio upward by about 0.05 at eight dichotomous
#' items per set, confirmed against an external TAM 2PL slope-group
#' anchor, while the corrected estimator is unbiased there. The linking
#' graph must connect all sets to a common scale, and each linking
#' person's response pattern must span a score range of at least 4 within
#' a set (at least four dichotomous items; six or more are recommended)
#' -- shorter patterns cannot support the correction and are refused.
#' The correction computes score distributions from the fitted
#' within-frame model, so it inherits violations of it: within-set
#' discrimination heterogeneity or guessing bias the recovered units
#' roughly in proportion, mistargeting by two logits adds a few per cent,
#' and person distributions concentrated far off-target (such as widely
#' separated modes) can bias the ratio substantially without warning --
#' inspect targeting before trusting units from such designs.
#'
#' The default hybrid standard errors combine the pairwise Godambe covariance,
#' a person bootstrap for set linking, and delta-method propagation. Each
#' linking replicate also draws the within-frame parameters from their joint
#' stage-one covariance: without that redraw the set-unit standard errors
#' understate by about 20\% and the unit tests reject a true null at 9-10\%;
#' with it they reject at 4.9\% over 1,200 simulated replicates, stable
#' across item counts, sample sizes, imbalance, and weak linking, and
#' matching a full-bootstrap benchmark. The corrected set-unit estimator
#' holds this calibration across designs: null size 3--5\% with 93--99\%
#' coverage over 5--15 items per set, unit ratios 1--2, partial credit
#' items, booklet missingness, pairwise-only person overlap, and person
#' skewness to 2.8 -- where it stays unbiased while a normal-population
#' MML anchor drifts. The estimator's fixed-design offset is small --
#' around one per cent on the unit ratio at eight dichotomous items per
#' set, decaying with set length -- and intervals remain calibrated to at
#' least five thousand linking persons at fine Monte Carlo resolution;
#' only near ten thousand does the offset approach the sampling error,
#' and there longer or polytomous sets, not more persons, buy further
#' accuracy. The person-group units are unaffected: estimated in the
#' within-frame conditional stage rather than through the linking
#' variance ratio, they show no detectable bias (0.000 over 60 replicates
#' at twelve items per set) with 95\% interval coverage, and they are the
#' more precise of the two: the conditional channel is close to twice as
#' precise as the linking channel at eight items per frame, narrowing to
#' about 15\% by twenty. A significant group unit beside a
#' non-significant set unit therefore reflects the design as much as the
#' data. With \code{se_method = "bootstrap"}, the complete
#' model is refitted to each person resample and all reported covariance comes
#' from the bootstrap distribution.
#'
#' Screen the items before interpreting a unit. Both units are ratios of
#' quantities estimated within frames, so an item that misbehaves in one
#' frame and not another enters the estimate as though it were a unit
#' difference. The two units are exposed differently, because their frames
#' are built differently.
#'
#' Person-group units are exposed to ordinary DIF. The same items are
#' answered by every group, so an item can sit differently in one group
#' than another: two items in twelve carrying a one-logit shift moved a
#' recovered ratio by 7 to 11\% in simulation, and four items whose
#' discrimination differed by half moved it by 10 to 12\%. Neither the
#' fitted object nor \code{\link{dif_anova}} can detect this. The fit holds
#' one location per item, shared across the frames it appears in and scaled
#' by the frame unit, so its per-frame estimates agree by construction; and
#' the grouping factor that defines the frames is constant within each
#' frame, so \code{dif_anova} refuses it. Testing the assumption means
#' stepping outside the model, which is what \code{\link{frame_invariance}}
#' does: it calibrates each frame separately, puts the locations on the
#' common scale, and compares them item by item.
#'
#' Item-set units are exposed to something else entirely. Sets partition
#' the items, so no item appears in two sets and DIF across sets is not
#' defined. What distorts a set unit is misfit CONCENTRATED IN ONE SET,
#' which the other set carries nothing to offset. This is the largest
#' effect measured anywhere for these estimators: with 8 items per set,
#' two over-discriminating items in one set moved a planted ratio of 1.40
#' to 1.17, two under-discriminating items moved it to 1.73, and four
#' over-discriminating items to 1.02 -- a real 40\% unit difference read
#' as none at all. The same misfit spread evenly across both sets cancels
#' almost exactly (1.41).
#'
#' Screen on the standardised fit residual, \code{fit_resid} in
#' \code{items}, and repair with \code{\link{drop_items}}. Of the fit
#' statistics compared in simulation it is the one that restores the ratio:
#' against a planted 1.40 read as 1.69 under two under-discriminating
#' items, dropping on the fit residual recovered 1.44 where the chi-square
#' item fit test reached only 1.64, and dropping the items actually planted
#' gave 1.43.
#'
#' Rank on it; do not threshold on it. A fixed cut such as
#' \code{abs(fit_resid) > 2} is a statement about detectability rather than
#' about magnitude, so it flags more items the more persons there are: with
#' two of eight items in a set discriminating twice as steeply, the
#' remaining items genuinely depart from the compromise the model settles
#' on, and that departure clears the cut in 0.8 of 14 sound items at 500
#' persons but 7.0 at 6,000. On the Rosenberg Self-Esteem data used by the
#' wording case study, at 6,000 respondents, the cut selects 7 of the 10
#' items and \code{\link{drop_items}} then refuses the drop for emptying a
#' set.
#'
#' Order the items by \code{abs(fit_resid)} instead, drop the largest
#' departure, and refit. Stop when the sets no longer differ in unit, which
#' is a question for the test in \code{efrm_vs_rasch$unit_tests} rather than
#' for the ratio: once there is no difference left, a further drop has
#' nothing to explain and is fitting noise. On the self-esteem data the
#' ranking puts Q8 first at 22.5 and Q4 second at 12.4 -- the same two items
#' a free-slope model fitted to the same respondents ranks lowest, and the
#' two orderings agree on the extreme three. Dropping Q8 alone moves the
#' ratio from 1.24 to 1.03 at p = 0.63; dropping Q4 as well does not improve
#' on that but returns a significant difference of 1.16, which is what
#' ignoring the stopping rule costs.
#'
#' Treat a badly fitting item as evidence about its set membership: an
#' item that shares no unit with its set, such as an ambivalently worded
#' item filed among the negatively worded ones, is better removed than
#' retained, and cannot be given a set of its own because a single item
#' carries no dispersion to estimate a unit from.
#'
#' The dichotomous model and the theory of frame-dependent units follow
#' Humphry (2005) and Humphry and Andrich (2008). The item-set linking step is
#' an error-corrected method-of-moments implementation of the variance-ratio
#' argument in Humphry (2005), rather than the likelihood proposed in section
#' 5.3 of that thesis. The polytomous, multigroup, and crossed-frame forms are
#' extensions implemented in this package. For these designs, the full
#' bootstrap gives the least conditional account of uncertainty.
#'
#' @param data Persons-by-items data (matrix or data frame, like
#'   \code{\link{rasch}}), plus a person-group column.
#' @param item_sets A named list mapping set names to item-column names, or
#'   a named character vector mapping item names to set names. Items not
#'   mentioned form their own set \code{"(rest)"} when a list is given.
#' @param groups Name of the person-group column in \code{data}, or a vector
#'   with one entry per person.
#'   Several columns define crossed group cells. Their units are returned in
#'   \code{phi_table}; \code{phi_factorial} and
#'   \code{phi_factorial_tests} contain the GLS factorial decomposition and
#'   omnibus Wald tests. Structurally unidentified units are refused. Very
#'   imprecise but identified units are retained with a warning.
#' @param id,factors,items,n_groups,adjust_N,na_codes As in
#'   \code{\link{rasch}}.
#' @param maxit,tol Outer iteration cap and convergence tolerance of the
#'   bilinear pairwise stage.
#' @param min_link_persons Minimum number of common persons required for a
#'   set pair to contribute to the unit linking.
#' @param se_method \code{"hybrid"} (sandwich + linking bootstrap + delta
#'   propagation; fast, default) or \code{"bootstrap"} (full person
#'   bootstrap of all stages).
#' @param boot_reps Bootstrap replicates; defaults to 300 for the linking
#'   bootstrap and 200 for the full bootstrap.
#' @return An object of classes \code{"rasch_efrm"} and \code{"rasch"}.
#'   Model-specific components include \code{frames}, \code{phi_table},
#'   \code{alpha_table}, \code{set_table}, common-unit item and threshold
#'   tables, group-specific \code{score_curves}, \code{efrm_vs_rasch}, and
#'   \code{linking}. See the extended frame of reference vignette for their
#'   interpretation.
#' @references
#' Andrich, D. (1982). An extension of the Rasch model for ratings providing
#' both location and dispersion parameters. Psychometrika, 47(1), 105--113.
#'
#' Andrich, D. and Luo, G. (2003). Conditional pairwise estimation in the
#' Rasch model for ordered response categories using principal components.
#' Journal of Applied Measurement, 4(3), 205--221.
#'
#' Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
#' Measuring in the Educational, Social and Health Sciences. Springer.
#'
#' Humphry, S. M. (2005). Maintaining a Common Arbitrary Unit in Social
#' Measurement. PhD thesis, Murdoch University.
#'
#' Humphry, S. M. (2010). Modeling the effects of person group factors on
#' discrimination. Educational and Psychological Measurement, 70(2),
#' 215--231.
#'
#' Humphry, S. M. (2012). Item set discrimination and the unit in the Rasch
#' model. Journal of Applied Measurement, 13(2), 165--180.
#'
#' Montuoro, P. and Humphry, S. M. (2024). Modeling the effect of reading
#' item clarity on item discrimination. Journal of Applied Measurement,
#' 24(3/4), 121--132.
#'
#' Humphry, S. M. and Andrich, D. (2008). Understanding the unit in the Rasch
#' model. Journal of Applied Measurement, 9(3), 249--264.
#' @seealso \code{\link{frame_invariance}}, which tests the item invariance
#'   this model assumes rather than imposing it, and \code{\link{drop_items}},
#'   which removes an item the test flags and refits. Also
#'   \code{\link{rasch}}, \code{\link{rasch_mfrm}},
#'   \code{\link{test_information}}, and \code{\link{simulate_efrm}}.
#' @examples
#' \donttest{
#' set.seed(1); Np <- 400
#' simP <- function(th, tau, r) { x <- 0:length(tau)
#'   p <- exp(r * (x * th - c(0, cumsum(tau)))); p / sum(p) }
#' grp <- rep(c("A", "B"), each = Np / 2)
#' phi <- c(A = 0.8, B = 1.25)
#' d <- seq(-1.5, 1.5, length.out = 10)
#' theta <- rnorm(Np)
#' X <- sapply(seq_along(d), function(i) sapply(seq_len(Np), function(n)
#'   sample(0:1, 1, prob = simP(theta[n], d[i], phi[grp[n]]))))
#' colnames(X) <- sprintf("I%02d", seq_along(d))
#' fit <- rasch_efrm(data.frame(X, grp = grp), item_sets = list(core = colnames(X)),
#'                   groups = "grp")
#' fit$phi_table
#' }
#' @export
rasch_efrm <- function(data, item_sets, groups, id = NULL, factors = NULL,
                       items = NULL, n_groups = NULL, adjust_N = NA,
                       na_codes = -1, maxit = 50, tol = 1e-7,
                       min_link_persons = 30,
                       se_method = c("hybrid", "bootstrap"),
                       boot_reps = NULL) {
  se_method <- match.arg(se_method)
  if (is.null(boot_reps)) boot_reps <- if (se_method == "hybrid") 300L else 200L
  # --- roles ----------------------------------------------------------------
  id_vec <- NULL; fac_df <- NULL; grp <- NULL; grp_name <- "group"
  grp_components <- NULL
  if (is.data.frame(data)) {
    nm <- names(data)
    if (is.character(groups) && length(groups) != nrow(data)) {
      # column name(s); a character vector of length nrow(data) is the
      # group values themselves and is handled below
      miss <- setdiff(groups, nm)
      if (length(miss))
        stop("group column(s) not found in the data: ",
             paste(miss, collapse = ", "))
      if (length(groups) == 1L) {
        grp <- data[[groups]]; grp_name <- groups
      } else {
        # SEVERAL frame-defining factors: the frames are their crossed
        # cells, and a factorial decomposition of the cell units is
        # reported in phi_factorial
        grp_components <- data[groups]
        grp <- interaction(grp_components, sep = ":", drop = TRUE)
        grp_name <- paste(groups, collapse = ":")
      }
    }
    if (is.character(id) && length(id) == 1L && id %in% nm) id_vec <- data[[id]]
    else if (!is.null(id) && length(id) == nrow(data)) id_vec <- id
    if (is.character(factors)) {
      missf <- setdiff(factors, nm)
      if (length(missf))
        stop("factor column(s) not found in the data: ",
             paste(missf, collapse = ", "))
      fac_df <- data[, factors, drop = FALSE]
    } else if (is.data.frame(factors) && nrow(factors) == nrow(data))
      fac_df <- factors
    drop_cols <- c(if (is.character(id)) id else NULL,
                   if (is.character(factors)) factors
                   else if (is.data.frame(factors)) intersect(names(factors), nm),
                   if (is.character(groups)) groups else NULL)
    item_cols <- if (!is.null(items)) intersect(items, nm) else setdiff(nm, drop_cols)
    X <- as.matrix(data[, item_cols, drop = FALSE])
  } else {
    X <- as.matrix(data)
  }
  if (is.null(grp)) {
    if (length(groups) != nrow(X))
      stop("'groups' must name a column of data or give one value per person")
    grp <- groups
  }
  grp <- factor(grp)
  if (nlevels(grp) < 1L) stop("no person groups found")
  if (is.null(id_vec)) id_vec <- seq_len(nrow(X))

  prep <- .prepare_X(X, na_codes = na_codes); X <- prep$X
  notes <- prep$notes
  m_item <- apply(X, 2, max, na.rm = TRUE); L <- ncol(X)

  # --- item sets --------------------------------------------------------------
  if (is.list(item_sets)) {
    if (is.null(names(item_sets)) || any(!nzchar(names(item_sets))))
      stop("item_sets must be a NAMED list (each element a set of item ",
           "names)")
    ov <- unlist(item_sets)
    if (anyDuplicated(ov))
      stop("item(s) assigned to more than one set: ",
           paste(unique(ov[duplicated(ov)]), collapse = ", "))
    unknown <- setdiff(ov, colnames(X))
    if (length(unknown))
      stop("item_sets name item(s) not in the data: ",
           paste(unknown, collapse = ", "))
    set_of <- rep(NA_character_, L); names(set_of) <- colnames(X)
    for (s in names(item_sets)) {
      hit <- intersect(item_sets[[s]], colnames(X))
      set_of[hit] <- s
    }
    if (anyNA(set_of)) set_of[is.na(set_of)] <- "(rest)"
  } else {
    if (is.null(names(item_sets))) stop("item_sets must be a named list or named vector")
    set_of <- as.character(item_sets)[match(colnames(X), names(item_sets))]
    if (anyNA(set_of))
      stop("item(s) missing from the item_sets map: ",
           paste(colnames(X)[is.na(set_of)], collapse = ", "))
    names(set_of) <- colnames(X)
  }
  sets_u <- sort(unique(set_of)); S <- length(sets_u)
  # STRING-SORTED group levels, matching .efrm_solve's internal ordering
  # exactly: interaction() orders levels first-factor-fastest, the solver
  # sorts labels as strings, and attaching the solver's phi positionally
  # to a differently ordered glevs swapped crossed cells (g2:N took
  # g1:S's unit in a 2-by-2, manufacturing a spurious factor effect)
  glevs <- sort(levels(grp)); G <- length(glevs)

  # --- virtual columns: item x group ------------------------------------------
  vmap <- expand.grid(item = colnames(X), group = glevs,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  vmap$set <- set_of[vmap$item]
  vmap$vkey <- paste(vmap$item, vmap$group, sep = ":")
  Xv <- matrix(NA_integer_, nrow(X), nrow(vmap),
               dimnames = list(NULL, vmap$vkey))
  for (v in seq_len(nrow(vmap))) {
    sel <- which(as.character(grp) == vmap$group[v])
    Xv[sel, v] <- X[sel, vmap$item[v]]
  }
  keep <- colSums(!is.na(Xv)) > 0L
  Xv <- Xv[, keep, drop = FALSE]; vmap <- vmap[keep, , drop = FALSE]
  rownames(vmap) <- NULL
  m_v <- m_item[vmap$item]
  thr_v <- threshold_index(m_v)

  # delta enumeration over underlying items (ordered by set, then item)
  iord <- order(match(set_of, sets_u), colnames(X))
  items_o <- colnames(X)[iord]
  thr_items <- threshold_index(m_item[items_o])
  Md <- nrow(thr_items)
  # map each virtual threshold row to its delta row
  drow <- vapply(seq_len(nrow(thr_v)), function(r) {
    it <- vmap$item[thr_v$item[r]]
    which(thr_items$item == match(it, items_o) & thr_items$k == thr_v$k[r])
  }, 1L)
  # per-set sum-zero blocks for dtilde
  set_of_drow <- set_of[items_o][thr_items$item]
  A_D <- matrix(0, Md, Md - S)
  cursor <- 0L
  for (s in sets_u) {
    rows <- which(set_of_drow == s); Ms <- length(rows)
    if (Ms < 2L) stop("set '", s, "' needs at least two thresholds")
    A_D[rows, cursor + seq_len(Ms - 1L)] <- rbind(diag(Ms - 1L), rep(-1, Ms - 1L))
    cursor <- cursor + Ms - 1L
  }

  # --- pairwise stage ----------------------------------------------------------
  pairs <- .efrm_filter_pairs(.pair_counts(Xv, m_v), vmap)
  if (!length(pairs)) stop("no informative within-frame item pairs")

  # phi-link check: two groups are joined only when they SHARE at least two
  # observed items of a common set. Sharing a set label alone is not enough:
  # the unit ratio phi_g/phi_h is identified through the common structural
  # thresholds of items both groups answered, and a spread comparison needs
  # at least two of them (disjoint item subsets of the same set leave the
  # ratio unidentified even though the old set-level check passed).
  has_data <- colSums(!is.na(Xv)) > 0L
  edges_g <- list()
  for (s in unique(vmap$set)) {
    by_grp <- lapply(glevs, function(g) {
      sel <- vmap$set == s & vmap$group == g & has_data
      unique(vmap$item[sel])
    })
    for (g1 in seq_along(glevs)) for (g2 in seq_len(g1 - 1L)) {
      if (length(intersect(by_grp[[g1]], by_grp[[g2]])) >= 2L)
        edges_g[[length(edges_g) + 1L]] <- c(g1, g2)
    }
  }
  comp_g <- .efrm_components(G, edges_g)
  if (length(unique(comp_g)) > 1L)
    stop("person groups are not linked by any common item set; relative units ",
         "(phi) are unidentified between: ",
         paste(tapply(glevs, comp_g, paste, collapse = "+"), collapse = " | "))

  sol <- .efrm_solve(Xv, thr_v, m_v, vmap, pairs, drow, A_D,
                     maxit = maxit, tol = tol)
  # STRUCTURAL non-identification (a flat direction of the information
  # along a unit) is an error: every common-unit quantity would silently
  # depend on the arbitrary point the optimiser stopped at. PRACTICAL
  # weakness -- an analytic SE above 5 log-units, i.e. the unit uncertain
  # beyond a factor of exp(5) ~ 150 -- keeps its estimate for sensitivity
  # work but warns loudly and is noted: the SE already says the data
  # carry essentially no unit information.
  bad_g <- sol$phi_unident | !is.finite(sol$se_log_phi)
  if (length(glevs) > 1L && any(bad_g))
    stop("the unit(s) of group(s) ", paste(glevs[bad_g], collapse = ", "),
         " are unidentified: the thresholds in their frames carry no ",
         "usable spread for phi to scale (information rank/conditioning ",
         "failure) -- refit without these groups, or include items whose ",
         "difficulties differ within the sets they answer")
  weak_g <- !bad_g & sol$se_log_phi > 5
  if (length(glevs) > 1L && any(weak_g)) {
    warning("the unit(s) of group(s) ",
            paste(glevs[weak_g], collapse = ", "), " are only weakly ",
            "identified (SE of log phi above 5): the estimates are kept ",
            "for sensitivity work, but the data carry essentially no ",
            "information about these units", call. = FALSE)
    notes <- c(notes, paste0(
      "weakly identified unit(s) for group(s) ",
      paste(glevs[weak_g], collapse = ", "),
      ": SE of log phi exceeds 5, so no unit inference is supportable"))
  }
  phi <- sol$phi; dtil <- sol$dtilde

  # equal-unit comparison on the same conditional information
  B0 <- A_D[drow, , drop = FALSE]
  bd0 <- qr.coef(qr(A_D), dtil); bd0[is.na(bd0)] <- 0
  glh0 <- .pcml_glh(drop(B0 %*% bd0), thr_v, pairs, m_v)
  for (it0 in 1:25) {
    gb <- drop(crossprod(B0, glh0$g)); Hb <- crossprod(B0, glh0$H %*% B0)
    step <- tryCatch(solve(Hb, gb), error = function(e)
      solve(Hb - diag(1e-8, nrow(Hb)), gb))
    lam <- 1; moved <- FALSE
    for (half in 1:30) {
      cand <- bd0 - lam * step
      g2 <- .pcml_glh(drop(B0 %*% cand), thr_v, pairs, m_v)
      if (is.finite(g2$ll) && g2$ll >= glh0$ll - 1e-12) {
        bd0 <- cand; glh0 <- g2; moved <- TRUE; break
      }
      lam <- lam / 2
    }
    if (!moved || max(abs(lam * step)) < 1e-7) break
  }

  # --- person-side linking (alpha, mu) ----------------------------------------
  # one builder for every linking path (point estimate, calibration-redraw
  # pool, outer bootstrap): person estimates u plus the correction moments
  # w = phi_w(R) and g = phi_g(R) per set (.person_link_moments), kept
  # coherent so every bootstrap replicate re-runs the corrected computation
  person_mats <- function(Xm, dt, phiv) {
    u <- w <- g <- matrix(NA_real_, nrow(Xm), S,
                          dimnames = list(NULL, sets_u))
    for (si in seq_len(S)) for (gl in glevs) {
      cols <- which(vmap$set == sets_u[si] & vmap$group == gl)
      if (!length(cols)) next
      # thresholds in dtilde units for these virtual columns
      tl <- lapply(cols, function(v) dt[drow[which(thr_v$item == v)]])
      Xs <- Xm[, cols, drop = FALSE]
      pe <- .person_estimates(Xs, tl, disc = phiv[gl])
      sel <- which(pe$n_items > 0 & !pe$extreme)
      u[sel, si] <- pe$theta[sel]
      lm <- .person_link_moments(Xs, tl, disc = phiv[gl])
      w[sel, si] <- lm$w[sel]; g[sel, si] <- lm$g[sel]
    }
    list(u = u, w = w, g = g)
  }
  if (S > 1L) {
    pm <- person_mats(Xv, dtil, phi)
    # calibration-noise propagation into the linking bootstrap: a
    # symmetric square root of each stage-1 covariance, used to redraw
    # (dtilde, log phi) per bootstrap replicate under the documented
    # item-side/person-side independence treatment
    mat_sqrt <- function(V) {
      if (is.null(V) || any(!is.finite(V))) return(NULL)
      ee <- eigen((V + t(V)) / 2, symmetric = TRUE)
      ee$vectors %*% (t(ee$vectors) * sqrt(pmax(ee$values, 0)))
    }
    K_dt <- length(dtil)
    # simulation-only escape hatches for the validation studies:
    # rasch.efrm_link_blockdiag = TRUE reproduces the marginal
    # (cross-covariance-free) draw for paired comparison with the joint
    # draw; rasch.efrm_link_draws overrides the parameter-draw pool size
    cj <- sol$cov_joint
    if (isTRUE(getOption("rasch.efrm_link_blockdiag", FALSE)) &&
        !is.null(cj)) {
      cj[seq_len(K_dt), K_dt + seq_along(phi)] <- 0
      cj[K_dt + seq_along(phi), seq_len(K_dt)] <- 0
    }
    L_j <- mat_sqrt(cj)
    regen <- if (!is.null(L_j)) function() {
      v <- drop(L_j %*% stats::rnorm(ncol(L_j)))
      person_mats(Xv, dtil + v[seq_len(K_dt)],
                  phi * exp(v[K_dt + seq_along(phi)]))
    } else NULL
    link <- .efrm_link_sets(pm$u, pm$w, pm$g, sets_u, min_link_persons,
                            boot_reps = boot_reps, regen = regen)
    alpha <- link$alpha; mu <- link$mu
  } else {
    alpha <- setNames(1, sets_u); mu <- setNames(0, sets_u)
    link <- list(alpha = alpha, se_log_alpha = setNames(0, sets_u), mu = mu,
                 cov_link = NULL, edges = data.frame())
  }

  # --- optional full person bootstrap of all stages ----------------------------
  boot <- NULL
  if (se_method == "bootstrap") {
    Npers <- nrow(Xv)
    boot_replicate <- function(idx) {
      Xb <- Xv[idx, , drop = FALSE]
      pb <- .efrm_filter_pairs(.pair_counts(Xb, m_v), vmap)
      if (!length(pb)) return(NULL)
      sb <- .efrm_solve(Xb, thr_v, m_v, vmap, pb, drow, A_D,
                        maxit = maxit, tol = tol)
      if (S > 1L) {
        pm_b <- person_mats(Xb, sb$dtilde, sb$phi)
        lb <- .efrm_link_sets(pm_b$u, pm_b$w, pm_b$g, sets_u,
                              min_link_persons, boot_reps = 0)
        ab <- lb$alpha; mb <- lb$mu
      } else { ab <- setNames(1, sets_u); mb <- setNames(0, sets_u) }
      db <- sb$dtilde / ab[set_of_drow] + mb[set_of_drow]
      c(log(sb$phi), log(ab), mb, db)
    }
    collect <- matrix(NA_real_, boot_reps, G + 2L * S + Md)
    for (r in seq_len(boot_reps)) {
      res <- tryCatch(boot_replicate(sample.int(Npers, Npers, replace = TRUE)),
                      error = function(e) NULL)
      if (!is.null(res)) collect[r, ] <- res
    }
    collect <- collect[stats::complete.cases(collect), , drop = FALSE]
    if (nrow(collect) < max(30, boot_reps / 2)) {
      warning("the full bootstrap failed in most replicates; ",
              "falling back to hybrid standard errors")
    } else boot <- collect
  }

  # --- assembly in arbitrary units ----------------------------------------------
  delta <- dtil / alpha[set_of_drow] + mu[set_of_drow]
  rho_v <- alpha[vmap$set] * phi[vmap$group]
  thr_v$tau <- delta[drow]

  # covariance of delta = dtilde/alpha + mu: the pairwise (conditional)
  # component plus the unit component propagated by the delta method; the
  # item-side and person-side information are treated as independent
  cov_delta <- sol$cov_dtilde / tcrossprod(alpha[set_of_drow])
  if (!is.null(boot)) {
    cov_delta <- stats::cov(boot[, G + 2L * S + seq_len(Md), drop = FALSE])
  } else if (S > 1L && !is.null(link$cov_link)) {
    Sidx <- match(set_of_drow, sets_u)
    Caa <- link$cov_link[seq_len(S), seq_len(S), drop = FALSE]
    Cam <- link$cov_link[seq_len(S), S + seq_len(S), drop = FALSE]
    Cmm <- link$cov_link[S + seq_len(S), S + seq_len(S), drop = FALSE]
    g1 <- -(dtil / alpha[set_of_drow])     # d delta / d log alpha
    T2 <- matrix(g1, Md, Md) * Cam[Sidx, Sidx, drop = FALSE]
    cov_delta <- cov_delta + tcrossprod(g1) * Caa[Sidx, Sidx, drop = FALSE] +
      T2 + t(T2) + Cmm[Sidx, Sidx, drop = FALSE]
  }
  cov_tau <- cov_delta[drow, drow, drop = FALSE]
  thr_v$se <- sqrt(pmax(diag(cov_tau), 0))
  thr_v$anchored <- FALSE
  # a virtual item threshold on a near-empty category is a boundary
  # artefact here too: flag it and report NA rather than a covariance-based
  # number, the same honesty rasch()/pcml()/rasch_mfrm() apply
  weak <- .pcml_weak_thresholds(Xv, m_v, thr_v, colnames(Xv))
  thr_v$weak <- weak$flag
  thr_v$se[weak$flag] <- NA_real_
  if (length(weak$notes)) notes <- c(notes, weak$notes)
  est <- list(model = "EFRM", thr = thr_v, cov_tau = cov_tau,
              loglik = sol$loglik, iterations = sol$iterations,
              converged = sol$converged, m = m_v, anchors = NULL,
              n_parameters = (Md - S) + (G - 1L) + 2L * (S - 1L))

  fac_all <- data.frame(g = as.character(grp), stringsAsFactors = FALSE)
  names(fac_all) <- grp_name
  if (!is.null(grp_components)) {
    gc_chr <- as.data.frame(lapply(grp_components, as.character),
                            stringsAsFactors = FALSE)
    fac_all <- cbind(fac_all, gc_chr)
  }
  if (!is.null(fac_df)) fac_all <- cbind(fac_all, fac_df)
  fit <- .assemble_fit("EFRM", Xv, est, id_vec, fac_all, n_groups, adjust_N,
                       notes, disc = rho_v)
  # every frame-defining factor (the crossed cell and its components) is
  # frame structure, not a testable DIF factor
  fit$frame_group <- c(grp_name,
                       if (!is.null(grp_components)) names(grp_components))

  # --- structural tables -----------------------------------------------------------
  se_lp <- if (!is.null(boot)) apply(boot[, seq_len(G), drop = FALSE], 2, sd)
           else unname(sol$se_log_phi)
  se_la <- if (!is.null(boot)) apply(boot[, G + seq_len(S), drop = FALSE], 2, sd)
           else unname(link$se_log_alpha)
  fit$phi_table <- data.frame(group = glevs, phi = unname(phi),
                              se_log_phi = se_lp)
  # factorial decomposition of the cell units: generalised least squares
  # of log phi_cell on sum-coded main effects (and the interaction when
  # every cell is observed), using the JOINT covariance of the cell
  # log-units -- from the bootstrap replicates when available, else the
  # solver's centred analytic covariance. The centring constraint
  # (sum log phi = 0) makes the cells dependent and the covariance
  # singular along the constant, so a spectral pseudo-inverse is used.
  # Coefficients are reported descriptively (estimate, SE); inference is
  # carried by MULTI-DEGREE-OF-FREEDOM Wald tests per term in
  # phi_factorial_tests (coefficient-wise z would not be an omnibus test
  # for factors with more than two levels).
  if (!is.null(grp_components)) {
    cell_lab <- strsplit(glevs, ":", fixed = TRUE)
    dd <- as.data.frame(do.call(rbind, cell_lab), stringsAsFactors = FALSE)
    names(dd) <- names(grp_components)
    for (cn in names(dd)) dd[[cn]] <- factor(dd[[cn]])
    estimable <- all(vapply(dd, nlevels, 1L) >= 2L) &&
      nrow(dd) > sum(vapply(dd, nlevels, 1L) - 1L) &&
      !anyNA(fit$phi_table$phi)
    if (estimable) {
      ctr <- stats::setNames(rep(list("contr.sum"), ncol(dd)), names(dd))
      full <- nrow(dd) >= prod(vapply(dd, nlevels, 1L))
      fml <- stats::as.formula(paste("~", paste(names(dd), collapse =
        if (full && ncol(dd) > 1L) " * " else " + ")))
      Xf <- stats::model.matrix(fml, dd, contrasts.arg = ctr)
      # the centring constraint (sum log phi = 0) makes the intercept
      # inestimable -- and it lies in the null space of the centred
      # covariance, so keeping it would make the GLS cross-product
      # singular. The sum-coded effects describe the centred cells fully.
      asg <- attr(Xf, "assign")
      Xf <- Xf[, asg != 0L, drop = FALSE]
      asg <- asg[asg != 0L]
      Sig <- if (!is.null(boot))
        stats::cov(boot[, seq_len(G), drop = FALSE])
      else sol$cov_log_phi
      eg <- eigen(Sig, symmetric = TRUE)
      pos <- eg$values > max(eg$values) * 1e-8
      Sinv <- eg$vectors[, pos, drop = FALSE] %*%
        (t(eg$vectors[, pos, drop = FALSE]) / eg$values[pos])
      XtS <- t(Xf) %*% Sinv
      V <- tryCatch(solve(XtS %*% Xf), error = function(e) NULL)
      if (!is.null(V)) {
        cf <- drop(V %*% (XtS %*% log(fit$phi_table$phi)))
        fit$phi_factorial <- data.frame(
          term = colnames(Xf), log_unit = cf,
          se = sqrt(pmax(diag(V), 0)), stringsAsFactors = FALSE)
        tls <- attr(stats::terms(fml), "term.labels")
        tests <- list()
        for (tno in seq_along(tls)) {
          ii <- which(asg == tno)
          if (!length(ii)) next
          Vt <- V[ii, ii, drop = FALSE]
          W <- tryCatch(drop(t(cf[ii]) %*% solve(Vt) %*% cf[ii]),
                        error = function(e) NA_real_)
          tests[[length(tests) + 1L]] <- data.frame(
            term = tls[tno], df = length(ii), wald = W,
            p = stats::pchisq(W, length(ii), lower.tail = FALSE),
            stringsAsFactors = FALSE)
        }
        fit$phi_factorial_tests <- do.call(rbind, tests)
      }
    }
  }
  fit$alpha_table <- data.frame(set = sets_u, alpha = unname(alpha),
                                se_log_alpha = se_la)
  fit$set_table <- data.frame(set = sets_u, mu = unname(mu),
                              alpha = unname(alpha),
                              n_items = as.integer(table(set_of)[sets_u]))
  # a distinct threshold is weak if any virtual threshold folded into it
  # rests on a near-empty category; carry that to the common-unit tables
  weak_d <- logical(length(delta))
  agg_w <- tapply(weak$flag, drow, any)
  weak_d[as.integer(names(agg_w))] <- as.logical(agg_w)
  se_arb <- sqrt(pmax(diag(cov_delta), 0)); se_arb[weak_d] <- NA_real_
  fit$thresholds_arbitrary <- data.frame(item = items_o[thr_items$item],
                                         set = set_of_drow,
                                         k = thr_items$k, delta = delta,
                                         se = se_arb, weak = weak_d)
  fit$item_arbitrary <- do.call(rbind, lapply(seq_along(items_o), function(i) {
    rows <- which(thr_items$item == i)
    weak_i <- any(weak_d[rows])
    data.frame(item = items_o[i], set = set_of[items_o[i]],
               location = mean(delta[rows]),
               se = if (weak_i) NA_real_ else
                 sqrt(max(mean(cov_delta[rows, rows, drop = FALSE]), 0)),
               weak = weak_i)
  }))
  rownames(fit$item_arbitrary) <- NULL

  # frame table with pooled fit
  fr <- unique(vmap[, c("set", "group")])
  fit$frames <- do.call(rbind, lapply(seq_len(nrow(fr)), function(j) {
    cols <- which(vmap$set == fr$set[j] & vmap$group == fr$group[j])
    npers <- sum(rowSums(!is.na(Xv[, cols, drop = FALSE])) > 0)
    gf <- .group_col_fit(fit$residuals, fit$moments, cols, disc = rho_v,
                         extreme = fit$person$extreme,
                         f_cell = fit$summary_stats$df_factor)
    data.frame(set = fr$set[j], group = fr$group[j],
               n_persons = npers, n_items = length(cols),
               alpha = unname(alpha[fr$set[j]]), phi = unname(phi[fr$group[j]]),
               rho = unname(alpha[fr$set[j]] * phi[fr$group[j]]),
               # with bootstrap replicates the SE of log rho comes from
               # the JOINT draws (log phi + log alpha per replicate), so
               # cross-stage dependence is captured; the analytic fallback
               # has no cross-stage covariance and treats the two stages
               # as uncorrelated, which the frames documentation states
               se_log_rho = if (!is.null(boot))
                 stats::sd(boot[, match(fr$group[j], glevs)] +
                           boot[, G + match(fr$set[j], sets_u)])
               else sqrt(
                 fit$alpha_table$se_log_alpha[match(fr$set[j], sets_u)]^2 +
                 fit$phi_table$se_log_phi[match(fr$group[j], glevs)]^2),
               origin = unname(mu[fr$set[j]]),
               infit_ms = gf$infit_ms, outfit_ms = gf$outfit_ms,
               fit_resid = gf$fit_resid, n_responses = gf$n)
  }))
  rownames(fit$frames) <- NULL

  # equal-unit comparison and score curves. The pairwise comparison is only
  # informative about the group units (phi): the within-frame likelihood is
  # invariant to the set units (alpha), which are identified person-side, so
  # their evidence is the Wald test on log alpha.
  wald_zero <- function(est, Sigma, term) {
    if (length(est) < 2L || is.null(Sigma) || any(!is.finite(Sigma)))
      return(NULL)
    ee <- eigen((Sigma + t(Sigma)) / 2, symmetric = TRUE)
    cut <- max(abs(ee$values)) * 1e-8
    use <- ee$values > cut
    if (!any(use)) return(NULL)
    Sinv <- ee$vectors[, use, drop = FALSE] %*%
      (t(ee$vectors[, use, drop = FALSE]) / ee$values[use])
    W <- drop(t(est) %*% Sinv %*% est)
    data.frame(term = term, df = sum(use), wald = W,
               p = stats::pchisq(W, sum(use), lower.tail = FALSE))
  }
  Sig_phi <- if (!is.null(boot))
    stats::cov(boot[, seq_len(G), drop = FALSE]) else sol$cov_log_phi
  Sig_alpha <- if (S > 1L) {
    if (!is.null(boot)) stats::cov(boot[, G + seq_len(S), drop = FALSE])
    else if (!is.null(link$cov_link))
      link$cov_link[seq_len(S), seq_len(S), drop = FALSE]
    else NULL
  } else NULL
  unit_omnibus <- do.call(rbind, Filter(Negate(is.null), list(
    if (G > 1L) wald_zero(log(phi), Sig_phi, "group units (phi)"),
    if (S > 1L) wald_zero(log(alpha), Sig_alpha, "set units (alpha)"))))

  ut <- rbind(
    if (G > 1L) data.frame(parameter = paste0("log phi[", glevs, "]"),
                           estimate = log(fit$phi_table$phi),
                           se = fit$phi_table$se_log_phi),
    if (S > 1L) data.frame(parameter = paste0("log alpha[", sets_u, "]"),
                           estimate = log(alpha),
                           se = fit$alpha_table$se_log_alpha))
  if (!is.null(ut)) {
    ut$z <- ut$estimate / ut$se
    ut$p <- 2 * pnorm(-abs(ut$z))
    ut$p_adj <- stats::p.adjust(ut$p, method = "holm")
    ut$significant <- ut$p_adj < 0.05
    rownames(ut) <- NULL
  }
  fit$unit_cov <- list(cov_dtilde = sol$cov_dtilde,
                       cov_log_phi = sol$cov_log_phi,
                       cov_joint = sol$cov_joint)
  fit$efrm_vs_rasch <- list(ll_efrm = sol$loglik, ll_equal = glh0$ll,
                            two_delta_ll = 2 * (sol$loglik - glh0$ll),
                            extra_parameters = G - 1L,
                            informative_for = if (G > 1L) "group units (phi)"
                              else "nothing: single group, set units are identified person-side",
                            unit_omnibus = unit_omnibus,
                            unit_tests = ut)
  grid <- seq(-6, 6, by = 0.1)
  fit$score_curves <- do.call(rbind, lapply(glevs, function(g) {
    r_i <- alpha[set_of] * phi[g]
    ew <- vapply(grid, function(th) sum(vapply(seq_len(L), function(i)
      r_i[i] * item_moments(th, delta[thr_items$item == match(colnames(X)[i], items_o)],
                            disc = r_i[i])$E, 0)), 0)
    info <- vapply(grid, function(th) sum(vapply(seq_len(L), function(i)
      r_i[i]^2 * item_moments(th, delta[thr_items$item == match(colnames(X)[i], items_o)],
                              disc = r_i[i])$V, 0)), 0)
    data.frame(group = g, theta = grid, expected_score = ew, sem = 1 / sqrt(info))
  }))
  fit$linking <- list(phi_edges = edges_g, alpha_edges = link$edges)
  fit$se_method <- if (!is.null(boot)) "bootstrap" else "hybrid"
  fit$boot_reps_used <- if (!is.null(boot)) nrow(boot) else NA_integer_
  fit$virtual_map <- vmap
  fit$set_of <- set_of
  fit <- .tag_tables(fit)
  class(fit) <- c("rasch_efrm", "rasch")
  fit
}

#' @export
print.rasch_efrm <- function(x, ...) {
  cat(sprintf("rasch extended frame of reference analysis: %d items in %d set(s) x %d group(s) = %d frames, %d persons\n",
              length(x$set_of), nrow(x$alpha_table), nrow(x$phi_table),
              nrow(x$frames), nrow(x$X)))
  cat(sprintf("Within-frame pairwise conditional ML: %s in %d iterations\n",
              if (x$est$converged) "converged" else "NOT converged",
              x$est$iterations))
  cat(sprintf("PSI %.3f, power of fit: %s\n", x$psi$PSI, x$power_of_fit))
  cat("\nPerson group units (phi):\n")
  print(x$phi_table, digits = 3, row.names = FALSE)
  cat("\nItem set units (alpha) and locations:\n")
  print(.fmt_df(merge(x$alpha_table, x$set_table[, c("set", "mu", "n_items")],
                      by = "set", sort = FALSE)), row.names = FALSE)
  cat(sprintf("\nEqual-unit comparison: 2(ll_EFRM - ll_equal) = %.3f with %d extra unit parameter(s)\n",
              x$efrm_vs_rasch$two_delta_ll, x$efrm_vs_rasch$extra_parameters))
  cat("(composite likelihood: descriptive; informative for ",
      x$efrm_vs_rasch$informative_for, ")\n", sep = "")
  if (!is.null(x$efrm_vs_rasch$unit_omnibus)) {
    cat("Omnibus Wald tests of equal units:\n")
    print(.fmt_df(x$efrm_vs_rasch$unit_omnibus), row.names = FALSE)
  }
  if (!is.null(x$efrm_vs_rasch$unit_tests)) {
    cat("Holm-adjusted exploratory unit contrasts (H0: unit = 1):\n")
    print(.fmt_df(x$efrm_vs_rasch$unit_tests), row.names = FALSE)
  }
  if (length(x$notes)) cat(sprintf("\nNotes: %s\n", paste(x$notes, collapse = "; ")))
  invisible(x)
}

#' Plot frame units
#'
#' Caterpillar plot of the frame units \code{rho_sg = alpha_s phi_g} on the
#' log scale, grouped by item set and coloured by person group, with 95 per
#' cent error bars; frames with pooled fit residuals beyond the band are
#' highlighted.
#'
#' @param fit A fitted object from \code{\link{rasch_efrm}}.
#' @param band Pooled fit residual band beyond which a frame is highlighted.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' \donttest{
#' # see ?rasch_efrm for a complete simulated example
#' }
#' @export
plot_frames <- function(fit, band = 2.5) {
  if (!inherits(fit, "rasch_efrm")) stop("plot_frames needs a rasch_efrm fit")
  fr <- fit$frames[order(fit$frames$set, fit$frames$rho), ]
  n <- nrow(fr)
  lr <- log(fr$rho)
  lo <- lr - 1.96 * fr$se_log_rho; hi <- lr + 1.96 * fr$se_log_rho
  labs <- paste0(fr$set, " \u00d7 ", fr$group)
  glev <- sort(unique(fr$group))
  colr <- .rr$pal[(match(fr$group, glev) - 1L) %% length(.rr$pal) + 1L]
  op <- par(mar = c(4.2, 9, 3.2, 1.5), mgp = c(2.5, 0.7, 0), tcl = -0.25,
            las = 1, col.axis = .rr$ink, col.lab = .rr$ink, col.main = .rr$ink,
            font.main = 2, cex.main = 1.15)
  on.exit(par(op))
  plot(NA, xlim = range(c(lo, hi, 0)) + c(-0.1, 0.1), ylim = c(0.5, n + 0.5),
       xlab = "log unit (log rho)", ylab = "", axes = FALSE, main = "")
  abline(h = seq_len(n), col = .rr$grid, lwd = 0.8)
  abline(v = 0, lty = 2, col = .rr$soft)
  axis(1, col = .rr$grid, col.ticks = .rr$soft)
  axis(2, at = seq_len(n), labels = labs, cex.axis = 0.75,
       col = .rr$grid, col.ticks = NA)
  misfit <- !is.na(fr$fit_resid) & abs(fr$fit_resid) > band
  segments(lo, seq_len(n), hi, seq_len(n), lwd = 2.2,
           col = ifelse(misfit, .rr$red, .rr$soft))
  points(lr, seq_len(n), pch = 21, cex = 1.5,
         bg = ifelse(misfit, .rr$red, colr), col = "white", lwd = 1.2)
  .rr_legend("bottomright", glev, pch = 21,
             pt.bg = .rr$pal[seq_along(glev)], col = "white", pt.cex = 1.2)
  if (any(misfit))
    mtext(sprintf("%d frame(s) with |pooled fit residual| > %.1f", sum(misfit), band),
          side = 3, line = 0.2, adj = 0, cex = 0.8, col = .rr$red)
  invisible(NULL)
}

#' Plot an item's characteristic curves across frames
#'
#' Plots the model expected-score curve for one item in each person group, with
#' observed class-interval means overlaid. Differences between the curves
#' reflect the fitted group units.
#'
#' @param fit A fitted object from \code{\link{rasch_efrm}}.
#' @param item Underlying item name.
#' @param n_groups Number of class intervals for the observed means.
#' @param grid Logit grid.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' \donttest{
#' # see ?rasch_efrm for a complete simulated example
#' }
#' @export
plot_icc_frames <- function(fit, item, n_groups = fit$n_groups,
                            grid = seq(-5, 5, 0.05)) {
  if (!inherits(fit, "rasch_efrm")) stop("plot_icc_frames needs a rasch_efrm fit")
  vm <- fit$virtual_map
  rows <- which(vm$item == item)
  if (!length(rows)) stop("no such item: ", item)
  thr <- fit$thresholds_arbitrary
  tau_i <- thr$delta[thr$item == item]
  mmax <- length(tau_i)
  op <- .rr_canvas(range(grid), c(0, mmax), "Person location (logits, common unit)",
                   "Expected score", item)
  on.exit(par(op))
  th_all <- fit$person$theta
  for (j in seq_along(rows)) {
    v <- rows[j]
    rho <- fit$disc[v]
    colr <- .rr$pal[(j - 1L) %% length(.rr$pal) + 1L]
    Ecurve <- vapply(grid, function(t) item_moments(t, tau_i, disc = rho)$E, 0)
    lines(grid, Ecurve, lwd = 2.6, col = colr)
    x <- fit$X[, v]; ok <- !is.na(th_all) & !is.na(x)
    if (sum(ok) >= 2 * n_groups) {
      ci <- cut(rank(th_all[ok], ties.method = "first"),
                min(n_groups, max(2, floor(sum(ok) / 15))), labels = FALSE)
      points(tapply(th_all[ok], ci, mean), tapply(x[ok], ci, mean),
             pch = 21, bg = colr, col = "white", cex = 1.3, lwd = 1)
    }
  }
  .rr_legend("topleft",
             sprintf("%s (rho %.3f)", vm$group[rows], fit$disc[rows]),
             lwd = 2.6, col = .rr$pal[seq_along(rows)])
  invisible(NULL)
}
