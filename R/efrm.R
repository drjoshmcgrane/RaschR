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
# Route 2 (person-side linking): alpha_s and set locations mu_s from persons
# common to two sets. In set s a person's natural coordinate is
# u_s = alpha_s * (theta - mu_s), hence u_b = r*u_a + c for a linked set pair.
# The response likelihood is integrated over masses estimated jointly with
# the link on a grid in u_a. This semiparametric mixing distribution can
# represent skewed, heavy-tailed and multimodal person populations without a
# prespecified normal shape.
# The pairwise log ratios and offsets are reconciled over the set-linking graph
# by weighted least squares. Truncated-score moments provide stable starting
# values and screen links whose score range is too weak to identify a scale.
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
# estimate of E[target(u)] over a broad class of person distributions. These
# moments now provide starting values and identification checks for the
# semiparametric likelihood link rather than the final set-unit estimate.
# Returns
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
    # categories) to provide a stable scale start and guard check. With only
    # 2 interior categories the score-map slope is too weak to support the
    # flexible distributional link reliably.
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

# Semiparametric likelihood link for two item sets. The common persons'
# distribution is represented by jointly estimated masses on a fixed grid in
# the first set's natural unit; the second set is evaluated at r*u + c. This
# retains the Rasch response model without prescribing a normal or unimodal
# person distribution. Response rows are compressed by
# observed-item pattern and weighted score before the EM iterations, so the
# link remains practical inside the person bootstrap.
.efrm_npml_pair <- function(Xm, vmap, tau_v, disc_v, sets_u, a, b, idx,
                            init_log_ratio, init_offset,
                            min_link_persons, grid_n = 61L) {
  ca <- which(vmap$set == sets_u[a]); cb <- which(vmap$set == sets_u[b])
  Xa <- Xm[idx, ca, drop = FALSE]; Xb <- Xm[idx, cb, drop = FALSE]
  oa <- !is.na(Xa); ob <- !is.na(Xb)
  keep <- rowSums(oa) > 0L & rowSums(ob) > 0L
  Xa <- Xa[keep, , drop = FALSE]; Xb <- Xb[keep, , drop = FALSE]
  oa <- oa[keep, , drop = FALSE]; ob <- ob[keep, , drop = FALSE]
  if (nrow(Xa) < min_link_persons) return(NULL)

  da <- disc_v[ca]; db <- disc_v[cb]
  score_a <- rowSums(sweep(Xa, 2L, da, "*"), na.rm = TRUE)
  score_b <- rowSums(sweep(Xb, 2L, db, "*"), na.rm = TRUE)
  pa <- apply(oa, 1L, paste0, collapse = "")
  pb <- apply(ob, 1L, paste0, collapse = "")
  key <- paste(pa, formatC(score_a, digits = 12L, format = "fg"),
               pb, formatC(score_b, digits = 12L, format = "fg"), sep = "\r")
  lev <- unique(key); grp <- match(key, lev)
  reps <- match(lev, key); count <- tabulate(grp, nbins = length(lev))
  oa <- oa[reps, , drop = FALSE]; ob <- ob[reps, , drop = FALSE]
  score_a <- score_a[reps]; score_b <- score_b[reps]
  # Group membership is observed and the person populations need not have the
  # same location, spread or shape. Keep one nonparametric margin per observed
  # group while estimating a common set transformation. A single pooled margin
  # would impose an avoidable population-distribution restriction whenever the
  # fitted group units differ.
  mix_group <- vapply(seq_len(nrow(oa)), function(i) {
    ja <- which(oa[i, ])
    jb <- which(ob[i, ])
    ga <- if (length(ja)) unique(vmap$group[ca[ja]]) else character(0)
    gb <- if (length(jb)) unique(vmap$group[cb[jb]]) else character(0)
    gg <- unique(c(ga, gb))
    if (length(gg) != 1L) stop("internal EFRM link mixes person groups")
    gg
  }, character(1L))
  mix_levels <- sort(unique(mix_group))
  mix_idx <- match(mix_group, mix_levels)

  # Log likelihood up to response-pattern constants, which cancel from the
  # posterior masses and do not depend on r or c.
  likelihood <- function(u, obs, score, taus, discs) {
    L <- outer(score, u)
    patterns <- apply(obs, 1L, paste0, collapse = "")
    for (pattern in unique(patterns)) {
      rows <- which(patterns == pattern)
      jj <- which(obs[rows[1L], ])
      den <- numeric(length(u))
      for (j in jj) {
        tt <- taus[[j]]; x <- 0:length(tt)
        lp <- discs[j] * (outer(u, x) -
          matrix(c(0, cumsum(tt)), nrow = length(u),
                 ncol = length(x), byrow = TRUE))
        mx <- apply(lp, 1L, max)
        den <- den + mx + log(rowSums(exp(lp - mx)))
      }
      L[rows, ] <- L[rows, , drop = FALSE] -
        matrix(den, nrow = length(rows), ncol = length(u), byrow = TRUE)
    }
    L
  }

  grid <- seq(-8, 8, length.out = as.integer(grid_n))
  La <- likelihood(grid, oa, score_a, tau_v[ca], da)
  par <- c(init_log_ratio, init_offset)
  if (any(!is.finite(par))) par <- c(0, 0)
  # Joint NPML by coordinate ascent. For a fixed link, the EM step estimates
  # the grid masses from both sets; for fixed masses, a bounded likelihood
  # step estimates the log scale ratio and offset. Using both sets in the
  # mixing-distribution step matters: estimating the masses from the first
  # set alone would turn the second into a plug-in prediction problem rather
  # than maximise the likelihood of the linked responses.
  add_weights <- function(L, logw)
    L + logw[mix_idx, , drop = FALSE]
  fit_weights <- function(L, logw, maxit = 100L) {
    conv <- FALSE; step <- Inf
    for (it in seq_len(maxit)) {
      A <- add_weights(L, logw); mx <- apply(A, 1L, max)
      post <- exp(A - mx); post <- post / rowSums(post)
      w <- matrix(NA_real_, nrow(logw), ncol(logw))
      for (h in seq_len(nrow(logw))) {
        sel <- mix_idx == h
        wh <- pmax(colSums(post[sel, , drop = FALSE] * count[sel]) /
                     sum(count[sel]), 1e-10)
        w[h, ] <- wh / sum(wh)
      }
      next_logw <- log(w)
      # Log changes are unstable for support points whose masses approach
      # zero. Convergence concerns the mixing distribution itself.
      step <- max(abs(w - exp(logw)))
      if (step < 1e-7) {
        logw <- next_logw; conv <- TRUE; break
      }
      logw <- next_logw
    }
    list(logw = logw, converged = conv, step = step)
  }
  logw <- matrix(-log(length(grid)), length(mix_levels), length(grid),
                 dimnames = list(mix_levels, NULL))
  converged <- FALSE
  op <- NULL
  last_step <- Inf
  for (outer in seq_len(20L)) {
    old <- par
    Lb <- likelihood(exp(par[1L]) * grid + par[2L], ob, score_b,
                     tau_v[cb], db)
    ew <- fit_weights(La + Lb, logw)
    logw <- ew$logw
    objective <- function(z) {
      Lbz <- likelihood(exp(z[1L]) * grid + z[2L], ob, score_b,
                        tau_v[cb], db)
      A <- add_weights(La + Lbz, logw); mx <- apply(A, 1L, max)
      -sum(count * (mx + log(rowSums(exp(A - mx)))))
    }
    op <- stats::optim(par, objective, method = "L-BFGS-B",
                       lower = c(log(0.1), -20), upper = c(log(10), 20),
                       control = list(maxit = 30L, factr = 1e7))
    par <- op$par
    last_step <- max(abs(par - old))
    if (op$convergence == 0L && last_step < 1e-3) {
      converged <- TRUE; break
    }
  }
  # Refresh the nuisance distribution at the reported link and evaluate the
  # corresponding joint log likelihood.
  Lb <- likelihood(exp(par[1L]) * grid + par[2L], ob, score_b,
                   tau_v[cb], db)
  ew <- fit_weights(La + Lb, logw)
  logw <- ew$logw; w <- exp(logw)
  A <- add_weights(La + Lb, logw); mx <- apply(A, 1L, max)
  ll <- sum(count * (mx + log(rowSums(exp(A - mx)))))
  converged <- converged || (!is.null(op) && op$convergence == 0L &&
                              last_step < 1e-3)
  edge_mass_by_group <- w[, 1L] + w[, ncol(w)]
  edge_mass <- max(edge_mass_by_group)
  if (isTRUE(getOption("rasch.efrm_link_debug", FALSE)))
    message("EFRM NPML: step=", signif(last_step, 4),
            ", mass step=", signif(ew$step, 4),
            ", optim=", op$convergence,
            ", max group edge mass=", signif(edge_mass, 4))
  if (!all(is.finite(par)) || !is.finite(ll) || edge_mass > 0.02) return(NULL)
  list(log_ratio = par[1L], offset = par[2L], n = sum(count),
       converged = converged, edge_mass = edge_mass, loglik = ll)
}

# Stage 2: alpha_s and set locations mu_s from persons common to set pairs,
# with semiparametric person-distribution linking. Score moments provide
# stable starting values and guard checks; the pair likelihood supplies the
# reported scale ratio and offset.
# Standard errors and the joint covariance of (log alpha, mu) come from a
# person bootstrap of the whole linking stage: each replicate resamples
# persons, refits the nuisance masses and link for each set pair, and
# re-solves the linking least squares. This captures the nonlinearity of the
# semiparametric link and its reconciliation over the linking graph.
.efrm_link_sets <- function(u_mat, w_mat, g_mat, sets_u, min_link_persons,
                            boot_reps = 300, regen = NULL, pair_link = NULL) {
  S <- ncol(u_mat)
  A <- rbind(diag(S - 1L), rep(-1, S - 1L))

  # one pass over a person index vector: per-edge corrected slopes and
  # offsets, then the two weighted least-squares solves. um/wm/gm default
  # to the point-estimate person matrices; the bootstrap can pass
  # regenerated ones (see below).
  link_once <- function(idx, um = u_mat, wm = w_mat, gm = g_mat,
                        hard = FALSE, pair_fun = pair_link) {
    edges <- list(); ls_est <- off_est <- off_n <- numeric(0)
    link_converged <- link_edge_mass <- link_loglik <- numeric(0)
    for (a in seq_len(S - 1)) for (b in (a + 1):S) {
      # Non-extreme scores are needed only for the corrected-moment start and
      # the weak-link screen. The likelihood itself must retain every common
      # person, including an extreme score in either set: dropping them would
      # condition on the realised responses without including that selection
      # in the likelihood.
      ok_start <- idx[is.finite(um[idx, a]) & is.finite(um[idx, b]) &
                        is.finite(wm[idx, a]) & is.finite(wm[idx, b]) &
                        is.finite(gm[idx, a]) & is.finite(gm[idx, b])]
      if (length(ok_start) < min_link_persons) next
      u1 <- um[ok_start, a]; u2 <- um[ok_start, b]
      d1 <- mean(gm[ok_start, a]); d2 <- mean(gm[ok_start, b])
      if (!is.finite(d1) || !is.finite(d2) || d1 < 0.1 || d2 < 0.1) {
        if (hard)
          stop("the score maps of sets '", sets_u[a], "' and '", sets_u[b],
               "' are too flat to link (degenerate map slope)")
        next
      }
      v1 <- (var(u1) - mean(wm[ok_start, a])) / d1^2
      v2 <- (var(u2) - mean(wm[ok_start, b])) / d2^2
      if (!is.finite(v1) || !is.finite(v2) || v1 <= 0 || v2 <= 0) {
        if (hard)
          stop("too little true person variance to link sets '", sets_u[a],
               "' and '", sets_u[b], "'")
        next
      }
      ls <- 0.5 * (log(v2) - log(v1))                # log(alpha_b / alpha_a)
      off <- mean(u2) - exp(ls) * mean(u1)
      np <- NULL
      if (!is.null(pair_fun)) {
        np <- pair_fun(a, b, idx, ls, off)
        if (is.null(np)) {
          if (hard)
            stop("the semiparametric likelihood link failed for sets '",
                 sets_u[a], "' and '", sets_u[b], "'")
          next
        }
        # Retain a usable point estimate with an explicit warning, but do not
        # admit a link that missed its numerical tolerance to a bootstrap
        # covariance calculation.
        if (!hard && !isTRUE(np$converged)) next
        ls <- np$log_ratio; off <- np$offset
      }
      edges[[length(edges) + 1L]] <- c(a, b)
      ls_est <- c(ls_est, ls)
      off_est <- c(off_est, off)  # alpha_b (mu_a - mu_b)
      off_n <- c(off_n, if (is.null(pair_fun)) length(ok_start) else np$n)
      link_converged <- c(link_converged,
                          if (is.null(np)) NA else np$converged)
      link_edge_mass <- c(link_edge_mass,
                          if (is.null(np)) NA else np$edge_mass)
      link_loglik <- c(link_loglik, if (is.null(np)) NA else np$loglik)
    }
    if (!length(edges)) {
      if (hard) stop("no set pairs share enough persons with informative ",
                     "score patterns to link the units: each person's ",
                     "pattern needs a score range of at least 4 within a ",
                     "set (at least four dichotomous items, or fewer ",
                 "polytomous ones) for semiparametric set linking")
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
    list(la = la, mu = mu, edges = edges, ls_est = ls_est, off_n = off_n,
         link_converged = link_converged, link_edge_mass = link_edge_mass,
         link_loglik = link_loglik)
  }

  N <- nrow(u_mat)
  point <- link_once(seq_len(N), hard = TRUE)
  point_converged <- !any(point$link_converged %in% FALSE)

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
  cov_alpha_phi <- NULL
  link_reps <- dtilde_reps <- NULL
  if (boot_reps > 0 && point_converged) {
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
    phi_reps <- NULL
    for (r in seq_len(boot_reps)) {
      mats <- if (n_draws == 0L) NULL else {
        di <- ((r - 1L) %% n_draws) + 1L
        if (is.null(draws[[di]])) draws[[di]] <- regen()
        draws[[di]]
      }
      b <- if (is.null(mats)) link_once(sample.int(N, N, replace = TRUE))
           else link_once(sample.int(N, N, replace = TRUE),
                          um = mats$u, wm = mats$w, gm = mats$g,
                          pair_fun = mats$pair_link %||% pair_link)
      if (!is.null(b)) reps[r, ] <- c(b$la, b$mu)
      if (!is.null(mats$log_phi)) {
        if (is.null(phi_reps))
          phi_reps <- matrix(NA_real_, boot_reps, length(mats$log_phi),
                             dimnames = list(NULL, names(mats$log_phi)))
        phi_reps[r, ] <- mats$log_phi
      }
      if (!is.null(mats$dtilde)) {
        if (is.null(dtilde_reps))
          dtilde_reps <- matrix(NA_real_, boot_reps, length(mats$dtilde))
        dtilde_reps[r, ] <- mats$dtilde
      }
    }
    complete <- stats::complete.cases(reps)
    reps_ok <- reps[complete, , drop = FALSE]
    # scale the requirement to what was asked for: a flat "< 30" rejects every
    # boot_reps below 30 without a single replicate having failed, and then
    # blames the design for it
    if (nrow(reps_ok) < max(2L, min(30L, boot_reps %/% 2L)))
      stop("the unit-linking bootstrap failed in most replicates; the linking ",
           "design is too weak for stable alpha estimation")
    cov_link <- stats::cov(reps_ok)
    link_reps <- reps_ok
    if (!is.null(dtilde_reps))
      dtilde_reps <- dtilde_reps[complete, , drop = FALSE]
    if (!is.null(phi_reps)) {
      joint_ok <- complete & stats::complete.cases(phi_reps)
      if (sum(joint_ok) >= 2L) {
        C <- stats::cov(cbind(reps[joint_ok, seq_len(S), drop = FALSE],
                              phi_reps[joint_ok, , drop = FALSE]))
        cov_alpha_phi <- C[seq_len(S), S + seq_len(ncol(phi_reps)),
                           drop = FALSE]
        dimnames(cov_alpha_phi) <- list(sets_u, colnames(phi_reps))
      }
    }
  }

  list(alpha = setNames(exp(point$la), sets_u),
       se_log_alpha = if (is.null(cov_link)) setNames(rep(NA_real_, S), sets_u)
         else setNames(sqrt(pmax(diag(cov_link)[seq_len(S)], 0)), sets_u),
       mu = setNames(point$mu, sets_u),
       cov_link = cov_link,
       cov_alpha_phi = cov_alpha_phi,
       link_reps = link_reps,
       dtilde_reps = dtilde_reps,
       edges = data.frame(set_a = sets_u[vapply(point$edges, `[`, 1L, 1)],
                          set_b = sets_u[vapply(point$edges, `[`, 1L, 2)],
                          n = point$off_n, log_slope = point$ls_est,
                          converged = point$link_converged,
                          edge_mass = point$link_edge_mass,
                          loglik = point$link_loglik))
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
#' Person-group units \eqn{\phi_g} are identified from common item thresholds
#' across groups. Item sets partition the items, so set units \eqn{\alpha_s}
#' are identified instead from persons observed in more than one set. The
#' set-linking graph and the group-by-set frame graph must each connect to a
#' common scale.
#'
#' Set units use a semiparametric likelihood for persons observed in each
#' linked pair of sets. For sets \eqn{a} and \eqn{b}, it maximises
#' \deqn{\prod_n\int P(X_{na}\mid u)P(X_{nb}\mid ru+c)\,dF_{g(n)}(u),}
#' where the masses of each observed group's \eqn{F_g}, the scale ratio
#' \eqn{r} and the offset \eqn{c} are estimated jointly on a fixed grid. This
#' avoids prescribing a normal or common person distribution across groups.
#' The conditional thresholds and group units are held fixed in this step;
#' only \eqn{r}, \eqn{c}, and the nuisance masses are estimated. The linked
#' parameters are then
#' \deqn{\delta_{ik}=\widetilde\delta_{ik}/\alpha_s+\mu_s,
#' \qquad \rho_{sg}=\alpha_s\phi_g.}
#' Score moments supply starting values and screen weak links. Response
#' patterns must span a score range of at least four within a set. Overlapping
#' item sets are not permitted. The public
#' convergence flag covers both estimation stages; \code{stage1_converged}
#' records the conditional stage separately.
#'
#' The hybrid covariance combines the pairwise Godambe covariance with a
#' person bootstrap for set linking. Each replicate jointly redraws the
#' within-frame thresholds and group units, then rebuilds the link. The joint
#' draws retain covariance among common-scale thresholds, set units and group
#' units. With \code{se_method = "bootstrap"}, the complete model is refitted
#' to each person resample.
#'
#' The \code{efrm_vs_rasch} component records the within-frame composite
#' log-likelihood comparison between group-dependent and equal group units.
#' This difference is descriptive and contains no information about set units,
#' which are identified at the linking stage. The accompanying Wald omnibus
#' tests provide inference for the group- and set-unit families.
#'
#' The model assumes that an item retains its location and discrimination
#' across the frames in which it appears, apart from the frame unit.
#' \code{\link{frame_invariance}} examines this assumption by separate frame
#' calibrations. Misfit concentrated within one item set can also distort its
#' estimated unit; inspect item fit and targeting before interpreting unit
#' differences. \code{\link{drop_items}} and \code{\link{resolve_frames}}
#' provide refitted sensitivity analyses.
#'
#' The dichotomous model follows Humphry (2005) and Humphry and Andrich
#' (2008). The polytomous, multigroup and crossed-frame forms are extensions
#' implemented in this package. The discrete nonparametric margin follows the
#' Rasch estimation approach of Follmann (1988); its use for linked item-set
#' units is an extension implemented here.
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
#' @param id Person identifier, either a column name or one value per row.
#'   EFRM data require one response row per person, so identifiers must be
#'   unique.
#' @param factors,items,n_groups,adjust_N,na_codes As in \code{\link{rasch}}.
#' @param maxit,tol Outer iteration cap and convergence tolerance of the
#'   bilinear pairwise stage.
#' @param min_link_persons Minimum number of common persons required for a
#'   set pair to contribute to the unit linking.
#' @param se_method \code{"hybrid"} (sandwich + linking bootstrap + delta
#'   propagation; fast, default) or \code{"bootstrap"} (full person
#'   bootstrap of all stages).
#' @param boot_reps Bootstrap replicates; defaults to 300 for the linking
#'   bootstrap and 200 for the full bootstrap. Use zero to omit unit
#'   uncertainty; otherwise at least 30 are required.
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
#' Follmann, D. (1988). Consistent estimation in the Rasch model based on
#' nonparametric margins. Psychometrika, 53, 553--562.
#' \doi{10.1007/BF02294407}
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
  .check_column_names(data)
  n_groups_requested <- n_groups
  se_method <- match.arg(se_method)
  if (is.null(boot_reps)) boot_reps <- if (se_method == "hybrid") 300L else 200L
  if (length(boot_reps) != 1L || !is.finite(boot_reps) || boot_reps < 0L ||
      boot_reps != floor(boot_reps))
    stop("boot_reps must be one non-negative whole number")
  boot_reps <- as.integer(boot_reps)
  if (boot_reps > 0L && boot_reps < 30L)
    stop("EFRM uncertainty needs either zero or at least 30 bootstrap replicates")
  # Simulation-only grid override used to verify that reported links are not
  # artefacts of the 61-point default. It is deliberately not a public model
  # option: changing the grid is a validation exercise, not an analyst choice.
  link_grid_n <- getOption("rasch.efrm_link_grid_n", 61L)
  if (length(link_grid_n) != 1L || !is.finite(link_grid_n) ||
      link_grid_n != floor(link_grid_n) || link_grid_n < 21L)
    stop("option rasch.efrm_link_grid_n must be one whole number of at least 21")
  link_grid_n <- as.integer(link_grid_n)
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
        grp <- .factor_cells(grp_components, sep = ":")
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
    # Frame-defining columns are already stored below as part of the frame
    # structure. Repeating them as ordinary factors creates duplicate names
    # and can make a later DIF or refit select the wrong column.
    if (!is.null(fac_df) && is.character(groups))
      fac_df <- fac_df[, !names(fac_df) %in% groups, drop = FALSE]
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
  present_id <- !is.na(id_vec)
  if (anyDuplicated(as.character(id_vec[present_id])))
    stop("rasch_efrm needs one response row per person; duplicate identifiers ",
         "would be treated as independent people by the set-link likelihood")
  # The crossed-cell column is internal metadata. Keep its readable name
  # unless it would collide with a component or an ordinary person factor.
  taken_factor_names <- unique(c(names(grp_components), names(fac_df)))
  if (grp_name %in% taken_factor_names) {
    grp_name <- ".frame_group"
    while (grp_name %in% taken_factor_names)
      grp_name <- paste0(grp_name, ".")
  }

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
  vmap$vkey <- as.character(.factor_cells(vmap[c("item", "group")],
                                           sep = ":"))
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
  if (!isTRUE(sol$converged))
    warning("EFRM estimation did NOT converge in ", sol$iterations,
            " iterations; increase maxit or inspect the frame design",
            call. = FALSE)
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
  make_pair_link <- function(Xm, dt, phiv) {
    tau_v_link <- lapply(seq_len(ncol(Xm)), function(v)
      dt[drow[which(thr_v$item == v)]])
    disc_v_link <- unname(phiv[vmap$group])
    function(a, b, idx, init_ls, init_off)
      .efrm_npml_pair(Xm, vmap, tau_v_link, disc_v_link, sets_u,
                      a, b, idx, init_ls, init_off, min_link_persons,
                      grid_n = link_grid_n)
  }
  if (S > 1L) {
    pm <- person_mats(Xv, dtil, phi)
    pair_link <- make_pair_link(Xv, dtil, phi)
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
      phi_draw <- phi * exp(v[K_dt + seq_along(phi)])
      dt_draw <- dtil + v[seq_len(K_dt)]
      c(person_mats(Xv, dt_draw, phi_draw),
        list(log_phi = log(phi_draw),
             dtilde = dt_draw,
             pair_link = make_pair_link(Xv, dt_draw, phi_draw)))
    } else NULL
    link <- .efrm_link_sets(pm$u, pm$w, pm$g, sets_u, min_link_persons,
                            boot_reps = boot_reps, regen = regen,
                            pair_link = pair_link)
    alpha <- link$alpha; mu <- link$mu
    if (any(link$edges$converged %in% FALSE)) {
      warning("one or more semiparametric set links stopped before the scale ",
              "step met its convergence tolerance; inspect fit$linking$alpha_edges",
              call. = FALSE)
      notes <- c(notes, paste(
        "one or more semiparametric set links stopped before the scale step",
        "met its convergence tolerance"))
    }
  } else {
    alpha <- setNames(1, sets_u); mu <- setNames(0, sets_u)
    link <- list(alpha = alpha, se_log_alpha = setNames(0, sets_u), mu = mu,
                 cov_link = NULL, cov_alpha_phi = NULL, edges = data.frame())
  }

  # --- optional full person bootstrap of all stages ----------------------------
  boot <- NULL
  if (se_method == "bootstrap" &&
      !any(link$edges$converged %in% FALSE)) {
    Npers <- nrow(Xv)
    boot_replicate <- function(idx) {
      Xb <- Xv[idx, , drop = FALSE]
      pb <- .efrm_filter_pairs(.pair_counts(Xb, m_v), vmap)
      if (!length(pb)) return(NULL)
      sb <- .efrm_solve(Xb, thr_v, m_v, vmap, pb, drow, A_D,
                        maxit = maxit, tol = tol)
      if (!isTRUE(sb$converged)) return(NULL)
      if (S > 1L) {
        pm_b <- person_mats(Xb, sb$dtilde, sb$phi)
        lb <- .efrm_link_sets(pm_b$u, pm_b$w, pm_b$g, sets_u,
                              min_link_persons, boot_reps = 0,
                              pair_link = make_pair_link(Xb, sb$dtilde, sb$phi))
        if (any(lb$edges$converged %in% FALSE)) return(NULL)
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

  # Covariance of delta = dtilde/alpha + mu. The full bootstrap supplies
  # delta directly. In the hybrid bootstrap, each calibration draw is used to
  # rebuild the set link; retaining that draw beside its linked alpha and mu
  # preserves their covariance. Adding the two variance components as if
  # independent would discard precisely this shared-calibration term.
  cov_delta <- sol$cov_dtilde / tcrossprod(alpha[set_of_drow])
  if (!is.null(boot)) {
    cov_delta <- stats::cov(boot[, G + 2L * S + seq_len(Md), drop = FALSE])
  } else if (S > 1L && !is.null(link$link_reps) &&
             !is.null(link$dtilde_reps)) {
    lr <- link$link_reps
    dr <- link$dtilde_reps
    delta_reps <- matrix(NA_real_, nrow(lr), Md)
    for (r in seq_len(nrow(lr))) {
      ar <- exp(lr[r, seq_len(S)])
      mr <- lr[r, S + seq_len(S)]
      delta_reps[r, ] <- dr[r, ] / ar[match(set_of_drow, sets_u)] +
        mr[match(set_of_drow, sets_u)]
    }
    cov_delta <- stats::cov(delta_reps)
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
  link_converged <- S == 1L || !any(link$edges$converged %in% FALSE)
  est <- list(model = "EFRM", thr = thr_v, cov_tau = cov_tau,
              loglik = sol$loglik, iterations = sol$iterations,
              converged = isTRUE(sol$converged) && link_converged,
              stage1_converged = sol$converged,
              m = m_v, anchors = NULL,
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
  # Several frame-specific copies of one item are calibration cells, not
  # additional administered items. Alpha and one raw-score conversion over
  # the expanded columns are therefore undefined. Retain both for the
  # one-frame reduction, where the fitted columns are the items and
  # the model is the ordinary Rasch model.
  expanded_cells <- any(table(vmap$item) > 1L)
  if (expanded_cells) {
    fit$alpha <- list(
      alpha = NA_real_, n = NA_integer_, applicable = FALSE,
      design_applicable = FALSE,
      reason = "not applicable when an item has several frame response cells")
    fit$score_table <- NULL
    fit$notes <- unique(c(fit$notes, paste(
      "a universal raw-score conversion is not defined across the expanded",
      "frame response cells; use score_curves and design-specific information")))
  } else fit$alpha$design_applicable <- TRUE
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
    # Recover component values from the recorded design rather than parsing
    # the composite label. Component levels may themselves contain colons.
    first <- match(glevs, as.character(grp))
    dd <- as.data.frame(lapply(grp_components, function(x)
      as.character(x[first])), stringsAsFactors = FALSE)
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
               # cross-stage dependence is captured. The hybrid calculation
               # uses the matching covariance from its joint calibration
               # redraws rather than adding the two marginal variances.
               se_log_rho = if (!is.null(boot))
                 stats::sd(boot[, match(fr$group[j], glevs)] +
                           boot[, G + match(fr$set[j], sets_u)])
               else {
                 is <- match(fr$set[j], sets_u)
                 ig <- match(fr$group[j], glevs)
                 cross <- if (is.null(link$cov_alpha_phi)) 0 else
                   link$cov_alpha_phi[is, ig]
                 sqrt(pmax(fit$alpha_table$se_log_alpha[is]^2 +
                             fit$phi_table$se_log_phi[ig]^2 + 2 * cross, 0))
               },
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
                       cov_log_alpha = Sig_alpha,
                       cov_joint = sol$cov_joint,
                       cov_log_alpha_phi = link$cov_alpha_phi)
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
  fit$linking <- list(
    phi_edges = edges_g,
    alpha_edges = link$edges,
    alpha_method = if (S > 1L)
      "finite-grid semiparametric maximum likelihood" else "not applicable",
    alpha_grid = if (S > 1L)
      c(lower = -8, upper = 8, points = link_grid_n) else NULL)
  fit$se_method <- if (!is.null(boot)) "bootstrap" else "hybrid"
  fit$boot_reps_used <- if (!is.null(boot)) nrow(boot) else NA_integer_
  fit$virtual_map <- vmap
  fit$set_of <- set_of
  fit$refit_spec <- list(
    groups = if (!is.null(grp_components)) names(grp_components) else grp_name,
    factors = setdiff(names(fit$factors), fit$frame_group),
    n_groups = n_groups_requested, adjust_N = adjust_N, na_codes = na_codes,
    maxit = maxit, tol = tol, min_link_persons = min_link_persons,
    se_method = se_method, boot_reps = boot_reps)
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
#' Plots the model expected-score curve for one item in each frame,
#' with observed class-interval means overlaid. Differences between the model
#' curves reflect the fitted frame units. A nominated non-frame person factor
#' separates the observed means within each frame for a DIF display.
#'
#' @param fit A fitted object from \code{\link{rasch_efrm}}.
#' @param item Underlying item name.
#' @param n_groups Number of class intervals for the observed means.
#' @param grid Logit grid.
#' @param group Optional person grouping vector, or one or more names of
#'   non-frame factors nominated in the fit. Several names define their
#'   factor-combination cells.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' \donttest{
#' # see ?rasch_efrm for a complete simulated example
#' }
#' @export
plot_icc_frames <- function(fit, item, n_groups = fit$n_groups,
                            grid = seq(-5, 5, 0.05), group = NULL) {
  if (!inherits(fit, "rasch_efrm")) stop("plot_icc_frames needs a rasch_efrm fit")
  vm <- fit$virtual_map
  rows <- which(vm$item == item)
  if (!length(rows)) stop("no such item: ", item)
  group_label <- NULL
  if (is.character(group) && length(group) < nrow(fit$X)) {
    bad <- intersect(group, fit$frame_group %||% character(0))
    if (length(bad))
      stop("frame-defining factors belong to the frame curves, not the DIF overlay: ",
           paste(bad, collapse = ", "))
    if (is.null(fit$factors) || !all(group %in% names(fit$factors)))
      stop("every named group must be a person factor in the fit")
    group_label <- paste(group, collapse = " x ")
    group <- if (length(group) == 1L) fit$factors[[group]] else
      .factor_cells(fit$factors[group], sep = ":")
  }
  if (!is.null(group) && length(group) != nrow(fit$X))
    stop("group must have one value per person")
  gf <- if (is.null(group)) NULL else droplevels(factor(group))
  if (!is.null(gf) && nlevels(gf) < 2L)
    stop("group must contain at least two observed levels")
  thr <- fit$thresholds_arbitrary
  tau_i <- thr$delta[thr$item == item]
  mmax <- length(tau_i)
  op <- .rr_canvas(range(grid), c(0, mmax),
                   "Person location (logits, common unit)", "Expected score",
                   if (is.null(group_label)) item else
                     sprintf("%s by %s", item, group_label))
  on.exit(par(op))
  th_all <- fit$person$theta
  frame_cols <- .rr$pal[(seq_along(rows) - 1L) %% length(.rr$pal) + 1L]
  group_pch <- if (is.null(gf)) 21 else
    rep(c(21:25, 0:20), length.out = nlevels(gf))
  for (j in seq_along(rows)) {
    v <- rows[j]
    rho <- fit$disc[v]
    colr <- frame_cols[j]
    Ecurve <- vapply(grid, function(t) item_moments(t, tau_i, disc = rho)$E, 0)
    lines(grid, Ecurve, lwd = 2.6, col = colr)
    x <- fit$X[, v]
    candidate <- !is.na(th_all) & !is.na(x) & !fit$person$extreme
    if (sum(candidate) >= 2 * n_groups) {
      ng <- min(n_groups, max(2, floor(sum(candidate) / 15)))
      ci <- .class_intervals(ifelse(is.na(x), NA_real_, th_all),
                             fit$person$extreme, ng)
      ok <- !is.na(ci)
      if (is.null(gf)) {
        points(tapply(th_all[ok], ci[ok], mean), tapply(x[ok], ci[ok], mean),
               pch = 21, bg = colr, col = "white", cex = 1.3, lwd = 1)
      } else {
        g <- gf[ok]
        for (k in seq_len(nlevels(gf))) {
          take <- g == levels(gf)[k]
          if (sum(take) < 2L) next
          points(tapply(th_all[ok][take], ci[ok][take], mean),
                 tapply(x[ok][take], ci[ok][take], mean),
                 pch = group_pch[k], bg = colr, col = colr,
                 cex = 1.15, lwd = 1)
        }
      }
    }
  }
  .rr_legend("topleft",
             sprintf("%s (rho %.3f)", vm$group[rows], fit$disc[rows]),
             lwd = 2.6, col = frame_cols)
  if (is.null(gf)) {
    .rr_legend("bottomright", c("Model", "Observed"),
               lwd = c(2.6, NA), pch = c(NA, 21),
               pt.bg = c(NA, .rr$ink), col = c(.rr$ink, "white"),
               pt.cex = 1.15)
  } else {
    .rr_legend("bottomright", levels(gf), pch = group_pch,
               pt.bg = "white", col = .rr$ink, pt.cex = 1.05)
  }
  invisible(NULL)
}
