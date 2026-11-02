#### Risk set regression calibration (paper Section 3.2) ####

# group events into validation risk sets: working back from the last event,
# merge risk sets until each holds more than min_pairs pairs (paper Section 5.1)
rrc_groups <- function(ev_time, pair_T, min_pairs) {
  np <- length(pair_T) - findInterval(ev_time, sort(pair_T), left.open = TRUE)
  E <- length(ev_time)
  grp <- integer(E); g <- 0L; open <- FALSE
  for (e in E:1) {
    if (!open) { g <- g + 1L; open <- TRUE }
    grp[e] <- g
    if (np[e] > min_pairs) open <- FALSE
  }
  max(grp) + 1L - grp
}

# Fit the error model in each validation risk set, build risk-set-specific x_hat,
# fit the survival model on x_hat and compute the sandwich variance (Supp. 7.1.3).
#   pairs: data.frame(id, l) of validated (subject, year) pairs
#   C, c:  n x T annual surrogate / true exposure; L: cum_lag_matrix(T)
rrc_fit <- function(pairs, C, c, time, ev, error, model, L, min_pairs, start = 0) {
  n <- nrow(C)
  pc <- c[cbind(pairs$id, pairs$l)]
  pC <- C[cbind(pairs$id, pairs$l)]
  if (error == "add") { py <- pc; px <- pC } else { py <- log(pc); px <- log(pC) }
  pT <- time[pairs$id]

  V <- matrix(FALSE, n, ncol(C)); V[cbind(pairs$id, pairs$l)] <- TRUE
  NV <- !V
  CV <- (c * V) %*% L                       # observed true exposure in validated years
  if (error == "add") {
    Nnv <- NV %*% L                         # number of non-validated years
    CSnv <- (C * NV) %*% L                  # surrogate sum over non-validated years
  }

  grp <- rrc_groups(ev$time, pT, min_pairs)
  G <- max(grp)
  p <- if (error == "add") 2 else 3
  zz <- numeric(length(ev$subj)); dz <- matrix(0, length(ev$subj), p)
  zc <- numeric(ev$n); dzc <- matrix(0, ev$n, p)
  fits <- vector("list", G)

  for (g in 1:G) {
    eg <- which(grp == g)
    S <- which(pT >= ev$time[eg[1]])        # nested risk sets: union = earliest one
    X <- cbind(1, px[S])
    XtX <- crossprod(X)
    a <- solve(XtX, crossprod(X, py[S]))
    res <- as.vector(py[S] - X %*% a)
    s2 <- mean(res^2)
    fits[[g]] <- list(S = S, X = X, res = res, s2 = s2,
                      B = if (error == "add") XtX else rbind(cbind(XtX, 0), c(0, 0, length(S))))

    il <- which(ev$eid %in% eg)
    sl <- ev$subj[il]; kl <- ev$k[ev$eid[il]]
    sc <- ev$id[eg]; kc <- ev$k[eg]
    if (error == "add") {
      zz[il] <- CV[cbind(sl, kl + 1)] + a[1] * Nnv[cbind(sl, kl + 1)] + a[2] * CSnv[cbind(sl, kl + 1)]
      dz[il, ] <- cbind(Nnv[cbind(sl, kl + 1)], CSnv[cbind(sl, kl + 1)])
      zc[eg] <- CV[cbind(sc, kc + 1)] + a[1] * Nnv[cbind(sc, kc + 1)] + a[2] * CSnv[cbind(sc, kc + 1)]
      dzc[eg, ] <- cbind(Nnv[cbind(sc, kc + 1)], CSnv[cbind(sc, kc + 1)])
    } else {
      # E(c | C) = exp(a0 + s2 / 2) * C^a1
      K <- exp(a[1] + s2 / 2)
      Cp <- C^a[2] * NV
      ks <- sort(unique(kc))
      P <- Cp %*% L[, ks + 1, drop = FALSE]
      Q <- (Cp * log(C)) %*% L[, ks + 1, drop = FALSE]
      cl <- match(kl, ks); cc <- match(kc, ks)
      Pl <- P[cbind(sl, cl)]; Pc <- P[cbind(sc, cc)]
      zz[il] <- CV[cbind(sl, kl + 1)] + K * Pl
      dz[il, ] <- cbind(K * Pl, K * Q[cbind(sl, cl)], K * Pl / 2)
      zc[eg] <- CV[cbind(sc, kc + 1)] + K * Pc
      dzc[eg, ] <- cbind(K * Pc, K * Q[cbind(sc, cc)], K * Pc / 2)
    }
  }

  fit <- pl_fit(zz, ev$eid, zc, model, start = start)
  b <- fit$beta; pt <- fit$parts

  # robust main-study score U~_i (eqs 7.37 / 7.38)
  if (model == "cox") { g_l <- zz; gc <- zc } else { g_l <- zz / pt$r; gc <- zc / (1 + b * zc) }
  term <- pt$r / pt$S0[ev$eid] * (g_l - pt$gbar[ev$eid])
  Ut <- numeric(n)
  Ut[ev$id] <- gc - pt$gbar
  rs <- rowsum(term, ev$subj)
  ri <- as.integer(rownames(rs))
  Ut[ri] <- Ut[ri] - rs[, 1]

  # U* = dU/d alpha for each group's parameters (eq 7.34)
  if (model == "cox") {
    bl <- pt$r * (1 + b * (zz - pt$gbar[ev$eid])) / pt$S0[ev$eid]
    contrib <- dzc - rowsum(bl * dz, ev$eid, reorder = FALSE)
  } else {
    be <- (pt$S0 - b * pt$Sz) / pt$S0^2
    contrib <- dzc / (1 + b * zc)^2 - be * rowsum(dz, ev$eid, reorder = FALSE)
  }
  Ustar <- rowsum(contrib, grp, reorder = FALSE)

  # validation contribution U* V_alpha U*' via subject-level scores (eq 7.35)
  s <- numeric(n)
  for (g in 1:G) {
    f <- fits[[g]]
    ag <- solve(t(f$B), Ustar[g, ])
    sc_pair <- f$X * f$res
    if (error == "mult") sc_pair <- cbind(sc_pair, f$res^2 - f$s2)
    v <- rowsum(as.vector(sc_pair %*% ag), pairs$id[f$S])
    vi <- as.integer(rownames(v))
    s[vi] <- s[vi] + v[, 1]
  }

  I <- pt$I
  list(beta = b, converged = fit$converged,
       se = sqrt((sum(Ut^2) + sum(s^2)) / I^2),   # full sandwich, eq 7.31
       se_beta_only = sqrt(sum(Ut^2) / I^2),       # ignores uncertainty in alpha
       se_model = fit$se, n_groups = G)
}
