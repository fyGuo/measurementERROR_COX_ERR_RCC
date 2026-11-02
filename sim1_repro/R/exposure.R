#### Annual exposure generation (paper Section 4.1) ####

# correlation matrix of annual surrogate exposure across T years
cor_matrix <- function(T, cor) {
  gap <- abs(outer(1:T, 1:T, "-"))
  switch(cor$type,
         ind = diag(T),
         ar1 = cor$rho^gap,
         dex = cor$rho^(gap^cor$theta))  # damped exponential, Munoz et al. 1992
}

# n x T draws from MVN(0, Sigma)
rmvn0 <- function(n, Sigma) {
  matrix(rnorm(n * ncol(Sigma)), n) %*% chol(Sigma)
}

# surrogate C and true c, both n x T (column l = calendar year l)
#   additive:        C ~ MVN(mu, Sigma),        c = g0 + g1 * C + e
#   multiplicative:  log C ~ MVN(mu, Sigma),    log c = g0 + g1 * log C + e
# Z0, e0: optional standard normal draws (n x T), reused during calibration
gen_exposure <- function(n, T, error, cor, ep, Z0 = NULL, e0 = NULL) {
  if (is.null(Z0)) Z0 <- rmvn0(n, cor_matrix(T, cor))
  if (is.null(e0)) e0 <- matrix(rnorm(n * T), n)
  Z <- ep$sd * Z0
  e <- ep$sd_e * e0
  if (error == "add") {
    C <- ep$mu + Z
    c <- ep$g0 + ep$g1 * C + e
  } else {
    logC <- ep$mu + Z
    C <- exp(logC)
    c <- exp(ep$g0 + ep$g1 * logC + e)
  }
  list(C = C, c = c)
}

# L[l, m] = 1 if l < m, so (annual %*% L)[, m] = cumulative exposure through year m - 1;
# column m therefore holds x_{m-1}, the exposure used for interval (m - 1, m]
cum_lag_matrix <- function(T) {
  L <- matrix(0, T, T)
  L[upper.tri(L)] <- 1
  L
}
