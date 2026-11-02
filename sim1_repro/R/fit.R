#### Partial likelihood fitting for Cox and ERR models (Breslow ties) ####

# Exposure enters as long vectors over (event, risk-set member):
#   zz: exposure of each risk-set member at that event, eid: event index of each entry
#   zc: exposure of the case at each event
# This lets exposure differ by risk set, which RRC needs.

# exposure from an n x T matrix X (column m = x_{m-1}) at each event
z_long <- function(X, ev) X[cbind(ev$subj, ev$k[ev$eid] + 1)]
z_case <- function(X, ev) X[cbind(ev$id, ev$k + 1)]

esum <- function(x, eid) as.vector(rowsum(x, eid, reorder = FALSE))

pl_parts <- function(b, zz, eid, zc, model) {
  if (model == "cox") {
    r <- exp(b * zz)
    S0 <- esum(r, eid); S1 <- esum(r * zz, eid); S2 <- esum(r * zz^2, eid)
    m1 <- S1 / S0
    list(ll = sum(b * zc - log(S0)), U = sum(zc - m1), I = sum(S2 / S0 - m1^2),
         r = r, S0 = S0, gbar = m1)
  } else {
    r <- 1 + b * zz; rc <- 1 + b * zc
    S0 <- esum(r, eid); Sz <- esum(zz, eid)
    list(ll = sum(log(rc) - log(S0)), U = sum(zc / rc - Sz / S0),
         I = sum(zc^2 / rc^2 - (Sz / S0)^2),
         r = r, S0 = S0, Sz = Sz, gbar = Sz / S0)
  }
}

# Maximise the partial likelihood with stats::nlminb (PORT trust-region optimiser,
# box constraints) using the analytic score and information from pl_parts().
# ERR is restricted to 1 + b z > 0 for all z.
pl_fit <- function(zz, eid, zc, model, start = 0, tol = 1e-10, maxit = 200) {
  lo <- -Inf; hi <- Inf
  if (model == "err") {
    if (any(zz > 0)) lo <- -1 / max(zz)
    if (any(zz < 0)) hi <- -1 / min(zz)
    hi <- min(hi, 50)                      # ERR likelihood can keep rising as beta -> Inf
    if (start <= lo) start <- lo / 2
    if (start >= hi) start <- hi / 2
  }
  # objective, gradient and hessian share one pl_parts() call per beta
  last_b <- NULL; last <- NULL
  parts <- function(b) {
    if (!identical(b, last_b)) { last <<- pl_parts(b, zz, eid, zc, model); last_b <<- b }
    last
  }
  opt <- nlminb(start,
                objective = function(b) { ll <- parts(b)$ll; if (is.finite(ll)) -ll else Inf },
                gradient  = function(b) -parts(b)$U,
                hessian   = function(b) as.matrix(parts(b)$I),
                lower = if (is.finite(lo)) lo + 1e-8 * abs(lo) else -Inf,
                upper = if (is.finite(hi)) hi - 1e-8 * abs(hi) else Inf,
                control = list(rel.tol = tol, iter.max = maxit, eval.max = 2 * maxit))
  b <- opt$par
  cur <- parts(b)
  converged <- opt$convergence == 0
  # an estimate stuck at the domain boundary is not a valid MLE
  if (is.finite(hi) && hi - b < 1e-3 * abs(hi)) converged <- FALSE
  if (is.finite(lo) && b - lo < 1e-3 * abs(lo)) converged <- FALSE
  list(beta = b, se = if (cur$I > 0) 1 / sqrt(cur$I) else NA_real_,
       converged = converged, parts = cur)
}
