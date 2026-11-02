#### Survival time generation (paper Section 4.1, Steps 1-6) ####

# relative rate r(x): Cox exp(beta x), ERR 1 + beta x (floored at 0, violations counted)
rel_rate <- function(x, model, beta) {
  if (model == "cox") exp(beta * x) else pmax(1 + beta * x, 0)
}

# x: n x T matrix, column m = cumulative true exposure x_{m-1} used for interval (m - 1, m]
# baseline hazard Weibull: Lambda0(t) = (nu * t)^theta
# returns observed time, event indicator and k = interval index (event in (k, k + 1])
gen_survival <- function(x, model, beta, nu, theta) {
  n <- nrow(x); T <- ncol(x)
  r <- rel_rate(x, model, beta)
  w <- (1:T)^theta - (0:(T - 1))^theta
  inc <- nu^theta * r * rep(w, each = n)
  H <- t(apply(inc, 1, cumsum))          # H[, j] = cumulative hazard at integer time j
  E <- rexp(n)                           # E = -log(1 - U)
  j <- rowSums(H <= E)                   # F_ij <= U_i < F_i,j+1
  event <- j < T
  Hj <- numeric(n)
  Hj[j > 0] <- H[cbind(which(j > 0), j[j > 0])]
  time <- rep(T, n)
  ie <- which(event)
  rj <- r[cbind(ie, j[ie] + 1)]
  time[ie] <- ((E[ie] - Hj[ie]) / (nu^theta * rj) + j[ie]^theta)^(1 / theta)
  list(time = time, D = as.integer(event), k = ifelse(event, j, NA_integer_),
       n_viol = if (model == "err") sum(1 + beta * x < 0) else 0L)
}

# event / risk-set structure in long format, events sorted by time
# risk set of event e = subjects with time >= t_e (all enter at time 0)
make_events <- function(time, D, k) {
  id <- which(D == 1)
  id <- id[order(time[id])]
  t_e <- time[id]
  risk <- lapply(t_e, function(t) which(time >= t))
  len <- lengths(risk)
  list(id = id, time = t_e, k = k[id], n = length(id),
       subj = unlist(risk), eid = rep(seq_along(id), len))
}
