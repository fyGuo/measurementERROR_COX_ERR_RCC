#### Annual exposure and validation sampling for Simulation 2 ####

# surrogate C and true c, both n x T (column l = year l); C first, then c from C
gen_exposure2 <- function(n, T, error, ep, rho_cs = RHO_CS) {
  # compound symmetry: shared subject effect plus independent yearly effect
  Z <- ep$sd * (sqrt(rho_cs) * rnorm(n) + sqrt(1 - rho_cs) * matrix(rnorm(n * T), n))
  e <- ep$sd_e * matrix(rnorm(n * T), n)
  if (error == "add") {
    C <- ep$mu + Z
    c <- ep$g0 + ep$g1 * C + e
  } else {
    C <- exp(ep$mu + Z)
    c <- exp(ep$g0 + ep$g1 * (ep$mu + Z) + e)
  }
  list(C = C, c = c)
}

# sample n_subj subjects and keep up to max_per_id of their follow-up years as (c, C) pairs
sample_validation2 <- function(time, n_subj, max_per_id, T) {
  vid <- sample.int(length(time), n_subj)
  nyr <- pmin(ceiling(time[vid]), T)
  keep <- pmin(nyr, max_per_id)
  l <- unlist(lapply(seq_along(vid), function(i) sort(sample.int(nyr[i], keep[i]))))
  data.frame(id = rep(vid, keep), l = l)
}
