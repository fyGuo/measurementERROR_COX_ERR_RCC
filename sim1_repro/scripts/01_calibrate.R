#### Calibrate exposure sd and baseline hazard scale nu per cell ####
# run from sim1_repro/: Rscript scripts/01_calibrate.R
#   sd: the cell's target (CELLS$target, CELLS$target_value in R/params.R)
#       sd_x45 / sd_X45: sd of true / surrogate cumulative exposure at year 45
#       naive_se: expected naive model-based SE of beta at n1 = N_MAIN
#   nu: expected cumulative incidence = TARGET_CIR
# cells are calibrated in parallel (future + furrr)

source("R/load.R")
library(future)
library(furrr)

n_big <- 20000
n_se <- 8           # naive_se target: mean SE over n_se cohorts of N_MAIN subjects
L <- cum_lag_matrix(T_STAR)
w <- (1:T_STAR)^THETA - (0:(T_STAR - 1))^THETA

calibrate_cell <- function(i) {
  cl <- CELLS[i, ]
  beta <- BETA[[cl$model]]
  set.seed(100 + i)                          # same draws for every trial sd
  Z0 <- rmvn0(n_big, cor_matrix(T_STAR, COR[[cl$cor]]))
  e0 <- matrix(rnorm(n_big * T_STAR), n_big)
  draw <- function(sd) gen_exposure(n_big, T_STAR, cl$error, COR[[cl$cor]],
                                    expo_params(cl$error, sd, cl$rho), Z0, e0)

  # expected CIR = mean(1 - exp(-nu^theta * A_i)), A_i = sum_m r(x_{m-1}) w_m
  calib_nu <- function(x) {
    A <- as.vector(rel_rate(x, cl$model, beta) %*% w)
    f_cir <- function(lnu) mean(1 - exp(-exp(lnu)^THETA * A)) - TARGET_CIR
    exp(uniroot(f_cir, c(-15, 5), tol = 1e-12)$root)
  }

  # mean naive SE over n_se cohorts of size N_MAIN (not one big cohort rescaled by
  # sqrt(n): with heavy-tailed exposure a bigger cohort holds more extreme subjects);
  # same event-time draws for every trial sd
  naive_se <- function(sd) {
    ex <- draw(sd)
    nu <- calib_nu(ex$c %*% L)
    set.seed(200 + i)
    mean(sapply(seq_len(n_se) - 1, function(b) {
      rows <- b * N_MAIN + seq_len(N_MAIN)
      x <- ex$c[rows, ] %*% L; X <- ex$C[rows, ] %*% L
      sv <- gen_survival(x, cl$model, beta, nu, THETA)
      ev <- make_events(sv$time, sv$D, sv$k)
      pl_fit(z_long(X, ev), ev$eid, z_case(X, ev), cl$model, start = beta / 2)$se
    }))
  }

  f_target <- switch(cl$target,
    sd_x45   = function(sd) sd((draw(sd)$c %*% L)[, 46]),
    sd_X45   = function(sd) sd((draw(sd)$C %*% L)[, 46]),
    naive_se = naive_se)
  bracket <- switch(cl$target,
    sd_x45   = c(0.01, if (cl$error == "add") 0.5 else 2.5),
    sd_X45   = c(0.01, if (cl$error == "add") 2 else 4),
    naive_se = if (cl$model == "err") c(1, 3) else c(0.3, 2.5))
  sd_hat <- uniroot(function(sd) f_target(sd) - cl$target_value, bracket, tol = 1e-6)$root
  ex <- draw(sd_hat)
  x <- ex$c %*% L; X <- ex$C %*% L
  nu <- calib_nu(x)

  # ERR domain: estimates are capped at -1 / min(exposure)
  data.frame(cell = cl$cell, model = cl$model, error = cl$error, cor = cl$cor,
             target = cl$target, target_value = cl$target_value, rho = cl$rho,
             sd = sd_hat, nu = nu,
             x_mean_k45 = mean(x[, 46]), x_sd_k45 = sd(x[, 46]), X_sd_k45 = sd(X[, 46]),
             x_max = max(x), x_min = min(x), X_min = min(X),
             cap_true = if (min(x) < 0) -1 / min(x) else Inf,
             cap_naive = if (min(X) < 0) -1 / min(X) else Inf,
             cor_cC = cor(as.vector(ex$c), as.vector(ex$C)))
}

plan(multisession, workers = min(nrow(CELLS), max(1, availableCores() - 1)))
# each cell sets its own seed, so future must not set or check seeds
out <- future_map(seq_len(nrow(CELLS)), calibrate_cell, .options = furrr_options(seed = NULL))
plan(sequential)

res <- do.call(rbind, out)
print(res, digits = 3)
write.csv(res, file.path(PROJ, "results", "calibration.csv"), row.names = FALSE)
