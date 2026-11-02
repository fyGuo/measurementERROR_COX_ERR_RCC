#### Calibrate the Weibull scale nu per cell to the target cumulative incidence ####
# run from sim2_repro/:  Rscript scripts/01_calibrate.R   (SIM2_SCENARIO=null for beta = 0)
# true x depends on rho(c, C) through the error e, so nu is calibrated per cell; the output
# also checks the realised corr(c, C) within a year and the level of x

source("R/load.R")
set.seed(2026110201)
n_cal <- 50000
L <- cum_lag_matrix(T_STAR)
w <- (1:T_STAR)^THETA - (0:(T_STAR - 1))^THETA

out <- do.call(rbind, lapply(CELLS$cell, function(k) {
  cl <- CELLS[CELLS$cell == k, ]
  ex <- gen_exposure2(n_cal, T_STAR, cl$error, expo_params(cl$error, cl$rho))
  x <- ex$c %*% L
  X <- ex$C %*% L
  # expected CIR = mean(1 - exp(-nu^theta * A_i)), A_i = sum_m r(x_{m-1}) w_m
  A <- as.vector(rel_rate(x, cl$model, BETA[[cl$model]]) %*% w)
  nu <- exp(uniroot(function(lnu) mean(1 - exp(-exp(lnu)^THETA * A)) - TARGET_CIR, c(-15, 5), tol = 1e-12)$root)
  data.frame(cell = k, model = cl$model, error = cl$error, rho = cl$rho, nu = nu,
             cor_cC = cor(ex$c[, 30], ex$C[, 30]), mean_c = mean(ex$c), mean_C = mean(ex$C),
             x_med = median(x), x45_mean = mean(x[, 46]), x45_sd = sd(x[, 46]), X45_sd = sd(X[, 46]),
             slope_x45 = unname(coef(lm(x[, 46] ~ X[, 46]))[2]),
             p_neg_C = mean(ex$C < 0), p_neg_c = mean(ex$c < 0), p_err_neg = mean(1 + BETA[["err"]] * x <= 0))
}))
print(format(out, digits = 3), row.names = FALSE)
write.csv(out, cal_file, row.names = FALSE)
