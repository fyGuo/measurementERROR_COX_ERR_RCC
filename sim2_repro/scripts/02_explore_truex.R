#### Explore exposure settings: Cox and ERR fitted to the true x, both true beta = 1 ####
# run from sim2_repro/:  Rscript scripts/02_explore_truex.R [n1] [cir] [reps] [grid]
# grid "base": level x spread x RHO_CS at the given cir; "focus": small or centred levels, cir 0.1-0.5
# For each exposure setting and generating model, fits both risk models to the true x and reports
# mean / median / SD of beta-hat. The aim is a setting where both fits are close to 1, which needs a
# small exposure level at event times (ERR -> about 1 / (1 - xbar) under Cox truth, Cox -> about
# 1 / (1 + xbar) under ERR truth) but enough spread in x for precision.

source("R/load.R")
library(future)
library(furrr)

args <- commandArgs(trailingOnly = TRUE)
n1 <- if (length(args) >= 1) as.integer(args[1]) else N_MAIN
cir <- if (length(args) >= 2) as.numeric(args[2]) else TARGET_CIR
reps <- if (length(args) >= 3) as.integer(args[3]) else 200
grid_name <- if (length(args) >= 4) args[4] else "base"
rho <- 0.6
L <- cum_lag_matrix(T_STAR)
w <- (1:T_STAR)^THETA - (0:(T_STAR - 1))^THETA

# exposure parameters for a setting (same forms as expo_params, with the level and spread free)
make_ep <- function(error, mean_C, spread, rho) {
  if (error == "add") {
    sd_C <- if (mean_C > 0 && spread < 5) spread * mean_C else spread / 1000; g1 <- 2 / 3
    list(mu = mean_C, sd = sd_C, g0 = 0, g1 = g1, sd_e = g1 * sd_C * sqrt(1 / rho^2 - 1))
  } else {
    s2 <- spread^2
    sd_e <- sqrt(log(((exp(s2) - 1) / rho^2 + 1) / exp(s2)))
    list(mu = log(mean_C) - s2 / 2, sd = spread, g1 = 1, sd_e = sd_e, g0 = log(2 / 3) - sd_e^2 / 2)
  }
}

calib_nu <- function(model, ep, error, rho_cs, cir) {
  set.seed(1)
  x <- gen_exposure2(20000, T_STAR, error, ep, rho_cs)$c %*% L
  A <- as.vector(pmax(rel_rate(x, model, 1), 0) %*% w)
  exp(uniroot(function(lnu) mean(1 - exp(-exp(lnu)^THETA * A)) - cir, c(-15, 5), tol = 1e-10)$root)
}

if (grid_name == "base") {
  grid <- expand.grid(model = c("cox", "err"), error = c("add", "mult"), mean_C = c(0.003, 0.006, 0.012),
                      spread_i = 1:3, rho_cs = c(0.5, 0.9), cir = cir, stringsAsFactors = FALSE)
  # spread: sd(C) / mean(C) for additive, sd(log C) for multiplicative
  grid$spread <- ifelse(grid$error == "add", c(0.5, 1, 2)[grid$spread_i], c(0.8, 1.2, 1.6)[grid$spread_i])
} else {
  # additive spread >= 5 is sd(C) x 1000 (needed for mean_C = 0, exposure centred at 0)
  ga <- expand.grid(model = c("cox", "err"), error = "add", mean_C = c(0, 0.003), spread = c(5, 8),
                    rho_cs = 0.9, cir = c(0.1, 0.3, 0.5), stringsAsFactors = FALSE)
  gm <- expand.grid(model = c("cox", "err"), error = "mult", mean_C = c(0.0015, 0.003), spread = c(0.5, 0.8),
                    rho_cs = 0.9, cir = c(0.1, 0.3, 0.5), stringsAsFactors = FALSE)
  grid <- rbind(ga, gm)
}

plan(multisession, workers = max(1, availableCores() - 1))
out <- do.call(rbind, lapply(seq_len(nrow(grid)), function(g) {
  s <- grid[g, ]
  ep <- make_ep(s$error, s$mean_C, s$spread, rho)
  nu <- calib_nu(s$model, ep, s$error, s$rho_cs, s$cir)
  r <- future_map(seq_len(reps), function(i) {
    set.seed(1000 * g + i)
    ex <- gen_exposure2(n1, T_STAR, s$error, ep, s$rho_cs)
    x <- ex$c %*% L
    sv <- gen_survival(x, s$model, 1, nu, THETA)
    ev <- make_events(sv$time, sv$D, sv$k)
    zz <- z_long(x, ev); zc <- z_case(x, ev)
    fc <- pl_fit(zz, ev$eid, zc, "cox", start = 1)
    fe <- pl_fit(zz, ev$eid, zc, "err", start = 1)
    c(cox = fc$beta, err = fe$beta, err_ok = fe$converged, xbar = mean(zc), sdx = sd(zc), D = ev$n,
      neg = mean(zz < 0))
  }, .options = furrr_options(seed = NULL))
  r <- do.call(rbind, r)
  ok <- r[, "err_ok"] == 1
  cbind(s[, c("model", "error", "mean_C", "spread", "rho_cs", "cir")],
        D = mean(r[, "D"]), xbar = mean(r[, "xbar"]), sd_x = mean(r[, "sdx"]), p_neg = mean(r[, "neg"]),
        cox_mean = mean(r[, "cox"]), cox_med = median(r[, "cox"]), cox_sd = sd(r[, "cox"]),
        err_mean = mean(r[ok, "err"]), err_med = median(r[, "err"]), err_sd = sd(r[ok, "err"]), err_nc = sum(!ok))
}))
plan(sequential)
options(width = 200)
print(format(out, digits = 3), row.names = FALSE)
write.csv(out, file.path(PROJ, "results", sprintf("explore_truex_%s_n%d_cir%g.csv", grid_name, n1, cir)), row.names = FALSE)
