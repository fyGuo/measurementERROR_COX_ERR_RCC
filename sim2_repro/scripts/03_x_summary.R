#### Distribution of the true cumulative exposure x per cell and scenario ####
# run from sim2_repro/:  Rscript scripts/03_x_summary.R
# x does not depend on beta, but the event times do; the x of cases at their event time and the
# risk-set mean x at event times (xbar) set how far apart the Cox and ERR fits are (about beta * xbar).
# One large data set (n = 20000) per cell and scenario -> results/x_summary.csv

source("R/load.R")
n_big <- 20000
L <- cum_lag_matrix(T_STAR)
q <- function(v) quantile(v, c(0.05, 0.25, 0.5, 0.75, 0.95))
scen <- list(paper = c(cox = 1, err = 1.109), sim1 = c(cox = log(1.5), err = 0.423), null = c(cox = 0, err = 0))

out <- do.call(rbind, lapply(names(scen), function(sc) {
  cal <- read.csv(file.path(PROJ, "results", if (sc == "paper") "calibration.csv" else sprintf("calibration_%s.csv", sc)))
  do.call(rbind, lapply(cal$cell, function(k) {
    cl <- CELLS[CELLS$cell == k, ]; nu <- cal$nu[cal$cell == k]; b <- scen[[sc]][[cl$model]]
    set.seed(5000 + k)
    x <- gen_exposure2(n_big, T_STAR, cl$error, expo_params(cl$error, cl$rho))$c %*% L
    sv <- gen_survival(x, cl$model, b, nu, THETA)
    ev <- make_events(sv$time, sv$D, sv$k)
    # person-years in follow-up: column m (interval m - 1 to m) while time > m - 1
    fu <- outer(sv$time, 0:(T_STAR - 1), ">")
    xpy <- x[fu]; xcase <- z_case(x, ev)
    xbar <- as.vector(rowsum(z_long(x, ev), ev$eid, reorder = FALSE)) / tabulate(ev$eid)
    data.frame(scenario = sc, cell = k, model = cl$model, error = cl$error, rho = cl$rho, beta = b,
               cir = mean(sv$D), event_age_med = median(ev$time),
               py_mean = mean(xpy), py_sd = sd(xpy), py_neg = mean(xpy < 0),
               t(setNames(q(xpy), paste0("py_q", c(5, 25, 50, 75, 95)))),
               x45_mean = mean(x[, 46]), x45_sd = sd(x[, 46]), x45_skew = mean(((x[, 46] - mean(x[, 46])) / sd(x[, 46]))^3),
               case_mean = mean(xcase), case_sd = sd(xcase),
               t(setNames(q(xcase), paste0("case_q", c(5, 25, 50, 75, 95)))),
               xbar_events = mean(xbar), beta_xbar = b * mean(xbar))
  }))
}))
options(width = 250)
print(format(out, digits = 3), row.names = FALSE)
write.csv(out, file.path(PROJ, "results", "x_summary.csv"), row.names = FALSE)
