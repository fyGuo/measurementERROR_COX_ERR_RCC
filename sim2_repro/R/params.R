#### Simulation 2 settings (paper Section 4.2, Table 2) ####

# design (paper Section 4.2 and Table 2 footnote), except the sample sizes: the paper's n1 = 2000 and
# n2 = 500 give about 200 events, too few for the ERR fit to settle even on the true x (README)
N_MAIN       <- 6000   # main study subjects (paper 2000)
T_STAR       <- 60     # years of follow-up
N_VAL_SUBJ   <- 2000   # internal validation subjects (paper 500)
MAX_PER_ID   <- 8      # pairs kept per validation subject (all years if fewer)
MIN_RS_PAIRS <- 20     # merge validation risk sets until each has > 20 pairs, as in Simulation 1
B_REPS       <- 2000

# true effects (Table 2), paired by the paper through median x = 0.204 (see EXPO)
# scenario, set by the environment variable SIM2_SCENARIO:
#   "paper" (default): the paper's beta above
#   "null": beta = 0 in both models, so both risk models are correctly specified; the hazard then
#           ignores x, Cox and ERR truth generate the same data and only the cox cells (1-6) are run
#   "sim1": Simulation 1's beta, log(1.5) (Cox) and 0.423 (ERR); at these smaller effects the Cox and ERR
#           fits on the true x differ by about 10% instead of about 40%
#   "beta1": beta = 1 in both models; the cox cells (1-6) are the same as in "paper"
SCENARIO <- Sys.getenv("SIM2_SCENARIO", "paper")
BETA <- switch(SCENARIO, null = c(cox = 0, err = 0), sim1 = c(cox = log(1.5), err = 0.423),
               beta1 = c(cox = 1, err = 1), c(cox = 1, err = 1.109))

# Weibull baseline hazard, as in Simulation 1; nu calibrated per (model, error) to TARGET_CIR
THETA      <- 4
TARGET_CIR <- 0.10

# surrogate C first, compound symmetry across years with correlation RHO_CS (paper: 0.5),
# then true c from the paper's error model, as in Simulation 1:
#   add:  C ~ MVN(mean_C, sd_C^2 R),          c = g0 + g1 C + e
#   mult: log C ~ MVN(mu, sd_log^2 R),        log c = g0 + g1 log C + e
# with e ~ N(0, sd_e^2) independent across subjects and years. Held fixed across rho: the
# distribution of C, the slope g1 = 2/3 and g0, so the naive bias is about -33% at every rho, as in
# the paper (mult: g1 = 1, E(c | C) proportional to C, with E(c) = mean_c). sd_e is solved from the
# paper's rho(c, C), the correlation of c and C within a year (on the original scale for both
# error models); for mult, g0 then keeps E(c) = mean_c.
# The paper's sd of 15 is not usable with beta = 1. mean_C = 0.012 puts the median x near the
# paper's 0.204, where 1 + beta x is not close to beta x, so ERR beta is identified.
# Additive: g0 = 0, so up to 27% of c (rho = 0.3) and 2% of C are negative and are dropped by the log
# calibration. Raising g0 to 0.024 made c almost always positive but raised the level of x relative to
# its spread, and ERR lost identification even on the true x (README, intercept test).
RHO_CS <- 0.5
EXPO <- list(add = list(mean_C = 0.012, sd_C = 0.006, g0 = 0, g1 = 2 / 3),
             mult = list(mean_C = 0.012, mean_c = 0.008, sd_log = 0.8))

expo_params <- function(error, rho) {
  e <- EXPO[[error]]
  if (error == "add") {
    # corr(c, C) = g1 sd_C / sqrt(g1^2 sd_C^2 + sd_e^2)
    list(mu = e$mean_C, sd = e$sd_C, g0 = e$g0, g1 = e$g1, sd_e = e$g1 * e$sd_C * sqrt(1 / rho^2 - 1))
  } else {
    # c = exp(g0) C exp(e): corr(c, C)^2 = (exp(s^2) - 1) / (exp(s^2) exp(sd_e^2) - 1), s = sd_log
    s2 <- e$sd_log^2
    sd_e <- sqrt(log(((exp(s2) - 1) / rho^2 + 1) / exp(s2)))
    list(mu = log(e$mean_C) - s2 / 2, sd = e$sd_log, g1 = 1, sd_e = sd_e,
         g0 = log(e$mean_c / e$mean_C) - sd_e^2 / 2)
  }
}

# Table 2 cells in paper order
CELLS <- expand.grid(rho = c(0.3, 0.6, 0.9), error = c("add", "mult"), model = c("cox", "err"),
                     stringsAsFactors = FALSE)[, 3:1]
CELLS$cell <- seq_len(nrow(CELLS))

# the six fits of Table 2: (risk model, calibration); NA = naive
FITS <- data.frame(name = c("naive_cox", "naive_err", "cox_A", "cox_M", "err_A", "err_M"),
                   model = c("cox", "err", "cox", "cox", "err", "err"),
                   cal = c(NA, NA, "add", "mult", "add", "mult"), stringsAsFactors = FALSE)

# baseline scale nu per cell, written by scripts/01_calibrate.R
cal_file <- file.path(PROJ, "results", if (SCENARIO == "paper") "calibration.csv" else sprintf("calibration_%s.csv", SCENARIO))
if (SCENARIO == "null") CELLS <- CELLS[CELLS$model == "cox", ]
if (file.exists(cal_file)) {
  CELLS <- merge(CELLS, read.csv(cal_file)[, c("cell", "nu")], by = "cell")
  CELLS <- CELLS[order(CELLS$cell), ]
}
