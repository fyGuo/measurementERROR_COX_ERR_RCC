#### Simulation 1 settings (paper Section 4.1, Table 1) ####

# design (Table 1 footnote)
N_MAIN       <- 2400   # main study subjects
T_STAR       <- 55     # years of follow-up
N_VAL_SUBJ   <- 800    # internal validation subjects
N_VAL_PAIR   <- 3000   # validated (c, C) pairs
MAX_PER_ID   <- 8      # max pairs per validation subject
MIN_RS_PAIRS <- 20     # merge validation risk sets until each has > 20 pairs (Section 5.1)
B_REPS       <- 2000

# true effects (Table 1); Table 6 would use log(5) and 1.905
BETA <- c(cox = log(1.5), err = 0.423)

# Weibull baseline hazard: shape fixed, scale nu calibrated per cell to TARGET_CIR
THETA      <- 4
TARGET_CIR <- 0.10     # paper: cumulative incidence 5-10%; upper end for more events

# correlation of annual surrogate exposure across years
# AR(1) and DEX values are the Uminer estimates in Uminer_Data_Feature_Exploration.R
COR <- list(ind = list(type = "ind"),
            ar1 = list(type = "ar1", rho = 0.62),
            dex = list(type = "dex", rho = 0.64, theta = 0.66))

# annual exposure (units 100 WLM) and error model; not reported in the paper (see README)
#   add:  C ~ MVN(mu, sd^2 R),      c = g0 + g1 C + N(0, sd_e^2)
#   mult: log C ~ MVN(mu, sd^2 R),  log c = g0 + g1 log C + N(0, sd_e^2)
# held fixed: mean annual C and c, slope g1, and correlation rho between c and C
# (on the log scale for mult; rho can be overridden per cell, see CELLS below);
# sd is calibrated per cell in scripts/01_calibrate.R to the cell's target.
# Both forms give E(c | C) slope 0.5 (mult: g1 = 1 so E(c | C) = (mean_c / mean_C) C).
# Low exposure level (mean x about 1.1 at year 45, beta x << 1) keeps ERR identifiable
# and keeps additive exposures above -1 / beta.
EXPO_BASE <- list(
  add  = list(mean_C = 0.05, mean_c = 0.025, g1 = 0.5, rho = 0.55),
  mult = list(mean_C = 0.05, mean_c = 0.025, g1 = 1.0, rho = 0.60)
)

expo_params <- function(error, sd, rho = EXPO_BASE[[error]]$rho) {
  b <- EXPO_BASE[[error]]
  sd_e <- b$g1 * sd * sqrt(1 / rho^2 - 1)
  if (error == "add") {
    list(mu = b$mean_C, sd = sd, g0 = b$mean_c - b$g1 * b$mean_C, g1 = b$g1, sd_e = sd_e)
  } else {
    mu <- log(b$mean_C) - sd^2 / 2
    list(mu = mu, sd = sd, g1 = b$g1, sd_e = sd_e,
         g0 = log(b$mean_c) - b$g1 * mu - (b$g1^2 * sd^2 + sd_e^2) / 2)
  }
}

# Table 1 cells in paper order ("resample from data" row not reproduced)
CELLS <- expand.grid(cor = c("ind", "ar1", "dex"), error = c("add", "mult"),
                     model = c("cox", "err"), stringsAsFactors = FALSE)[, 3:1]
CELLS$cell <- seq_len(nrow(CELLS))

# per-cell calibration target for the exposure sd (scripts/01_calibrate.R)
#   sd_x45:   sd of true cumulative exposure at year 45
#   sd_X45:   sd of surrogate cumulative exposure at year 45
#   naive_se: expected naive model-based SE of beta at n1 = 2400
# Cox additive: sd_X45 = 1.17 gives the paper's naive SE (about 0.055).
# Cox multiplicative and ERR multiplicative: matched to the paper's naive SE
#   directly; the sd_X45 route overshoots with lognormal C (its tail carries extra
#   information). rho = 0.8 keeps the error sd, and so the tail of true x, moderate.
# ERR additive: the aim is to show RRC removes the measurement-error bias, not to match
#   the paper's SE. Normal C is kept; too little spread leaves ERR estimates skewed,
#   too much makes some x negative and caps beta (1 + beta x > 0). A grid over
#   rho in {0.55, 0.7, 0.85} and sd_X45 in {0.6, 0.75, 0.9, 1.2} gave sd_X45 = 0.9,
#   rho = 0.85 as the best setting for all three structures (fewest negative x).
CELLS$target <- ifelse(CELLS$error == "add", "sd_X45", "naive_se")
CELLS$target_value <- c(1.17, 1.17, 1.17, 0.059, 0.059, 0.061,
                        0.9, 0.9, 0.9, 0.083, 0.092, 0.091)
CELLS$rho <- ifelse(CELLS$error == "mult", 0.8,
                    ifelse(CELLS$model == "err", 0.85, EXPO_BASE$add$rho))
# ERR multiplicative independent: at rho = 0.8 the error sd (1.51) was large enough for deaths to
#   deplete the highest true exposures, which the log-linear calibration cannot follow, leaving an
#   RRC bias of about -20% at any n. rho = 0.95 keeps the naive SE target with a smaller error sd.
CELLS$rho[CELLS$cell == 10] <- 0.95

# Resulting settings (from results/calibration.csv, 2026-10-30; earlier: calibration_v2.csv). sd is the sd of C
# (additive) or of log C (multiplicative); sd_e is the error sd.
#   cell model error cor  C        target   value  rho   sd      sd_e      nu
#    1   cox   add   ind  normal   sd_X45   1.170  0.55  0.1764  0.1339    0.00908
#    2   cox   add   ar1  normal   sd_X45   1.170  0.55  0.0869  0.0660    0.00917
#    3   cox   add   dex  normal   sd_X45   1.170  0.55  0.0627  0.0476    0.00920
#    4   cox   mult  ind  lognorm  naive_se 0.059  0.80  1.3709  1.0282    0.00918
#    5   cox   mult  ar1  lognorm  naive_se 0.059  0.80  1.1223  0.8417    0.00921
#    6   cox   mult  dex  lognorm  naive_se 0.061  0.80  0.9226  0.6920    0.00922
#    7   err   add   ind  normal   sd_X45   0.900  0.85  0.1350  0.0418    0.00943
#    8   err   add   ar1  normal   sd_X45   0.900  0.85  0.0668  0.0207    0.00943
#    9   err   add   dex  normal   sd_X45   0.900  0.85  0.0482  0.0149    0.00943
#   10   err   mult  ind  lognorm  naive_se 0.083  0.95  2.1885  0.7193    0.00945
#   11   err   mult  ar1  lognorm  naive_se 0.092  0.80  1.8059  1.3544    0.00946
#   12   err   mult  dex  lognorm  naive_se 0.091  0.80  1.7200  1.2900    0.00945

# calibrated exposure sd and baseline scale nu per cell (written by scripts/01_calibrate.R)
cal_file <- file.path(PROJ, "results", "calibration.csv")
if (file.exists(cal_file)) {
  CELLS <- merge(CELLS, read.csv(cal_file)[, c("cell", "sd", "nu")], by = "cell")
}
