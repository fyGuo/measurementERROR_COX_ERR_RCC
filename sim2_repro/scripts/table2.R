#### Summarise Simulation 2: mean and median estimated beta (and relative bias) of each fit per cell ####
# run from sim2_repro/:  Rscript scripts/table2.R <tag>
# <fit>_b, <fit>_bmed: mean over converged fits and median over all fits of beta-hat;
# relative bias of each fit is against the true beta of the risk model it fits (Cox 1, ERR 1.109);
# drop_pct is the share of validation pairs the log calibration drops (c or C <= 0)

source("R/load.R")
args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args) >= 1) args[1] else "pilot"

tab <- do.call(rbind, lapply(CELLS$cell, function(k) {
  f <- file.path(PROJ, "results", tag, sprintf("cell%02d.rds", k))
  if (!file.exists(f)) return(NULL)
  d <- readRDS(f)
  if ("error" %in% names(d)) d <- d[is.na(d$error), ]
  cl <- CELLS[CELLS$cell == k, ]
  row <- data.frame(cell = k, model = cl$model, error = cl$error, rho = cl$rho, B = nrow(d),
                    events = mean(d$n_events), nonpos = mean(d$n_nonpos) / (N_MAIN * T_STAR),
                    drop_pct = if ("n_pairs_pos" %in% names(d)) 100 * mean(1 - d$n_pairs_pos / d$n_pairs) else NA)
  # fits on the true x first, as references for the risk-model misspecification
  fits <- rbind(data.frame(name = c("true_cox", "true_err"), model = c("cox", "err")), FITS[, c("name", "model")])
  for (i in seq_len(nrow(fits))) {
    nm <- fits$name[i]; b <- BETA[[fits$model[i]]]; ok <- d[[paste0("conv_", nm)]] %in% TRUE
    row[[paste0(nm, "_b")]] <- mean(d[[paste0("b_", nm)]][ok])
    row[[paste0(nm, "_bmed")]] <- median(d[[paste0("b_", nm)]])
    # relative bias is undefined at beta = 0
    row[[nm]] <- if (b != 0) 100 * (mean(d[[paste0("b_", nm)]][ok]) / b - 1) else NA
    row[[paste0(nm, "_med")]] <- if (b != 0) 100 * (median(d[[paste0("b_", nm)]]) / b - 1) else NA
    row[[paste0(nm, "_nc")]] <- sum(!ok)
  }
  row
}))
options(width = 220)
print(format(tab, digits = 3), row.names = FALSE)
write.csv(tab, file.path(PROJ, "results", sprintf("table2_%s.csv", tag)), row.names = FALSE)
