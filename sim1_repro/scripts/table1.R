#### Build Table 1 from saved replicates, next to the paper's values ####
# run from sim1_repro/:  Rscript scripts/table1.R <tag>

source("R/load.R")

args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args) >= 1) args[1] else "pilot"
files <- list.files(file.path(PROJ, "results", tag), pattern = "^cell\\d+\\.rds$", full.names = TRUE)
raw <- lapply(files, readRDS)
cols <- unique(unlist(lapply(raw, names)))
raw <- do.call(rbind, lapply(raw, function(x) { for (v in setdiff(cols, names(x))) x[[v]] <- NA; x[cols] }))
if (!"error" %in% names(raw)) raw$error <- NA
cat("failed replicates:", sum(!is.na(raw$error)), "\n")
if (any(!is.na(raw$error))) print(table(raw$cell[!is.na(raw$error)], raw$error[!is.na(raw$error)]))

summ <- function(b, se, beta) {
  c(relbias = 100 * (mean(b) - beta) / beta, se = mean(se), se_sd = mean(se) / sd(b),
    cover = 100 * mean(abs(b - beta) <= qnorm(0.975) * se))
}

tab <- do.call(rbind, lapply(split(raw, raw$cell), function(d) {
  ok <- is.na(d$error) & d$conv %in% TRUE
  d <- d[ok, ]
  cl <- CELLS[CELLS$cell == d$cell[1], ]
  n1 <- if ("n1" %in% names(d)) d$n1[1] else N_MAIN
  tr <- summ(d$b_true, d$se_true, d$beta[1])
  nv <- summ(d$b_naive, d$se_naive, d$beta[1])
  rr <- summ(d$b_rrc, d$se_rrc, d$beta[1])
  data.frame(cell = cl$cell, model = cl$model, error = cl$error, cor = cl$cor, n1 = n1, B = nrow(d),
             cir = round(100 * mean(d$n_events) / n1, 1),
             true_bias = tr[["relbias"]],
             naive_bias = nv[["relbias"]], rrc_bias = rr[["relbias"]],
             naive_se = nv[["se"]], rrc_se = rr[["se"]],
             naive_se_sd = nv[["se_sd"]], rrc_se_sd = rr[["se_sd"]],
             naive_cov = nv[["cover"]], rrc_cov = rr[["cover"]])
}))

# paper Table 1 (rows without "resample from data")
paper <- data.frame(
  cell = 1:12,
  p_naive_bias = c(-46.97, -52.40, -52.15, -48.22, -52.69, -52.74, -49.34, -53.41, -53.98, -50.12, -55.56, -56.03),
  p_rrc_bias   = c(-3.32, -4.55, -4.80, -3.64, -4.90, -5.04, -4.78, -6.14, -6.99, -4.96, -6.86, -7.80),
  p_naive_se   = c(0.054, 0.056, 0.055, 0.059, 0.059, 0.061, 0.077, 0.086, 0.089, 0.083, 0.092, 0.091),
  p_rrc_se     = c(0.093, 0.105, 0.105, 0.096, 0.111, 0.112, 0.106, 0.111, 0.110, 0.110, 0.114, 0.113),
  p_rrc_cov    = c(94.80, 94.30, 94.25, 94.55, 94.25, 94.15, 94.40, 94.25, 94.10, 94.20, 94.05, 94.00))
tab <- merge(tab, paper, by = "cell")

options(width = 200)
print(format(tab, digits = 3), row.names = FALSE)
write.csv(tab, file.path(PROJ, "results", sprintf("table1_%s.csv", tag)), row.names = FALSE)
