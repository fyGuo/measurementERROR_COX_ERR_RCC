#### Format results/table1_<tag>.csv in the layout of the paper's Table 1 (markdown) ####
# run from sim1_repro/ after scripts/table1.R:  Rscript scripts/table1_format.R <tag>

source("R/load.R")

args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args) >= 1) args[1] else "full"
tab <- read.csv(file.path(PROJ, "results", sprintf("table1_%s.csv", tag)))

# median relative bias, as a robust companion to the mean (ERR estimates are right-skewed)
med <- do.call(rbind, lapply(sprintf("%s/cell%02d.rds", file.path(PROJ, "results", tag), tab$cell), function(f) {
  d <- readRDS(f); d <- d[d$conv %in% TRUE, ]; b <- d$beta[1]
  data.frame(cell = d$cell[1], naive_med = 100 * (median(d$b_naive) - b) / b,
             rrc_med = 100 * (median(d$b_rrc) - b) / b)
}))
tab <- merge(tab, med, by = "cell")

lab <- c(add = "Additive", mult = "Multiplicative", ind = "Independent", ar1 = "AR(1)", dex = "DEX",
         cox = "Cox", err = "ERR")
f <- function(x, d) formatC(x, format = "f", digits = d)
pair <- function(ours, paper, d) sprintf("%s (%s)", f(ours, d), f(paper, d))

out <- c(sprintf("Table 1 reproduction (tag `%s`). Each entry: ours (paper). Resample-from-data rows not reproduced.", tag), "",
  "| Model | Error | Generation of C | B | Naive rel. bias % | RRC rel. bias % | Naive SE | RRC SE | RRC 95% CP | RRC SE/SD | Median rel. bias %: naive / RRC |",
  "|---|---|---|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(tab))) {
  r <- tab[i, ]
  out <- c(out, sprintf("| %s | %s | %s | %d | %s | %s | %s | %s | %s | %s | %s / %s |",
    lab[r$model], lab[r$error], lab[r$cor], r$B,
    pair(r$naive_bias, r$p_naive_bias, 1), pair(r$rrc_bias, r$p_rrc_bias, 1),
    pair(r$naive_se, r$p_naive_se, 3), pair(r$rrc_se, r$p_rrc_se, 3),
    pair(r$rrc_cov, r$p_rrc_cov, 1), f(r$rrc_se_sd, 2), f(r$naive_med, 1), f(r$rrc_med, 1)))
}
writeLines(out, file.path(PROJ, "results", sprintf("table1_%s.md", tag)))
cat(out, sep = "\n")
