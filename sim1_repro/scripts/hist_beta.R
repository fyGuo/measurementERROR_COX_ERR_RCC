#### Histograms of beta-hat per Table 1 cell against a normal curve ####
# run from sim1_repro/:  Rscript scripts/hist_beta.R <tag>
# writes results/hist_beta_<tag>_true.png and _rrc.png; uses every replicate, converged or not

source("R/load.R")
args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args) >= 1) args[1] else "full_v2"
FILL <- "#2a78d6"; EDGE <- "#fcfcfb"; INK <- "#0b0b0b"; MUTED <- "#52514e"; REF <- "#C2410C"
lab <- c(cox = "Cox", err = "ERR", add = "Additive", mult = "Multiplicative",
         ind = "independent", ar1 = "AR(1)", dex = "DEX")
skew <- function(x) mean(((x - mean(x)) / sd(x))^3)

panel <- function(cl, col) {
  d <- readRDS(file.path(PROJ, "results", tag, sprintf("cell%02d.rds", cl$cell)))
  b <- d[[col]]; beta <- d$beta[1]
  rng <- range(b, beta); br <- seq(rng[1], rng[2], length.out = 41)
  h <- hist(b, breaks = br, plot = FALSE)
  xs <- seq(rng[1], rng[2], length.out = 300); ns <- dnorm(xs, mean(b), sd(b))
  plot(h, freq = FALSE, col = FILL, border = EDGE, main = "", xlab = "", ylab = "", axes = FALSE,
       xlim = rng, ylim = c(0, max(h$density, ns) * 1.12))
  axis(1, col = MUTED, col.axis = MUTED, cex.axis = 0.9, tcl = -0.3)
  lines(xs, ns, col = INK, lwd = 1.6)
  abline(v = beta, col = REF, lwd = 1.6, lty = 2)
  title(main = sprintf("%s, %s, %s", lab[[cl$model]], lab[[cl$error]], lab[[cl$cor]]),
        adj = 0, cex.main = 1, font.main = 2, col.main = INK, line = 0.9)
  mtext(sprintf("skew %.2f   mean %.3f   median %.3f", skew(b), mean(b), median(b)),
        side = 3, adj = 0, line = -0.1, cex = 0.72, col = MUTED)
}

fig <- function(col, name, ttl) {
  png(file.path(PROJ, "results", sprintf("hist_beta_%s_%s.png", tag, name)), 1800, 2000, res = 170)
  par(mfrow = c(4, 3), mar = c(2.6, 1, 2.6, 1), oma = c(0, 0, 3.2, 0), family = "sans")
  for (k in CELLS$cell) panel(CELLS[CELLS$cell == k, ], col)
  mtext(ttl, side = 3, outer = TRUE, adj = 0.01, line = 1.7, font = 2, cex = 1.05, col = INK)
  mtext("Bars: all 2,000 replicates.  Black curve: normal with the same mean and SD.  Orange dashed line: true beta.",
        side = 3, outer = TRUE, adj = 0.01, line = 0.5, cex = 0.72, col = MUTED)
  dev.off()
}
fig("b_true", "true", "Estimator using true exposure x (no measurement error)")
fig("b_rrc", "rrc", "RRC-corrected estimator")
