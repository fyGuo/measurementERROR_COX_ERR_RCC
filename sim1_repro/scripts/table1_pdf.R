#### Render results/table1_<tag>.csv as a LaTeX table and compile to PDF ####
# run from sim1_repro/ after scripts/table1.R:  Rscript scripts/table1_pdf.R <tag>
# needs pdflatex (booktabs, multirow, pdflscape, threeparttable)

source("R/load.R")

args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args) >= 1) args[1] else "full"
tab <- read.csv(file.path(PROJ, "results", sprintf("table1_%s.csv", tag)))
tab <- tab[order(tab$cell), ]

# median relative bias, as a robust companion to the mean (ERR estimates are right-skewed)
tab$naive_med <- tab$rrc_med <- NA
for (i in seq_len(nrow(tab))) {
  d <- readRDS(file.path(PROJ, "results", tag, sprintf("cell%02d.rds", tab$cell[i])))
  d <- d[d$conv %in% TRUE, ]; b <- d$beta[1]
  tab$naive_med[i] <- 100 * (median(d$b_naive) - b) / b
  tab$rrc_med[i] <- 100 * (median(d$b_rrc) - b) / b
}

lab <- c(add = "Additive", mult = "Multiplicative", ind = "Independent", ar1 = "AR(1)", dex = "DEX",
         cox = "Cox", err = "ERR")
f <- function(x, d) sub("^-", "$-$", formatC(x, format = "f", digits = d))   # math minus
pp <- function(x, d) sprintf("{\\color{gray}(%s)}", f(x, d))   # paper value

rows <- character(0)
for (i in seq_len(nrow(tab))) {
  r <- tab[i, ]
  first_model <- (i - 1) %% 6 == 0
  first_error <- (i - 1) %% 3 == 0
  rows <- c(rows, paste0(
    if (first_model) sprintf("\\multirow{6}{*}{%s}", lab[r$model]) else "", " & ",
    if (first_error) sprintf("\\multirow{3}{*}{%s}", lab[r$error]) else "", " & ",
    lab[r$cor], " & ", r$B, " & ",
    f(r$naive_bias, 1), " ", pp(r$p_naive_bias, 1), " & ",
    f(r$rrc_bias, 1), " ", pp(r$p_rrc_bias, 1), " & ",
    f(r$naive_se, 3), " ", pp(r$p_naive_se, 3), " & ",
    f(r$rrc_se, 3), " ", pp(r$p_rrc_se, 3), " & ",
    f(r$rrc_cov, 1), " ", pp(r$p_rrc_cov, 1), " & ",
    f(r$rrc_se_sd, 2), " & ", f(r$naive_med, 1), " & ", f(r$rrc_med, 1), " \\\\"))
  if (i %% 6 == 0 && i < nrow(tab)) rows <- c(rows, "\\midrule")
  else if (i %% 3 == 0 && i < nrow(tab)) rows <- c(rows, "\\cmidrule(l){2-12}")
}

tex <- c(
  "\\documentclass[10pt]{article}",
  "\\usepackage[margin=0.6in,landscape]{geometry}",
  "\\usepackage{booktabs,multirow,threeparttable,xcolor}", "\\usepackage[labelformat=empty,justification=raggedright,singlelinecheck=false]{caption}",
  "\\pagestyle{empty}",
  "\\begin{document}",
  "\\begin{table}[ht]\\centering\\small",
  "\\begin{threeparttable}",
  sprintf("\\caption*{\\textbf{Table 1 (reproduction).} Simulation 1: relative bias, estimated SE and coverage of the naive and RRC estimators of $\\beta$ (B = 2000 replicates per cell).}"),
  "\\begin{tabular}{lll r rr rr r r rr}",
  "\\toprule",
  " & & & & \\multicolumn{2}{c}{Relative bias (\\%)} & \\multicolumn{2}{c}{Estimated SE} & 95\\% CP & SE/SD & \\multicolumn{2}{c}{Median rel.\\ bias (\\%)} \\\\",
  "\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}\\cmidrule(lr){9-9}\\cmidrule(lr){10-10}\\cmidrule(l){11-12}",
  "Model & Error & Generation of $C$ & $B_{\\mathrm{conv}}$ & Naive & RRC & Naive & RRC & RRC & RRC & Naive & RRC \\\\",
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\footnotesize",
  "\\item Entries are our estimates with the paper's Table 1 values in {\\color{gray}(gray)}. Rows with $C$ resampled from the Uminer data are not reproduced.",
  "\\item $n_1 = 2400$, $T^* = 55$, $n_2 = 800$ validation subjects with 3000 pairs (at most 8 per subject); true $\\beta = \\log 1.5$ (Cox) and $0.423$ (ERR); cumulative incidence about 10\\%.",
  "\\item $B_{\\mathrm{conv}}$: replicates in which all three fits (true, naive, RRC) converged; ERR fits ending on the parameter boundary are excluded.",
  "\\item Relative bias $= 100(\\bar{\\hat\\beta} - \\beta)/\\beta$; estimated SE is the mean model-based SE (naive) or sandwich SE (RRC); SE/SD is mean SE over the empirical SD of $\\hat\\beta$; 95\\% CP is coverage of the Wald interval. Median relative bias uses the median of $\\hat\\beta$ in place of the mean.",
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}",
  "\\end{document}")

out_dir <- file.path(PROJ, "results")
tex_file <- file.path(out_dir, sprintf("table1_%s.tex", tag))
writeLines(tex, tex_file)
build <- file.path(tempdir(), "latex"); dir.create(build, showWarnings = FALSE)
status <- system2("pdflatex", c("-interaction=nonstopmode", "-halt-on-error",
                                paste0("-output-directory=", build), tex_file), stdout = FALSE)
if (status != 0) stop("pdflatex failed; see log in ", build)
file.copy(file.path(build, sprintf("table1_%s.pdf", tag)), out_dir, overwrite = TRUE)
cat("wrote", file.path(out_dir, sprintf("table1_%s.pdf", tag)), "\n")
