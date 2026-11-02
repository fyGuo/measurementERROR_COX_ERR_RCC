#### Render results/table2_<tag>.csv in the layout of the paper's Table 2 and compile to PDF ####
# entries are estimated beta (mean, median); the paper's relative biases are converted to beta
# run from sim2_repro/ after scripts/table2.R:  Rscript scripts/table2_pdf.R <tag>
# needs pdflatex (booktabs, threeparttable, xcolor)

source("R/load.R")

args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args) >= 1) args[1] else "pilot"
tab <- read.csv(file.path(PROJ, "results", sprintf("table2_%s.csv", tag)))
tab <- tab[order(tab$cell), ]

# paper's Table 2, cells in CELLS order, columns in FITS order
paper <- matrix(c(
  -33.5, -28.3,  -4.4,  -5.7,  5.9,  6.6,
  -32.8, -27.9,  -4.0,  -5.5,  5.6,  6.1,
  -31.4, -27.2,  -3.8,  -5.4,  5.3,  5.7,
  -34.7, -29.0,  -6.1,  -5.7,  4.8,  4.3,
  -33.9, -28.4,  -5.6,  -5.3,  4.5,  4.1,
  -32.6, -27.8,  -5.4,  -5.0,  4.0,  3.7,
  -42.1, -35.8, -14.9, -15.2, -5.5, -6.6,
  -40.6, -33.0, -14.5, -14.9, -5.1, -6.4,
  -38.9, -31.2, -14.2, -14.7, -4.6, -6.2,
  -43.0, -36.2, -16.8, -16.0, -6.7, -6.3,
  -41.7, -33.4, -15.8, -15.6, -5.8, -6.1,
  -40.4, -31.9, -14.7, -15.2, -4.8, -5.7), ncol = 6, byrow = TRUE)

lab <- c(add = "Additive", mult = "Multiplicative", cox = "Cox", err = "ERR")
null <- SCENARIO == "null"
paper_vals <- SCENARIO == "paper"   # the paper's values only apply at its own beta
bfmt <- function(b) formatC(b, format = "fg", digits = 4)
truth <- if (null) c(cox = "Cox and ERR true $\\beta = 0$") else
  c(cox = sprintf("Cox true $\\beta = %s$", bfmt(BETA[["cox"]])), err = sprintf("ERR true $\\beta = %s$", bfmt(BETA[["err"]])))
true_lab <- if (null) "true: 0 for both models" else sprintf("true: Cox %s, ERR %s", bfmt(BETA[["cox"]]), bfmt(BETA[["err"]]))
f <- function(x, d = 1) ifelse(is.na(x), "--", sub("^-", "$-$", formatC(x, format = "f", digits = d)))
pp <- function(x) sprintf("{\\color{gray}(%s)}", f(x, 2))   # paper value
# correctly specified fit: RRC with the generating risk model and error model
# (at beta = 0 both risk models are correct, so both RRC fits with the right error model)
correct <- function(r) {
  m <- if (null) c("cox", "err") else r$model
  paste0(m, "_", substr(toupper(r$error), 1, 1))
}
# columns: fits on the true x (references, no paper value), then the paper's six fits
cols <- c("true_cox", "true_err", FITS$name)
fit_model <- c("cox", "err", FITS$model)
nk <- length(cols) + 2   # rho, fits, dropped pairs

panel <- function(cell_fun, last_fun) {
  rows <- character(0)
  for (i in seq_len(nrow(tab))) {
    r <- tab[i, ]
    if ((i - 1) %% 3 == 0) {
      if (i > 1) rows <- c(rows, "\\midrule")
      rows <- c(rows, sprintf("\\multicolumn{%d}{l}{\\textbf{%s}\\quad %s, %s independent error} \\\\",
                              nk, if (null) lab[r$error] else paste0(lab[r$model], "+", lab[r$error]),
                              truth[r$model], tolower(lab[r$error])))
    }
    cells <- vapply(seq_along(cols), function(j) {
      s <- cell_fun(r, i, j)
      if (cols[j] %in% correct(r)) sprintf("\\boldmath\\textbf{%s}", s) else s
    }, "")
    rows <- c(rows, paste0("\\quad ", r$rho, " & ", paste(cells, collapse = " & "), " & ", last_fun(r), " \\\\"))
  }
  rows
}

# panel A: mean over converged fits, paper value in gray
rows_a <- panel(function(r, i, j) {
  s <- f(r[[paste0(cols[j], "_b")]], 2)
  if (j > 2 && paper_vals) paste(s, pp(BETA[[fit_model[j]]] * (1 + paper[i, j - 2] / 100))) else s
}, function(r) "")
# panel B: median over all fits, with the number of non-converged fits in brackets
rows_b <- panel(function(r, i, j) {
  nc <- r[[paste0(cols[j], "_nc")]]
  paste0(f(r[[paste0(cols[j], "_bmed")]], 2), if (nc > 0) sprintf(" [%d]", nc) else "")
}, function(r) f(r$drop_pct))

head_rows <- c(
  "\\toprule",
  sprintf(" & \\multicolumn{%d}{c}{Estimated $\\beta$ (%s)} & Pairs \\\\", length(cols), true_lab),
  sprintf("\\cmidrule(lr){2-%d}", length(cols) + 1),
  " & \\multicolumn{2}{c}{True $x$} & \\multicolumn{2}{c}{Naive} & \\multicolumn{4}{c}{RRC corrected} & dropped \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-9}",
  "$\\rho(c, C)$ & Cox & ERR & Cox & ERR & Cox.A & Cox.M & ERR.A & ERR.M & (\\%) \\\\",
  "\\midrule")

x_med <- if (file.exists(cal_file)) {
  cal <- read.csv(cal_file); sprintf("%.2f (additive) and %.2f (multiplicative)",
                                     median(cal$x_med[cal$error == "add"]), median(cal$x_med[cal$error == "mult"]))
} else "--"

tex <- c(
  "\\documentclass[10pt]{article}",
  "\\usepackage[margin=0.5in,landscape]{geometry}",
  "\\usepackage{booktabs,threeparttable,xcolor}",
  "\\usepackage[labelformat=empty,justification=raggedright,singlelinecheck=false]{caption}",
  "\\pagestyle{empty}",
  "\\begin{document}",
  "\\begin{table}[ht]\\centering\\footnotesize",
  "\\begin{threeparttable}",
  sprintf("\\caption*{\\textbf{Table 2 (%s, %s, $B = %d$ per cell).} Estimated exposure effects of cumulative total exposure %s (bold cells with correct model specification).}",
          switch(SCENARIO, null = "null effect", paper = "reproduction", sprintf("scenario %s", SCENARIO)), gsub("_", "\\\\_", tag), min(tab$B),
          if (null) "when the true effect is zero, so both risk models are correctly specified" else "with model misspecification"),
  sprintf("\\begin{tabular}{l %s r}", paste(rep("r", length(cols)), collapse = "")),
  head_rows,
  sprintf("\\multicolumn{%d}{l}{\\textit{A. Mean $\\hat\\beta$ over converged fits%s}} \\\\[2pt]", nk,
          if (!paper_vals) "" else "; the paper's values in {\\color{gray}(gray)}, converted from its relative bias"),
  rows_a,
  "\\midrule[\\heavyrulewidth]",
  sprintf("\\multicolumn{%d}{l}{\\textit{B. Median $\\hat\\beta$ over all fits; [number of fits not converged or stopped at the cap $\\beta = 50$]}} \\\\[2pt]", nk),
  rows_b,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\footnotesize",
  if (null) "\\item Cox model $\\lambda(t) = \\lambda_0(t)\\exp(\\beta x(t))$ and ERR model $\\lambda(t) = \\lambda_0(t)(1 + \\beta x(t))$, both with true $\\beta = 0$: the hazard does not depend on $x$, so Cox and ERR truth generate the same data and every fit should estimate 0. True $x$: the same risk model fitted to the true cumulative exposure."
  else sprintf("\\item Cox model $\\lambda(t) = \\lambda_0(t)\\exp(\\beta x(t))$ with true $\\beta_{\\mathrm{Cox}} = %s$; ERR model $\\lambda(t) = \\lambda_0(t)(1 + \\beta x(t))$ with true $\\beta_{\\mathrm{ERR}} = %s$. Each fit is compared with the true $\\beta$ of the model it fits%s. True $x$: the same risk model fitted to the true cumulative exposure, separating risk-model misspecification and small-sample bias from measurement error.",
               bfmt(BETA[["cox"]]), bfmt(BETA[["err"]]),
               if (paper_vals) "; the paper's value is $\\beta(1 + \\mathrm{bias}/100)$" else "; the paper's values (for $\\beta = 1$ and 1.109) are not shown"),
  sprintf("\\item $n_1 = %d$, $n_2 = %d$ (at most 8 pairs each), $T^* = %d$, cumulative incidence about 10\\%%. Surrogate $C$ has compound symmetry ($\\rho = 0.5$) across years, $E(C) = 0.012$; true $c = \\gamma_1 C + \\epsilon$ with $\\gamma_1 = 2/3$ (additive) or $\\log c = \\gamma_0 + \\log C + \\epsilon$ with $E(c) = 0.008$ (multiplicative), $\\sigma_\\epsilon$ set by $\\rho(c, C)$; median $x$ %s.", N_MAIN, N_VAL_SUBJ, T_STAR, x_med),
  "\\item Pairs dropped: share of validation pairs with $c \\le 0$ or $C \\le 0$, which the log (.M) calibration cannot use. -- : no fit converged, or not recorded. A median of 50.00 means at least half of the fits stopped at the cap $\\beta = 50$.",
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}",
  "\\end{document}")

out_dir <- file.path(PROJ, "results")
tex_file <- file.path(out_dir, sprintf("table2_%s.tex", tag))
writeLines(tex, tex_file)
build <- file.path(tempdir(), "latex"); dir.create(build, showWarnings = FALSE)
status <- system2("pdflatex", c("-interaction=nonstopmode", "-halt-on-error",
                                paste0("-output-directory=", build), tex_file), stdout = FALSE)
if (status != 0) stop("pdflatex failed; see log in ", build)
file.copy(file.path(build, sprintf("table2_%s.pdf", tag)), out_dir, overwrite = TRUE)
