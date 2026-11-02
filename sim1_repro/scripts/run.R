#### Run replicates for all Table 1 cells in parallel (future + furrr) ####
# run from sim1_repro/:  Rscript scripts/run.R <tag> <B> [cells] [n1] [n2]
#   e.g. Rscript scripts/run.R pilot 200      Rscript scripts/run.R full 2000
#        Rscript scripts/run.R full_v5 2000 all 6000 2000   (validated pairs scale with n2)
# saves results/<tag>/cell<k>.rds (one row per replicate); finished cells are skipped

source("R/load.R")
library(future)
library(furrr)

args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args) >= 1) args[1] else "pilot"
B <- if (length(args) >= 2) as.integer(args[2]) else 200
cells <- if (length(args) >= 3 && args[3] != "all") as.integer(strsplit(args[3], ",")[[1]]) else CELLS$cell
n1 <- if (length(args) >= 4) as.integer(args[4]) else N_MAIN
n2 <- if (length(args) >= 5) as.integer(args[5]) else N_VAL_SUBJ
n_pair <- round(N_VAL_PAIR * n2 / N_VAL_SUBJ)
n_workers <- max(1, availableCores() - 1)

plan(multisession, workers = n_workers)
# one_rep() seeds itself from (cell, rep), so future must not set or check seeds:
# seed = NULL keeps the default RNG and the same results as a sequential run
opts <- furrr_options(seed = NULL)

out_dir <- file.path(PROJ, "results", tag)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (k in cells) {
  f <- file.path(out_dir, sprintf("cell%02d.rds", k))
  if (file.exists(f)) { cat("skip cell", k, "\n"); next }
  cl <- CELLS[CELLS$cell == k, ]
  t0 <- Sys.time()
  rows <- future_map(1:B, function(r) {
    tryCatch(one_rep(cl, r, n1 = n1, n2 = n2, n_pair = n_pair), error = function(e) data.frame(cell = k, rep = r, error = conditionMessage(e)))
  }, .options = opts)
  cols <- unique(unlist(lapply(rows, names)))
  res <- do.call(rbind, lapply(rows, function(x) { for (v in setdiff(cols, names(x))) x[[v]] <- NA; x[cols] }))
  saveRDS(res, f)
  cat(sprintf("cell %2d (%s/%s/%s): %d reps in %.1f min\n", k, cl$model, cl$error, cl$cor,
              B, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

plan(sequential)
