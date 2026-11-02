#### Run Simulation 2 replicates in parallel (future + furrr) ####
# run from sim2_repro/:  Rscript scripts/run.R <tag> <B> [cells]
# saves results/<tag>/cell<k>.rds (one row per replicate); finished cells are skipped

source("R/load.R")
library(future)
library(furrr)

args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args) >= 1) args[1] else "pilot"
B <- if (length(args) >= 2) as.integer(args[2]) else 48
cells <- if (length(args) >= 3 && args[3] != "all") as.integer(strsplit(args[3], ",")[[1]]) else CELLS$cell

plan(multisession, workers = max(1, availableCores() - 1))
# one_rep() seeds itself from (cell, rep), so future must not set seeds
opts <- furrr_options(seed = NULL)

out_dir <- file.path(PROJ, "results", tag)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (k in cells) {
  f <- file.path(out_dir, sprintf("cell%02d.rds", k))
  if (file.exists(f)) { cat("skip cell", k, "\n"); next }
  cl <- CELLS[CELLS$cell == k, ]
  t0 <- Sys.time()
  rows <- future_map(1:B, function(r) {
    tryCatch(one_rep(cl, r), error = function(e) data.frame(cell = k, rep = r, error = conditionMessage(e)))
  }, .options = opts)
  cols <- unique(unlist(lapply(rows, names)))
  res <- do.call(rbind, lapply(rows, function(x) { for (v in setdiff(cols, names(x))) x[[v]] <- NA; x[cols] }))
  saveRDS(res, f)
  cat(sprintf("cell %2d (%s/%s/rho %.1f): %d reps in %.1f min\n", k, cl$model, cl$error, cl$rho,
              B, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

plan(sequential)
