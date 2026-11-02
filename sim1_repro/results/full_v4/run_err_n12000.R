setwd("/Users/fuyuguo/yueyang_measurement_error/sim1_repro"); source("R/load.R")
library(future); library(furrr)
SP <- "/private/tmp/claude-501/-Users-fuyuguo-yueyang-measurement-error-sim1-repro/5d137b51-f4b9-4f36-9e01-10438228c5c3/scratchpad"
L <- cum_lag_matrix(T_STAR); n <- 12000; n2 <- 4000; np <- 15000; B <- 96; W <- availableCores()
one <- function(cl, r) {
  set.seed(888000 + 1000 * cl$cell + r); beta <- BETA[[cl$model]]; md <- cl$model
  ex <- gen_exposure(n, T_STAR, cl$error, COR[[cl$cor]], expo_params(cl$error, cl$sd, cl$rho))
  x <- ex$c %*% L
  sv <- gen_survival(x, "err", beta, cl$nu, THETA); ev <- make_events(sv$time, sv$D, sv$k)
  pairs <- sample_validation(sv$time, n2, np, MAX_PER_ID, T_STAR)
  ft <- pl_fit(z_long(x, ev), ev$eid, z_case(x, ev), md, start = beta)
  X <- ex$C %*% L
  fn <- pl_fit(z_long(X, ev), ev$eid, z_case(X, ev), md, start = beta / 2)
  fr <- rrc_fit(pairs, ex$C, ex$c, sv$time, ev, cl$error, md, L, MIN_RS_PAIRS, start = beta)
  out <- data.frame(cell = cl$cell, r = r, beta = beta, n_ev = ev$n, g = fr$n_groups, b_true = ft$beta, b_naive = fn$beta, se_naive = fn$se, b_rrc = fr$beta,
                    se_rrc = fr$se, se_true = ft$se, ct = ft$converged, cn = fn$converged, cr = fr$converged)
  rm(ex, x, X, sv, ev, pairs, ft, fn, fr); gc(FALSE); out
}
plan(multisession, workers = W); f <- file.path(SP, "rrc_v3.rds"); res <- if (file.exists(f)) readRDS(f) else NULL; t0 <- Sys.time()
for (k in c(10, 7, 8, 9, 11, 12)) {
  cl <- CELLS[CELLS$cell == k, ]
  for (st in seq(1, B, by = W)) {
    if (!is.null(res) && sum(res$cell == k) >= min(st + W - 1, B)) next
    rs <- st:min(st + W - 1, B)
    res <- rbind(res, do.call(rbind, future_map(rs, function(r) tryCatch(one(cl, r), error = function(e) NULL), .options = furrr_options(seed = NULL))))
    saveRDS(res, file.path(SP, "rrc_v3.rds"))
    cat(k, max(rs), round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
  }
}
plan(sequential); cat("DONE\n")
