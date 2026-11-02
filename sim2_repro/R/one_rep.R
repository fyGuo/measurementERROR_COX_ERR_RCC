#### One Simulation 2 replicate: naive and RRC fits under both risk and error models ####

one_rep <- function(cell, rep, base_seed = 20261102) {
  seed <- base_seed + 100000 * cell$cell + rep
  set.seed(seed)
  L <- cum_lag_matrix(T_STAR)

  ex <- gen_exposure2(N_MAIN, T_STAR, cell$error, expo_params(cell$error, cell$rho))
  x <- ex$c %*% L
  X <- ex$C %*% L
  sv <- gen_survival(x, cell$model, BETA[[cell$model]], cell$nu, THETA)
  ev <- make_events(sv$time, sv$D, sv$k)
  pairs <- sample_validation2(sv$time, N_VAL_SUBJ, MAX_PER_ID, T_STAR)

  # the log calibration needs positive values: under additive data it uses only pairs with
  # c, C > 0, and a non-positive C is floored so its calibrated c is about 0
  pos <- ex$c[cbind(pairs$id, pairs$l)] > 0 & ex$C[cbind(pairs$id, pairs$l)] > 0
  C_pos <- pmax(ex$C, 1e-8)

  out <- data.frame(cell = cell$cell, rep = rep, seed = seed, n_events = ev$n, n_pairs = nrow(pairs),
                    n_nonpos = sum(ex$C <= 0) + sum(ex$c <= 0), n_pairs_pos = sum(pos))
  for (m in c("cox", "err")) {
    f <- pl_fit(z_long(x, ev), ev$eid, z_case(x, ev), m, start = BETA[[m]])
    out[[paste0("b_true_", m)]] <- f$beta; out[[paste0("conv_true_", m)]] <- f$converged
  }
  for (k in seq_len(nrow(FITS))) {
    fm <- FITS[k, ]
    f <- if (is.na(fm$cal)) {
      pl_fit(z_long(X, ev), ev$eid, z_case(X, ev), fm$model, start = BETA[[fm$model]] / 2)
    } else if (fm$cal == "add") {
      rrc_fit(pairs, ex$C, ex$c, sv$time, ev, "add", fm$model, L, MIN_RS_PAIRS, start = BETA[[fm$model]])
    } else {
      rrc_fit(pairs[pos, ], C_pos, ex$c, sv$time, ev, "mult", fm$model, L, MIN_RS_PAIRS, start = BETA[[fm$model]])
    }
    out[[paste0("b_", fm$name)]] <- f$beta
    out[[paste0("se_", fm$name)]] <- f$se
    out[[paste0("conv_", fm$name)]] <- f$converged
  }
  out
}
