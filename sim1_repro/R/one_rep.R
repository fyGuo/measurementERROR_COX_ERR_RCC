#### One simulation replicate ####

# cell: one row of CELLS (model, error, cor, cell, nu)
# returns one row of estimates; the dataset itself is not kept and can be
# regenerated exactly from (cell, rep)
one_rep <- function(cell, rep, base_seed = 20261028, n1 = N_MAIN, n2 = N_VAL_SUBJ, n_pair = N_VAL_PAIR) {
  seed <- base_seed + 100000 * cell$cell + rep
  set.seed(seed)
  beta <- BETA[[cell$model]]
  L <- cum_lag_matrix(T_STAR)

  ex <- gen_exposure(n1, T_STAR, cell$error, COR[[cell$cor]], expo_params(cell$error, cell$sd, cell$rho))
  x <- ex$c %*% L                         # true cumulative exposure, column m = x_{m-1}
  X <- ex$C %*% L                         # surrogate cumulative exposure
  sv <- gen_survival(x, cell$model, beta, cell$nu, THETA)
  ev <- make_events(sv$time, sv$D, sv$k)
  pairs <- sample_validation(sv$time, n2, n_pair, MAX_PER_ID, T_STAR)

  f_true  <- pl_fit(z_long(x, ev), ev$eid, z_case(x, ev), cell$model, start = beta)
  f_naive <- pl_fit(z_long(X, ev), ev$eid, z_case(X, ev), cell$model, start = beta / 2)
  f_rrc   <- rrc_fit(pairs, ex$C, ex$c, sv$time, ev, cell$error, cell$model, L,
                     MIN_RS_PAIRS, start = beta)

  data.frame(cell = cell$cell, rep = rep, seed = seed, beta = beta, n1 = n1, n2 = n2,
             n_events = ev$n, n_viol = sv$n_viol, n_groups = f_rrc$n_groups,
             b_true = f_true$beta, se_true = f_true$se,
             b_naive = f_naive$beta, se_naive = f_naive$se,
             b_rrc = f_rrc$beta, se_rrc = f_rrc$se,
             se_rrc_beta_only = f_rrc$se_beta_only, se_rrc_model = f_rrc$se_model,
             conv = f_true$converged & f_naive$converged & f_rrc$converged)
}
