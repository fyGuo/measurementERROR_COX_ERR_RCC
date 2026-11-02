#### Component checks: run from sim1_repro/ with Rscript tests/test_components.R ####

source("R/load.R")
library(survival)

ok <- function(name, cond) cat(sprintf("[%s] %s\n", if (cond) "PASS" else "FAIL", name))

set.seed(1)
cell <- data.frame(model = "cox", error = "add", cor = "ar1", cell = 1, nu = 0.012)
L <- cum_lag_matrix(T_STAR)
ep <- expo_params("add", 0.1)
ex <- gen_exposure(N_MAIN, T_STAR, "add", COR$ar1, ep)
x <- ex$c %*% L; X <- ex$C %*% L

#### 1. Exposure moments ####
ok("additive C mean / sd", abs(mean(ex$C) - ep$mu) < 0.01 && abs(sd(ex$C) - ep$sd) < 0.01)
ok("AR(1) lag-1 correlation", abs(cor(ex$C[, 10], ex$C[, 11]) - 0.62) < 0.05)
ok("error model slope recovered", abs(coef(lm(as.vector(ex$c) ~ as.vector(ex$C)))[2] - ep$g1) < 0.01)

#### 2. Cox engine equals coxph (Breslow) ####
sv <- gen_survival(x, "cox", BETA[["cox"]], 0.012, THETA)
ev <- make_events(sv$time, sv$D, sv$k)
cat("events:", ev$n, "\n")
f <- pl_fit(z_long(x, ev), ev$eid, z_case(x, ev), "cox")
# counting-process data: subject i, interval (k, min(k + 1, time_i)] carries x_k
nk <- ceiling(sv$time)
long <- data.frame(id = rep(seq_len(N_MAIN), nk), k = sequence(nk) - 1)
long$t_in <- long$k
long$t_out <- pmin(long$k + 1, sv$time[long$id])
long$d <- as.integer(sv$D[long$id] == 1 & long$t_out == sv$time[long$id])
long$x <- x[cbind(long$id, long$k + 1)]
cf <- coxph(Surv(t_in, t_out, d) ~ x, data = long, ties = "breslow")
ok("Cox beta matches coxph", abs(f$beta - coef(cf)) < 1e-6)
ok("Cox SE matches coxph", abs(f$se - sqrt(vcov(cf)[1])) < 1e-6)

#### 3. ERR engine reaches the maximum of the partial likelihood ####
sv2 <- gen_survival(x, "err", BETA[["err"]], 0.014, THETA)
ev2 <- make_events(sv2$time, sv2$D, sv2$k)
zz <- z_long(x, ev2); zc <- z_case(x, ev2)
fe <- pl_fit(zz, ev2$eid, zc, "err")
ll <- function(b) pl_parts(b, zz, ev2$eid, zc, "err")$ll
hi <- if (min(zz) < 0) -1 / min(zz) else 20
op <- optimize(ll, c(-1 / max(zz) + 1e-6, hi - 1e-6), maximum = TRUE, tol = 1e-10)
ok("ERR beta matches direct maximisation", abs(fe$beta - op$maximum) < 1e-4)

#### 4. RRC with one forced risk set equals ordinary regression calibration ####
pairs <- sample_validation(sv$time, N_VAL_SUBJ, N_VAL_PAIR, MAX_PER_ID, T_STAR)
ok("validation: 3000 pairs, <= 8 per subject",
   nrow(pairs) == N_VAL_PAIR && max(table(pairs$id)) <= MAX_PER_ID)
r1 <- rrc_fit(pairs, ex$C, ex$c, sv$time, ev, "add", "cox", L, min_pairs = 1e9)
a <- coef(lm(ex$c[cbind(pairs$id, pairs$l)] ~ ex$C[cbind(pairs$id, pairs$l)]))
V <- matrix(FALSE, N_MAIN, T_STAR); V[cbind(pairs$id, pairs$l)] <- TRUE
xhat <- (ex$c * V + (a[1] + a[2] * ex$C) * !V) %*% L
f_rc <- pl_fit(z_long(xhat, ev), ev$eid, z_case(xhat, ev), "cox")
ok("one-group RRC equals RC", r1$n_groups == 1 && abs(r1$beta - f_rc$beta) < 1e-6)

#### 5. Sandwich: with alpha uncertainty the SE grows ####
r2 <- rrc_fit(pairs, ex$C, ex$c, sv$time, ev, "add", "cox", L, MIN_RS_PAIRS)
cat(sprintf("RRC: groups %d, beta %.3f, se %.3f, se(beta only) %.3f, se(model) %.3f\n",
            r2$n_groups, r2$beta, r2$se, r2$se_beta_only, r2$se_model))
ok("sandwich SE > SE ignoring alpha", r2$se > r2$se_beta_only)
