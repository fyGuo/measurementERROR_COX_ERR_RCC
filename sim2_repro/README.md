# Simulation 2 reproduction (paper Section 4.2, Table 2)

Simulation 2 of Yang and others: naive and regression-calibration (RRC) estimates of the effect of cumulative
exposure when the risk model (Cox or excess relative risk, ERR) and the measurement error model (additive or
multiplicative) may be misspecified. The estimators are Simulation 1's, reused unchanged from `../sim1_repro/R/`.

**Status (2026-11-02).** The design is settled. The chosen scenario is `SIM2_SCENARIO=beta1` (true beta = 1
for both Cox and ERR) at **n1 = 6000, n2 = 2000, 10% incidence**. No full run has been done at these sizes;
all pilots used n1 = 2000, n2 = 500. A 22-replicate timing run at the new sizes is in `results/beta1_timing`
(see Part A.4). The full run is meant to be done on a faster Mac.

The README has five parts: A how to run (including on another machine), B the current design, C the
scenarios, D the files, E the history of what was tried and decided, with the result folder of each step.

---

## Part A — How to run

### A.1 Setup on a new machine

1. Clone the repository: `git clone https://github.com/fyGuo/measurementERROR_COX_ERR_RCC.git`.
   `sim2_repro/` needs `sim1_repro/` next to it (`R/load.R` sources `../sim1_repro/R/`), so keep the
   repository layout.
2. R 4.5 (developed on R 4.5.3) with `future` (1.70.0) and `furrr` (0.4.0):
   `install.packages(c("future", "furrr"))`. Nothing else is needed for the simulation.
3. For the PDF tables only: `pdflatex` with `booktabs`, `threeparttable`, `xcolor` and `caption`
   (all in a standard TeX Live or MacTeX install).
4. Run every command from `sim2_repro/`.

### A.2 Commands

Every script reads the scenario from the environment variable `SIM2_SCENARIO` (default `paper`), so prefix
each command with it. For the chosen setting:

```sh
cd sim2_repro
export SIM2_SCENARIO=beta1
Rscript scripts/01_calibrate.R              # optional: results/calibration_beta1.csv is already in the repo
Rscript scripts/run.R beta1_full 2000       # the full run; results/beta1_full/cellNN.rds
Rscript scripts/table2.R beta1_full         # results/table2_beta1_full.csv
Rscript scripts/table2_pdf.R beta1_full     # results/table2_beta1_full.pdf (needs pdflatex)
```

- `run.R <tag> <B> [cells]`: runs B replicates per cell in parallel and saves one file per cell,
  `results/<tag>/cellNN.rds`, when that cell finishes. A cell whose file exists is skipped, so an interrupted
  run resumes where it stopped (the cell in progress restarts). A subset of cells: `run.R beta1_full 2000 7,8,9`.
- Parallelism: `future::plan(multisession)` with all cores but one (`availableCores() - 1`).
- Reproducibility: each replicate seeds itself with `20261102 + 100000 * cell + rep` (`R/one_rep.R`), so
  results do not depend on the number of workers or the machine, given the same R version.
- Long runs: run in the background and keep the Mac awake, e.g.
  `caffeinate -i nohup Rscript scripts/run.R beta1_full 2000 > run_beta1_full.log 2>&1 &`.
  `run.R` prints one line per finished cell with its time.

### A.3 Bringing results back

After the run, on the fast Mac: `git add sim2_repro/results/beta1_full sim2_repro/results/table2_beta1_full.*`,
commit and push; then `git pull` here. The `.rds` files are tracked on purpose (they are the results), and the
whole repository is about 40 MB.

### A.4 Run time

Cost grows with n1 times the number of events (each fit loops over every risk-set member at every event).

| Setting | Time per cell (this Mac, 11 workers) | 2000 replicates x 12 cells |
|---|---|---|
| n1 = 2000, n2 = 500 (pilots) | about 0.3 min per 48 replicates | about 2.5 hours |
| n1 = 6000, n2 = 2000 (current) | see `results/beta1_timing` (22 replicates) and the note below | see below |

Timing run at the current sizes (`SIM2_SCENARIO=beta1 Rscript scripts/run.R beta1_timing 22`): pending; the
per-cell times will be added here.

---

## Part B — Current design

All settings are in `R/params.R`; generation is in `R/exposure2.R`; one replicate is `R/one_rep.R`.

### B.1 Sizes and survival

| Quantity | Value | Paper |
|---|---|---|
| Main study n1 | 6000 | 2000 |
| Validation subjects n2 | 2000, at most 8 (c, C) pairs each from their follow-up years | 500 |
| Follow-up T* | 60 years | 60 |
| Replicates B | 2000 | 2000 |
| Cumulative incidence | 10% | about 10% |
| Baseline hazard | Weibull, Lambda0(t) = (nu t)^theta, theta = 4; nu calibrated per cell | not given |
| Validation risk sets | merged until each has at least 20 pairs, as in Simulation 1 | |

The sizes differ from the paper because at n1 = 2000 (about 200 events) the ERR fit does not settle even on
the true x (Part E.8).

### B.2 Exposure: surrogate C first, then true c (the paper's error model)

Annual exposures for subject i, year j = 1..60. Cumulative exposure x(t) sums the true c over the years
before t (one-year lag, `cum_lag_matrix`); X(t) does the same with C.

- C has compound symmetry across years with correlation 0.5 (paper): a shared subject effect plus an
  independent yearly effect.
- **Additive:** C ~ MVN(0.012, 0.006^2 R), c = g0 + g1 C + e, with g0 = 0, g1 = 2/3.
- **Multiplicative:** log C ~ MVN(mu, 0.8^2 R) with E(C) = 0.012, log c = g0 + log C + e (g1 = 1), with g0
  set so that E(c) = 0.008.
- e ~ N(0, sd_e^2), independent across subjects and years. sd_e is solved in closed form from the paper's
  rho(c, C) in {0.3, 0.6, 0.9}, the correlation of c and C within a year, on the original scale for both error
  models (`expo_params`).
- Held fixed across rho: the distribution of C and E(c | C) = (2/3) C. So the naive estimate is attenuated by
  about 2/3 at every rho, as the paper's flat naive column implies; rho only changes the residual variance.
- Exposure level: median x about 0.2 (paper 0.204). The paper's sd of 15 cannot be used with beta = 1.
- Consequences: additive C is negative in 2% of years and c in 4-27% (more at rho = 0.3); x has mean 0.23,
  and x at year 45 has mean 0.36 and SD 0.13-0.16 (additive, symmetric) or 0.22-0.28 (multiplicative,
  skewness 2-5). Full summary: `results/x_summary.csv` (`scripts/03_x_summary.R`).

### B.3 Fits per replicate

On the main study (`R/one_rep.R`), each with Simulation 1's partial-likelihood code (`pl_fit`, `rrc_fit`):

| Column | Fit |
|---|---|
| True x: Cox, ERR | each risk model fitted to the true x (reference: what each model estimates without measurement error) |
| Naive: Cox, ERR | fitted to the surrogate X |
| Cox.A, ERR.A | RRC with the additive calibration c = a0 + a1 C |
| Cox.M, ERR.M | RRC with the multiplicative calibration log c = a0 + a1 log C, E(c given C) = exp(a0 + s2 / 2) C^a1 |

- ERR is restricted to 1 + beta x > 0 and capped at beta = 50; a fit at a bound counts as not converged.
- The log calibration can only use pairs with c, C > 0; under additive error the other pairs are dropped
  (`n_pairs_pos`, reported as "Pairs dropped") and C <= 0 is floored at 1e-8 for prediction.

### B.4 The table

`scripts/table2.R` writes, per cell and fit, the mean beta-hat over converged fits (`<fit>_b`), the median over
all fits (`<fit>_bmed`), the number not converged (`<fit>_nc`) and, for the paper scenario, relative bias.
`scripts/table2_pdf.R` lays this out like the paper's Table 2: panel A mean beta-hat (with the paper's values,
converted to beta, in gray for the `paper` scenario only), panel B median beta-hat with non-convergence counts
and the share of pairs dropped. Bold marks the RRC fit with the correct risk and error model. Results are
reported as estimated beta, not relative bias.

---

## Part C — Scenarios (`SIM2_SCENARIO`)

| Value | True beta (Cox, ERR) | Cells | Calibration file | Note |
|---|---|---|---|---|
| `paper` (default) | 1, 1.109 | 12 | `calibration.csv` | the paper's pairing of the two effects at x = 0.204 |
| `beta1` (**chosen**) | 1, 1 | 12 | `calibration_beta1.csv` | cells 1-6 identical to `paper` |
| `sim1` | log(1.5), 0.423 | 12 | `calibration_sim1.csv` | Simulation 1's effects |
| `null` | 0, 0 | 1-6 | `calibration_null.csv` | hazard ignores x, so Cox and ERR truth coincide |

Cells, in the paper's order: 1-3 Cox + additive, 4-6 Cox + multiplicative, 7-9 ERR + additive,
10-12 ERR + multiplicative, each at rho(c, C) = 0.3, 0.6, 0.9. The calibration (nu) depends on the scenario's
beta and on the exposure design, not on n1, so the stored files are valid at any sample size.

---

## Part D — Files

| File | Purpose |
|---|---|
| `R/load.R` | sources Simulation 1's `exposure.R` (`cum_lag_matrix`), `survival.R`, `fit.R`, `rrc.R`, then this folder's files |
| `R/params.R` | sizes, scenario and beta, exposure settings (`EXPO`, `expo_params`), the 12 cells, the 6 fits, per-cell nu |
| `R/exposure2.R` | `gen_exposure2` (C first, then c) and `sample_validation2` |
| `R/one_rep.R` | one replicate: generation, survival, all fits |
| `scripts/01_calibrate.R` | nu per cell for 10% incidence, plus checks of corr(c, C), slope of E(x given X), level of x |
| `scripts/run.R` | parallel replicates -> `results/<tag>/cellNN.rds` |
| `scripts/table2.R` | summary per cell and fit -> `results/table2_<tag>.csv` |
| `scripts/table2_pdf.R` | Table 2 layout -> `results/table2_<tag>.pdf` and `.tex` |
| `scripts/02_explore_truex.R` | exploration: Cox and ERR on the true x over a grid of exposure settings (Part E.6) |
| `scripts/03_x_summary.R` | distribution of x per cell and scenario -> `results/x_summary.csv` (Part E.9) |

Result folders (all pilots at n1 = 2000, n2 = 500 unless noted):

| Folder / table | Scenario and design | Replicates |
|---|---|---|
| `pilot_v1`, `table2_pilot_v1.*` | paper; c generated first, then C (rejected, E.1) | 48 |
| `pilot_v2`, `table2_pilot_v2.*` | paper; C first, additive sd(C) = 0.012 (E.2) | 48 |
| `pilot`, `table2_pilot.*` | paper; current design (additive cells rerun with sd(C) = 0.006, multiplicative copied from `pilot_v2`) | 48 |
| `pilot_v4`, `table2_pilot_v4.*` | paper; additive g0 = 0.024 (rejected, E.4) | 48 |
| `null_pilot`, `table2_null_pilot.*` | null | 48 |
| `sim1_pilot`, `table2_sim1_pilot.*` | sim1 | 48 |
| `beta1_pilot`, `table2_beta1_pilot.*` | beta1 (cells 1-6 copied from `pilot`) | 48 |
| `beta1_timing`, `table2_beta1_timing.*` | beta1 at **n1 = 6000, n2 = 2000** (timing) | 22 |
| `explore_truex_*.csv` | true-x fits over exposure grids (E.6) | 200 per setting |
| `calibration_v1.csv`, `calibration_v2.csv` | nu for `pilot_v1`, `pilot_v2` | |

Pilots before the "Pairs dropped" column was added do not record it (shown as --).

---

## Part E — History (2026-11-02)

**E.1 First design: true c first, then C = c + u (classical error).** 48-replicate pilot (`pilot_v1`). RRC
overcorrected even when both models were right (Cox.A on Cox + additive: +394, +146, +23% at rho 0.3, 0.6,
0.9). Reason: c had compound symmetry, so summing years averages out the independent error but not each
subject's level; per-year calibration E(c_j given C_j) has slope rho^2 and overcorrects the sum. Also ERR was
not identified: x was about 2.5 at death ages, 1 + beta x is then about beta x, and most ERR fits stopped at
the cap. The paper's naive bias is nearly flat in rho (about -33%), which classical error cannot give but a
model c = g0 + g1 C + e with fixed g1 does. **Decision:** generate C first and c from C, the paper's error
model, as in Simulation 1, and lower the exposure level to median x about 0.2.

**E.2 C first, then c (`pilot_v2`, additive sd(C) = 0.012).** Cox results matched the paper: naive Cox about
-30 to -37% (paper -33%), Cox.A on Cox + additive -3.5, 4.5, 3.7% (paper about -4%), Cox.M on
Cox + multiplicative 6.7, -4.9, 1.8%. Under additive error about half of all ERR fits, including fits on the
true x, stopped at the bound 1 + beta x > 0, because normal C with sd equal to its mean gives negative x.
**Decision:** additive sd(C) = 0.006 (now 2% negative C; boundary hits at most 3 of 48).

**E.3 Why ERR.M is unstable under additive error.** The log calibration drops pairs with c or C <= 0 (up to
27% at rho 0.3), which flattens the fitted slope (a1 = 0.28 instead of 1 at rho 0.3) and the calibrated x
(spread over level 0.15 instead of 0.44). Cox is unaffected by a common shift in x, but ERR depends on the
level of x, so its likelihood goes flat (SD of beta-hat up to 20). ERR.M under multiplicative error is fine.

**E.4 Intercept test: additive g0 = 0.024 (`pilot_v4`).** Lifting c so that only 0.8% is negative fixed Cox.M
(-22.8, -5.2, 1.2% instead of +18, +24, +11%), but ERR lost identification even on the true x (SD about 20,
9-38 of 48 fits at the cap): a common shift raises the level of x relative to its spread. With additive normal
error at rho = 0.3, c can be mostly positive only if E(c) is about 2.4 times its yearly SD, which caps the
spread of x at about 10% of its level. **Decision:** g0 = 0; report ERR.M under additive error as error-model
misspecification, with the share of pairs dropped.

**E.5 Why ERR fits under Cox truth are far from 1.109.** 1.109 = (e^0.204 - 1) / 0.204 matches exp(beta x) and
1 + beta x only at x = 0.204. Within a risk set the ERR log-rate slope is beta / (1 + beta xbar), so the ERR
fit tends to about 1 / (1 - xbar) under Cox truth and the Cox fit to about 1 / (1 + xbar) under ERR truth;
xbar, the mean x of those at risk at event times, is about 0.38. The multiplicative right tail makes the gap
larger. RRC recovers the fit on the true x, so the remaining gap is risk-model misspecification. The naive ERR
fit is pulled down by about 2/3, which partly cancels it. Under ERR + multiplicative the naive ERR fit is about
2/3 of the true-x fit at every rho (0.71, 0.70, 0.66); the trend across rho comes from the true-x fit and
Monte Carlo noise, not from measurement error.

**E.6 Exploration: can the Cox and ERR fits on the true x both be near 1?** (`scripts/02_explore_truex.R`,
`results/explore_truex_*.csv`, both true beta = 1.) Both near 1 needs xbar near 0. At about 200 events no
setting works (where xbar is small, the SD of beta-hat is 0.7-2). Additive with E(C) = 0, sd(C) = 0.005,
compound symmetry 0.9 and 30% incidence gives Cox 1.00 / ERR 0.96 (Cox truth) and 1.01 / 0.99 (ERR truth),
SD about 0.25, but half of all exposures are negative. Multiplicative cannot get there (about 25,000 events
needed). **Decision:** not adopted; keep the paper's design.

**E.7 Other effect sizes.** Only beta x matters: rescaling x by 1/k and beta by k gives the same data, and
exp(beta x) and 1 + beta x differ at second order, about (beta x)^2 / 2.
- `null` (beta = 0, `null_pilot`): every fit, naive included, is near 0 (means -0.11 to 0.06, except additive
  ERR means 0.15-0.35 from right skew). Non-differential error cannot create an effect, so this scenario only
  tests that the methods add no bias and, if wanted, their type I error. The 2000-replicate run was stopped.
- `sim1` (beta = log 1.5 and 0.423, `sim1_pilot`): beta xbar falls from about 0.38 to 0.15, so the Cox fits
  agree closely, but ERR's SD (0.4-0.6) is now as large as beta and its mean runs high even under ERR truth.
- `beta1` (beta = 1 for both, `beta1_pilot`): **chosen.** Under ERR truth the true-x ERR fit centres near 1
  (medians 0.90-1.20; means 1.09-1.53 at n1 = 2000) and ERR.A / ERR.M recover it; the Cox fit gives 0.55-0.80
  (misspecification, about 1 / (1 + 0.38) = 0.72). Under Cox truth ERR gives about 1.6-3.0.

**E.8 Sample size.** Under ERR truth with beta = 1, the ERR fit on the true x is consistent but needs events:

| Cell | Incidence | n1 | Events | Mean | Median | SD |
|---|---|---|---|---|---|---|
| ERR + additive, rho 0.6 | 10% | 2000 | 200 | 1.22 | 0.95 | 1.08 |
| | 10% | 5000 | 490 | 1.03 | 1.01 | 0.65 |
| | 30% | 5000 | 1490 | 0.93 | 0.95 | 0.29 |
| ERR + multiplicative, rho 0.6 | 10% | 2000 | 200 | 1.14 | 1.05 | 0.63 |
| | 10% | 5000 | 500 | 1.02 | 1.03 | 0.43 |
| | 10% | 10000 | 1000 | 0.98 | 0.94 | 0.28 |

(44 replicates each.) With beta = 1.109 the same check gave mean 1.13 and median 1.11 at n1 = 20000 (cell 11).
**Decision:** n1 = 6000, n2 = 2000, incidence 10%, in every scenario (about 600 events).

**E.9 Distribution of x** (`results/x_summary.csv`). x is the same in all scenarios (only beta changes); the
event times and the x of cases depend on beta. Person-years: mean 0.23, median 0.20 (additive) or 0.16-0.17
(multiplicative), 95th percentile 0.55-0.70. Cases at event time: median about 0.40, 95th percentile 0.72
(additive) and up to 1.33 (multiplicative, rho 0.3, paper scenario). xbar at events is about 0.38 everywhere,
so beta xbar is about 0.38-0.43 (paper, beta1), 0.15 (sim1) and 0 (null).
