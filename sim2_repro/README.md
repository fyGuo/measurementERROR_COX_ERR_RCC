# Simulation 2 reproduction (paper Section 4.2, Table 2)

Status (2026-11-02): settings changed to n1 = 6000 and n2 = 2000 (paper: 2000 and 500), incidence 10%,
in every scenario, so that the ERR fit has about 600 events. Chosen scenario: `SIM2_SCENARIO=beta1` (Cox and
ERR true beta = 1). No simulation has been run at the new sizes; all pilots below used n1 = 2000, n2 = 500.
Calibrations (nu) do not depend on n1 and stay valid. Earlier pilots: `pilot_v1` (c first), `pilot_v2`
(additive sd(C) = 0.012), `pilot_v4` (additive g0 = 0.024); `pilot` is the paper scenario.

## Layout

- `R/load.R` sources Simulation 1's estimators unchanged from `../sim1_repro/R/`
  (`exposure.R` for `cum_lag_matrix`, `survival.R`, `fit.R`, `rrc.R`), then this folder's files.
- `R/params.R`: sizes, true effects, exposure settings, the 12 cells and the 6 fits.
- `R/exposure2.R`: data generation (true c first, then C) and validation sampling.
- `R/one_rep.R`: one replicate: true-x, naive and four RRC fits (Cox/ERR x additive/multiplicative calibration).
- `scripts/01_calibrate.R`: Weibull scale nu per cell for 10% cumulative incidence, plus checks of corr(c, C) and x -> `results/calibration.csv`.
- `scripts/run.R <tag> <B> [cells]`: parallel run (future + furrr) -> `results/<tag>/cell<k>.rds`.
- `scripts/table2.R <tag>`: relative bias (mean and median) per fit and cell -> `results/table2_<tag>.csv`.
- `scripts/02_explore_truex.R [n1] [cir] [reps] [base|focus]`: Cox and ERR fitted to the true x over a grid of
  exposure settings (both true beta = 1) -> `results/explore_truex_*.csv`. Exploration only.
- `scripts/03_x_summary.R`: distribution of the true x (person-years, year 45, cases at event time, risk-set
  mean at events) per cell and scenario -> `results/x_summary.csv`.
- `scripts/table2_pdf.R <tag>`: mean and median estimated beta in the paper's Table 2 layout (paper values converted to beta, in gray) -> `results/table2_<tag>.pdf`.

## Current settings

n1 = 6000, T* = 60, n2 = 2000 with up to 8 pairs each (paper and pilots: n1 = 2000, n2 = 500, about 4000 pairs), at least 20 pairs per validation
risk set, theta = 4, target incidence 10%. Cox beta = 1, ERR beta = 1.109. Surrogate C has compound symmetry
(rho 0.5) across years: additive C ~ N(0.012, 0.006^2), multiplicative log C normal with E(C) = 0.012 and
sd(log C) = 0.8. True c = g0 + g1 C + e (log c = g0 + log C + e), with E(c) = 0.008, so E(c | C) has slope 2/3
at every rho; sd(e) is solved in closed form from rho(c, C) in {0.3, 0.6, 0.9} (original scale for both error
models). Median x is about 0.205 (additive; paper 0.204) and 0.17 (multiplicative). Calibration confirms corr(c, C) and slope 0.667 in all cells.

## Pilot findings (`results/pilot`, `results/table2_pilot.pdf`, 48 reps per cell)

Earlier pilots: `pilot_v1` (c first, then C) and `pilot_v2` (C first, additive sd(C) = 0.012).
In `pilot`, additive sd(C) = 0.006: negative C 2%, negative c 4-27%; ERR boundary hits fell from about half
of the fits to at most 3 of 48. The multiplicative cells are copied from `pilot_v2` (their settings did not change).

1. Cox fits on Cox data match the paper: naive Cox about -30 to -38% (paper -33%); Cox.A on Cox+additive
   -4.6, 5.4, 6.6%; Cox.M on Cox+multiplicative 6.7, -4.9, 1.8% (paper about -5%).
2. ERR fits on Cox data are risk-model misspecification, not measurement error: the ERR fit on the true x is
   already +44 to +115% (additive) and +107 to +169% (multiplicative), and ERR.A / ERR.M land near it. 1.109
   matches exp(beta x) only at x = 0.204, and (e^x - 1) / x grows with x, so high-x subjects pull the slope up.
3. Precision is the main problem: with about 200 events and sd(x45) 0.13 (additive) or 0.23-0.27
   (multiplicative), the true-x Cox fit has SD 0.14-0.6 and the true-x ERR fit has SD 0.7-2.0 (beta = 1.1).
   ERR estimates are right-skewed, so their mean relative bias is +8 to +30% even on the true x.
4. ERR.M under additive error stays unstable (SD up to 20): the log calibration drops negative c.

Intercept test (2026-11-02, additive cells in `pilot_v4`): additive
g0 = 0.024 (E(c) = 0.032, negative c 0.8% at rho 0.3, median x 0.93). Cox.M on Cox+additive improves
(-22.8, -5.2, 1.2% from +18, +24, +11%), but ERR loses identification: even the true-x ERR fit has SD about 20
and 9-38 of 48 fits stop at the cap. A shift common to all subjects cancels in Cox risk sets but raises the
level of x relative to its spread, and ERR depends on that ratio. With additive normal error at rho = 0.3, c
can only be mostly positive if E(c) is about 2.4 times the per-year sd of c, which caps the coefficient of
variation of x near 0.1. Rejected; params.R keeps g0 = 0 (results in `pilot_v4`).

Decision: keep the paper's design (n1, incidence, median x about 0.2) and report as is. The table adds fits on
the true x (Cox and ERR) as references, medians with non-convergence counts, and the share of validation
pairs the log calibration drops (`n_pairs_pos`, recorded from the full run on).

## Exploration: can the Cox and ERR fits on the true x both be near 1? (2026-11-02, not adopted)

Within a risk set the ERR log-rate slope is beta / (1 + beta xbar), so the misspecified fit tends to
1 / (1 - xbar) (ERR under Cox truth) or 1 / (1 + xbar) (Cox under ERR truth), xbar = exposure at event times.
Both near 1 needs xbar near 0. `results/explore_truex_base_*.csv` and `explore_truex_focus_*.csv`:
- At the paper's size (about 200 events) no setting works: where xbar is small enough, the SD of beta-hat is 0.7-2.
- Additive: E(C) = 0, sd(C) = 0.005, RHO_CS = 0.9 and 30% incidence gives Cox 1.00 / ERR 0.96 (Cox truth) and
  Cox 1.01 / ERR 0.99 (ERR truth), SD about 0.25, but half of all exposures are negative.
- Multiplicative: positive, skewed x leaves the misspecified fit at 0.8-1.3 even at 50% incidence; about
  25,000 events would be needed.
The ERR fit on the true x under ERR truth is consistent (cell 11 at n1 = 20000: mean 1.13, median 1.11); at
n1 = 2000 it is right-skewed (cell 8: mean 1.52, median 1.07, SD 1.6), so medians are the better summary.
Decision: keep the paper's design.

## Null scenario: true beta = 0 (2026-11-02)

`SIM2_SCENARIO=null` sets beta = 0 in both models, so both risk models are correctly specified; the hazard
ignores x, so Cox and ERR truth give the same data and only cells 1-6 are run (calibration in
`results/calibration_null.csv`). Same exposure design as above. Prefix every script with the variable, e.g.
`SIM2_SCENARIO=null Rscript scripts/run.R null_full 2000`, then `table2.R null_full` and
`table2_pdf.R null_full` (no paper values in this table). Pilot (`null_pilot`, 48 reps): all mean and median
estimates within about 0.1 of 0 except additive ERR means of 0.15-0.35 (right skew, Monte Carlo SE about 0.1).
The 2000-replicate run `null_full` (about 75 minutes) was started and stopped before any cell finished; not needed for now.

## Simulation 1 effects: Cox beta = log(1.5), ERR beta = 0.423 (2026-11-02)

`SIM2_SCENARIO=sim1`, all 12 cells, calibration in `results/calibration_sim1.csv`, same exposure design. With
xbar about 0.36 at event times the misspecified fits on the true x should differ from the other model's
beta by about 10% (ERR about 0.47 under Cox truth, Cox about 0.37 under ERR truth) instead of about 40%. The
paper's values are not shown in this table, since its relative biases are for beta = 1 and 1.109.
Pilot: `sim1_pilot` (48 reps).

## beta = 1 in both models (2026-11-02, chosen setting)

`SIM2_SCENARIO=beta1`: the paper design (table2_pilot settings) with ERR true beta = 1 instead of 1.109.
Cells 1-6 (Cox truth) are identical to the "paper" scenario (same seeds and nu), so `beta1_pilot` reuses
`pilot` cells 1-6 and only cells 7-12 were rerun. Calibration in `results/calibration_beta1.csv`.
