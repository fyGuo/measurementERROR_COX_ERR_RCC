# Reproduction of Simulation 1 (Yang, Spiegelman, Wang — Table 1)

This folder re-implements **Simulation 1** of the paper (Section 4.1, Table 1) from the paper's text.
None of the original simulation scripts are reused (they do not match the paper; see the notes in
the parent folder's conversation history).

The README is organised as a walkthrough. Part A tells you how to run things. Part B explains the
setting **one step at a time**, in the order the data are generated and analysed. Each step says
what the paper says, what we do, the exact formula and values, where it is in the code, and how to
check it yourself. Part C lists the open decision.

---

## Part A — How to run

All commands run from this folder (`sim1_repro/`). R 4.5 with base packages, `future` + `furrr`
(parallelism: `plan(multisession)`, all cores but one) and `survival` (used only in the tests).

| Step | Command | What it does | Output |
|---|---|---|---|
| 1 | `Rscript tests/test_components.R` | 9 component checks; all should print `PASS` | screen |
| 2 | `Rscript scripts/01_calibrate.R` | per cell: exposure sd and baseline hazard scale nu (Steps 2 and 5) | `results/calibration.csv` |
| 3 | `Rscript scripts/run.R pilot 200` | 200 replicates per cell | `results/pilot/cellNN.rds` |
| 4 | `Rscript scripts/table1.R pilot` | Table 1 metrics next to the paper's values | screen + `results/table1_pilot.csv` |
| 5 | `Rscript scripts/run.R full 2000` | full run, about 1 hour | `results/full/cellNN.rds` |
| 6 | `Rscript scripts/table1.R full` | final table | `results/table1_full.csv` |

Notes
- `run.R <tag> <B> [cells]` skips cells whose `.rds` already exists. To redo a cell, delete its file.
  Run a subset with e.g. `Rscript scripts/run.R pilot 200 1,7`.
- If you change anything in `R/params.R` that affects exposure or hazard, rerun step 2 **and** delete
  old results before rerunning step 3.
- Simulated datasets are **not saved**, only one row of estimates per replicate. Any dataset can be
  rebuilt exactly (see "Inspect one replicate" below).

### Inspect one replicate interactively

```r
setwd("sim1_repro"); source("R/load.R")
cell <- CELLS[CELLS$cell == 7, ]          # e.g. ERR + additive + independent
one_rep(cell, rep = 1)                     # same numbers as row 1 of results/<tag>/cell07.rds

# or step through the pieces by hand (this is exactly what one_rep() does):
set.seed(20261028 + 100000 * cell$cell + 1)
L  <- cum_lag_matrix(T_STAR)
ex <- gen_exposure(N_MAIN, T_STAR, cell$error, COR[[cell$cor]], expo_params(cell$error, cell$sd))
x  <- ex$c %*% L; X <- ex$C %*% L          # true / surrogate cumulative exposure
sv <- gen_survival(x, cell$model, BETA[[cell$model]], cell$nu, THETA)
ev <- make_events(sv$time, sv$D, sv$k)
pairs <- sample_validation(sv$time, N_VAL_SUBJ, N_VAL_PAIR, MAX_PER_ID, T_STAR)
```

### Files

| File | Content | Walkthrough step |
|---|---|---|
| `R/params.R` | every constant and the 12 cells | all |
| `R/exposure.R` | correlation matrices, exposure generation, cumulative-sum matrix | 1–3 |
| `R/survival.R` | event-time generation, risk sets | 4–5 |
| `R/validation.R` | validation sampling | 6 |
| `R/fit.R` | Cox / ERR partial likelihood | 7 |
| `R/rrc.R` | risk set regression calibration + sandwich variance | 8–9 |
| `R/one_rep.R` | one replicate end to end | 10 |
| `scripts/01_calibrate.R` | calibration of sd and nu | 2, 5 |
| `scripts/run.R`, `scripts/table1.R` | replicate loop, summary | 10–11 |
| `tests/test_components.R` | checks | — |

---

## Part B — The setting, step by step

### Overview: the 12 cells

Table 1 crosses 2 survival models × 2 error models × 4 ways to generate C. We reproduce 3 of the
4 generation methods, so **12 cells**. The "resample from data" rows are left blank (the Uminer data
are not available).

| cell | model | error | correlation of C | cell | model | error | correlation of C |
|---|---|---|---|---|---|---|---|
| 1 | Cox | additive | independent | 7 | ERR | additive | independent |
| 2 | Cox | additive | AR(1) | 8 | ERR | additive | AR(1) |
| 3 | Cox | additive | DEX | 9 | ERR | additive | DEX |
| 4 | Cox | multiplicative | independent | 10 | ERR | multiplicative | independent |
| 5 | Cox | multiplicative | AR(1) | 11 | ERR | multiplicative | AR(1) |
| 6 | Cox | multiplicative | DEX | 12 | ERR | multiplicative | DEX |

Fixed design (paper, Table 1 footnote), in `R/params.R`:

| Quantity | Value | Name in code |
|---|---|---|
| main study subjects n1 | 2400 | `N_MAIN` |
| follow-up years T* | 55 | `T_STAR` |
| validation subjects n2 | 800 | `N_VAL_SUBJ` |
| validation pairs | 3000, at most 8 per subject | `N_VAL_PAIR`, `MAX_PER_ID` |
| replicates B | 2000 | `B_REPS` |
| true β (Cox) | log(1.5) = 0.405 | `BETA["cox"]` |
| true β (ERR) | 0.423 | `BETA["err"]` |

Notation: subject i = 1..2400, calendar year l = 1..55. **C** = surrogate (error-prone) annual
exposure, **c** = true annual exposure, **X** and **x** = their cumulative sums. Units: 100 WLM.

---

### Step 1 — Correlation of annual exposure across years

**Paper.** The 55 annual surrogate exposures of one subject are correlated; Table 1 uses an
independent, an AR(1) and a damped exponential (DEX) structure, with parameters "chosen to
approximate" the Uminer data (values not reported).

**We do.** Correlation between years l and l' with gap d = |l − l'|:

| Structure | corr(d) | Values | Source of values |
|---|---|---|---|
| independent | 1 if d = 0, else 0 | — | — |
| AR(1) | ρ^d | ρ = 0.62 | Uminer AR(1) fit in the original `Uminer_Data_Feature_Exploration.R` |
| DEX | ρ^(d^θ) | ρ = 0.64, θ = 0.66 | same file |

DEX decays much more slowly than AR(1) at long gaps (θ < 1), so cumulative exposure varies more
between subjects under DEX than under AR(1), and more under AR(1) than independent.

**Code.** `cor_matrix()` in `R/exposure.R`; values in `COR` in `R/params.R`.

**Check.** `cor_matrix(5, COR$ar1)`; the test "AR(1) lag-1 correlation" checks generated data.

---

### Step 2 — Surrogate exposure C and true exposure c

**Paper.** Generate C first, then c from C (a regression-calibration / Berkson-type model):
- additive: C ~ MVN(μ_C, Σ_C), c = γ0 + γ1·C + ε
- multiplicative: C ~ log-MVN(μ_C, Σ_C), log c = γ0 + γ1·log C + ε

with values chosen to approximate Uminer (not reported).

**We do.** For each subject, draw a 55-vector Z ~ MVN(0, R) with R from Step 1, and ε ~ N(0, sd_e²)
independent across years:

| | additive | multiplicative |
|---|---|---|
| surrogate | C = μ + sd·Z | log C = μ + sd·Z |
| true | c = γ0 + γ1·C + ε | log c = γ0 + γ1·log C + ε |

The parameters are not free; they are derived (function `expo_params()`) from these fixed targets
in `EXPO_BASE`:

| Target | additive | multiplicative | Why |
|---|---|---|---|
| mean annual C | 0.05 | 0.05 | sets the exposure level (see Step 5) |
| mean annual c | 0.025 | 0.025 | so E(c \| C) has slope 0.5 in both models |
| slope γ1 | 0.5 | 1.0 | additive: E(c\|C) = 0.5·C; multiplicative with γ1 = 1: E(c\|C) = 0.5·C too, so both error models attenuate by the same amount (≈ −50%, like the paper's naive bias) |
| correlation ρ | cor(c, C) = 0.55 | cor(log c, log C) = 0.60 | moderate error; sets sd_e |

Derived formulas:
- additive: μ = mean_C, γ0 = mean_c − γ1·mean_C (= 0), sd_e = γ1·sd·√(1/ρ² − 1)
- multiplicative: μ = log(mean_C) − sd²/2, sd_e = γ1·sd·√(1/ρ² − 1),
  γ0 = log(mean_c) − γ1·μ − (γ1²·sd² + sd_e²)/2

**sd is calibrated per cell** (`scripts/01_calibrate.R`) to the cell's target in `CELLS$target` /
`CELLS$target_value` (`R/params.R`). The targets differ by model because the constraints differ:

| cells | target | why |
|---|---|---|
| Cox additive (1–3) | sd of surrogate cumulative X at year 45 = 1.17 | naive Cox SE ≈ 1/√(events·Var X); 1.17 gives the paper's naive SE (≈ 0.055) with ~240 events |
| Cox multiplicative (4–6) | mean naive model-based SE at n1 = 2400 = paper's value (0.059 / 0.059 / 0.061) | the sd(X) route overshoots here: lognormal tails carry more information than Var X suggests |
| ERR additive (7–9) | sd of surrogate cumulative X at year 45 = 0.9, with ρ = 0.85 | chosen from a grid (ρ ∈ {0.55, 0.7, 0.85}, sd(X₄₅) ∈ {0.6, 0.75, 0.9, 1.2}) to show RRC removing the measurement-error bias; less spread leaves ERR estimates skewed, more makes some x negative and caps β |
| ERR multiplicative (10–12) | mean naive model-based SE = paper's value (0.083 / 0.092 / 0.091), with ρ = 0.8 | ERR needs cumulative exposure with a large spread relative to its mean |

The naive-SE target is computed as the mean SE over 8 simulated cohorts of exactly 2400 subjects
(not one large cohort rescaled by √n, which is wrong for heavy-tailed exposure). Matching sd or SE
across correlation structures also keeps the cells comparable: without it DEX cells would have about
3× the exposure spread of independent cells.

For the multiplicative cells (Cox and ERR) the error correlation is raised to cor(log c, log C) = 0.8
(`CELLS$rho`). At this spread, ρ = 0.6 makes sd_e so large that true x gets an extreme tail (relative
risks up to 10^47 in a 20,000-subject draw) and the Cox fit fails. Attenuation is set by
E(c | C) = 0.5·C, not by ρ, so the naive bias stays at about −50%.

Calibrated values (`results/calibration.csv`, used by `results/full_v3`; earlier ones are kept in
`results/calibration_v2.csv` for `results/full_v2` and `results/calibration_v1.csv` for `results/full`):

| cells | ρ | sd (additive: of C; mult: of log C) | sd_e |
|---|---|---|---|
| Cox additive ind / AR(1) / DEX | 0.55 | 0.176 / 0.087 / 0.063 | 0.134 / 0.066 / 0.048 |
| Cox multiplicative ind / AR(1) / DEX | 0.80 | 1.37 / 1.12 / 0.92 | 1.03 / 0.84 / 0.69 |
| ERR additive ind / AR(1) / DEX | 0.85 | 0.135 / 0.067 / 0.048 | 0.042 / 0.021 / 0.015 |
| ERR multiplicative ind / AR(1) / DEX | 0.95 / 0.80 / 0.80 | 2.19 / 1.81 / 1.72 | 0.72 / 1.35 / 1.29 |

In the Cox additive cells 5–15% of subjects have negative cumulative x. This is harmless for Cox
but is why this design cannot be reused for ERR.

**Code.** `gen_exposure()` in `R/exposure.R`; `EXPO_BASE`, `expo_params()`, `CELLS$target`, `CELLS$rho` in `R/params.R`.

**Check.** Test 1 ("error model slope recovered") fits lm(c ~ C) and recovers γ1.

---

### Step 3 — Cumulative exposure and timing

**Paper.** x_ij = Σ_{k ≤ j} c_ik and X_ij = Σ_{k ≤ j} C_ik; the hazard in year interval (k, k+1] uses
x_ik, with x_i0 = 0.

**We do.** Exactly that. Exposure through year k drives risk in (k, k+1]: a one-year lag, and nobody
is exposed during (0, 1]. In code, cumulative exposure is stored as an n × 55 matrix whose
**column m holds x_{m−1}**, built as `c %*% L` with L[l, m] = 1 if l < m.

This alignment is used both when generating events (Step 4) and when fitting (Step 7). The old
scripts in the parent folder misaligned these by one year.

**Code.** `cum_lag_matrix()` in `R/exposure.R`.

---

### Step 4 — Survival times

**Paper (Section 4.1, Steps 1–6).** Weibull(ν, θ) baseline, cumulative incidence
F_ij = 1 − exp{−ν^θ Σ_{k=0}^{j−1} r(x_ik)[(k+1)^θ − k^θ]}, draw U, find j with F_ij ≤ U < F_i,j+1,
invert within the interval. r(x) = exp(βx) for Cox, 1 + βx for ERR.

**We do.** The same, written with cumulative hazard:
1. H_ij = ν^θ Σ_{k=0}^{j−1} r(x_ik)·[(k+1)^θ − k^θ], the cumulative hazard at integer time j.
2. Draw E_i = −log(1 − U_i) ~ Exp(1).
3. j_i = number of integer times with H_ij ≤ E_i, so the event falls in (j_i, j_i + 1].
4. T0_i = [ (E_i − H_ij) / (ν^θ·r(x_ij)) + j_i^θ ]^{1/θ}, which is the paper's formula.
5. If j_i = 55, there is no event during follow-up and the subject is censored at T* = 55.
6. No other censoring (paper: "we assume no censoring effect").

Event times are continuous, so there are no ties. Subjects with an event in the first year are
kept (x = 0 there).

For ERR, r = 1 + βx must be positive. If it is negative for some person-year, the rate is set to 0
and counted in `n_viol`. With the current settings, no violations occurred.

**Code.** `gen_survival()` and `rel_rate()` in `R/survival.R`.

---

### Step 5 — Baseline hazard and the exposure level

**Paper.** (ν, θ) chosen so cumulative incidence is 5–10%.

**We do.**
- θ = 4 (`THETA`): events concentrate in later follow-up years.
- ν is calibrated per cell so the expected cumulative incidence is **10%** (`TARGET_CIR`), the
  top of the paper's range, which gives about 240 events. The expected cumulative incidence is
  mean_i[1 − exp(−ν^θ A_i)] with A_i = Σ_k r(x_ik)[(k+1)^θ − k^θ], solved for ν on 20,000
  simulated subjects. Result: ν ≈ 0.0092 for Cox cells and 0.0094 for ERR cells.

**Why the exposure level is low (mean c = 0.025 per year, mean x ≈ 1.1 at year 45).** We arrived at
this level after two failed attempts, and it drives the ERR results:
- In the ERR model, the relative rate between two subjects is (1 + βx_i)/(1 + βx_j). If βx is large
  compared with 1, this is close to x_i/x_j for any large β. The likelihood then keeps rising as β
  grows and the estimate runs away. Our first attempt (mean x ≈ 3.4) gave estimates of 10^13.
- For the additive cells, x is normal and can be negative. The ERR estimate cannot exceed
  −1/min(x), because 1 + βx must stay positive for everyone. A larger spread makes min(x) more
  negative and caps β below its true value. Our first attempt capped β at 0.24, below the true 0.423.
- So for additive exposure there is a ceiling on how much spread we can have. With mean x ≈ 1.1 and
  sd 0.5 at year 45, the cap is about 0.7–1.1 (see `cap_true` / `cap_naive` in
  `results/calibration.csv`).

The paper's median cumulative exposure was 0.204 (Supplement 7.2), also low, which is consistent
with this.

**Code.** `scripts/01_calibrate.R`; `THETA`, `TARGET_CIR` in `R/params.R`.

---

### Step 6 — Internal validation study

**Paper.** Randomly sample n2 = 800 subjects from the main study, then 3000 (c, C) pairs with no more
than 8 per subject.

**We do.**
1. Sample 800 of the 2400 subjects.
2. Pool every (subject, year) these subjects were followed: years 1..ceiling(T_i).
3. Shuffle the pool and keep pairs in order, skipping a subject once it has 8 pairs, until 3000
   are kept.

In these pairs both c and C are observed. Everywhere else only C is observed.

**Code.** `sample_validation()` in `R/validation.R`.

**Check.** Test "validation: 3000 pairs, <= 8 per subject".

---

### Step 7 — Fitting Cox and ERR models

**Paper.** Partial likelihood, Breslow ties; ERR by maximising the ERR partial likelihood (Section 2.2).

**We do.** One engine for both models. At each event time t_e, the risk set is every subject with
T_i ≥ t_e (everyone enters at time 0). The exposure used is the value for the interval containing
t_e. Log partial likelihood:

ℓ(β) = Σ_e [ log r(z_case,e) − log Σ_{j ∈ risk set e} r(z_je) ]

It is maximised by `stats::nlminb` (PORT trust-region optimiser with box constraints), given the
analytic score and observed information. For ERR, β is restricted so that 1 + βz > 0 for every z
used, and is capped at 50. A fit is marked not converged, and excluded from the table, if nlminb
reports non-convergence or the estimate ends on the boundary.

Exposure is passed per (event, risk-set member), not per subject. This lets the RRC exposure differ
between risk sets (Step 8). Three fits per replicate:

| Fit | Exposure | Purpose |
|---|---|---|
| true | x (true cumulative) | reference: the best we could do without error |
| naive | X (surrogate cumulative) | the uncorrected estimate (paper's "Naive") |
| RRC | x̂ from Step 8 | the corrected estimate (paper's "RRC") |

Each cell fits only its own model: Cox cells fit Cox, ERR cells fit ERR.

**SE.** For true and naive, the model-based SE is 1/√(observed information). For RRC, the sandwich
SE (Step 9).

**Code.** `pl_parts()`, `pl_fit()` in `R/fit.R`.

**Check.** Tests "Cox beta/SE matches coxph" (to 1e-6) and "ERR beta matches direct maximisation".

---

### Step 8 — Risk set regression calibration (RRC)

**Paper (Section 3.2).** For each event time t_k: take the validation pairs still at risk (T ≥ t_k),
fit the error model c | C there to get α̂(t_k), predict c for everyone, and refit the survival model
with the risk-set-specific x̂. Merge sparse late risk sets.

**We do.**
1. Sort events by time.
2. **Validation risk set** for event time t: all pairs from validation subjects with T_i ≥ t. These
   sets are nested and shrink over time.
3. **Merging**: working back from the last event, merge risk sets until each has more than 20 pairs.
   This is the rule the paper used for Uminer (Section 5.1). With 10% incidence nearly all subjects
   survive, so merging rarely happens, and there are about as many groups as events (≈ 240).
4. In each group, fit the error model by least squares. This equals GEE with an independence
   working correlation; clustering by subject is handled in the variance (Step 9).
   - additive: c = α0 + α1·C
   - multiplicative: log c = α0 + α1·log C, plus σ̂² = mean squared residual
5. Predict annual true exposure for every non-validated subject-year:
   - additive: ĉ = α̂0 + α̂1·C
   - multiplicative: ĉ = exp(α̂0 + σ̂²/2)·C^α̂1 (the log-normal mean; paper Section 2.3)
6. Cumulative x̂_i(t_e) = sum over years ≤ k_e of (observed c if that year is a validation pair,
   else ĉ). The **whole** exposure history is recalibrated with the α̂ of the event's own risk set.
7. Fit the survival model with x̂ (Step 7).

**Code.** `rrc_groups()`, `rrc_fit()` in `R/rrc.R`.

**Check.** Test "one-group RRC equals RC": with one forced group, RRC equals ordinary regression
calibration from `lm()`.

---

### Step 9 — Variance of the RRC estimate

**Paper (Supplement 7.1.3).** Sandwich: Var(β̂_RRC) = I⁻¹ [H_β + U*·V_α·U*ᵀ] I⁻¹.

**We do.** The same formula, written with sums instead of averages:
- I: observed information of the RRC fit.
- H_β = Σ_i Ũ_i², where Ũ_i is the robust (Lin–Wei) score of subject i (eqs 7.37/7.38).
- U* = ∂U/∂α: how the score changes with each group's error-model parameters (eq 7.34). For
  multiplicative error, σ² is included as a parameter because it enters x̂.
- V_α = B⁻¹ M B⁻¹: the GEE sandwich for all groups' α jointly, with scores summed within
  validation subject, so pairs from the same subject are treated as correlated.
- It is computed without forming the large matrix: s_i = Σ_g (U*_g B_g⁻¹)·u_ig, and the second term
  equals Σ_i s_i².

Also saved, for comparison: `se_rrc_beta_only` (ignores uncertainty in α) and `se_rrc_model`
(model-based).

**Code.** End of `rrc_fit()` in `R/rrc.R`.

**Check.** Test "sandwich SE > SE ignoring alpha".

---

### Step 10 — One replicate and seeds

`one_rep(cell, rep)` (`R/one_rep.R`) runs Steps 2–9 and returns one row:
`n_events`, `n_viol`, `n_groups`, then β̂ and SE for the true, naive and RRC fits, and `conv`
(all three fits converged).

Seed = 20261028 + 100000·cell + rep, so every replicate is reproducible and independent of how
the work is split across cores.

---

### Step 11 — Table 1 metrics

`scripts/table1.R`, over converged replicates of each cell:

| Column | Definition |
|---|---|
| relative bias (%) | 100·(mean β̂ − β)/β |
| estimated SE | mean of the replicate SEs |
| SE/SD | mean SE ÷ empirical SD of β̂ (≈ 1 if the SE is right) |
| 95% coverage (%) | share of replicates with \|β̂ − β\| ≤ 1.96·SE |

Also shown: `true_bias` (bias of the true-x fit, which shows finite-sample bias that has nothing
to do with measurement error) and `cir` (cumulative incidence, %). The paper's values are printed
alongside (`p_*` columns).

---

## Part C — Current status and open decision

**Cox cells, revised design** (pilot B = 200, `results/table1_pilot_v2.csv`; `results/full` used the
old design in `results/calibration_v1.csv`):
- Additive (1–3): naive SE 0.055–0.056 (paper 0.054–0.056), naive bias −51% to −52%, RRC bias −1% to
  −4%, RRC SE 0.108–0.112 (paper 0.093–0.105), SE/SD 0.96–1.07, RRC coverage 96–97.5%.
- Multiplicative AR(1)/DEX (5–6): naive SE 0.060 (paper 0.059–0.061), SE/SD ≈ 1, RRC coverage 93.5–94%.
  RRC bias is −2.3% / +4.4% (paper about −5%).
- Multiplicative independent (4): naive SE 0.065 (target 0.059), SE/SD 0.88, RRC bias +9.7%. Reaching
  the paper's precision from independent years needs very heavy-tailed annual exposure, and a few
  extreme subjects then dominate the fit. This is a limitation of the lognormal-independent design.

**ERR cells, revised design** (pilot B = 200, `results/table1_pilot_v2.csv`):
- True-x fits are now close to unbiased in every cell (−8% to +7%; old design +18% to +21% in the
  multiplicative cells), so the remaining RRC behaviour reflects measurement error, not a lack of
  information.
- Additive (7–9): naive bias −47% to −54%, RRC bias −7% to +6%, which tracks the true-x fit's bias
  in each cell (RRC removes the measurement-error part). RRC SE/SD 1.05–1.17, coverage 92–93%,
  6% / 2.5% / 1% of replicates not converged (fits on the 1 + βx > 0 boundary). SEs (0.30–0.33)
  are about 3× the paper's: additive normal exposure cannot carry more ERR information without
  negative x.
- Multiplicative (10–12): naive SE matches the paper (0.085–0.092), naive bias −47% to −56%,
  RRC bias −10% / −3% / +8%, RRC SE about 0.18 (paper 0.11), coverage 85–92%. The independent cell
  needs heavy-tailed annual exposure, as in Cox cell 4.

Paper choices this reproduction cannot recover: the exposure means and variances, the error-model
coefficients, the correlation parameters actually used, ν and θ, and the "resample from data" rows.
