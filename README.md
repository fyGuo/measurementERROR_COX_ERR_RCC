# Measurement error in cumulative exposure: Cox, ERR and regression calibration

Reproduction of the simulations in Yang and others, *Adjustment for measurement error in cumulative exposure
variable in the presence of model misspecification* (`Yue_Yang_measurement_Error.pdf`).

| Folder | Contents |
|---|---|
| `sim1_repro/` | Simulation 1 (Section 4.1, Table 1): RRC under correct models. See its README. |
| `sim2_repro/` | Simulation 2 (Section 4.2, Table 2): risk-model and error-model misspecification. Uses `sim1_repro/R/`, so keep both folders side by side. See its README for setup, how to run on another machine, the design and the history of decisions. |
| `Base_Function_R_Code/` | the authors' original R code (reference only; not used by the reproductions) |
| `Latency_Visualization_Code/`, `RRC_Macro_SAS_Code/`, `RRC_Macro_Donna_Website/`, `Reference Code and Package/` | reference material: SAS RRC macros, latency plots, the rERR package |

Requirements: R 4.5 with `future` and `furrr`; `pdflatex` for the PDF tables. Result files (`.rds`, `.csv`,
`.pdf`) are tracked in git so the tables can be reproduced without rerunning.
