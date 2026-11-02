#### Load Simulation 2 code; run scripts from the project root (sim2_repro/) ####
# survival generation, partial likelihood and RRC estimators are Simulation 1's, unchanged

if (!exists("PROJ")) PROJ <- getwd()
SIM1 <- file.path(dirname(PROJ), "sim1_repro", "R")
for (f in c("exposure.R", "survival.R", "fit.R", "rrc.R")) source(file.path(SIM1, f))
for (f in c("params.R", "exposure2.R", "one_rep.R")) source(file.path(PROJ, "R", f))
