#### Load all project code; run scripts from the project root (sim1_repro/) ####

if (!exists("PROJ")) PROJ <- getwd()
for (f in c("params.R", "exposure.R", "survival.R", "fit.R", "validation.R", "rrc.R", "one_rep.R")) {
  source(file.path(PROJ, "R", f))
}
