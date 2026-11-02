setwd("~/OneDrive - Harvard University/JUL2015-FASRC/analysis")

library(haven)

crude_exp<-read_sas("crude_exp.sas7bdat")
crude_exp_1<-read_sas("crude_exp_1.sas7bdat")
mv_exp<-read_sas("mv_exp.sas7bdat")
mv_exp_1<-read_sas("~/OneDrive - Harvard University/JUL2015-FASRC/rrc/mv_exp.sas7bdat")
all.equal(mv_exp,mv_exp_1)
#[1] "Component “predexpxm”: Mean relative difference: 0.1467942"
#[2] "Component “i”: Mean relative difference: 0.6"

valid_2009_short<-read_sas("~/OneDrive - Harvard University/JUL2015-FASRC/data/valid_2009_short.sas7bdat")

library(survival)

model.coxph.rrc.adjusted<-coxph(Surv(time,CASE)~(predexpxm+SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5)+strata(age_tv),data=mv_exp)
summary(model.coxph.rrc.adjusted)
#predexpxm       SMOKE      native    hispanic       COMP5       COMP6 
#1.38326702  1.63554065 -0.88283684 -1.03435642  0.51738618  0.71297533 
#yy2         yy3         yy4         yy5 
#0.12275396  0.20207523 -0.08915978 -0.53429846

# surtotxm==cumexp100
# The following cumulative sum procedure is not necessary !!!
# mv_exp_cumsum<-transform(mv_exp,predcumexp100=ave(predexpxm,id,FUN=cumsum))

model.coxph.rrc.adjusted_1<-coxph(Surv(time,CASE)~(predexpxm+SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5)+strata(age_tv),data=mv_exp_1)
summary(model.coxph.rrc.adjusted_1)
#predexpxm       SMOKE      native    hispanic       COMP5       COMP6 
#1.30513512  1.68191613 -0.89941655 -1.01603984  0.47703233  0.63468918 
#yy2         yy3         yy4         yy5 
#0.16660939  0.26196295 -0.02973221 -0.46103043







