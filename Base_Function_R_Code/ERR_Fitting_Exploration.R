#### Load Uminer Data ####

setwd("~/Dropbox/Uminers/uminers/uminers/2015/data")

library(haven)
new_long<-read_sas("new_long.sas7bdat")
# 79750*26
# long format
valid_2009<-read_sas("valid_2009.sas7bdat")
# 2833*97
# 862 unique IDs
ca_new<-read_sas("ca_new.sas7bdat")
# 2009*166
# wide format: 1 id per row
exp_vs<-read.table("exp_vs.dat")
# 2833*4 
# unclear what the columns are yet
View(ca_new)
colnames(ca_new)
#[1] "ID"           "FYR"          "FYRNM"        "ETHNIC"       "BIRTH_YY"    
#[6] "SMOKE"        "COAL"         "NONCOAL"      "FAGE"         "EXPSTRT"     
#[11] "EXPEND"       "YREXP75"      "YREXP76"      "YREXP77"      "YREXP78"     
#[16] "YREXP79"      "YREXP80"      "YREXP81"      "YREXP68"      "YREXP69"     
#[21] "YREXP70"      "YREXP71"      "YREXP72"      "YREXP73"      "YREXP74"     
#[26] "YREXP61"      "YREXP62"      "YREXP63"      "YREXP64"      "YREXP65"     
#[31] "YREXP66"      "YREXP67"      "YREXP60"      "YREXP82"      "YREXP59"     
#[36] "YREXP84"      "YREXP85"      "YREXP58"      "YREXP83"      "YREXP56"     
#[41] "YREXP57"      "YREXP52"      "YREXP53"      "YREXP54"      "YREXP55"     
#[46] "YREXP50"      "YREXP51"      "YREXP49"      "YREXP48"      "YREXP47"     
#[51] "YREXP42"      "YREXP43"      "YREXP44"      "YREXP45"      "YREXP46"     
#[56] "CUMEXP42"     "CUMEXP43"     "CUMEXP44"     "CUMEXP45"     "CUMEXP46"    
#[61] "CUMEXP47"     "CUMEXP48"     "CUMEXP49"     "CUMEXP50"     "CUMEXP51"    
#[66] "CUMEXP52"     "CUMEXP53"     "CUMEXP54"     "CUMEXP55"     "CUMEXP56"    
#[71] "CUMEXP57"     "CUMEXP58"     "CUMEXP59"     "CUMEXP60"     "CUMEXP61"    
#[76] "CUMEXP62"     "CUMEXP63"     "CUMEXP64"     "CUMEXP65"     "CUMEXP66"    
#[81] "CUMEXP67"     "CUMEXP68"     "CUMEXP69"     "CUMEXP70"     "CUMEXP71"    
#[86] "CUMEXP72"     "CUMEXP73"     "CUMEXP74"     "CUMEXP75"     "CUMEXP76"    
#[91] "CUMEXP77"     "CUMEXP78"     "CUMEXP79"     "CUMEXP80"     "CUMEXP81"    
#[96] "CUMEXP82"     "CUMEXP83"     "CUMEXP84"     "CUMEXP85"     "YREXP86"     
#[101] "CUMEXP86"     "YREXP87"      "CUMEXP87"     "YREXP88"      "CUMEXP88"    
#[106] "YREXP89"      "CUMEXP89"     "YREXP90"      "CUMEXP90"     "YREXP91"     
#[111] "CUMEXP91"     "YREXP92"      "CUMEXP92"     "YREXP93"      "CUMEXP93"    
#[116] "AGE_END"      "AGE_STRT"     "MINETIME"     "AGE"          "LUNG_CA"     
#[121] "NATIVE"       "HISPANIC"     "COMP5"        "COMP6"        "COMP12"      
#[126] "DATEBIRTH"    "byear"        "first_age"    "first_age2"   "first_age3"  
#[131] "EXPTOT93"     "year41"       "vital_status" "DLA"          "COD"         
#[136] "VER"          "year"         "DATEDIED"     "can"          "all_death"   
#[141] "all_death_c"  "cohort1"      "cohort2"      "cohort3"      "age_new"     
#[146] "cumexp94"     "cumexp95"     "cumexp96"     "cumexp97"     "cumexp98"    
#[151] "cumexp99"     "cumexp100"    "cumexp101"    "cumexp102"    "cumexp103"   
#[156] "cumexp104"    "cumexp105"    "cumexp106"    "cumexp107"    "cumexp108"   
#[161] "cumexp109"    "cumexp110"    "cumexp111"    "cumexp112"    "CALYEAR"     
#[166] "CALYEAR2" 


#### Exploration: linERR Package ####

library(linERR)
#model.err<-fit.linERR(Surv(AGE_STRT,AGE,CASE)~SMOKE+ETHNIC+COMP5+COMP6|CUMEXP,data=new_long,beta=NULL,lag=0)

agedf<-data.frame(matrix(0,nrow=nrow(ca_new),ncol=112-42+1))
colnames(agedf)<-paste0("AGE",42:112)
for(j in 1:(112-42+1)){
  for(i in 1:nrow(ca_new)){
    if (j+42-1<ca_new$BIRTH_YY[i]){
      agedf[i,j]<-0
    }
    else if (j+42-1>=ca_new$BIRTH_YY[i]&&j+42-1<=ca_new$CALYEAR[i]){
      agedf[i,j]<-j+42-1-ca_new$BIRTH_YY[i]
    }
    else if (j+42-1>ca_new$CALYEAR[i]){
      agedf[i,j]<-0
    }
  }
}
#View(agedf)
ca_new<-cbind(ca_new,agedf)

model.err<-fit.linERR(Surv(AGE_STRT,AGE_END,LUNG_CA)~SMOKE+NATIVE+HISPANIC+COMP5+COMP6|
                        CUMEXP42+CUMEXP43+CUMEXP44+CUMEXP45+CUMEXP46+CUMEXP47+CUMEXP48+CUMEXP49+
                        CUMEXP50+CUMEXP51+CUMEXP52+CUMEXP53+CUMEXP54+CUMEXP55+CUMEXP56+CUMEXP57+CUMEXP58+CUMEXP59+
                        CUMEXP60+CUMEXP61+CUMEXP62+CUMEXP63+CUMEXP64+CUMEXP65+CUMEXP66+CUMEXP67+CUMEXP68+CUMEXP69+
                        CUMEXP70+CUMEXP71+CUMEXP72+CUMEXP73+CUMEXP74+CUMEXP75+CUMEXP76+CUMEXP77+CUMEXP78+CUMEXP79+
                        CUMEXP80+CUMEXP81+CUMEXP82+CUMEXP83+CUMEXP84+CUMEXP85+CUMEXP86+CUMEXP87+CUMEXP88+CUMEXP89+
                        CUMEXP90+CUMEXP91+CUMEXP92+CUMEXP93+cumexp94+cumexp95+cumexp96+cumexp97+cumexp98+cumexp99+
                        cumexp100+cumexp101+cumexp102+cumexp103+cumexp104+cumexp105+cumexp106+cumexp107+cumexp108+cumexp109+
                        cumexp110+cumexp111+cumexp112,
                      data=ca_new, ages=ca_new[,167:237],beta=NULL,lag=0)
#Error in `$<-.data.frame`(`*tmp*`, "id", value = c(1, 2, 3, 4, 5, 6, 7,  : 
#                                                     replacement has 2009 rows, data has 1717
#2009-1717
as.data.frame(colSums(is.na(ca_new)))
# NAs from 1-292 in CUMEXP63-CUMEXP93, no other NAs except in DATEDIED for non-cases 1159

# Fix the NAs in ca_new 
for(i in 63:93){
  ca_new[,colnames(ca_new)==paste0("CUMEXP",i)]<-ca_new[,colnames(ca_new)==paste0("CUMEXP",i-1)]+ca_new[,colnames(ca_new)==paste0("YREXP",i)]
}

model.err<-fit.linERR(Surv(AGE_STRT,AGE_END,LUNG_CA)~SMOKE+NATIVE+HISPANIC+COMP5+COMP6|
                        YREXP42+YREXP43+YREXP44+YREXP45+YREXP46+YREXP47+YREXP48+YREXP49+
                        YREXP50+YREXP51+YREXP52+YREXP53+YREXP54+YREXP55+YREXP56+YREXP57+YREXP58+YREXP59+
                        YREXP60+YREXP61+YREXP62+YREXP63+YREXP64+YREXP65+YREXP66+YREXP67+YREXP68+YREXP69+
                        YREXP70+YREXP71+YREXP72+YREXP73+YREXP74+YREXP75+YREXP76+YREXP77+YREXP78+YREXP79+
                        YREXP80+YREXP81+YREXP82+YREXP83+YREXP84+YREXP85+YREXP86+YREXP87+YREXP88+YREXP89+
                        YREXP90+YREXP91+YREXP92+YREXP93,
                      data=ca_new, ages=ca_new[,167:218],beta=NULL,lag=0)
# Parameter Summary Table:
#   Estimate  Std. Error  Test Stat.    p-value
# dose             0.1159973   0.2935579  0.39514273 0.69273756
# ca_new$SMOKE     1.8316799   1.0237075  1.78926106 0.07357278
# ca_new$NATIVE    0.5078312   1.0132107  0.50120985 0.61622344
# ca_new$HISPANIC -2.0286230   1.0305083 -1.96856545 0.04900301
# ca_new$COMP5    13.0183530 899.3532635  0.01447524 0.98845083
# ca_new$COMP6    13.1507964 899.3532410  0.01462250 0.98833335
# 
# AIC:  290.3557 
# Deviance:  278.3557 
# Informative risk sets:  28

model.err.2<-fit.linERR(Surv(AGE_STRT,AGE_END,LUNG_CA)~SMOKE+NATIVE+HISPANIC+COMP5+COMP6|
                        CUMEXP42+CUMEXP43+CUMEXP44+CUMEXP45+CUMEXP46+CUMEXP47+CUMEXP48+CUMEXP49+
                        CUMEXP50+CUMEXP51+CUMEXP52+CUMEXP53+CUMEXP54+CUMEXP55+CUMEXP56+CUMEXP57+CUMEXP58+CUMEXP59+
                        CUMEXP60+CUMEXP61+CUMEXP62+CUMEXP63+CUMEXP64+CUMEXP65+CUMEXP66+CUMEXP67+CUMEXP68+CUMEXP69+
                        CUMEXP70+CUMEXP71+CUMEXP72+CUMEXP73+CUMEXP74+CUMEXP75+CUMEXP76+CUMEXP77+CUMEXP78+CUMEXP79+
                        CUMEXP80+CUMEXP81+CUMEXP82+CUMEXP83+CUMEXP84+CUMEXP85+CUMEXP86+CUMEXP87+CUMEXP88+CUMEXP89+
                        CUMEXP90+CUMEXP91+CUMEXP92+CUMEXP93,
                      data=ca_new, ages=ca_new[,167:218],beta=NULL,lag=0)
# Parameter Summary Table:
#   Estimate Std. Error Test Stat.       p-value
# dose             81.9313272  2.7849887  29.418909 3.147130e-190
# ca_new$SMOKE      1.7467662  1.0218738   1.709376  8.738138e-02
# ca_new$NATIVE     0.5552339  0.9990048   0.555787  5.783565e-01
# ca_new$HISPANIC  -2.1321159  1.0303026  -2.069407  3.850787e-02
# ca_new$COMP5    246.9858853         NA         NA            NA
# ca_new$COMP6    247.1392102         NA         NA            NA
# 
# AIC:  270.7456 
# Deviance:  258.7456 
# Informative risk sets:  28 
# Warning message:
#   In summary.fit.linERR(model.err.2) : Hessian is not positive definite


col_order<-c(56:101,103,105,107,109,111,113,115,146:164)
model.err<-fit.linERR(Surv(AGE_STRT,AGE_END,LUNG_CA)~SMOKE|
                        CUMEXP42+CUMEXP43+CUMEXP44,
                      ages=,
                      data=ca_new,beta=NULL,lag=0)

ca_new$age

fit.1 <- fit.linERR(Surv(entryage, exitage, leu)~sex|dose1+dose2+dose3+dose4+dose5+dose6+
                      dose7+dose8+dose9+dose10+dose11+dose12+dose13+dose14+dose15+dose16+
                      dose17+dose18+dose19+dose20+dose21+dose22+dose23+dose24+dose25+dose26+
                      dose27+dose28+dose29+dose30+dose31+dose32, data=cohort1, beta=NULL, 
                    ages=cohort1[, 7:38], lag=2)




#### Exploration: Uminer Data with coxph() ####

setwd("~/OneDrive - Harvard University/JUL2015-FASRC/data")

library(haven)
final_data<-read_sas("final_data.sas7bdat")
# 88843 rows, 40 columns

summary(final_data)
# 438 NAs for EXPSTRT, EXPEND, AGESTRT
# 1 NA for CEXP
# 42783 NAs for LASTE

final_data$YR_LAST<-as.numeric(final_data$year)-1900

final_data_1<-final_data[which((is.na(final_data$EXPSTRT)==FALSE&final_data$YR_TV>=final_data$EXPSTRT&final_data$YR_TV<=final_data$YR_LAST)|(is.na(final_data$EXPSTRT)==TRUE&final_data$CUMEXP!=0)),]
# 88405 rows, 41 columns
# manual double check
# View(final_data[which(is.na(final_data$EXPSTRT)==TRUE),])

final_data_1$surf<-ifelse(final_data_1$type=="e"&final_data_1$YREXP!=0,1,0)
final_data_1$white<-ifelse(final_data_1$ETHNIC==1,1,0)
final_data_1$hispanic<-ifelse(final_data_1$ETHNIC==2,1,0)
final_data_1$native<-ifelse(final_data_1$ETHNIC==3,1,0)
final_data_1$time<-ifelse(final_data_1$CASE==0,1,1-(12-final_data_1$mnth)/12)
final_data_1$yy1<-ifelse(final_data_1$YR_TV<70,1,0)
final_data_1$yy2<-ifelse(final_data_1$YR_TV>=70&final_data_1$YR_TV<80,1,0)
final_data_1$yy3<-ifelse(final_data_1$YR_TV>=80&final_data_1$YR_TV<90,1,0)
final_data_1$yy4<-ifelse(final_data_1$YR_TV>=90&final_data_1$YR_TV<100,1,0)
final_data_1$yy5<-ifelse(final_data_1$YR_TV>=100,1,0)
final_data_1$cumexp100<-final_data_1$CUMEXP/100
# final_data_1: 88405*51

final_data_1$AGE_OUT<-final_data_1$AGE_TV+final_data_1$time

library(survival)
# these two COXPH models output match with draft results
model.coxph.naive<-coxph(Surv(time,CASE)~(cumexp100)+strata(AGE_TV),data=final_data_1)
summary(model.coxph.naive)
model.coxph.adjusted<-coxph(Surv(time,CASE)~(cumexp100+yy2+yy3+yy4+yy5+SMOKE+native+hispanic+COMP5+COMP6)+strata(AGE_TV),data=final_data_1)
summary(model.coxph.adjusted)


glmfit<-glm(CASE~cumexp100+SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5,family="poisson",data=final_data_1)
glmfit$coefficients
#(Intercept)   cumexp100       SMOKE      native    hispanic       COMP5 
#-10.7849190   0.4386409   1.6945203  -1.2991482  -1.3036259   0.5397272 
#COMP6         yy2         yy3         yy4         yy5 
#0.8467359   1.1771163   2.0819081   2.4633279   2.6870119
glmfit.background<-glm(CASE~SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5,family="poisson",data=final_data_1)
glmfit.background$coefficients


#### Exploration: likelihood-maximization of ERR model with Poisson count assumption ####

# ERR model with Poisson count assumptions
S <- sum(final_data_1$CASE) # total number of cases
bcvr <- cbind(rep(1, nrow(final_data_1)), final_data_1$SMOKE,
              final_data_1$native,final_data_1$hispanic,final_data_1$COMP5,
              final_data_1$COMP6,final_data_1$yy2,final_data_1$yy3,final_data_1$yy4,final_data_1$yy5) # covariance matrix
M <- ncol(bcvr) # number of covariates at the background risk term
Sk <- as.numeric(final_data_1$CASE%*%bcvr) # inner product of cases with covariates
doses <- cbind(final_data_1$cumexp100) # covariates of the RR term
Nd <- ncol(doses) # number of covariates at the RR term
# log-likelihood
llk <- function(b) sum(log(dpois(final_data_1$CASE, 1*exp(bcvr%*%b[1:M])*(1 + doses%*%b[(M + 1):(M + Nd)]))))
# gradient of the log-likelihood with respect the model parameters
grd <- function(t){ 
  pyEi <- as.numeric(1*exp(bcvr%*%t[1:M]))
  ERRi <- as.numeric(1 + doses%*%t[(M + 1):(M + Nd)])
  c( Sk - ((ERRi*pyEi)%*%bcvr), (final_data_1$CASE/ERRi)%*%doses - (pyEi)%*%doses )
}
# definition of constrains for ERR parameters: beta > -1/max(D)
A <- matrix(c(rep(0, M), 1), ncol=M + 1)
B <- 1/max(final_data_1$cumexp100)
constraints <- list(ineqA = A, ineqB = B) # constrains, RR term must be positive
# model fitting
library(maxLik)
err.fit <- maxLik(logLik = llk, grad = grd, start = c(as.numeric(glmfit.background$coefficients), .1), constraints = constraints)
summary(err.fit)
# Maximum Likelihood estimation
# BFGS maximization, 46 iterations
# Return code 0: successful convergence 
# Log-Likelihood: -867.4137 
# 11  free parameters
# Estimates:
#   Estimate Std. error t value  Pr(> t)    
# [1,] -11.5125     0.8444 -13.634  < 2e-16 ***
#   [2,]   1.5794     0.3901   4.049 5.15e-05 ***
#   [3,]  -1.0017     0.7181  -1.395 0.163050    
# [4,]  -1.2087     0.2358  -5.125 2.97e-07 ***
#   [5,]   0.1290     0.4892   0.264 0.791958    
# [6,]   0.4595     0.4690   0.980 0.327169    
# [7,]   1.2891     0.6362   2.026 0.042749 *  
#   [8,]   2.2324     0.6033   3.700 0.000215 ***
#   [9,]   2.6637     0.5993   4.445 8.80e-06 ***
#   [10,]   2.8961     0.5991   4.834 1.34e-06 ***
#   [11,]   5.0210     2.0214   2.484 0.012995 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Warning: constrained likelihood estimation. Inference is probably wrong
# Constrained optimization based on constrOptim 
# 1  outer iterations, barrier value -0.0004706798 

err.fit2 <- maxLik(logLik = llk, start = c(as.numeric(glmfit.background$coefficients), 3), constraints = constraints, control=list(iterlim=5000))
# change the starting values from 5 to 3
# both converges after 2500+ iterations
# yield different results for 11th parameter: 5.06353 to 4.18253

#### Exploration: linERR package continued with 1 exposure variable only ####


# ERR model fit using linERR package
library(linERR)
final_data_1$AGE_OUT<-final_data_1$AGE_TV+final_data_1$time
model.err.adjusted<-fit.linERR(Surv(AGE_TV,AGE_OUT,CASE)~SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5|cumexp100,
                               data=final_data_1, ages=final_data_1[,24],beta=NULL,lag=0)
# so far, doesn't work because the function requires more than one doses column

model.err.adjusted<-fit.linERR1(Surv(AGE_TV,AGE_OUT,CASE)~SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5|cumexp100,
                               data=final_data_1, ages=final_data_1[,24],beta=NULL,lag=0)

fit.linERR1<-function (formula, beta = NULL, data, ages, lag = 0) 
{
  Call <- match.call()
  if (missing(data)) 
    data <- environment(formula)
  mf <- match.call()
  mm <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, mm)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- suppressWarnings(eval(mf, parent.frame()))
  mt <- attr(mf, "terms")
  attr(mt, "predvars") <- gsub(" ", "", attr(mt, "predvars"), 
                               fixed = TRUE)
  loglin.part <- strsplit(strsplit(x = as.character(attr(mt, 
                                                         "predvars")[3]), split = "|", fixed = TRUE)[[1]][1], 
                          split = "+", fixed = TRUE)
  lin.part <- strsplit(strsplit(x = as.character(attr(mt, "predvars")[3]), 
                                split = "|", fixed = TRUE)[[1]], split = "+", fixed = TRUE)
  for (i in 1:length(lin.part)) {
    if (length(grep("strata", lin.part[[i]])) != 0) 
      lin.part[[i]] <- lin.part[[i]][-grep("strata", lin.part[[i]], 
                                           fixed = TRUE)]
    lin.part[[i]] <- gsub("\n", "", lin.part[[i]])
  }
  lin.part1 <- lin.part[[2]]
  max.exp <- length(lin.part1)
  doses <- matrix(nrow = dim(data)[1], ncol = max.exp)
  doses[, 1] <- cbind(eval(parse(text = paste0(Call$data, "$", 
                                               lin.part1[1]))))
  # Would require max.exp>=2 to work
  #for (i in 2:max.exp) {
  #  doses[, i] <- cbind(eval(parse(text = paste0(Call$data, 
  #                                               "$", lin.part1[i]))))
  #}
  data$cumd <- apply(doses, 1, sum)
  end <- round(eval(parse(text = paste0(Call$data, "$", strsplit(x = attr(mt, 
                                                                          "predvars")[[2]], ",", fixed = TRUE)[[1]][2]))), 8)
  status <- eval(parse(text = paste0(Call$data, "$", gsub(")", 
                                                          "", x = strsplit(x = attr(mt, "predvars")[[2]], ",", 
                                                                           fixed = TRUE)[[1]][3]))))
  data[, ncol(data) + 1] <- end
  data[, ncol(data) + 1] <- status
  data <- data[, c(1:4, ncol(data) - 1, ncol(data), (5:ncol(data))[-c(ncol(data) - 
                                                                        1, ncol(data))])]
  data[, ncol(data)] <- NULL
  data[, ncol(data)] <- NULL
  colnames(data)[5] <- "end"
  colnames(data)[6] <- "status"
  if (length(grep("strata", formula)) == 0) 
    model <- coxph(eval(parse(text = attr(mt, "predvars")[[2]])) ~ 
                     cumd, data = data)
  if (length(grep("strata", formula)) > 0) 
    model <- coxph(eval(parse(text = attr(mt, "predvars")[[2]])) ~ 
                     cumd + strata(eval(parse(text = formula[[3]][[3]][[3]][[2]]))), 
                   data = data)
  rsets <- as.data.frame(coxph.detail(model, riskmat = TRUE)$riskmat)
  rsets$id <- seq(1, dim(data)[1], 1)
  id <- rsets$id
  ages <- as.data.frame(cbind(id, ages))
  doses <- as.data.frame(cbind(id, doses))
  data_2 <- as.list(data)
  data_2 <- lapply(data_2, as.numeric)
  rsets_2 <- as.list(rsets)
  rsets_2 <- lapply(rsets_2, as.numeric)
  ages_2 <- as.list(ages)
  ages_2 <- lapply(ages_2, as.numeric)
  doses_2 <- as.list(doses)
  doses_2 <- lapply(doses_2, as.numeric)
  covariates1 <- NULL
  suppressWarnings(if (is.na(as.integer(loglin.part[[1]]))) {
    covariates1 <- sapply(paste0(Call$data, "$", loglin.part[[1]]), 
                          function(x) eval(parse(text = x)))
  })
  covars1_2 <- ifelse(loglin.part[[1]] == 1, NA, apply(covariates1, 
                                                       2, as.list))
  if (all(is.na(covars1_2))) 
    covars1_2 <- NULL
  covars1_2 <- lapply(covars1_2, as.numeric)
  colnames(rsets)[1:(length(rsets_2) - 1)] <- round(as.numeric(colnames(rsets)[1:(length(rsets_2) - 
                                                                                    1)]), 8)
  names(rsets_2)[1:(length(rsets_2) - 1)] <- round(as.numeric(colnames(rsets)[1:(length(rsets_2) - 
                                                                                   1)]), 8)
  failtimes <- as.list(colnames(rsets)[1:(length(rsets_2) - 
                                            1)])
  failtimes <- lapply(failtimes, as.numeric)
  suppressWarnings(if (is.null(beta)) {
    if (!is.na(as.integer(loglin.part[[1]])) & all(as.integer(loglin.part[[1]])) == 
        1 & length(lin.part) <= 2) {
      beta <- 0.1
    }
    if (!is.na(as.integer(loglin.part[[1]])) & all(as.integer(loglin.part[[1]])) == 
        1 & length(lin.part) > 2) {
      beta <- rep(0.1, length(lin.part) - 1)
    }
    if (is.na(as.integer(loglin.part[[1]])) & length(lin.part) == 
        2) {
      if (length(grep("strata", formula)) == 0) 
        beta <- c(0.1, coxph(eval(parse(text = attr(mt, 
                                                    "predvars")[[2]])) ~ as.matrix(covariates1), 
                             data = data)$coef)
      if (length(grep("strata", formula)) != 0) 
        beta <- c(0.1, coxph(eval(parse(text = attr(mt, 
                                                    "predvars")[[2]])) ~ as.matrix(covariates1) + 
                               strata(eval(parse(text = formula[[3]][[3]][[3]][[2]]))), 
                             data = data)$coef)
    }
    if (is.na(as.integer(loglin.part[[1]])) & length(lin.part) > 
        2) {
      if (length(grep("strata", formula)) == 0) 
        beta <- c(0.1, coxph(eval(parse(text = attr(mt, 
                                                    "predvars")[[2]])) ~ as.matrix(covariates1), 
                             data = data)$coef, rep(0.1, length(lin.part) - 
                                                      2))
      if (length(grep("strata", formula)) != 0) 
        beta <- c(0.1, coxph(eval(parse(text = attr(mt, 
                                                    "predvars")[[2]])) ~ as.matrix(covariates1) + 
                               strata(eval(parse(text = formula[[3]][[3]][[3]][[2]]))), 
                             data = data)$coef, rep(0.1, length(lin.part) - 
                                                      2))
    }
  })
  add.arguments <- function(f, n) {
    t = paste("arg <- alist(", paste(sapply(1:n, function(i) paste("x", 
                                                                   i, "=", sep = "")), collapse = ","), ")", sep = "")
    formals(f) <- eval(parse(text = t))
    f
  }
  p.est <- function() {
    beta3 <- vector()
    for (i in 1:length(beta)) {
      beta3[i] <- eval(parse(text = paste0("x", i)))
    }
    beta_2 <- as.list(beta3)
    beta_2 <- lapply(beta_2, as.numeric)
    suppressWarnings(if (loglin.part[[1]] == 1) {
      ncovs1 <- 0
    }
    else {
      ncovs1 <- length(loglin.part[[1]])
    })
    res1 <- .Call("llhood", beta_2, data_2, rsets_2, ages_2, 
                  doses_2, covars1_2, as.integer(length(doses_2[[1]])), 
                  as.integer(length(rsets_2) - 1), as.integer(length(data_2)), 
                  failtimes, as.integer(ncovs1), as.integer(max.exp), 
                  as.numeric(lag))
    return(-sum(res1 * is.finite(res1), na.rm = T))
  }
  p.est <- add.arguments(p.est, length(beta))
  beta <- as.list(beta)
  names(beta) <- paste0("x", seq(1, length(beta)))
  llim <- c(-1/max(data$cumd), rep(-Inf, sum(length(covars1_2), 
                                             length(lin.part) - 2)))
  res <- mle(p.est, start = beta, method = "L-BFGS-B", lower = llim, 
             upper = Inf)
  vcov <- vcov(res)
  aic <- AIC(res)
  attr(res, "lowb") <- llim
  attr(res, "beta") <- beta
  attr(res, "max.exp") <- max.exp
  attr(res, "covariates1") <- covariates1
  attr(res, "data_2") <- data_2
  attr(res, "rsets_2") <- rsets_2
  attr(res, "doses_2") <- doses_2
  attr(res, "ages_2") <- ages_2
  attr(res, "call")$start <- beta
  attr(res, "call")$lower <- llim
  attr(res, "vcov") <- vcov
  attr(res, "aic") <- aic
  attr(res, "Call") <- Call
  attr(res, "llike") <- -attr(res, "details")$value
  attr(res, "deviance") <- 2 * attr(res, "details")$value
  class(res) <- c("fit.linERR")
  return(res)
}


#### Exploration: likelihood-maximization of ERR model using optim() ####


# checkout distribution of 132 cases
final_data_cases<-final_data_1[which(final_data_1$CASE==1),c(1,11,18,24,26,39,52)]
summary(final_data_cases)
length(unique(final_data_cases$AGE_TV))
# 39 unique AGE_TV among 132 cases
length(unique(final_data_cases$AGE_OUT))
# 113 unique AGE_OUT among 132 cases
View(final_data_cases)
# 1/2/3 ties at each distinct death ages/times 

# ERR model fit using partial likelihood approach
S <- sum(final_data_1$CASE) # total number of cases
#bcvr <- cbind(rep(1, nrow(final_data_1)), final_data_1$SMOKE,
#              final_data_1$native,final_data_1$hispanic,final_data_1$COMP5,
#              final_data_1$COMP6,final_data_1$yy2,final_data_1$yy3,final_data_1$yy4,final_data_1$yy5) # covariance matrix
bcvr <- cbind(final_data_1$SMOKE,
              final_data_1$native,final_data_1$hispanic,final_data_1$COMP5,
              final_data_1$COMP6,final_data_1$yy2,final_data_1$yy3,final_data_1$yy4,final_data_1$yy5) # covariance matrix
M <- ncol(bcvr) # number of covariates at the background risk term, no intercept
Sk <- as.numeric(final_data_1$CASE%*%bcvr) # inner product of cases with covariates
doses <- cbind(final_data_1$cumexp100) # covariates of the RR term
Nd <- ncol(doses) # number of covariates at the RR term

# partial log-likelihood term for each unique failure times
d_times<-unique(final_data_cases$AGE_OUT)
# define failure rate lambda(beta,i) function
#fr<-function(b,i){
#  rate<-exp(bcvr[i,]%*%b[1:M])*(1+doses[i,]%*%b[(M+1):(M+Nd)])
#  return(rate)
#}
#b<-c(2,-1,-1,1,1,1,1,-1,-1,0.5)

llk_err<-function(b){
  L<-0
  for(i in 1:length(d_times)){
    # create death set: may have ties, up to 3 ties in each Death set
    D_set<-which(final_data_1$CASE==1&final_data_1$AGE_OUT==d_times[i])
    D_tienum<-length(D_set)
    D_age<-unique(final_data_1$AGE_TV[D_set])
    D_time<-unique(final_data_1$AGE_OUT[D_set])
    # create risk set: stratified by AGE_TV, survival past death time
    S_set<-which(final_data_1$AGE_TV==D_age&final_data_1$AGE_OUT>=D_time)
    # write log likelihood term adjusting for ties using Efron method
    tie_adjust<-mean(exp(bcvr[D_set,]%*%b[1:M])*(1+as.matrix(doses[D_set,])%*%b[(M + 1):(M + Nd)]))
    denom_adjust<-0
    for(r in 1:D_tienum){
      denom_adjust<-denom_adjust+log(sum(exp(bcvr[S_set,]%*%b[1:M])*(1+as.matrix(doses[S_set,])%*%b[(M + 1):(M + Nd)]))-(r-1)*tie_adjust)
    }
    L<-L+sum(bcvr[D_set,]%*%b[1:M])+sum(log(1+as.matrix(doses[D_set,])%*%b[(M + 1):(M + Nd)]))-denom_adjust
  }
  return(L)
}
# definition of constrains for ERR parameters: beta > -1/max(D)
A <- matrix(c(rep(0, M), 1), ncol=M + 1)
B <- 1/max(final_data_1$cumexp100)
constraints <- list(ineqA = A, ineqB = B) # constrains, RR term must be positive

model.err.pllk<-maxLik(logLik=llk_err,start=c(as.numeric(model.coxph.adjusted$coefficients[2:10]),1),constraints=constraints,control=list(iterlim=5000))
summary(model.err.pllk)
# converges after 2293 iterations
# Warning messages:
# 1: In sqrt(diag(vc)) : NaNs produced
# 2: In sqrt(diag(vc)) : NaNs produced
# std error estimates for yy2-yy5 are NAs
# inference using this function is wrong anyway, so not a problem
model.err.pllk$estimate
#[1]  1.8004681 -0.6880004 -0.9124938  0.3533800  0.5079917  0.3712569  0.4656554  0.1326030
#[9] -0.2874889  0.8415461
# Likelihood value=848.3824
model.coxph.adjusted$coefficients
#cumexp100        yy2        yy3        yy4        yy5      SMOKE     native 
#0.2976242  0.1397019  0.1048294 -0.2646445 -0.7092440  1.6713391 -0.8613983 
#hispanic      COMP5      COMP6 
#-0.9537376  0.3015765  0.4851583

model.err.pllk2<-maxLik(logLik=llk_err,start=c(1,as.numeric(model.coxph.adjusted$coefficients[2:10]),1),constraints=constraints,control=list(iterlim=5000),method="BFGS")
#[1]  0.9999990  1.6737207 -0.8078639 -0.9455561  0.2581093  0.3988092  0.1833801  0.1684174
#[9] -0.1634893 -0.5786441  0.7509628
# Likelihood value=848.0932

llk_err_neg<-function(b){
  L<-0
  for(i in 1:length(d_times)){
    # create death set: may have ties, up to 3 ties in each Death set
    D_set<-which(final_data_1$CASE==1&final_data_1$AGE_OUT==d_times[i])
    D_tienum<-length(D_set)
    D_age<-unique(final_data_1$AGE_TV[D_set])
    D_time<-unique(final_data_1$AGE_OUT[D_set])
    # create risk set: stratified by AGE_TV, survival past death time
    S_set<-which(final_data_1$AGE_TV==D_age&final_data_1$AGE_OUT>=D_time)
    # write log likelihood term adjusting for ties using Efron method
    tie_adjust<-mean(exp(bcvr[D_set,]%*%b[1:M])*(1+as.matrix(doses[D_set,])%*%b[(M + 1):(M + Nd)]))
    denom_adjust<-0
    for(r in 1:D_tienum){
      denom_adjust<-denom_adjust+log(sum(exp(bcvr[S_set,]%*%b[1:M])*(1+as.matrix(doses[S_set,])%*%b[(M + 1):(M + Nd)]))-(r-1)*tie_adjust)
    }
    L<-L+sum(bcvr[D_set,]%*%b[1:M])+sum(log(1+as.matrix(doses[D_set,])%*%b[(M + 1):(M + Nd)]))-denom_adjust
  }
  return(-L)
}
model.err.optim<-optim(par=c(as.numeric(glmfit.background$coefficients[2:10]),0.1),fn=llk_err_neg,control=list(maxit=5000))
model.err.optim$convergence
# 0 indicates successful convergence
model.err.optim$par
# model.coxph.adjusted$coefficients[2:10] as initial values
#[1]  1.6917439 -0.4440592 -0.9432146  0.3840316  0.5205391  0.4910298  0.4827413  0.1607169
#[9] -0.2783559  0.7470721
# Likelihood value=848.3997
# Unstable results depending on initial values!!!
# glmfit.background$coefficients[2:10] as initial values
#$par
#[1]  1.7532493 -0.3235582 -0.9782685  0.2307318  0.3768119  1.2806187
#[7]  1.2365792  0.9264990  0.5350139  0.7219235
#$value
#[1] 849.676

model.err.optim2<-optim(par=c(1,as.numeric(model.coxph.adjusted$coefficients[2:10]),1),fn=llk_err_neg,control=list(maxit=5000),method="BFGS")
#[1]  1.0000000  1.6738417 -0.8080762 -0.9455624  0.2583927  0.3990654  0.1833369  0.1684278
#[9] -0.1635049 -0.5786871  0.7508694
# Likelihood value=848.0932

#### Exploration: likelihood-maximization of RR model ####


# RR model instead of ERR model, not successful
# beta_10*cumexp100 instead of (1+beta_10*cumexp100)
llk_rr<-function(b){
  L<-0
  for(i in 1:length(d_times)){
    # create death set: may have ties, up to 3 ties in each Death set
    D_set<-which(final_data_1$CASE==1&final_data_1$AGE_OUT==d_times[i])
    D_tienum<-length(D_set)
    D_age<-unique(final_data_1$AGE_TV[D_set])
    D_time<-unique(final_data_1$AGE_OUT[D_set])
    # create risk set: stratified by AGE_TV, survival past death time
    S_set<-which(final_data_1$AGE_TV==D_age&final_data_1$AGE_OUT>=D_time)
    # write log likelihood term adjusting for ties using Efron method
    tie_adjust<-mean(exp(bcvr[D_set,]%*%b[1:M])*(as.matrix(doses[D_set,])%*%b[(M + 1):(M + Nd)]))
    denom_adjust<-0
    for(r in 1:D_tienum){
      denom_adjust<-denom_adjust+log(sum(exp(bcvr[S_set,]%*%b[1:M])*(as.matrix(doses[S_set,])%*%b[(M + 1):(M + Nd)]))-(r-1)*tie_adjust)
    }
    L<-L+sum(bcvr[D_set,]%*%b[1:M])+sum(log(as.matrix(doses[D_set,])%*%b[(M + 1):(M + Nd)]))-denom_adjust
  }
  return(L)
}
# definition of constrains for RR parameters: beta > 0
A <- matrix(c(rep(0, M), 1), ncol=M + 1)
B <- 0
constraints <- list(ineqA = A, ineqB = B) # constrains, RR term must be positive

model.rr.pllk<-maxLik(logLik=llk_rr,start=c(as.numeric(glmfit.background$coefficients[2:10]),0.1),constraints=constraints)
summary(model.rr.pllk)

#### Exploration: likelihood-maximization of Cox model using optim() VS maxLik() ####

# check the maxLik function by writing cox-ph likelihood
llk_coxph<-function(b){
  L<-0
  for(i in 1:length(d_times)){
    # create death set: may have ties, up to 3 ties in each Death set
    D_set<-which(final_data_1$CASE==1&final_data_1$AGE_OUT==d_times[i])
    D_tienum<-length(D_set)
    D_age<-unique(final_data_1$AGE_TV[D_set])
    D_time<-unique(final_data_1$AGE_OUT[D_set])
    # create risk set: stratified by AGE_TV, survival past death time
    S_set<-which(final_data_1$AGE_TV==D_age&final_data_1$AGE_OUT>=D_time)
    # write log likelihood term adjusting for ties using Efron method
    tie_adjust<-mean(exp(bcvr[D_set,]%*%b[1:M]+as.matrix(doses[D_set,])%*%b[(M + 1):(M + Nd)]))
    denom_adjust<-0
    for(r in 1:D_tienum){
      denom_adjust<-denom_adjust+log(sum(exp(bcvr[S_set,]%*%b[1:M]+as.matrix(doses[S_set,])%*%b[(M + 1):(M + Nd)]))-(r-1)*tie_adjust)
    }
    L<-L+sum(bcvr[D_set,]%*%b[1:M])+sum(as.matrix(doses[D_set,])%*%b[(M + 1):(M + Nd)])-denom_adjust
  }
  return(L)
}
llk_coxph_neg<-function(b){
  L<-0
  for(i in 1:length(d_times)){
    # create death set: may have ties, up to 3 ties in each Death set
    D_set<-which(final_data_1$CASE==1&final_data_1$AGE_OUT==d_times[i])
    D_tienum<-length(D_set)
    D_age<-unique(final_data_1$AGE_TV[D_set])
    D_time<-unique(final_data_1$AGE_OUT[D_set])
    # create risk set: stratified by AGE_TV, survival past death time
    S_set<-which(final_data_1$AGE_TV==D_age&final_data_1$AGE_OUT>=D_time)
    # write log likelihood term adjusting for ties using Efron method
    tie_adjust<-mean(exp(bcvr[D_set,]%*%b[1:M]+as.matrix(doses[D_set,])%*%b[(M + 1):(M + Nd)]))
    denom_adjust<-0
    for(r in 1:D_tienum){
      denom_adjust<-denom_adjust+log(sum(exp(bcvr[S_set,]%*%b[1:M]+as.matrix(doses[S_set,])%*%b[(M + 1):(M + Nd)]))-(r-1)*tie_adjust)
    }
    L<-L+sum(bcvr[D_set,]%*%b[1:M])+sum(as.matrix(doses[D_set,])%*%b[(M + 1):(M + Nd)])-denom_adjust
  }
  return(-L)
}
model.coxph.optim<-optim(par=c(2,-1,-1,1,1,1,1,1,-1,0.2),fn=llk_coxph_neg,control=list(maxit=5000))
# Default is NM method
#[1]  1.72379606 -0.63024428 -0.98177563  0.57077837  0.74114755  0.80665403  0.74114661
#[8]  0.36324588 -0.08490365  0.29579247
# Likelihood value=848.5624
model.coxph.optim2<-optim(par=c(1,2,-1,-1,1,1,1,1,1,-1,0.2),fn=llk_coxph_neg,control=list(maxit=5000),method="BFGS")
# Need to specify an additional starting value
# BFGS method much closer to coxph output
# [1]  1.0000000  1.6711032 -0.8638338 -0.9538940  0.3015149  0.4851061  0.1401379  0.1052908
# [9] -0.2642538 -0.7088733  0.2976138

model.coxph.adjusted$coefficients
#cumexp100        yy2        yy3        yy4        yy5      SMOKE     native 
#0.2976242  0.1397019  0.1048294 -0.2646445 -0.7092440  1.6713391 -0.8613983 
#hispanic      COMP5      COMP6 
#-0.9537376  0.3015765  0.4851583

# definition of constrains for Cox PH RR parameters: none
#A <- matrix(c(rep(0, M), 1), ncol=M + 1)
#B <- 0
#constraints <- list(ineqA = A, ineqB = B) # constrains, RR term must be positive
model.coxph.pllk<-maxLik(logLik=llk_coxph,start=c(2,-1,-1,1,1,1,1,1,-1,0.2),control=list(iterlim=5000),method="BFGS")
# if no constraints, maxLik defaults to NR method, produces NA in gradient
# should specify other methods, choose "BFGS" based on stability
# [1]  1.6714001 -0.8607765 -0.9536891  0.3017926  0.4853464  0.1398611
# [7]  0.1049535 -0.2645066 -0.7091322  0.2976267
# NaNs produced, no std estimates for yy2-yy5

#### Exploration: rERR package ####

# https://github.com/fbr600/rERR
# Github rERR package trial
install.packages("devtools")
devtools::install_github("fbr600/rERR")

library(rERR)
formula.rERR<-Surv(AGE_STRT,AGE,CASE)~loglin(SMOKE,native,hispanic,COMP5,COMP6,yy2,yy3,yy4,yy5)+lin(dose_cum)
model.rERR<-f_fit_linERR_ef(formula.rERR,data=final_data_1,id_name="ID",dose_name="YREXP",time_name="AGE_TV",covars_names=c("SMOKE","native","hispanic","COMP5","COMP6","yy2","yy3","yy4","yy5","AGE_TV"),lag=0,exclusion_done=T)
#Error in nrows_cases_[i] <- dt$row[which(dt$id == id_cases[i])] : 
#  replacement has length zero

dt1<-f_to_event_table_ef_all(formula=formula.rERR,data=final_data_1,id_name="ID",dose_name="YREXP",
                                          time_name="AGE_TV",covars_names=c("SMOKE","native","hispanic","COMP5","COMP6","yy2","yy3","yy4","yy5","AGE_TV"))
model.rERR<-f_fit_linERR_all(formula.rERR,data=dt1,id_name="ID",time_name="AGE_TV",lag=0)
#Error in nrows_cases_[i] <- dt$row[which(dt$id == id_cases[i])] : 
#  replacement has length zero



#### What we most likely will use ####

# Next Step: Simulation based on real data
# Set beta_10=1
# Calculate lambda(age)=exp(beta1*SMOKE+...+beta9*yy5)*(1+beta10*cumexp100)
# Calculate f(age)=lambda(age)*exp(-sum lambda(age from entry to current-1?))
# Or calculate f(age) = lambda(age)*(1-sum f(age from entry to current-1))?
# Normalize f_1(age)=f(age)/sum f(age from entry to end) to real pdf


setwd("~/OneDrive - Harvard University/JUL2015-FASRC/analysis/")

library(haven)
crude_exp<-read_sas("crude_exp.sas7bdat")
# 88405 rows, 57 columns
#[1] "id"        "age_tv"    "predexpxm" "surtotxm"  "BIRTH_YY"  "LUNG_CA"  
#[7] "can"       "EXPSTRT"   "EXPEND"    "SMOKE"     "ETHNIC"    "age_new"  
#[13] "year"      "mnth"      "AGE_STRT"  "FYR"       "COMP5"     "COMP6"    
#[19] "COMP12"    "_NAME_"    "CUMEXP"    "YR_TV"     "YREXP"     "comp_m"   
#[25] "yrn"       "AGE"       "CEXP"      "CASE"      "y1"        "y2"       
#[31] "y3"        "y4"        "y5"        "y6"        "y7"        "LASTE"    
#[37] "type"      "type1"     "type2"     "type3"     "time"      "etype"    
#[43] "yr_last0"  "yr_last"   "surf"      "white"     "hispanic"  "native"   
#[49] "cumexp100" "cumexp500" "yy1"       "yy2"       "yy3"       "yy4"      
#[55] "yy5"       "casenew"   "i" 

crude_exp_1<-read_sas("crude_exp_1.sas7bdat")
mv_exp<-read_sas("mv_exp.sas7bdat")
vs<-read_sas("~/OneDrive - Harvard University/JUL2015-FASRC/data/valid_2009_short.sas7bdat")


library(survival)
# COXPH RRC crude match in point estimate, not in SE
# COXPH RRC adjusted doesn't match in both with draft results
model.coxph.rrc.crude<-coxph(Surv(time,CASE)~(predexpxm)+strata(age_tv),data=crude_exp)
summary(model.coxph.rrc.crude)
# 1.8165    6.1502   0.2261
model.coxph.rrc.crude_2<-coxph(Surv(time,CASE)~(predexpxm)+strata(age_tv),data=mv_exp)
summary(model.coxph.rrc.crude_2)
# 1.8084    6.1006   0.2265
model.coxph.rrc.adjusted<-coxph(Surv(time,CASE)~(predexpxm+SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5)+strata(age_tv),data=mv_exp)
summary(model.coxph.rrc.adjusted)
# 1.38327   3.98791  0.27532


mv_exp$age_out<-mv_exp$age_tv+mv_exp$time

glmfit<-glm(CASE~predexpxm+SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5,family="poisson",data=mv_exp)
glmfit$coefficients
#(Intercept)   predexpxm       SMOKE      native    hispanic       COMP5 
#-11.0714043   2.6304324   1.6249479  -1.2743255  -1.4272818   0.7349994 
#COMP6         yy2         yy3         yy4         yy5 
#1.1625064   0.9693951   1.9190164   2.3693374   2.5890539
glmfit.background<-glm(CASE~SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5,family="poisson",data=mv_exp)
glmfit.background$coefficients


# checkout distribution of 132 cases
mv_exp_cases<-mv_exp[which(mv_exp$CASE==1),c("id","age_tv","age_out","time","predexpxm","cumexp100","CASE","YR_TV")]
summary(mv_exp_cases)
length(unique(mv_exp_cases$age_tv))
# 39 unique age_tv among 132 cases
length(unique(mv_exp_cases$age_out))
# 113 unique age_out among 132 cases
View(mv_exp_cases)
# 1/2/3 ties at each distinct death ages/times 

# ERR model fit using partial likelihood approach
S <- sum(mv_exp$CASE) # total number of cases
#bcvr <- cbind(rep(1, nrow(mv_exp)), mv_exp$SMOKE,
#              mv_exp$native,mv_exp$hispanic,mv_exp$COMP5,
#              mv_exp$COMP6,mv_exp$yy2,mv_exp$yy3,mv_exp$yy4,mv_exp$yy5) # covariance matrix
bcvr <- cbind(mv_exp$SMOKE,mv_exp$native,mv_exp$hispanic,
              mv_exp$COMP5,mv_exp$COMP6,
              mv_exp$yy2,mv_exp$yy3,mv_exp$yy4,mv_exp$yy5) # covariance matrix
M <- ncol(bcvr) # number of covariates at the background risk term, no intercept
Sk <- as.numeric(mv_exp$CASE%*%bcvr) # inner product of cases with covariates
doses <- cbind(mv_exp$predexpxm) # covariates of the RR term
Nd <- ncol(doses) # number of covariates at the RR term

# partial log-likelihood term for each unique failure times
d_times<-unique(mv_exp_cases$age_out)
# define failure rate lambda(beta,i) function
#fr<-function(b,i){
#  rate<-exp(bcvr[i,]%*%b[1:M])*(1+doses[i,]%*%b[(M+1):(M+Nd)])
#  return(rate)
#}
#b<-c(2,-1,-1,1,1,1,1,-1,-1,0.5)

llk_rrc_err<-function(b){
  L<-0
  for(i in 1:length(d_times)){
    # create death set: may have ties, up to 3 ties in each Death set
    D_set<-which(mv_exp$CASE==1&mv_exp$age_out==d_times[i])
    D_tienum<-length(D_set)
    D_age<-unique(mv_exp$age_tv[D_set])
    D_time<-unique(mv_exp$age_out[D_set])
    # create risk set: stratified by AGE_TV, survival past death time
    S_set<-which(mv_exp$age_tv==D_age&mv_exp$age_out>=D_time)
    # write log likelihood term adjusting for ties using Efron method
    tie_adjust<-mean(exp(bcvr[D_set,]%*%b[1:M])*(1+as.matrix(doses[D_set,])%*%b[(M + 1):(M + Nd)]))
    denom_adjust<-0
    for(r in 1:D_tienum){
      denom_adjust<-denom_adjust+log(sum(exp(bcvr[S_set,]%*%b[1:M])*(1+as.matrix(doses[S_set,])%*%b[(M + 1):(M + Nd)]))-(r-1)*tie_adjust)
    }
    L<-L+sum(bcvr[D_set,]%*%b[1:M])+sum(log(1+as.matrix(doses[D_set,])%*%b[(M + 1):(M + Nd)]))-denom_adjust
  }
  return(L)
}
# definition of constrains for ERR parameters: beta > -1/max(D)
A <- matrix(c(rep(0, M), 1), ncol=M + 1)
B <- 1/max(mv_exp$predexpxm)
constraints <- list(ineqA = A, ineqB = B) # constrains, RR term must be positive

library(maxLik)
model.rrc.err.pllk<-maxLik(logLik=llk_rrc_err,start=c(as.numeric(model.coxph.rrc.adjusted$coefficients[2:10]),1),constraints=constraints,control=list(iterlim=5000))
summary(model.rrc.err.pllk)
# convergence after 803 iterations
# NaN produced in SE estimates for yy2-4
model.rrc.err.pllk$estimate
#[1]  1.65112510 -0.83529339 -1.04232597  0.59709344  0.71388323  0.34886506
#[7]  0.36396283  0.05321453 -0.40036750  2.75318073
# Likelihood value=850.1801
model.rrc.err.pllk_1<-maxLik(logLik=llk_rrc_err,start=c(as.numeric(glmfit.background$coefficients[2:10]),1),constraints=constraints,control=list(iterlim=5000))
#[1]  1.8545281 -0.6950775 -1.1030993  0.3188755  0.4529549  1.2507354
#[7]  1.3061497  1.0440657  0.5583048  3.1817539
# Likelihood value=851.8329
# convergence after 1893 iterations

# definition of constrains for ERR parameters: beta > -1/max(D)
A <- matrix(c(rep(0, M+1), 1), ncol=M + 2)
B <- 1/max(mv_exp$predexpxm)
constraints <- list(ineqA = A, ineqB = B)
model.rrc.err.pllk2<-maxLik(logLik=llk_rrc_err,start=c(1,as.numeric(model.coxph.rrc.adjusted$coefficients[2:10]),1),constraints=constraints,control=list(iterlim=5000),method="BFGS")
summary(model.rrc.err.pllk2)
# convergence after 37 iterations
# Std all Inf
model.rrc.err.pllk2$estimate
#[1]  1.6440122 -0.8923009 -1.0417609  0.5603437  0.6826110  0.1637983
#[7]  0.1773501 -0.1384405 -0.5875772  2.8331062  1.0000000
# Likelihood value=850.1238
model.rrc.err.pllk2_1<-maxLik(logLik=llk_rrc_err,start=c(1,as.numeric(glmfit.background$coefficients[2:10]),1),constraints=constraints,control=list(iterlim=5000),method="BFGS")
#[1]  1.6439931 -0.8925359 -1.0417179  0.5604210  0.6826383  0.1636410
#[7]  0.1771389 -0.1386723 -0.5878630  2.8323667  1.0000000
# Likelihood value=850.1238
# convergence after 34 iterations

llk_rrc_err_neg<-function(b){
  L<-0
  for(i in 1:length(d_times)){
    # create death set: may have ties, up to 3 ties in each Death set
    D_set<-which(mv_exp$CASE==1&mv_exp$age_out==d_times[i])
    D_tienum<-length(D_set)
    D_age<-unique(mv_exp$age_tv[D_set])
    D_time<-unique(mv_exp$age_out[D_set])
    # create risk set: stratified by AGE_TV, survival past death time
    S_set<-which(mv_exp$age_tv==D_age&mv_exp$age_out>=D_time)
    # write log likelihood term adjusting for ties using Efron method
    tie_adjust<-mean(exp(bcvr[D_set,]%*%b[1:M])*(1+as.matrix(doses[D_set,])%*%b[(M + 1):(M + Nd)]))
    denom_adjust<-0
    for(r in 1:D_tienum){
      denom_adjust<-denom_adjust+log(sum(exp(bcvr[S_set,]%*%b[1:M])*(1+as.matrix(doses[S_set,])%*%b[(M + 1):(M + Nd)]))-(r-1)*tie_adjust)
    }
    L<-L+sum(bcvr[D_set,]%*%b[1:M])+sum(log(1+as.matrix(doses[D_set,])%*%b[(M + 1):(M + Nd)]))-denom_adjust
  }
  return(-L)
}
model.rrc.err.optim<-optim(par=c(as.numeric(model.coxph.rrc.adjusted$coefficients[2:10]),1),fn=llk_rrc_err_neg,control=list(maxit=5000))
#[1]  1.65112510 -0.83529339 -1.04232597  0.59709344  0.71388323  0.34886506
#[7]  0.36396283  0.05321453 -0.40036750  2.75318073
#Likelihood value=850.1801
model.rrc.err.optim2<-optim(par=c(1,as.numeric(model.coxph.rrc.adjusted$coefficients[2:10]),1),fn=llk_rrc_err_neg,control=list(maxit=5000),method="BFGS")
#[1]  1.6440124 -0.8923007 -1.0417609  0.5603435  0.6826109  0.1637967
#[7]  0.1773484 -0.1384422 -0.5875788  2.8331063  1.0000000
#Likelihood value=850.1238

model.coxph.rrc.adjusted$coefficients
#predexpxm       SMOKE      native    hispanic       COMP5       COMP6 
#1.38326702  1.63554065 -0.88283684 -1.03435642  0.51738618  0.71297533 
#yy2         yy3         yy4         yy5 
#0.12275396  0.20207523 -0.08915978 -0.53429846

