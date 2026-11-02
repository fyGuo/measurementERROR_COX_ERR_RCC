#### Explore data features of Uminer main study and validation study data ####
# to conduct simulation mimicking the nature of available Uminer data

# Estimate features of YREXP distribution from main study
library(haven)

final_data<-read_sas("~/OneDrive - Harvard University/JUL2015-FASRC/data/final_data.sas7bdat")

final_data$YR_LAST<-as.numeric(final_data$year)-1900
final_data_1<-final_data[which((is.na(final_data$EXPSTRT)==FALSE&final_data$YR_TV>=final_data$EXPSTRT&final_data$YR_TV<=final_data$YR_LAST)|(is.na(final_data$EXPSTRT)==TRUE&final_data$CUMEXP!=0)),]
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
final_data_1$AGE_OUT<-final_data_1$AGE_TV+final_data_1$time


valid_data<-read_sas("~/OneDrive - Harvard University/JUL2015-FASRC/data/valid_2009_short.sas7bdat")
# VALUE is error-prone YREXP from final_data
# EVALUE is the true YREXP 
valid_data$VALUE100<-valid_data$VALUE/100
valid_data$EVALUE100<-valid_data$EVALUE/100
valid_data$yy1<-ifelse(valid_data$YEAR<70,1,0)
valid_data$yy2<-ifelse(valid_data$YEAR>=70&valid_data$YEAR<80,1,0)
valid_data$yy3<-ifelse(valid_data$YEAR>=80&valid_data$YEAR<90,1,0)
valid_data$yy4<-ifelse(valid_data$YEAR>=90&valid_data$YEAR<100,1,0)
valid_data$yy5<-ifelse(valid_data$YEAR>=100,1,0)

#
unique_death_years<-sort(unique(round(final_data_1$AGE_TV[final_data_1$CASE==1],digits=0)))
RS_size<-rep(NA,length(unique_death_years))
RS_size_1<-rep(NA,length(unique_death_years))

for(i in 1:length(unique_death_years)){
  RS_size[i]<-sum(valid_data$AGE_TV==unique_death_years[i])
  RS_size_1[i]<-sum(valid_data$AGE_TV>=unique_death_years[i])
}
#Based on RS_size, 19 with an additional RS=1 at the end

unique_death_months<-sort(unique(final_data_1$AGE_OUT[final_data_1$CASE==1]))
RS_size_2<-rep(NA,length(unique_death_months))

for(i in 1:length(unique_death_months)){
  RS_size_2[i]<-sum(valid_data$AGE_TV>=unique_death_months[i]&valid_data$AGE_TV<1+unique_death_months[i])
}
#Based on RS_size_2, 49 with an additional RS=1 at the end



# Create variable RETIRETIME, RETIREYEAR, RETIREAGE
# the number of years (row number) after which a miner retired from the mining industry
# match to the calendar year and the miner's age
# selected by continuous YREXP=0 after this year and age
# investigate the distribution of RETIRETIME, RETIREYEAR, RETIREAGE

library(dplyr)

final_data_1<-final_data_1%>%
  group_by(ID)%>%
  mutate(RETIRETIME=cumsum(rle(YREXP)$lengths)[length(rle(YREXP)$lengths)-1]+1,
         RETIREYEAR=YR_TV[RETIRETIME],
         RETIREAGE=AGE_TV[RETIRETIME])
summary(final_data_1$RETIRETIME)
hist(final_data_1$RETIRETIME)

final_data_retire<-final_data_1%>%
  select(ID,BIRTH_YY,EXPSTRT,EXPEND,AGE_STRT,AGE,RETIRETIME,RETIREYEAR,RETIREAGE)%>%
  group_by(ID)%>%
  distinct()
summary(final_data_retire$RETIRETIME)
hist(final_data_retire$RETIRETIME)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.000   4.000   8.000   9.863  14.000  30.000
summary(final_data_retire$RETIREYEAR)
hist(final_data_retire$RETIREYEAR)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#59.00   70.00   77.00   75.55   82.00   86.00
summary(final_data_retire$RETIREAGE)
hist(final_data_retire$RETIREAGE)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#14.00   29.00   35.00   36.88   43.00   79.00
summary(final_data_retire$AGE_STRT)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#11.00   22.00   26.00   28.02   32.00   65.00
hist(final_data_retire$AGE_STRT)

final_data_beforeretire<-final_data_1[1,]
for(i in 1:nrow(final_data_retire)){
  id<-final_data_retire$ID[i]
  retiretime<-final_data_retire$RETIRETIME[i]
  final_data_beforeretire<-rbind(final_data_beforeretire,final_data_1[which(final_data_1$ID==id)[1:(retiretime-1)],])
}
final_data_beforeretire<-final_data_beforeretire[-1,]
# number of rows = 20625, which matches
# sum(final_data_retire$RETIRETIME)-nrow(final_data_retire)
# [1] 20625


# Calculate the average intraclass correlation rho_{CS} for YREXP
# Convert into correlation rho_{AR} for YREXP
# Generate YREXP before RETIRETIME as YREXP ~ MVN (mu_C, Sigma_C), for t*=60
library(nlme)
model.AR1<-gls(YREXP~1,data=final_data_beforeretire,
               correlation=corAR1(form=~1|ID),method="REML")
#mean(final_data_beforeretire$YREXP)
#[1] 7.912407
YREXP.mean<-as.numeric(model.AR1$coefficients)
# 7.882844
YREXP.sd<-model.AR1$sigma
# 17.06189
YREXP.arcor<-as.numeric(coef(model.AR1$modelStruct$corStruct,unconstrained=FALSE))
# 0.6155869

# Goal: obtain the correct covariance matrix from GLS summary output
# R GLS summary output does not give covariance matrix directly
# Instead it gives correlation, SD and weights which require additional calculation
# There is a built-in function getVarCov used with gls to obtain the covariance matrix
# However it has a bug and we recommend another way
# We provide you with the following getVarCov.fix function and you can use it directly
# It will give the correct covariance matrix to be used in Homeworks

#Do NOT change this function, just run it so we define a new function
getVarCov.fix <- function (obj, individual = 1, ...) {
  S <- corMatrix(obj$modelStruct$corStruct)[[individual]]
  if (!is.null(obj$modelStruct$varStruct)) {
    ## ind <- obj$groups == individual                                      # index from dataset
    dimensions.Ri <- sapply(corMatrix(obj$modelStruct$corStruct), nrow)     # dimentions from corMatrix()
    ind <- rep(seq_along(dimensions.Ri), dimensions.Ri) == individual       # index from dimensions
    vw <- 1/varWeights(obj$modelStruct$varStruct)[ind]
  }
  else vw <- rep(1, nrow(S))
  vars <- (obj$sigma * vw)^2
  result <- t(S * sqrt(vars)) * sqrt(vars)
  class(result) <- c("marginal", "VarCov")
  attr(result, "group.levels") <- names(obj$groups)
  result
}

#Can change model input here, do not change individual input
getVarCov.fix(model.AR1, individual = 1) 


# Create variable AGE_START
# uniform discrete according to empirical AGE_START

# Generate sim-VALUE|sim-eVALUE
# Estimate VALUE=alpha0+alpha1*EVALUE+e in validation data
# alpha0, alpha1, var(e)

# 
error.model.LR<-lm(EVALUE~VALUE,data=valid_data)
error.model.CS<-gls(EVALUE~VALUE,data=valid_data,
                    correlation=corCompSymm(form=~1|ID),method="REML")
error.model.AR1<-gls(EVALUE~VALUE,data=valid_data,
                    correlation=corAR1(form=~1|ID),method="REML")
anova(error.model.CS,error.model.AR1)
AIC(error.model.LR)
#[1] 15309.18

# choose error.model.AR1

# Correlation structure: damped exponential covariance structure
# Reference: Munoz Carey Rosner et al, Biometrics 1992
install.packages("RDEC", repos="http://R-Forge.R-project.org")



x_median<-median(final_data_1$cumexp100)
x_mean<-mean(final_data_1$cumexp100)
beta_cox<-c(log(1.5),log(5),log(10))
alpha_err_median1<-(exp(beta_cox)-1)/(1-x_median*(exp(beta_cox)-1))
alpha_err_mean1<-(exp(beta_cox)-1)/(1-x_mean*(exp(beta_cox)-1))
alpha_err_median2<-(exp(beta_cox*x_median)-1)/x_median
alpha_err_mean2<-(exp(beta_cox*x_mean)-1)/x_mean

gap<-c(0:55)
cor_dex<-0.64^(gap^0.66)
cor_ar1<-0.62^gap
cor_ar1-cor_dex



