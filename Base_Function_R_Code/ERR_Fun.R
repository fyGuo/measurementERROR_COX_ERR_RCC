#### ERR Fitting Function ####

ERR<-function(data,outcome,covariates,exposure,times,strata,start){
  # case/death vector
  case <- as.matrix(data[,outcome])
  # covariance matrix
  bcvr <- as.matrix(data[,covariates])
  # number of covariates at the background risk term, no intercept
  if(is.null(covariates)==TRUE){
    M <- 0
  }
  else{
    M <- ncol(bcvr)
  }
  # covariates of the RR term
  doses <- as.matrix(data[,exposure]) 
  # number of covariates at the RR term
  Nd <- ncol(doses) 
  # unique failure times
  timein<-times[1]
  timeout<-times[2]
  d_times<-unique(as.matrix(data[which(case==1),timeout]))
  
  # define log likelihood function
  llk_err<-function(b){
    L<-0
    for(i in 1:length(d_times)){
      # create death set: may have ties
      D_set<-which(as.matrix(data[,outcome])==1&as.matrix(data[,timeout])==d_times[i])
      D_tienum<-length(D_set)
      D_age<-as.numeric(unique(as.matrix(data[D_set,timein])))
      D_time<-as.numeric(unique(as.matrix(data[D_set,timeout])))
      # create risk set: stratified by AGE_TV, survival past death time
      S_set<-which(as.matrix(data[,strata])==D_age&as.matrix(data[,timeout])>=D_time)
      # write log likelihood term adjusting for ties using Breslow method
      if(is.null(covariates)==TRUE){
        L<-L+sum(log(1+as.matrix(doses[D_set,])%*%b[(M + 1):(M + Nd)]))-D_tienum*log(sum((1+as.matrix(doses[S_set,])%*%b[(M + 1):(M + Nd)])))
      }
      else{
        L<-L+sum(bcvr[D_set,]%*%b[1:M])+sum(log(1+as.matrix(doses[D_set,])%*%b[(M + 1):(M + Nd)]))-D_tienum*log(sum(exp(bcvr[S_set,]%*%b[1:M])*(1+as.matrix(doses[S_set,])%*%b[(M + 1):(M + Nd)])))
      }
    }
    return(L)
  }
  
  # definition of constrains for ERR parameters: beta > -1/max(D)
  A <- matrix(c(rep(0, M), 1), ncol=M + 1)
  B <- 1/max(doses)
  constraints <- list(ineqA = A, ineqB = B) # constrains, RR term must be positive
  
  # optimize log likelihood function
  library(maxLik)
  model.err<-maxLik(logLik=llk_err,start=start,constraints=constraints,control=list(iterlim=5000),method="BFGS")
  
  return(model.err)
}

#### Trial on Uminer Data ####

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

library(survival)
model.coxph.adjusted<-coxph(Surv(time,CASE)~(cumexp100+SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5)+strata(AGE_TV),data=final_data_1)
model.coxph.adjusted.b<-coxph(Surv(time,CASE)~(cumexp100+SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5)+strata(AGE_TV),data=final_data_1,ties="breslow")


model.err.adjusted<-ERR_breslow(data=final_data_1,outcome=c("CASE"),
                        covariates=c("SMOKE","native","hispanic","COMP5","COMP6","yy2","yy3","yy4","yy5"),
                        exposure=c("cumexp100"),times=c("AGE_TV","AGE_OUT"),strata=c("AGE_TV"),
                        start=c(as.numeric(model.coxph.adjusted.b$coefficients[2:10]),1))

summary(model.err.adjusted)
0.7504

model.err.unadjusted<-ERR_breslow(data=final_data_1,outcome=c("CASE"),
                                covariates=c(),
                                exposure=c("cumexp100"),times=c("AGE_TV","AGE_OUT"),strata=c("AGE_TV"),
                                start=c(1))

#### ERR Sandwich Covariance Matrix Function ####

ERR_sandwichcov<-function(data,outcome,covariates,exposure,times,strata,beta){
  # case/death vector
  case <- as.matrix(data[,outcome])
  # covariance matrix
  bcvr <- as.matrix(data[,covariates])
  # number of covariates at the background risk term, no intercept
  if(is.null(covariates)==TRUE){
    M <- 0
  }
  else{
    M <- ncol(bcvr)
  }
  # covariates of the RR term
  doses <- as.matrix(data[,exposure]) 
  # number of covariates at the RR term
  Nd <- ncol(doses) 
  # unique failure times
  timein<-times[1]
  timeout<-times[2]
  d_times<-unique(as.matrix(data[which(case==1),timeout]))
  
  b<-beta
  #b<-as.numeric(model.err.adjusted$estimate)
  if(is.null(covariates)==TRUE){
    # calculate risk function: 88405*1
    risk<-(1+doses%*%b[(M + 1):(M + Nd)])
    # calculate doses adjustment weight
    adjdoses<-doses/(1+doses%*%b[(M + 1):(M + Nd)])
    # combine with covariates matrix: 88405*10
    design<-matrix(adjdoses,nrow=nrow(data),ncol=Nd)
  }
  else{
    # calculate risk function: 88405*1
    risk<-exp(bcvr%*%b[1:M])*(1+doses%*%b[(M + 1):(M + Nd)])
    # calculate doses adjustment weight
    adjdoses<-doses/(1+doses%*%b[(M + 1):(M + Nd)])
    # combine with covariates matrix: 88405*10
    design<-cbind(bcvr,adjdoses)
  }
  
  # initialize first derivative G squared matrix (M+Nd)*(M+Nd)
  G2<-matrix(0,nrow=M+Nd,ncol=M+Nd)
  #G2<-0
  # initialize second derivative H matrix (M+Nd)*(M+Nd)
  H<-matrix(0,nrow=M+Nd,ncol=M+Nd)
  #H<-0
  
  #for(i in 1:length(d_times)){
  #  # create death set: may have ties
  #  D_set<-which(as.matrix(data[,outcome])==1&as.matrix(data[,timeout])==d_times[i])
  #  D_tienum<-length(D_set)
  #  D_age<-as.numeric(unique(as.matrix(data[D_set,timein])))
  #  D_time<-as.numeric(unique(as.matrix(data[D_set,timeout])))
  #  # create risk set: stratified by AGE_TV, survival past death time
  #  S_set<-which(as.matrix(data[,strata])==D_age&as.matrix(data[,timeout])>=D_time)
  #  # Gi first derivative vector (M+Nd)*1
  #  Gi<-colSums(matrix(design[D_set,],nrow=D_tienum,ncol=M+Nd))-D_tienum*sum(design[S_set,]*risk[S_set])/sum(risk[S_set])
  #  # Gi-squared first derivative matrix (M+Nd)*(M+Nd)
  #  Gi2<-Gi^2
  #  G2<-G2+Gi2
  #  # Hi second derivative matrix (M+Nd)*(M+Nd)
  #  Hi<--sum(adjdoses[D_set]^2)+D_tienum*(sum(adjdoses[S_set]*risk[S_set])^2)/(sum(risk[S_set])^2)
  #  H<-H+Hi
  #}
  
  # loop through unique death times
  for(i in 1:length(d_times)){
    # create death set: may have ties
    D_set<-which(as.matrix(data[,outcome])==1&as.matrix(data[,timeout])==d_times[i])
    D_tienum<-length(D_set)
    D_age<-as.numeric(unique(as.matrix(data[D_set,timein])))
    D_time<-as.numeric(unique(as.matrix(data[D_set,timeout])))
    # create risk set: stratified by AGE_TV, survival past death time
    S_set<-which(as.matrix(data[,strata])==D_age&as.matrix(data[,timeout])>=D_time)
    # Gi first derivative vector (M+Nd)*1
    Gi<-t(t(colSums(matrix(design[D_set,],nrow=D_tienum,ncol=M+Nd))))-D_tienum*as.matrix(colSums(design[S_set,]*risk[S_set]))/sum(risk[S_set])
    # Gi-squared first derivative matrix (M+Nd)*(M+Nd)
    Gi2<-Gi%*%t(Gi)
    G2<-G2+Gi2
    # Hi second derivative matrix (M+Nd)*(M+Nd)
    Hi<--D_tienum*(sum(risk[S_set])*t(design[S_set,]*risk[S_set])%*%design[S_set,]-as.matrix(colSums(design[S_set,]*risk[S_set]))%*%colSums(design[S_set,]*risk[S_set]))/(sum(risk[S_set])^2)
    Hi[M+Nd,M+Nd]<--sum(adjdoses[D_set]^2)+D_tienum*(sum(adjdoses[S_set]*risk[S_set])^2)/(sum(risk[S_set])^2)
    H<-H+Hi
  }
  
  # covaraince matrix
  COV_sandwich1<-solve(-H)%*%G2%*%solve(-H)
  COV_sandwich2<-solve(G2)
  return(list("COV_strict"=COV_sandwich1,"COV_approx"=COV_sandwich2))
}

model.err.adjusted.cov<-ERR_sandwichcov(data=final_data_1,outcome=c("CASE"),
                        covariates=c("SMOKE","native","hispanic","COMP5","COMP6","yy2","yy3","yy4","yy5"),
                        exposure=c("cumexp100"),times=c("AGE_TV","AGE_OUT"),strata=c("AGE_TV"),
                        beta=as.numeric(model.err.adjusted$estimate))


model.err.naive<-ERR(data=final_data_1,outcome=c("CASE"),
                        covariates=NULL,
                        exposure=c("cumexp100"),times=c("AGE_TV","AGE_OUT"),strata=c("AGE_TV"),
                        start=c(1))
summary(model.err.naive)
# 1.6085
# Likelihood = 880.4163
model.err.naive.cov<-ERR_sandwichcov(data=final_data_1,outcome=c("CASE"),
                     covariates=NULL,
                     exposure=c("cumexp100"),times=c("AGE_TV","AGE_OUT"),strata=c("AGE_TV"),
                     beta=as.numeric(model.err.naive$estimate))
# the function doesn't work in this case
# run line by line and obtain result
# -H=3.273648, G=3.256194
# Cov_sandwich_strict=0.303841
# Cov_sandwich_approximate=0.307107




mv_exp<-read_sas("~/OneDrive - Harvard University/JUL2015-FASRC/analysis/mv_exp.sas7bdat")
model.coxph.rrc.adjusted<-coxph(Surv(time,CASE)~(predexpxm+SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5)+strata(age_tv),data=mv_exp)

model.coxph.rrc.adjusted.b<-coxph(Surv(time,CASE)~(predexpxm+SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5)+strata(age_tv),data=mv_exp,ties="breslow")
model.coxph.rrc.adjusted.b$coefficients
#predexpxm       SMOKE      native    hispanic       COMP5       COMP6         yy2         yy3         yy4 
#1.38250767  1.63518121 -0.88252998 -1.03405730  0.51738098  0.71250150  0.12274896  0.20187410 -0.08966401 
#yy5 
#-0.53436790

mv_exp$age_out<-mv_exp$age_tv+mv_exp$time

model.err.rrc.adjusted<-ERR(data=mv_exp,outcome=c("CASE"),
                        covariates=c("SMOKE","native","hispanic","COMP5","COMP6","yy2","yy3","yy4","yy5"),
                        exposure=c("predexpxm"),times=c("age_tv","age_out"),strata=c("age_tv"),
                        start=c(as.numeric(model.coxph.adjusted.b$coefficients[2:10]),2))

summary(model.err.rrc.adjusted)
model.err.rrc.adjusted.cov<-ERR_sandwichcov(data=mv_exp,outcome=c("CASE"),
                                            covariates=c("SMOKE","native","hispanic","COMP5","COMP6","yy2","yy3","yy4","yy5"),
                                            exposure=c("predexpxm"),times=c("age_tv","age_out"),strata=c("age_tv"),
                                            beta=c(1.644,-0.893,-1.042,0.56,0.683,0.164,0.177,-0.139,-0.588,2.832))
#     SMOKE     native   hispanic      COMP5      COMP6        yy2        yy3        yy4        yy5 
#0.19528403 0.52976243 0.05918627 0.23173766 0.22124143 0.41052025 0.34154175 0.36015884 0.37563174 
#predexpxm 
#1.58675140 

#     SMOKE     native   hispanic      COMP5      COMP6        yy2        yy3        yy4        yy5 
#0.12473502 0.52340626 0.05861613 0.24231063 0.22043096 0.43513190 0.45577376 0.46272597 0.48496670 
#predexpxm 
#1.56541669 

#### CoxPH Sandwich Covariance Matrix Function ####

CoxPH.sandwichcov<-function(data,outcome,covariates,exposure,times,strata,beta){
  # case/death vector
  case <- as.matrix(data[,outcome])
  # covariance matrix
  bcvr <- as.matrix(data[,covariates])
  # number of covariates at the background risk term, no intercept
  M <- ncol(bcvr) 
  # covariates of the RR term
  doses <- as.matrix(data[,exposure]) 
  # number of covariates at the RR term
  Nd <- ncol(doses) 
  # unique failure times
  timein<-times[1]
  timeout<-times[2]
  d_times<-unique(as.matrix(data[which(case==1),timeout]))
  
  b<-as.numeric(model.coxph.adjusted.b$coefficients)
  # combine with covariates matrix: 88405*10
  design<-cbind(bcvr,doses)
  # calculate risk function: 88405*1
  risk<-exp(design%*%b)

  
  # initialize first derivative G squared matrix (M+Nd)*(M+Nd)
  G2<-matrix(0,nrow=M+Nd,ncol=M+Nd)
  # initialize second derivative H matrix (M+Nd)*(M+Nd)
  H<-matrix(0,nrow=M+Nd,ncol=M+Nd)
  
  # loop through unique death times
  for(i in 1:length(d_times)){
    # create death set: may have ties
    D_set<-which(as.matrix(data[,outcome])==1&as.matrix(data[,timeout])==d_times[i])
    D_tienum<-length(D_set)
    D_age<-as.numeric(unique(as.matrix(data[D_set,timein])))
    D_time<-as.numeric(unique(as.matrix(data[D_set,timeout])))
    # create risk set: stratified by AGE_TV, survival past death time
    S_set<-which(as.matrix(data[,strata])==D_age&as.matrix(data[,timeout])>=D_time)
    # Gi first derivative vector (M+Nd)*1
    Gi<-t(t(colSums(matrix(design[D_set,],nrow=D_tienum,ncol=M+Nd))))-D_tienum*as.matrix(colSums(design[S_set,]*risk[S_set]))/sum(risk[S_set])
    # Gi-squared first derivative matrix (M+Nd)*(M+Nd)
    Gi2<-Gi%*%t(Gi)
    G2<-G2+Gi2
    # Hi second derivative matrix (M+Nd)*(M+Nd)
    Hi<--D_tienum*(sum(risk[S_set])*t(design[S_set,]*risk[S_set])%*%design[S_set,]-as.matrix(colSums(design[S_set,]*risk[S_set]))%*%colSums(design[S_set,]*risk[S_set]))/(sum(risk[S_set])^2)
    H<-H+Hi
  }
  
  # covaraince matrix
  COV_sandwich1.coxph<-solve(-H)%*%G2%*%solve(-H)
  COV_sandwich2.coxph<-solve(G2)
}

# doesn't match model.coxph.adjusted.b
# explore further 
