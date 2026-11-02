#### Load functions ERR_breslow(), ERR_sandwichcov(), gen_data() ####

library(maxLik)
library(MASS)
library(sampling)
#library(dplyr)
library(survival)


ERR_breslow<-function(data,outcome,covariates,exposure,times,strata,start){
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
  d_times<-which(case==1)
  
  # define log likelihood function
  llk_err<-function(b){
    L<-0
    for(i in 1:length(d_times)){
      # case i in d_times
      D_set<-as.numeric(d_times[i])
      D_age<-as.numeric(data[D_set,timein])
      D_time<-as.numeric(data[D_set,timeout])
      # create risk set: stratified by AGE_TV, survival past death time
      if(is.null(strata)==TRUE){
        S_set<-which(as.matrix(data[,timeout])>=D_time)
      }
      else{
        S_set<-which(as.matrix(data[,strata])==D_age&as.matrix(data[,timeout])>=D_time)
      }
      # write log likelihood term adjusting for ties using Breslow method
      if(is.null(covariates)==TRUE){
        L<-L+log(1+as.matrix(doses[D_set,])%*%b[(M + 1):(M + Nd)])-log(sum((1+as.matrix(doses[S_set,])%*%b[(M + 1):(M + Nd)])))
      }
      else{
        L<-L+bcvr[D_set,]%*%b[1:M]+log(1+as.matrix(doses[D_set,])%*%b[(M + 1):(M + Nd)])-log(sum(exp(bcvr[S_set,]%*%b[1:M])*(1+as.matrix(doses[S_set,])%*%b[(M + 1):(M + Nd)])))
      }
    }
    return(L)
  }
  
  # definition of constrains for ERR parameters: 1+beta*doses>0, beta > -1/max(D)
  A <- matrix(c(rep(0, M), 1), ncol=M + 1)
  B <- 1/max(doses)
  constraints <- list(ineqA = A, ineqB = B) # constrains, RR term must be positive
  
  # optimize log likelihood function
  library(maxLik)
  model.err<-maxLik(logLik=llk_err,start=start,constraints=constraints,control=list(iterlim=5000),method="BFGS")
  
  return(model.err)
}

ERR_sandwichcov<-function(data,outcome,covariates,exposure,times,strata,beta){
  data=data.main
  outcome=c("D")
  covariates=c()
  exposure=c("x")
  times=c("time","time_out")
  strata=c("time")
  beta=model.err.true$estimate
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
  
  #for(i in 1:length(d_times)){
  #  # create death set: may have ties
  #  D_set<-which(as.matrix(data[,outcome])==1&as.matrix(data[,timeout])==d_times[i])
  #  D_tienum<-length(D_set)
  #  D_age<-as.numeric(unique(as.matrix(data[D_set,timein])))
  #  D_time<-as.numeric(unique(as.matrix(data[D_set,timeout])))
  #  # create risk set: stratified by AGE_TV, survival past death time
  #  S_set<-which(as.matrix(data[,strata])==D_age&as.matrix(data[,timeout])>=D_time)
  #  # Gi first derivative vector (M+Nd)*1
  #  Gi<-t(t(colSums(matrix(design[D_set,],nrow=D_tienum,ncol=M+Nd))))-D_tienum*as.matrix(sum(design[S_set,]*risk[S_set]))/sum(risk[S_set])
  #  # Gi-squared first derivative matrix (M+Nd)*(M+Nd)
  #  Gi2<-Gi%*%t(Gi)
  #  G2<-G2+Gi2
  #  # Hi second derivative matrix (M+Nd)*(M+Nd)
  #  Hi<--D_tienum*(sum(risk[S_set])*t(design[S_set,]*risk[S_set])%*%design[S_set,]-as.matrix(sum(design[S_set,]*risk[S_set]))%*%sum(design[S_set,]*risk[S_set]))/(sum(risk[S_set])^2)
  #  Hi[M+Nd,M+Nd]<--sum(adjdoses[D_set]^2)+D_tienum*(sum(adjdoses[S_set]*risk[S_set])^2)/(sum(risk[S_set])^2)
  #  H<-H+Hi
  #}
  
  # covaraince matrix
  COV_sandwich1<-solve(-H)%*%G2%*%solve(-H)
  COV_sandwich2<-solve(G2)
  return(list("COV_strict"=COV_sandwich1,"COV_approx"=COV_sandwich2))
}

# In final_data, mean(yrexp100)=0.01840276, var(yrexp100)=0.007875351, C
# In valid_data, mean(value100)=0.01893011, var(value100)=0.0004642924, c

## Generate AR(1) correlation matrix
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}

## Generate data
gen_data=function(n,m,rho,beta,t_star,nu,theta){
  # n: main study sample size
  # m: validation study sample size
  # alpha: measurement error model coef
  # beta: ERR effect model coef
  # t_star: administrative study period
  # nu,theta: parameters of baseline hazard weibull distribution
  
  #n<-2000
  #m<-500
  #alpha<-c(0.001,0.9)
  #beta<-2
  #t_star<-20
  #nu<-0.03
  #theta<-4
  
  # generate error-prone yearly exposure C_{ij}, i indicates subject 1-n, j indicates time 1-t_star
  # AR(1) correlation matrix
  Sigma_C<-0.002*ar1_cor(n=t_star,rho=0.996)
  # independent correlation matrix
  #Sigma_C<-0.002*diag(1,t_star)
  C<-mvrnorm(n,mu=rep(0.2,t_star),Sigma_C)
  # generate true yearly exposure c_{ij}, i, j as above
  #c<-alpha[1]+alpha[2]*C+mvrnorm(n,rep(0,t_star),diag(rep(1/17*(1-0.7^2)*var(as.numeric(C)),t_star)))
  alpha1<-0.01893011-mean(C)*rho/sqrt(17)
  alpha2<-rho/sqrt(17)
  c<-alpha1+alpha2*C+mvrnorm(n,rep(0,t_star),(1-rho^2)/17*Sigma_C)
  # create data frame with columns: id, t, C, c
  data.generate<-data.frame(id=rep(c(1:n),t_star),t=rep(c(1:t_star),each=n), C=as.vector(C),c=as.vector(c))
  # calculate cumulative total error-prone exposure X
  data.generate$X<-ave(data.generate$C, data.generate$id, FUN=cumsum)
  # calculate cumulative total true exposure x
  data.generate$x<-ave(data.generate$c, data.generate$id, FUN=cumsum)
  # generate starting age variable agestart, changing age variable time
  agestart<-sample(c(10:60),size=n,replace=TRUE)
  data.generate$age<-rep(agestart,t_star)+data.generate$t-1
  # generate censoring time V~exp(0.01)
  # ignore censorying for now
  #V<-rexp(n,rate=0.01)
  #V[V>t_star]=t_star
  # calculate cumulative incidence function F_{ij}, i, j as above
  #nu<-0.03
  #theta<-4
  # nu decreases, event rate increases, very sensitive, recommend change at 0.001 scale
  # theta increases, event rate increases, not too sensitive, can change at 1 scale
  # beta increases, event rate increases
  data.generate$cumrisk<-ave((1+beta*data.generate$x)*((data.generate$t+1)^theta-data.generate$t^theta), data.generate$id, FUN=cumsum)
  data.generate$ciF<-1-exp(-nu^theta*data.generate$cumrisk)
  1-aggregate(ciF~t,data.generate,mean)[t_star-1,]
  #event rate = 1-aggregate(ciF~t,data.generate,mean)[t_star-1,]
  #summary(data.generate$ciF)
  # generate outcome time T0, follow-up time T=min(T0,V,t_star) and event indicator D=I(T0<V)
  U<-runif(n)
  T0<-rep(NA,n)
  T<-rep(NA,n)
  for(i in 1:n){
    data.temp<-data.generate[which(data.generate$id==i),]
    data.temp$ciF.align<-c(0,data.temp$ciF[1:(t_star-1)])
    data.temp$cumrisk.align<-c(0,data.temp$cumrisk[1:(t_star-1)])
    if(U[i]>max(data.temp$ciF.align)){
      T0[i]<-2*t_star
    }
    else{
      index<-max(which(data.temp$ciF.align<=U[i]))
      T0[i]<-round((data.temp$t[index]^theta-1/(1+beta*data.generate$x[index])*(nu^(-theta)*log(1-U[i])+data.temp$cumrisk.align[index]))^(1/theta),3)
    }
    #T[i]<-min(T0[i],V[i],t_star)
    T[i]<-min(T0[i],t_star)
    data.temp$fT<-T[i]
    #data.temp$D<-c(rep(0,ceiling(T[i])-1),ifelse(T0[i]<=min(V[i],t_star),1,0),rep(NA,t_star-ceiling(T[i])))
    data.temp$D<-c(rep(0,ceiling(T[i])-1),ifelse(T0[i]<=t_star,1,0),rep(NA,t_star-ceiling(T[i])))
    if(i==1){
      data.main<-data.temp[1:ceiling(T[i]),]
    }
    else{
      data.main<-rbind(data.main,data.temp[1:ceiling(T[i]),])
    }
  }
  # reformat data into counting process structure
  data.main$t_in<-data.main$t-1
  data.main$t_out<-data.main$t_in+ifelse(data.main$D==1,data.main$fT-data.main$t+1,1)
  data.main$t_dur<-ifelse(data.main$D==1,data.main$fT-data.main$t+1,1)
  data.main$age_out<-data.main$age+ifelse(data.main$D==1,data.main$fT-data.main$t+1,1)
  #sum(data.main$D)/n
  # data.main consists of all observations before event, censoring or end of study
  
  # random sample m subjects into a validation data
  # would like to control no more than 3 measurements per subject
  # for now, take all measurements into validation data
  id.main<-unique(data.main$id)
  id.val<-sample(id.main,size=m,replace=FALSE)
  data.val<-data.main[data.main$id %in% id.val,]
  
  list(data.main=data.main,data.val=data.val,alpha=c(alpha1,alpha2),beta=beta)
}

#### Function for one-run simulation with ERR ####

simulation_one_run<-function(n,m,rho,beta,t_star,nu,theta){
  data.full<-gen_data(n,m,rho,beta,t_star,nu,theta)
  data.main<-data.full$data.main
  data.val<-data.full$data.val
  D.rate<-sum(data.main$D)/n
  
  model.err.true<-ERR_breslow(data=data.main,outcome=c("D"),covariates=c(),
                              exposure=c("x"),times=c("t_in","t_out"),strata=c(),
                              start=c(1))
  
  model.err.true.strata<-ERR_breslow(data=data.main,outcome=c("D"),covariates=c(),
                              exposure=c("x"),times=c("t_in","t_out"),strata=c("t_in"),
                              start=c(1))
  
  model.coxph.true<-coxph(Surv(t_in,t_out,D)~(x),data=data.main,ties="breslow")
  # Surv(t_in,t_out,D)~x, Surv(t_dur,D)~x+strata(t_in) are equivalent and correct way of fitting
  
  model.err.ep<-ERR_breslow(data=data.main,outcome=c("D"),covariates=c(),
                              exposure=c("X"),times=c("t_in","t_out"),strata=c(),
                              start=c(1))
  model.err.ep.strata<-ERR_breslow(data=data.main,outcome=c("D"),covariates=c(),
                            exposure=c("X"),times=c("t_in","t_out"),strata=c("t_in"),
                            start=c(1))
  # RRC using simple linear regression
  model.me<-lm(c~C,data=data.val)
  data.main$c_rrc<-predict(model.me,data.frame(C=data.main$C))
  data.main$x_rrc<-ave(data.main$c_rrc, data.main$id, FUN=cumsum)
  
  model.err.rrc<-ERR_breslow(data=data.main,outcome=c("D"),covariates=c(),
                             exposure=c("x_rrc"),times=c("t_in","t_out"),strata=c(),
                             start=c(1))
  model.err.rrc.strata<-ERR_breslow(data=data.main,outcome=c("D"),covariates=c(),
                             exposure=c("x_rrc"),times=c("t_in","t_out"),strata=c("t_in"),
                             start=c(1))
  
  list(data.full=data.full,
       betas.err=c(model.err.true$estimate,model.err.ep$estimate,model.err.rrc$estimate,model.coxph.true$coefficients),
       betas.err.strata=c(model.err.true.strata$estimate,model.err.ep.strata$estimate,model.err.rrc.strata$estimate))
}


#### Trial: ERR simulation with varying parameters ####

results1<-simulation_one_run(n=2000,m=500,
                             rho=0.6,beta=2,t_star=20,
                             nu=0.03,theta=4)
results1$betas.err
results1$betas.err.strata
sum(results1$data.full$data.main$D)/2000
# if ERR and CoxPH fit with strata t_in
#
#
# if ERR and CoxPH fit without strata
#-0.2490080 -0.2258608 -0.2497863  0.2142383 
#[1] 0.4245

#-1.8552508 -0.2183034 -2.1404469  1.4332035
#2.0611213 0.2451506 1.5706521
#0.176
results2<-simulation_one_run(n=2000,m=500,
                             rho=0.6,beta=2,t_star=20,
                             nu=0.025,theta=5)
results2$betas.err
results2$betas.err.strata
sum(results2$data.full$data.main$D)/2000
# if ERR and CoxPH fit with strata t_in
#
#
# if ERR fit without strata, this is the correct way
# since data generation is done without strata
# changes in nu and theta won't affect ERR beta estimates, will affect CoxPH beta estimates
#-0.2417310 -0.2203973 -0.2436167  0.5068278 
#[1] 0.2965

#-1.7714548 -0.2040467 -1.9660203 -1.6508859
#-0.9101391 -0.1535187 -1.4137116
#0.046

results3<-simulation_one_run(n=2000,m=500,
                             rho=0.6,beta=10,t_star=20,
                             nu=0.03,theta=4)
results3$betas.err
results3$betas.err.strata
sum(results3$data.full$data.main$D)/2000
#-1.8372257 -0.2134408 -2.1002984  1.7544715
#4.529310 3.042741 7.575624
#0.3385

results5<-simulation_one_run(n=2000,m=500,
                             rho=0.3,beta=2,t_star=20,
                             nu=0.03,theta=4)
results5$betas.err
results5$betas.err.strata
sum(results5$data.full$data.main$D)/2000
#-1.8540573 -0.2115331 -2.3076302 -0.6359180
#-0.53940186 -0.04419746 -0.57075357
#0.1765

set.seed(211)
results.sim1<-matrix(NA,nrow=12,ncol=8)
parameter.sim1<-matrix(NA,nrow=12,ncol=3)
counter<-1
for(i in 1:3){
  rho=0.3*(4-i)
  for(j in 1:2){
    beta=2*5^(j-1)
    for(k in 1:2){
      theta=5-k
      results.temp<-simulation_one_run(n=2000,m=500,
                                       rho=rho,beta=beta,t_star=20,
                                       nu=0.03,theta=theta)
      results.sim1[counter,]<-c(as.numeric(results.temp$betas.err),as.numeric(results.temp$betas.err.strata),sum(results.temp$data.full$data.main$D)/2000)
      parameter.sim1[counter,]<-c(rho,beta,theta)
      print(counter)
      counter<-counter+1
    }
  }
}
results.sim1<-as.data.frame(results.sim1)
colnames(results.sim1)<-c("err_true","err_ep","err_rrc","cox_true","err_true_strata","err_ep_strata","err_rrc_strata","prevalence")
View(results.sim1)
parameter.sim1<-as.data.frame(parameter.sim1)
colnames(parameter.sim1)<-c("rho","beta","theta")



results4<-matrix(NA,nrow=15,ncol=3)
for(i in 1:15){
  data.full<-gen_data(n=2000,m=500,
                      rho=0.6,beta=2,t_star=20,
                      nu=0.03,theta=4)
  data.main<-data.full$data.main
  data.val<-data.full$data.val
  
  # RRC using simple linear regression
  model.me<-lm(c~C,data=data.val)
  data.main$c_rrc<-predict(model.me,data.frame(C=data.main$C))
  data.main$x_rrc<-ave(data.main$c_rrc, data.main$id, FUN=cumsum)
  
  model.err.true<-ERR_breslow(data=data.main,outcome=c("D"),covariates=c(),
                              exposure=c("x"),times=c("t_in","t_out"),strata=c(),
                              start=c(1))
  model.err.ep<-ERR_breslow(data=data.main,outcome=c("D"),covariates=c(),
                              exposure=c("X"),times=c("t_in","t_out"),strata=c(),
                              start=c(1))
  model.err.rrc<-ERR_breslow(data=data.main,outcome=c("D"),covariates=c(),
                              exposure=c("x_rrc"),times=c("t_in","t_out"),strata=c(),
                              start=c(1))
  
  results4[i,]<-c(model.err.true$estimate,model.err.ep$estimate,model.err.rrc$estimate)
}
results4<-as.data.frame(results4)
colnames(results4)<-c("err_true","err_ep","err_rrc")
colMeans(results4)
apply(results4,2,var)





err.swcov<-ERR_sandwichcov(data=data.main,outcome=c("D"),
                           covariates=c(),
                           exposure=c("x"),times=c("t_in","t_out"),strata=c("t_in"),
                           beta=model.err.true$estimate)

results<-c()
for(j in 1:100){
  data.full<-gen_data(n=2000,
                      m=500,
                      alpha=c(0.001,0.9),
                      beta=0.5,
                      t_star=20)
  data.main<-data.full$data.main
  data.main$time_out<-data.main$time+ifelse(data.main$D==1,data.main$fT-data.main$t+1,1)
  data.main$time_dur<-ifelse(data.main$D==1,data.main$fT-data.main$t+1,1)
  model.err.true<-ERR(data=data.main,outcome=c("D"),
                      covariates=c(),
                      exposure=c("x"),times=c("time","time_out"),strata=c("time"),
                      start=c(0.5))
  results<-c(results,model.err.true$estimate)
}

mean(results)
[1] -0.1761454
var(results)
[1] 4.162685e-05
sd(results)
[1] 0.006451887

COV_sandwich1
0.002650386
COV_sandwich2
0.004157547

results_2<-c()
for(j in 1:20){
  data.full<-gen_data(n=2000,
                      m=500,
                      alpha=c(0.001,0.9),
                      beta=0.5,
                      t_star=20)
  data.main<-data.full$data.main
  data.main$time_out<-data.main$time+ifelse(data.main$D==1,data.main$fT-data.main$t+1,1)
  data.main$time_dur<-ifelse(data.main$D==1,data.main$fT-data.main$t+1,1)
  model.err.ep<-ERR(data=data.main,outcome=c("D"),
                      covariates=c(),
                      exposure=c("X"),times=c("time","time_out"),strata=c("time"),
                      start=c(0.5))
  results_2<-c(results_2,model.err.ep$estimate)
}

0.5
true x--ERR
mean(results)
-0.1761454
ep X--ERR
mean(results_2)
-0.1598
mean(results_3)
-0.17523




