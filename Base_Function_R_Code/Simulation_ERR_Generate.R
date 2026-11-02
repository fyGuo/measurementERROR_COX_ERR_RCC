#### Simulate/Generate Survival Data According to ERR ####

library(MASS)
library(survival)
library(dplyr)

gen_data_err=function(n,m,rho,beta,t_star,nu,theta){
  # n: main study sample size
  # m: validation study sample size
  # alpha: measurement error model coef
  # beta: ERR effect model coef
  # t_star: administrative study period
  # nu,theta: parameters of baseline hazard function weibull distribution
  
  #n<-2000
  #m<-500
  #rho<-0.6
  #beta<-2
  #t_star<-20
  #nu<-1
  #theta<-1
  
  # generate true yearly exposure C_{ij}, i indicates subject 1-n, j indicates time 1-t_star
  # independent correlation matrix
  Sigma_c<-2*diag(1,t_star)
  c<-mvrnorm(n,mu=rep(5,t_star),Sigma_c)
  
  # generate error-prone yearly exposure c_{ij}, i, j as above
  alpha1<-7.5-5*rho*sqrt(3/2)
  alpha2<-rho*sqrt(3/2)
  C<-alpha1+alpha2*c+mvrnorm(n,rep(0,t_star),(1-rho^2)*(3/2)*Sigma_c)
  
  # create data frame with columns: id 1-2000, t 1-20, C, c
  data.generate<-data.frame(id=rep(c(1:n),t_star),t=rep(c(1:t_star),each=n), C=as.vector(C),c=as.vector(c))
  
  # calculate cumulative total true exposure x
  data.generate$x<-ave(data.generate$c, data.generate$id, FUN=cumsum)
  # calculate cumulative total error-prone exposure X
  data.generate$X<-ave(data.generate$C, data.generate$id, FUN=cumsum)
  
  # generate censoring time V~exp(0.01)
  # ignore censorying for now
  #V<-rexp(n,rate=0.01)
  #V[V>t_star]=t_star
  
  # calculate cumulative incidence function F_{ij}, i, j as above
  # cumrisk sum from 0-1,0-2,...,0-20, need to re-align later
  #data.generate$cumrisk<-ave(exp(beta*data.generate$x)*((data.generate$t+1)^theta-data.generate$t^theta), data.generate$id, FUN=cumsum)
  #data.generate$ciF<-1-exp(-nu^theta*data.generate$cumrisk)
  #event rate = 1-aggregate(ciF~t,data.generate,mean)[t_star-1,]
  
  # generate outcome time T0, follow-up time T=min(T0,t_star) and event indicator D=I(T0<t_star)
  U<-runif(n)
  T0<-rep(NA,n)
  T<-rep(NA,n)
  
  #nu<-0.003
  #theta<-6
  
  data.generate<-data.generate %>%
    group_by(id) %>%
    mutate(t.align=lag(t,default=0),x.align=lag(x,default=0))
  data.generate<-as.data.frame(data.generate)
  data.generate$cumrisk.align<-ave((1+beta*data.generate$x.align)*(data.generate$t^theta-data.generate$t.align^theta), data.generate$id, FUN=cumsum)
  data.generate$ciF.align<-1-exp(-nu^theta*data.generate$cumrisk.align)
  
  for(i in 1:n){
    # select 20 rows of data for subject i
    data.temp<-data.generate[which(data.generate$id==i),]
    # calculate cumulative incidence function F_{ij}, i,j as above
    # Note that cumrisk sum from 0 to t_star-1, need to re-align t and x variable
    # Set t0=0, x(0)=0
    # Variables with .align goes from 0 to t_star-1, otherwise from 0 to t_star
    #data.temp$t.align<-c(0, data.temp$t[1:(t_star-1)])
    #data.temp$x.align<-c(0, data.temp$x[1:(t_star-1)])
    #data.temp$cumrisk.align<-cumsum((1+beta*data.temp$x.align)*((data.temp$t)^theta-data.temp$t.align^theta))
    #data.temp$ciF.align<-1-exp(-nu^theta*data.temp$cumrisk.align)
    # re-align cumrisk, ciF, now sum from 0-0,0-1,...,0-19
    #data.temp$ciF.align<-c(1-exp(-nu^theta),data.temp$ciF[1:(t_star-1)])
    #data.temp$cumrisk.align<-c(0,data.temp$cumrisk[1:(t_star-1)])
    
    # calculate event time
    if(U[i]>max(data.temp$ciF.align)){
      T0[i]<-2*t_star
    }
    else{
      # find j such that F_{ij}<=U_{i}<F_{i,j+1}
      index<-max(which(data.temp$ciF.align<=U[i]))
      # calculate T0_{i}, euqation (A.4), round to 3 digits
      T0[i]<-round((data.temp$t[index]^theta-1/(1+beta*data.temp$x[index])*(nu^(-theta)*log(1-U[i])+data.temp$cumrisk.align[index]))^(1/theta),3)
    }
    
    # define fT failure time variable
    T[i]<-min(T0[i],t_star)
    data.temp$fT<-T[i]
    
    # define event indicator: D = I(T0_{i}<=t_star)
    data.temp$D<-c(rep(0,ceiling(T[i])-1),ifelse(T0[i]<=t_star,1,0),rep(NA,t_star-ceiling(T[i])))
    
    # save newly added variables for subject i in a new dataset: data.main
    # for subject i with event at period j, only include j entries, discard all entries after event
    if(i==1){
      data.main<-data.temp[1:ceiling(T[i]),]
    }
    else{
      data.main<-rbind(data.main,data.temp[1:ceiling(T[i]),])
    }
  }
  
  #sum(data.main$D)/2000
  
  # create variables to fit
  data.main$t_in<-data.main$t-1
  data.main$t_out<-data.main$t_in+ifelse(data.main$D==1,data.main$fT-data.main$t+1,1)
  data.main$t_dur<-ifelse(data.main$D==1,data.main$fT-data.main$t+1,1)
  # data.main consists of all observations before event, censoring or end of study
  
  # random sample m subjects into a validation data
  # would like to control no more than 3 measurements per subject
  # for now, take all measurements into validation data
  id.main<-unique(data.main$id)
  id.val<-sample(id.main,size=m,replace=FALSE)
  data.val<-data.main[data.main$id %in% id.val,]
  
  list(data.main=data.main,alpha=c(alpha1,alpha2),beta=beta)
}

# run the data generation function
output.err<-gen_data_err(n=3000,m=500,
                     rho=0.6,beta=1,t_star=40,
                     nu=0.03,theta=6)

mainstudy.err<-output.err$data.main
#View(mainstudy)
# check prevalence
sum(mainstudy.err$D)/3000

model.coxph.err<-coxph(Surv(t_in,t_out,D)~x,data=mainstudy.err,ties="breslow")
model.coxph.err$coefficients

#### Fit ERR Model using Breslow for ties ####


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
      # create risk set: stratified by AGE_TV, entry before timein, survival past death timeout
      # in our simulated dataset, since year always increases by 1 year,
      # stratified by AGE_TV does not affect S_set creation
      if(is.null(strata)==TRUE){
        S_set<-which(as.matrix(data[,timein])<=D_age&as.matrix(data[,timeout])>=D_time)
      }
      else{
        S_set<-which(as.matrix(data[,strata])==D_age&as.matrix(data[,timein])<=D_age&as.matrix(data[,timeout])>=D_time)
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
  
  # definition of constrains for ERR parameters: 1+beta*doses>0, beta > -1/min(D)? 
  if(M==0){
    A<-as.matrix(1)
  }
  else{
    A <- matrix(c(rep(0, M), 1), ncol=M + 1)
  }
  B <- 1/min(doses)
  constraints <- list(ineqA = A, ineqB = B) # constrains, RR term must be positive
  
  # optimize log likelihood function
  library(maxLik)
  model.err<-maxLik(logLik=llk_err,start=start,constraints=constraints,control=list(iterlim=5000),method="BFGS")
  
  return(model.err)
}


#### Fit ERR Model naively not considering ties ####

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
      S_set<-which(as.matrix(data[,timein])<=D_age&as.matrix(data[,timeout])>=D_time)
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

model.err.true.1<-ERR_breslow(data=mainstudy.err,outcome=c("D"),covariates=c(),
                              exposure=c("x"),times=c("t_in","t_out"),strata=c(),
                              start=c(2))
model.err.true.1$estimate
# with n=5000,m=500,rho=0.6,beta=2,t_star=20,nu=0.03,theta=6
# takes about 5 minutes to run once
# prevalence=0.5184
# 1.456084



model.err.true.2<-ERR_breslow(data=mainstudy.err,outcome=c("D"),covariates=c(),
                              exposure=c("x"),times=c("t_in","t_out"),strata=c("t_in"),
                              start=c(1))
model.err.true.2$estimate

model.err.true.3<-ERR(data=mainstudy.err,outcome=c("D"),covariates=c(),
                              exposure=c("x"),times=c("t_in","t_out"),strata=c(),
                              start=c(1))
model.err.true.3$estimate
# with n=5000,m=500,rho=0.6,beta=2,t_star=20,nu=0.03,theta=6
# takes about 5 minutes to run once
# prevalence=0.5184
# 1.453942


# repetitively run above process, not efficient, run about 10min
results<-matrix(NA,nrow=20,ncol=4)
for(i in 1:20){
  print(i)
  output.err<-gen_data_err(n=2000,m=500,
                           rho=0.6,beta=2,t_star=20,
                           nu=0.03,theta=6)
  
  mainstudy.err<-output.err$data.main
  
  model.err.true<-ERR_breslow(data=mainstudy.err,outcome=c("D"),covariates=c(),
                              exposure=c("x"),times=c("t_in","t_out"),strata=c(),
                              start=c(2))
  model.err.ep<-ERR_breslow(data=mainstudy.err,outcome=c("D"),covariates=c(),
                            exposure=c("X"),times=c("t_in","t_out"),strata=c(),
                            start=c(2))
  model.coxph.err<-coxph(Surv(t_in,t_out,D)~x,data=mainstudy.err,ties="breslow")
  
  results[i,]<-c(model.err.true$estimate,model.err.ep$estimate,model.coxph.err$coefficients,sum(mainstudy.err$D)/2000)
}
#7:36-51
results<-as.data.frame(results)
colnames(results)<-c("err_true","err_ep","coxph_err_true","prevalence")
View(results)

colMeans(results)
#err_true         err_ep coxph_err_true     prevalence 
#2.8511735      1.4496689      0.1066048      0.5298750
#err_true         err_ep coxph_err_true     prevalence 
#2.6217303      0.3928985      0.1098473      0.5237750
apply(results,2,var)
#err_true         err_ep coxph_err_true     prevalence 
#1.494870e+01   8.894607e+00   2.891871e-03   4.839145e-05
#err_true         err_ep coxph_err_true     prevalence 
#7.6351082076   0.8058764953   0.0025756135   0.0001150914
# estimates range too large



output.err<-gen_data_err(n=2000,m=500,
                         rho=0.6,beta=2,t_star=20,
                         nu=0.03,theta=6)

mainstudy.err<-output.err$data.main
#View(mainstudy)
# check prevalence
sum(mainstudy.err$D)/2000

output.err2<-gen_data_err(n=2000,m=500,
                         rho=0.6,beta=2,t_star=20,
                         nu=0.035,theta=6)

mainstudy.err2<-output.err2$data.main
#View(mainstudy)
# check prevalence
sum(mainstudy.err2$D)/2000

model.err.true.1<-ERR_breslow(data=mainstudy.err,outcome=c("D"),covariates=c(),
                              exposure=c("x"),times=c("t_in","t_out"),strata=c(),
                              start=c(2))
model.err.true.1$estimate

model.err.true.2<-ERR_breslow(data=mainstudy.err2,outcome=c("D"),covariates=c(),
                              exposure=c("x"),times=c("t_in","t_out"),strata=c(),
                              start=c(2))
model.err.true.2$estimate



library(coxed)

hazard.err<-function(t){
  
}

simdata<-sim.survdata(N=1000,T=20,xvars=1,cernsor=0,beta=c(1),num.data.frames=1,type="tvc")



