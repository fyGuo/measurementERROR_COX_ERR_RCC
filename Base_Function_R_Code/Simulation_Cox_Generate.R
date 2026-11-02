#### Simulate/Generate Survival Data According to Cox ####

library(MASS)
library(survival)
library(dplyr)

gen_data_cox=function(n,m,rho,beta,t_star,nu,theta){
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
  #nu<-0.05
  #theta<-6
  
  # generate true yearly exposure C_{ij}, i indicates subject 1-n, j indicates time 1-t_star
  # independent correlation matrix
  Sigma_c<-0.0225*diag(1,t_star)
  c<-mvrnorm(n,mu=rep(0.45,t_star),Sigma_c)
  
  # generate error-prone yearly exposure c_{ij}, i, j as above
  alpha1<-0.5-0.45*rho*sqrt(0.04/0.0225)
  alpha2<-rho*sqrt(0.04/0.0225)
  C<-alpha1+alpha2*c+mvrnorm(n,rep(0,t_star),(1-rho^2)*(0.04/0.0225)*Sigma_c)
  
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
  
  data.generate<-data.generate %>%
    group_by(id) %>%
    mutate(t.align=lag(t,default=0),x.align=lag(x,default=0))
  data.generate<-as.data.frame(data.generate)
  data.generate$cumrisk.align<-ave(exp(beta*data.generate$x.align)*(data.generate$t^theta-data.generate$t.align^theta), data.generate$id, FUN=cumsum)
  data.generate$ciF.align<-1-exp(-nu^theta*data.generate$cumrisk.align)
  
  #nu<-0.003
  #theta<-6
  
  for(i in 1:n){
    # select 20 rows of data for subject i
    data.temp<-data.generate[which(data.generate$id==i),]
    # calculate cumulative incidence function F_{ij}, i,j as above
    # Note that cumrisk sum from 0 to t_star-1, need to re-align t and x variable
    # Set t0=0, x(0)=0
    # Variables with .align goes from 0 to t_star-1, otherwise from 0 to t_star
    #data.temp$t.align<-c(0, data.temp$t[1:(t_star-1)])
    #data.temp$x.align<-c(0, data.temp$x[1:(t_star-1)])
    #data.temp$cumrisk.align<-cumsum(exp(beta*data.temp$x.align)*((data.temp$t)^theta-data.temp$t.align^theta))
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
      T0[i]<-round((data.temp$t[index]^theta-exp(-beta*data.temp$x[index])*(nu^(-theta)*log(1-U[i])+data.temp$cumrisk.align[index]))^(1/theta),3)
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
output<-gen_data_cox(n=5000,m=500,
             rho=0.6,beta=log(2),t_star=20,
             nu=0.015,theta=4)

mainstudy<-output$data.main
#View(mainstudy)
# check prevalence
sum(mainstudy$D)/5000

# cox proportional hazard model on true x, should return beta=2
model.coxph.true.1<-coxph(Surv(t_in,t_out,D)~x,data=mainstudy,ties="breslow")
model.coxph.true.1$coefficients

# equivalent ways of writing the model
model.coxph.true.2<-coxph(Surv(t_dur,D)~x+strata(t_in),data=mainstudy,ties="breslow")
model.coxph.true.2$coefficients

# cox proportional hazard model on error-prone x
model.coxph.ep<-coxph(Surv(t_in,t_out,D)~X,data=mainstudy,ties="breslow")
model.coxph.ep$coefficients


#### Fit CoxPH Model using Breslow for ties ####

Cox_breslow<-function(data,outcome,covariates,exposure,times,strata,start){
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
        L<-L+as.matrix(doses[D_set,])%*%b[(M + 1):(M + Nd)]-log(sum(exp(as.matrix(doses[S_set,])%*%b[(M + 1):(M + Nd)])))
      }
      else{
        L<-L+bcvr[D_set,]%*%b[1:M]+log(1+as.matrix(doses[D_set,])%*%b[(M + 1):(M + Nd)])-log(sum(exp(bcvr[S_set,]%*%b[1:M])*(1+as.matrix(doses[S_set,])%*%b[(M + 1):(M + Nd)])))
      }
    }
    return(L)
  }
  
  # definition of constrains for ERR parameters: 1+beta*doses>0, beta > -1/max(D)
  A <- matrix(c(rep(0, M), 1), ncol=M + 1)
  B <- 0
  constraints <- list(ineqA = A, ineqB = B) # constrains, RR term must be positive
  
  # optimize log likelihood function
  library(maxLik)
  model.cox<-maxLik(logLik=llk_err,start=start,constraints=constraints,control=list(iterlim=5000),method="BFGS")
  
  return(model.cox)
}

model.coxph.true.3<-Cox_breslow(data=mainstudy,outcome=c("D"),covariates=c(),
                                exposure=c("x"),times=c("t_in","t_out"),strata=c(),
                                start=c(2))
model.coxph.true.3$estimate


# repetitively run above process, not efficient, run about 10min
set.seed(211)
results1<-matrix(NA,nrow=20,ncol=5)
for(i in 1:20){
  output<-gen_data_cox(n=2000,m=500,
                       rho=0.6,beta=2,t_star=20,
                       nu=0.003,theta=6)
  
  mainstudy<-output$data.main
  
  model.coxph.true.1<-coxph(Surv(t_in,t_out,D)~x,data=mainstudy,ties="breslow")
  model.coxph.ep<-coxph(Surv(t_in,t_out,D)~X,data=mainstudy,ties="breslow")
  
  results1[i,]<-c(model.coxph.true.1$coefficients,as.numeric(model.coxph.true.1$var),model.coxph.ep$coefficients,as.numeric(model.coxph.ep$var),sum(mainstudy$D)/2000)
}
results1<-as.data.frame(results1)
colnames(results1)<-c("coxph_true","var_true","coxph_ep","var_ep","prevalence")
View(results1)

# repetitive mean of beta coef
colMeans(results1)
#coxph_true    var_true    coxph_ep      var_ep  prevalence 
#1.849083240 0.003993760 0.653743247 0.001602436 0.444725000
# repetitive var of beta coef
apply(results1,2,var)
#coxph_true     var_true     coxph_ep       var_ep   prevalence 
#4.020295e-03 5.384054e-08 8.865805e-04 8.099404e-09 1.556441e-04

results2<-matrix(NA,nrow=20,ncol=5)
for(i in 1:20){
  output<-gen_data_cox(n=2000,m=500,
                       rho=0.6,beta=1.5,t_star=20,
                       nu=0.003,theta=6)
  
  mainstudy<-output$data.main
  
  model.coxph.true.1<-coxph(Surv(t_in,t_out,D)~x,data=mainstudy,ties="breslow")
  model.coxph.ep<-coxph(Surv(t_in,t_out,D)~X,data=mainstudy,ties="breslow")
  
  results2[i,]<-c(model.coxph.true.1$coefficients,as.numeric(model.coxph.true.1$var),model.coxph.ep$coefficients,as.numeric(model.coxph.ep$var),sum(mainstudy$D)/2000)
}
results2<-as.data.frame(results2)
colnames(results2)<-c("coxph_true","var_true","coxph_ep","var_ep","prevalence")
View(results2)

# repetitive mean of beta coef
colMeans(results2)
#coxph_true   var_true   coxph_ep     var_ep prevalence 
#1.42142717 0.09487731 0.59773382 0.05335872 0.01275000 
# repetitive var of beta coef
apply(results2,2,var)
#coxph_true     var_true     coxph_ep       var_ep   prevalence 
#1.022740e-01 4.806771e-04 6.801638e-02 1.391016e-04 8.223684e-06


results3<-matrix(NA,nrow=20,ncol=5)
for(i in 1:20){
  output<-gen_data_cox(n=2000,m=500,
                       rho=0.6,beta=2,t_star=20,
                       nu=0.0035,theta=6)
  
  mainstudy<-output$data.main
  
  model.coxph.true.1<-coxph(Surv(t_in,t_out,D)~x,data=mainstudy,ties="breslow")
  model.coxph.ep<-coxph(Surv(t_in,t_out,D)~X,data=mainstudy,ties="breslow")
  
  results3[i,]<-c(model.coxph.true.1$coefficients,as.numeric(model.coxph.true.1$var),model.coxph.ep$coefficients,as.numeric(model.coxph.ep$var),sum(mainstudy$D)/2000)
}
results3<-as.data.frame(results3)
colnames(results3)<-c("coxph_true","var_true","coxph_ep","var_ep","prevalence")
View(results3)

# repetitive mean of beta coef
colMeans(results3)
#coxph_true    var_true    coxph_ep      var_ep  prevalence 
#1.820064844 0.002964649 0.608634788 0.001123870 0.653925000 
# repetitive var of beta coef
apply(results3,2,var)
#coxph_true     var_true     coxph_ep       var_ep   prevalence 
#2.906872e-03 1.356406e-08 9.395859e-04 2.890774e-09 4.316513e-05

results4<-matrix(NA,nrow=20,ncol=5)
for(i in 1:20){
  output<-gen_data_cox(n=2000,m=500,
                       rho=0.6,beta=2,t_star=20,
                       nu=0.003,theta=7)
  
  mainstudy<-output$data.main
  
  model.coxph.true.1<-coxph(Surv(t_in,t_out,D)~x,data=mainstudy,ties="breslow")
  model.coxph.ep<-coxph(Surv(t_in,t_out,D)~X,data=mainstudy,ties="breslow")
  
  results4[i,]<-c(model.coxph.true.1$coefficients,as.numeric(model.coxph.true.1$var),model.coxph.ep$coefficients,as.numeric(model.coxph.ep$var),sum(mainstudy$D)/2000)
}
results4<-as.data.frame(results4)
colnames(results4)<-c("coxph_true","var_true","coxph_ep","var_ep","prevalence")
View(results4)

# repetitive mean of beta coef
colMeans(results4)
#coxph_true   var_true   coxph_ep     var_ep prevalence 
#1.91940074 0.01998575 0.83298963 0.01045882 0.06507500
# repetitive var of beta coef
apply(results4,2,var)
#  coxph_true     var_true     coxph_ep       var_ep   prevalence 
#1.486549e-02 6.194093e-06 1.553858e-02 1.779359e-06 3.266513e-05


results5<-matrix(NA,nrow=100,ncol=5)
for(i in 1:100){
  output<-gen_data_cox(n=2000,m=500,
                       rho=0.6,beta=2,t_star=20,
                       nu=0.003,theta=6)
  
  mainstudy<-output$data.main
  
  model.coxph.true.1<-coxph(Surv(t_in,t_out,D)~x,data=mainstudy,ties="breslow")
  model.coxph.ep<-coxph(Surv(t_in,t_out,D)~X,data=mainstudy,ties="breslow")
  
  results5[i,]<-c(model.coxph.true.1$coefficients,as.numeric(model.coxph.true.1$var),model.coxph.ep$coefficients,as.numeric(model.coxph.ep$var),sum(mainstudy$D)/2000)
}
results5<-as.data.frame(results5)
colnames(results5)<-c("coxph_true","var_true","coxph_ep","var_ep","prevalence")
View(results5)

# repetitive mean of beta coef
colMeans(results5)
#coxph_true    var_true    coxph_ep      var_ep  prevalence 
#1.839026946 0.003974402 0.658451392 0.001604115 0.446420000
# repetitive var of beta coef
apply(results5,2,var)
#coxph_true     var_true     coxph_ep       var_ep   prevalence 
#4.437123e-03 4.045413e-08 1.699465e-03 7.659852e-09 1.277259e-04


results6<-matrix(NA,nrow=50,ncol=5)
for(i in 1:50){
  output<-gen_data_cox(n=5000,m=500,
                       rho=0.6,beta=1,t_star=20,
                       nu=0.007,theta=4)
  
  mainstudy<-output$data.main
  
  model.coxph.true.1<-coxph(Surv(t_in,t_out,D)~x,data=mainstudy,ties="breslow")
  model.coxph.ep<-coxph(Surv(t_in,t_out,D)~X,data=mainstudy,ties="breslow")
  
  results6[i,]<-c(model.coxph.true.1$coefficients,as.numeric(model.coxph.true.1$var),model.coxph.ep$coefficients,as.numeric(model.coxph.ep$var),sum(mainstudy$D)/5000)
}
results6<-as.data.frame(results6)
colnames(results6)<-c("coxph_true","var_true","coxph_ep","var_ep","prevalence")
View(results6)

# repetitive mean of beta coef
colMeans(results6)
#coxph_true     var_true     coxph_ep       var_ep   prevalence 
#0.9402509288 0.0009843811 0.3909790741 0.0005049063 0.5600400000
apply(results6,2,var)
#coxph_true     var_true     coxph_ep       var_ep   prevalence 
#9.409503e-04 1.195171e-09 4.612664e-04 1.720947e-10 5.936163e-05

results7<-matrix(NA,nrow=50,ncol=5)
for(i in 1:50){
  output<-gen_data_cox(n=5000,m=500,
                       rho=0.6,beta=1,t_star=20,
                       nu=0.007,theta=4)
  
  mainstudy<-output$data.main
  
  model.coxph.true.1<-coxph(Surv(t_in,t_out,D)~x,data=mainstudy,ties="breslow")
  model.coxph.ep<-coxph(Surv(t_in,t_out,D)~X,data=mainstudy,ties="breslow")
  
  results7[i,]<-c(model.coxph.true.1$coefficients,as.numeric(model.coxph.true.1$var),model.coxph.ep$coefficients,as.numeric(model.coxph.ep$var),sum(mainstudy$D)/5000)
}
results7<-as.data.frame(results7)
colnames(results7)<-c("coxph_true","var_true","coxph_ep","var_ep","prevalence")
View(results7)

# repetitive mean of beta coef
colMeans(results7)

apply(results7,2,var)


