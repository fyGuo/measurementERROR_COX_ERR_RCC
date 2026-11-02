#### Load functions gen_data_cox() ####

library(maxLik)
library(MASS)
library(sampling)
library(survival)

gen_data_cox=function(n,m,rho,beta,t_star,nu,theta){
  # n: main study sample size
  # m: validation study sample size
  # alpha: measurement error model coef
  # beta: cox effect model coef
  # t_star: administrative study period
  # nu,theta: parameters of baseline hazard weibull distribution
  
  n<-2000
  m<-500
  alpha<-c(0.001,0.9)
  beta<-2
  t_star<-20
  nu<-0.03
  theta<-4
  
  # generate error-prone yearly exposure C_{ij}, i indicates subject 1-n, j indicates time 1-t_star
  # AR(1) correlation matrix
  #Sigma_C<-0.002*ar1_cor(n=t_star,rho=0.996)
  # independent correlation matrix
  Sigma_C<-0.002*diag(1,t_star)
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
  # beta increases, event rate increas es
  data.generate$cumrisk<-ave(exp(beta*data.generate$x)*((data.generate$t+1)^theta-data.generate$t^theta), data.generate$id, FUN=cumsum)
  data.generate$ciF<-1-exp(-nu^theta*data.generate$cumrisk)
  #1-aggregate(ciF~t,data.generate,mean)[t_star-1,]
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
      T0[i]<-round((data.temp$t[index]^theta-exp(-beta*data.generate$x[index])*(nu^(-theta)*log(1-U[i])+data.temp$cumrisk.align[index]))^(1/theta),3)
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


gen_data_cox=function(n,m,rho,beta,t_star,nu,theta){
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
  
  # generate true yearly exposure C_{ij}, i indicates subject 1-n, j indicates time 1-t_star
  # AR(1) correlation matrix
  #Sigma_C<-0.002*ar1_cor(n=t_star,rho=0.996)
  # independent correlation matrix
  #Sigma_c<-0.0005*diag(1,t_star)
  #c<-mvrnorm(n,mu=rep(0.2,t_star),Sigma_c)
  Sigma_c<-0.3*diag(1,t_star)
  c<-mvrnorm(n,mu=rep(5,t_star),Sigma_c)
  # generate error-prone yearly exposure c_{ij}, i, j as above
  #c<-alpha[1]+alpha[2]*C+mvrnorm(n,rep(0,t_star),diag(rep(1/17*(1-0.7^2)*var(as.numeric(C)),t_star)))
  #alpha1<-0.01840276-mean(c)*rho*sqrt(17)
  #alpha2<-rho*sqrt(17)
  alpha1<-0
  alpha2<-1
  C<-c+mvrnorm(n,rep(0,t_star),(1-rho^2)*2*Sigma_c)
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
  # beta increases, event rate increas es
  data.generate$cumrisk<-ave(exp(beta*data.generate$x)*((data.generate$t+1)^theta-data.generate$t^theta), data.generate$id, FUN=cumsum)
  data.generate$ciF<-1-exp(-nu^theta*data.generate$cumrisk)
  #1-aggregate(ciF~t,data.generate,mean)[t_star-1,]
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
      T0[i]<-round((data.temp$t[index]^theta-exp(-beta*data.generate$x[index])*(nu^(-theta)*log(1-U[i])+data.temp$cumrisk.align[index]))^(1/theta),3)
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

#### Function for one-run simulation with Cox ####

simulation_one_run_cox<-function(n,m,rho,beta,t_star,nu,theta){
  data.full<-gen_data_cox(n,m,rho,beta,t_star,nu,theta)
  data.main<-data.full$data.main
  data.val<-data.full$data.val
  D.rate<-sum(data.main$D)/n
  
  model.err.true<-ERR_breslow(data=data.main,outcome=c("D"),covariates=c(),
                              exposure=c("x"),times=c("t_in","t_out"),strata=c(),
                              start=c(1))
  
  model.coxph.true<-coxph(Surv(t_in,t_out,D)~(x),data=data.main,ties="breslow")
  # Surv(t_in,t_out,D)~x, Surv(t_dur,D)~x+strata(t_in) are equivalent and correct way of fitting
  
  model.err.ep<-ERR_breslow(data=data.main,outcome=c("D"),covariates=c(),
                            exposure=c("X"),times=c("t_in","t_out"),strata=c(),
                            start=c(1))
  model.coxph.ep<-coxph(Surv(t_in,t_out,D)~(X),data=data.main,ties="breslow")
  
  # RRC using simple linear regression
  model.me<-lm(c~C,data=data.val)
  data.main$c_rrc<-predict(model.me,data.frame(C=data.main$C))
  data.main$x_rrc<-ave(data.main$c_rrc, data.main$id, FUN=cumsum)
  
  model.err.rrc<-ERR_breslow(data=data.main,outcome=c("D"),covariates=c(),
                             exposure=c("x_rrc"),times=c("t_in","t_out"),strata=c(),
                             start=c(1))
  model.coxph.rrc<-coxph(Surv(t_in,t_out,D)~(x_rrc),data=data.main,ties="breslow")
  
  
  list(data.full=data.full,
       betas.err=c(model.err.true$estimate,model.err.ep$estimate,model.err.rrc$estimate),
       betas.cox=c(model.coxph.true$coefficients,model.coxph.ep$coefficients,model.coxph.rrc$coefficients))
}

#### Trial: Cox simulation with varying parameters ####

results1<-simulation_one_run_cox(n=2000,m=500,
                             rho=0.6,beta=2,t_star=20,
                             nu=0.03,theta=4)
results1$betas.err
results1$betas.cox
sum(results1$data.full$data.main$D)/2000

results2.sim2<-results2
results2<-matrix(NA,nrow=10,ncol=3)
for(i in 1:10){
  data.full<-gen_data_cox(n=2000,m=500,
                          rho=0.6,beta=2,t_star=20,
                          nu=0.03,theta=4)
  data.main<-data.full$data.main
  data.val<-data.full$data.val
  
  # RRC using simple linear regression
  model.me<-lm(c~C,data=data.val)
  data.main$c_rrc<-predict(model.me,data.frame(C=data.main$C))
  data.main$x_rrc<-ave(data.main$c_rrc, data.main$id, FUN=cumsum)
  
  model.coxph.true<-coxph(Surv(t_in,t_out,D)~(x),data=data.main,ties="breslow")
  model.coxph.ep<-coxph(Surv(t_in,t_out,D)~(X),data=data.main,ties="breslow")
  model.coxph.rrc<-coxph(Surv(t_in,t_out,D)~(x_rrc),data=data.main,ties="breslow")
  
  results2[i,]<-c(model.coxph.true$coefficients,model.coxph.ep$coefficients,model.coxph.rrc$coefficients)
}

results2<-as.data.frame(results2[1:10,])
colnames(results2)<-c("cox_true","cox_ep","cox_rrc")
colMeans(results2)
#cox_true    cox_ep   cox_rrc 
#1.3843195 0.1710627 1.1775734 
apply(results2,MARGIN=2,FUN=var)
#1.21312 0.06579977 3.107211


