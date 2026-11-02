#### Simulate/Generate Survival Data According to Cox ####

library(MASS)
library(dplyr)


gen_data_cox<-function(n,t_star,mu_c,var_c,mu_C,var_C,rho,lambda,nu,theta,beta){
  
  # n: number of subjects/unique id
  # t_star: maximum length of follow-up, end of study years
  
  # mu_c, var_c: mean and variance of true exposure variable
  # mu_C, var_C: mean and variance of error-prone exposure variable
  # rho: desired correlation between true and error-prone exposure variables
  
  # lambda: censoring rate ~ exp(lambda)
  # nu, theta: parameters of baseline hazard Weibull distribution
  
  # beta: desired true exposure effect
  
  # test value reference: xiaomei paper appendix 
  #n<-1000
  #t_star<-50
  #mu_c<-0.45
  #var_c<-0.0225
  #mu_C<-0.5
  #var_C<-0.04
  #rho<-0.6
  #lambda<-0.01
  #nu<-0.02
  #theta<-6
  #beta<-0.5
  
  # generate true yearly exposure c_{ij} ~ MVN(mu_c, Sigma_c)
  # i indicates subject ID 1-n
  # j indicates year order 1-t_star
  # start with independent correlation matrix between true yearly exposure
  mu_c_vec<-rep(mu_c,t_star)
  Sigma_c<-var_c*diag(1,t_star)
  c<-mvrnorm(n,mu_c_vec,Sigma_c)
  
  # generate error-prone yearly exposure C_{ij}
  # i and j as above
  # assuming the error model C_{ij}|c_{ij} to be additive
  alpha_1<-mu_C-rho*sqrt(var_C/var_c)*mu_c
  alpha_2<-rho*sqrt(var_C/var_c)
  e<-mvrnorm(n,rep(0,t_star),var_C/var_c*(1-rho^2)*Sigma_c)
  C<-alpha_1+alpha_2*c+e
  
  # create data frame with the following columns
  # ID: 1-n, t: 1-t_star, c, C
  data.generate<-data.frame(id=rep(c(1:n),t_star),t=rep(c(1:t_star),each=n),c=as.vector(c),C=as.vector(C))
  
  # calculate and create cumulative total exposure by year x_{ij}, X_{ij}
  data.generate$x<-ave(data.generate$c, data.generate$id, FUN=cumsum)
  data.generate$X<-ave(data.generate$C, data.generate$id, FUN=cumsum)
  
  # generate survival time T0
  # create new columns t.align and x.align, which corresponds to 0,t1,...,t_{t_star-1} and 0,x1,...,x_{t_star-1} 
  # to make calculation of cumulative risk and cumulative incidence function easier
  data.generate<-data.generate %>%
    group_by(id) %>%
    mutate(t.align=lag(t,default=0),x.align=lag(x,default=0)) %>%
    arrange(id)
  data.generate<-as.data.frame(data.generate)
  
  # calculate cumulative incidence function F_{ij} at each observation time
  data.generate$cumrisk<-ave(exp(beta*data.generate$x.align)*(data.generate$t^theta-data.generate$t.align^theta), data.generate$id, FUN=cumsum)
  data.generate$ciF<-1-exp(-nu^theta*data.generate$cumrisk)
  
  # generate random U(0,1) comparison variable for each id U_{i}
  U<-runif(n)
  data.generate$U_vec<-rep(U,each=t_star)
  data.generate$U_comp<-(data.generate$ciF<=data.generate$U_vec)
  
  # select the row for each id which satisfies F_{ij}<=U_{i}<=F_{i,j+1}
  data.incidence<-data.generate %>%
    group_by(id) %>%
    filter(U_comp==TRUE) %>%
    slice_max(order_by=t)
  
  # Note: the number of subjects in data.incidence could be fewer than n
  # F_{i1}=1-exp(-nu^theta)
  # there is a F_{i1}=1-exp(-nu^theta) chance that U_{i}<F_{i1}
  # indicating that subject i should not have survived past entering time 1
  # therefore these subjects should be removed from the simulated data
  # try to avoid combinations of nu and theta when F_{i1} is far away from 0 and too many subjects will be removed
  
  # calculate survival time T0 for each id
  # for those id that F_{i,t_star}<U_{i}, set T0 to be a big constant 2*t_star>t_star
  data.incidence$T0<-(data.incidence$t^theta-exp(-beta*data.incidence$x)*(nu^(-theta)*log(1-data.incidence$U_vec)+data.incidence$cumrisk))^(1/theta)
  data.incidence$T0[which(data.incidence$t==t_star)]<-2*t_star
  
  n_new<-length(unique(data.incidence$id))
  # generate censoring time V~exp(lambda)
  # if lambda=0, i.e. no censoring, set V=2*t_star
  if(lambda==0){
    V<-rep(2*t_star,n_new)
  } else {
    V<-rexp(n_new,rate=lambda)
  }
  data.incidence$V<-V
  
  # follow-up time fT=min(T0,V,t_star)
  data.incidence$fT<-pmin(data.incidence$T0,data.incidence$V,t_star)
  # disease status by id Di=I(T0<=min(V,t_star))
  data.incidence$Di<-as.numeric(data.incidence$T0<=pmin(data.incidence$V,t_star))
  
  # add columns T0,V,fT,Di to data.generate
  # subset data.generate to only those rows where t<=T (survive,not censored,before end of study)
  data.incidence<-as.data.frame(data.incidence)
  data.incidence.subset<-data.incidence[,c("id","T0","V","fT","Di")]
  data.generate<-data.generate %>% 
    filter(id %in% unique(data.incidence$id)) %>%
    full_join(data.incidence.subset,by="id") %>%
    filter(t<=fT)
  # create survival time by id and t: time_{ij}
  # =1, if survived or at end of study
  # =fT-t, if censored or became a case
  data.generate$time<-pmin(data.generate$fT-data.generate$t,1)
  data.generate$time[which(data.generate$fT==t_star)]<-1
  # create disease status by id and t: D_{ij}
  # =1, if became a case during [t,t+1]
  # =0, if survived, censored or at end of study
  data.generate$Dij<-data.generate$Di
  data.generate$Dij[which(data.generate$time==1)]<-0

  # create time_in, time_out columns 
  data.generate$t_in<-data.generate$t
  data.generate$t_out<-data.generate$t+data.generate$time
  
  return(data.generate)
  # which contains the following columns
  # id, t: id and measurement times
  # c: true yearly exposure, C: error-prone yearly exposure
  # x: true cumulative exposure, X: error-prone cumulative exposure
  # t.align, x.align: t,x with 1 year lag to calculate cumrisk, ciF more easily
  # cumrisk, ciF: cumulative incidence function
  # U_vec, U_comp: generate U(0,1) variable and make comparison with ciF
  # T0, V, fT, Di: simualted survival time, simulated censoring time, observed failure time, survival status of subject i
  # time, Dij, t_in, t_out: counting process data structure feature, time between 0-1, survival status from t-t+1
}


#### Test: Varying nu ####

# test out the generating function with varying nu
# holding other parameters constant
# nu decreases, cumulative incidence rate decreases;
# nu increases, 1-exp(-nu^theta) increases, sample size decreases

library(survival)

nu_list<-c(0.005,0.01,0.1,0.8)

for(j in 1:4){
  #nu<-nu_list[j]
  B<-50
  results<-matrix(rep(NA,B*4),nrow=B,ncol=4)
  for(i in 1:B){
  data.generate<-gen_data_cox(n<-1000,t_star<-50,
                              mu_c<-0.45,var_c<-0.0225,
                              mu_C<-0.5,var_C<-0.04,
                              rho<-0.6,lambda<-0,
                              nu<-nu_list[j],theta<-6,
                              beta<-0.5)
  results[i,1]<-sum(data.generate$Dij)/length(unique(data.generate$id))
  model.coxph.true<-coxph(Surv(t_in,t_out,Dij)~x,data=data.generate,ties="breslow")
  #model.coxph<-coxph(Surv(time,Dij)~x+strata(t),data=data.generate,ties="breslow")
  results[i,2]<-as.numeric(model.coxph.true$coefficients)
  model.coxph.ep<-coxph(Surv(t_in,t_out,Dij)~X,data=data.generate,ties="breslow")
  results[i,3]<-as.numeric(model.coxph.ep$coefficients)
  results[i,4]<-length(unique(data.generate$id))
  }
  colnames(results)<-c("cir","beta.true","beta.ep","sample size")
  results<-as.data.frame(results)
  
  print(colMeans(results))
  #apply(results,2,var)

}

#cir    beta.true      beta.ep  sample size 
#0.9856800    0.5020860    0.1975874 1000.0000000 
#cir    beta.true      beta.ep  sample size 
#1.0000000    0.4970803    0.2066906 1000.0000000 
#cir    beta.true      beta.ep  sample size 
#1.0000000    0.5066423    0.2250255 1000.0000000 
#cir   beta.true     beta.ep sample size 
#1.0000000   0.5155710   0.2637954 766.1000000

#set CIR: 20%-50%


