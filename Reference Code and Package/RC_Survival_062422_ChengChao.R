rm(list=ls())
library("spatstat")
library(survival)
library(nleqslv)
library(sandwich)
library(magic)

gen_time=function(theta=8,v,LP) {
  U = runif(length(LP))
  1/v*(-log(1-U)*exp(-LP))^(1/theta)
}

####### generate data
gen_data=function(n=10000,n2=250,rho=0.25,corma=0.4,beta123=c(log(1.5),log(1.5),log(1.5))) {
  # n=5000
  # n2=200
  # rho=0.8
  # corma=0.6
  # beta123=c(log(1.5),log(1.5),log(1.5))
  n1=n-n2
  # generate exposure (A and As), mediator (M), and confounder (W)
  AMW.mat=MASS::mvrnorm(n=n,mu=c(0,0,0),Sigma=matrix(c(1,corma,0.2,corma,1,0.2,0.2,0.2,1),ncol=3,nrow=3,byrow=T))
  A=AMW.mat[,1]
  M=AMW.mat[,2]
  W=AMW.mat[,3]
  As=A+rnorm(n,mean=0,sd=sqrt(1/rho^2 - 1))
  # mediator model parameter
  alpha = c(0,c(solve(matrix(c(1,0.2,0.2,1),ncol=2,nrow=2,byrow=T)) %*% c(corma,0.2)))
  sigma2.alpha = 1 - c(t(c(corma,0.2)) %*% solve(matrix(c(1,0.2,0.2,1),byrow=T,ncol=2)) %*% c(corma,0.2))
  # generate outcome (Time and Event)
  U=runif(n)
  LP = beta123[1]*A + beta123[2]*M + beta123[3]*A*M + log(1.2)*W
  Time.tilde = gen_time(theta=6,v=0.013,LP=LP)
  C = rexp(n,rate=0.01)
  C[C>50] = 50
  Event = as.numeric(Time.tilde<C)
  Time = ifelse(Time.tilde<C,Time.tilde,C)
  Time = round(Time,1)
  # outcome model parameter
  beta=c(beta123,log(1.2)) # coefficient
  # main study / validation study
  I = rep(0,n); I[(n1+1):n]=1
  Event[which(I==1)]=NA; #Time[which(I==1)]=NA;M[which(I==1)]=NA 
  list(data=data.frame(A=A,M=M,Time=Time,Event=Event,As=As,W=W,I=I),
       alpha=alpha,
       sigma2.alpha=sigma2.alpha,
       beta=beta)
}

get_true_effect.f=function(data, a0=0, a1=1, cw=0,true=TRUE) {
  #t=0.5;a0=0;a1=1;cw=0
  alpha=data$alpha
  sigma2 = data$sigma2.alpha
  beta = data$beta
  # mediation effect on log hazard ratio scale under rare outcome
  NIE.la = (beta[2]*alpha[2] + beta[3]*alpha[2]*a1)*(a1-a0)
  NDE.la = (beta[1]+beta[3]*(alpha[1]+alpha[2]*a0+alpha[3]*cw+beta[2]*sigma2))*(a1-a0) + 0.5*beta[3]^2*sigma2*(a1^2-a0^2)
  TE.la = NIE.la + NDE.la
  MP.la=NIE.la/TE.la
  return(c(NIE =NIE.la, NDE = NDE.la,
              TE=TE.la,MP=MP.la))
}

# gold standard approach
# assume we can observe the true exposure among all individuals
gold.f=function(Outcome.formula,Mediator.formula,data,Iname,Tname,Dname,is.gold=1) {
  if (is.gold==0) {
    data$A = data$As
  }
  data.main = data[data[,Iname]==0,]
  data.main = data.main[order(data.main[,Tname]),]
  data=rbind(data.main,data[data[,Iname]==1,])
  n=dim(data)[1]
  n1 = dim(data.main)[1]
  n2=n-n1
  # mediator model
  alpha.lm = lm(Mediator.formula,data)
  alpha = alpha.lm$coefficient
  s2alpha = mean((alpha.lm$residuals)^2)
  Xm = model.matrix(Mediator.formula,data)
  M = data[,all.vars(Mediator.formula)[1]]
  # estimating equation for alpha
  # type 1, bread; type 2, meat
  ee_alpha = function(alpha,type=1) {
    alpha.coef = alpha[-length(alpha)]
    s2alpha= alpha[length(alpha)]
    resi=M - c(Xm %*% alpha.coef)
    out1=Xm * resi
    out2= s2alpha - resi^2
    out = cbind(out1,out2)
    if (type==1) return(apply(out,2,sum))
    if (type==2) return(out)
  }
  alpha = c(alpha,s2alpha)
  Ua = ee_alpha(alpha=alpha,type=2)
  Ia = numDeriv::jacobian(ee_alpha, x=alpha, method="simple", type=1)/n
  cov_a = (solve(Ia) %*% (crossprod(Ua)/n) %*% t(solve(Ia)))/n
  #outcome model
  beta = coxph(Outcome.formula,data=data.main,ties="breslow")$coefficients
  Xt = model.matrix(Outcome.formula,data.main)[,-1]
  Time = data.main[,Tname]
  D = data.main[,Dname]
  # create some quantities for all unique time (for ties)
  unique_t<-unique(Time)
  index_t1 = num_t <- c()
  for (i in (1:length(unique_t))) {
    index=which(Time==unique_t[i])
    num_t[i] = length(index)
    index_t1[i] = index[1]
  }
  # estimating equation for beta
  ee_beta = function(beta,type=1) {
    exp_LP = exp(c(Xt %*% beta))
    E_no = rev(cumsum(rev(exp_LP)))
    E_deno = apply(apply(apply(Xt * exp_LP,2,rev),2,cumsum),2,rev)
    E = E_deno/E_no
    # This is for the ties based on Breslow's approximation
    E = matrix(unlist(rep(as.data.frame(t(E[index_t1,])),num_t)),ncol=length(beta),byrow=T)
    out = (Xt - E) * (D %*% t(rep(1,dim(Xt)[2])))
    if (type==1) return(apply(out,2,sum))
    if (type==2) return(out)
  }
  beta=nleqslv(x=beta,ee_beta)$x
  Ub = ee_beta(beta=beta,type=2)
  Ib = numDeriv::jacobian(ee_beta, x=beta, method="simple", type=1)/n1
  cov_b = (solve(Ib) %*% (crossprod(Ub)/n1) %*% t(solve(Ib)))/n1
  return(list(Point=c(alpha,beta),COV=magic::adiag(cov_a,cov_b)))
}


orc1.f=function(Outcome.formula,Mediator.formula,ME1.formula,
                data,Iname,Tname,Dname,Aname,Asname,Mname,Wname) {
  data.main = data[data[,Iname]==0,]
  data.main = data.main[order(data.main[,Tname]),]
  data.va   = data[data[,Iname]==1,]
  data=rbind(data.main,data.va)
  n=dim(data)[1]
  n1 = dim(data.main)[1]
  n2=n-n1
  # measurement error model
  gamma.lm = lm(ME1.formula,data=data.va)
  gamma.coef = gamma.lm$coefficient
  s2gamma = mean((gamma.lm$residuals)^2)
  Xe = model.matrix(gamma.lm,data.va)
  Ae = data.va[,all.vars(ME1.formula)[1]]
  # estimating equation for alpha
  # type 1, bread; type 2, meat
  ee_gamma = function(gamma,type=1) {
    gamma.coef = gamma[-length(gamma)]
    s2gamma= gamma[length(gamma)]
    resi=Ae - c(Xe %*% gamma.coef)
    out1=Xe * resi
    out2= s2gamma - resi^2
    out3 = cbind(out1,out2)
    out=matrix(0,ncol=dim(out3)[2],nrow=n)
    out[(n1+1):n,] = out3
    if (type==1) return(apply(out,2,sum))
    if (type==2) return(out)
  }
  gamma = c(gamma.coef,s2gamma)
  Ug = ee_gamma(gamma=gamma,type=2)
  Ig = numDeriv::jacobian(ee_gamma, x=gamma, method="simple", type=1)
  #cov_g = (solve(Ig) %*% (crossprod(Ug)/n2) %*% t(solve(Ig)))/n2
  # mediator model
  Xm = model.matrix(Mediator.formula,data)
  Xm[1:n1,Aname] = c(model.matrix(ME1.formula,data.main) %*% gamma[-length(gamma)])
  M = data[,all.vars(Mediator.formula)[1]]
  alpha.lm = lm(M~Xm-1)
  alpha.coef = alpha.lm$coefficients
  s2alpha = (sum((alpha.lm$residuals)^2) - n1*(alpha.coef[2]^2)*gamma[length(gamma)])/n
  alpha = c(alpha.coef,s2alpha)
  # Outcome model
  Xt0 = as.matrix(data.main)
  Xt0[,Aname] = (alpha[2]*(Xt0[,Mname]-alpha[1]- c(Xt0[,Wname,drop=F] %*% alpha[-c(1:2,length(alpha))]))*gamma[length(gamma)]
                 +c(model.matrix(ME1.formula,as.data.frame(Xt0)) %*% gamma[-length(gamma)])*alpha[length(alpha)])/(
                   alpha[2]^2*gamma[length(gamma)]+alpha[length(alpha)])
  Xt0=as.data.frame(Xt0)
  beta = coxph(Outcome.formula,data=Xt0,ties="breslow")$coefficients
  # now we calculate the sandwich variance
  param = c(alpha,beta,gamma)
  param_sep = c(0,cumsum(c(length(alpha),length(beta),length(gamma))))
  # joint estimating equation for alpha and beta
  # type 1, bread; type 2, meat
  ee_alpha_beta = function(param, type=1) {
    alpha = param[(param_sep[1]+1):param_sep[2]]
    beta = param[(param_sep[2]+1):param_sep[3]]
    gamma = param[(param_sep[3]+1):param_sep[4]]
    # ee for mediator model
    Xm[1:n1,Aname] = c(model.matrix(ME1.formula,data.main) %*% gamma[-length(gamma)])
    M = data[,all.vars(Mediator.formula)[1]]
    alpha.coef = alpha[-length(alpha)]
    s2alpha= alpha[length(alpha)]
    resi=M - c(Xm %*% alpha.coef)
    out1=Xm * resi
    offset.term = c(rep((alpha.coef[2]^2)*gamma[length(gamma)],n1),rep(0,n2))
    out2= s2alpha+offset.term - resi^2
    out = cbind(out1,out2)
    # ee for outcome model
    Xt0 = as.matrix(data.main)
    Xt0[,Aname] = (alpha[2]*(Xt0[,Mname]-alpha[1]- c(Xt0[,Wname,drop=F] %*% alpha[-c(1:2,length(alpha))]))*gamma[length(gamma)]
                   +c(model.matrix(ME1.formula,as.data.frame(Xt0)) %*% gamma[-length(gamma)])*alpha[length(alpha)])/(
                     alpha[2]^2*gamma[length(gamma)]+alpha[length(alpha)])
    Xt0=as.data.frame(Xt0)
    Xt = model.matrix(Outcome.formula,Xt0)[,-1]
    Time = data.main[,Tname]
    D = data.main[,Dname]
    unique_t<-unique(Time)
    index_t1 = num_t <- c()
    for (i in (1:length(unique_t))) {
      index=which(Time==unique_t[i])
      num_t[i] = length(index)
      index_t1[i] = index[1]
    }
    exp_LP = exp(c(Xt %*% beta))
    E_no = rev(cumsum(rev(exp_LP)))
    E_deno = apply(apply(apply(Xt * exp_LP,2,rev),2,cumsum),2,rev)
    E = E_deno/E_no
    # This is for the ties based on Breslow's approximation
    E = matrix(unlist(rep(as.data.frame(t(E[index_t1,])),num_t)),ncol=length(beta),byrow=T)
    out3=rbind((Xt - E) * (D %*% t(rep(1,dim(Xt)[2]))),matrix(0,ncol=length(beta),nrow=n2))
    out = cbind(out,out3)
    if (type==1) return(apply(out,2,sum))
    if (type==2) return(out)
  }
  Uab = ee_alpha_beta(param=param,type=2)
  Iabg = numDeriv::jacobian(ee_alpha_beta, x=param, method="simple", type=1)/n
  Iabg_ab = Iabg[,(param_sep[1]+1):param_sep[3]]
  Iabg_g  = Iabg[,(param_sep[3]+1):param_sep[4]]
  Meat = crossprod(Uab + t(Iabg_g %*% solve(Ig) %*% t(Ug)))
  cov_ab = (solve(Iabg_ab) %*% (Meat/n) %*% t(solve(Iabg_ab)))/n
  return(list(Point=c(alpha,beta),COV=cov_ab))
}

orc2.f=function(Outcome.formula,Mediator.formula,ME1.formula,ME2.formula,data,
                Iname,Tname,Dname,Aname,Asname,Mname,Wname) {
  data.main = data[which(data[,Iname]==0),]
  data.main = data.main[order(data.main[,Tname]),]
  data.va   = data[which(data[,Iname]==1),]
  data=rbind(data.main,data.va)
  n=dim(data)[1]
  n1 = dim(data.main)[1]
  n2=n-n1
  ##########################
  # measurement error models
  ##########################
  ## measurement error model 1
  gamma1.lm = lm(ME1.formula,data=data.va)
  gamma1.coef = gamma1.lm$coefficient
  s2gamma1 = mean((gamma1.lm$residuals)^2)
  ## measurement error model 2
  gamma2.lm = lm(ME2.formula,data=data.va)
  gamma2.coef = gamma2.lm$coefficient
  s2gamma2 = mean((gamma2.lm$residuals)^2)
  # design matrices
  Xe1 = model.matrix(gamma1.lm,data.va)
  Xe2 = model.matrix(gamma2.lm,data.va)
  Ae = data.va[,all.vars(ME1.formula)[1]]
  ### all parameters in ME models
  gamma = c(gamma1.coef,s2gamma1,gamma2.coef,s2gamma2)
  gamma_sep = c(0,cumsum(c(length(gamma1.coef)+1,length(gamma2.coef)+1)))
  gamma1 = gamma[(gamma_sep[1]+1):gamma_sep[2]]
  gamma2 = gamma[(gamma_sep[2]+1):gamma_sep[3]]
  # estimating equation for gamma
  # type 1, bread; type 2, meat
  ee_gamma = function(gamma,type=1) {
    gamma1 = gamma[(gamma_sep[1]+1):gamma_sep[2]]
    gamma2 = gamma[(gamma_sep[2]+1):gamma_sep[3]]
    gamma1.coef = gamma1[-length(gamma1)]
    s2gamma1= gamma1[length(gamma1)]
    gamma2.coef = gamma2[-length(gamma2)]
    s2gamma2= gamma2[length(gamma2)]
    # ee for gamma1
    resi=Ae - c(Xe1 %*% gamma1.coef)
    out1=Xe1 * resi
    out2= s2gamma1 - resi^2
    out3 = cbind(out1,out2)
    out4=matrix(0,ncol=dim(out3)[2],nrow=n)
    out4[(n1+1):n,] = out3
    # ee for gamma2
    resi=Ae - c(Xe2 %*% gamma2.coef)
    out5=Xe2 * resi
    out6= s2gamma2 - resi^2
    out7 = cbind(out5,out6)
    out8=matrix(0,ncol=dim(out7)[2],nrow=n)
    out8[(n1+1):n,] = out7
    out=cbind(out4,out8)
    if (type==1) return(apply(out,2,sum))
    if (type==2) return(out)
  }
  Ug = ee_gamma(gamma=gamma,type=2)
  Ig = numDeriv::jacobian(ee_gamma, x=gamma, method="simple", type=1)
  # mediator model
  Xm = model.matrix(Mediator.formula,data)
  Xm[1:n1,Aname] = c(model.matrix(ME1.formula,data.main) %*% gamma1[-length(gamma1)])
  M = data[,all.vars(Mediator.formula)[1]]
  alpha.lm = lm(M~Xm-1)
  alpha.coef = alpha.lm$coefficients
  s2alpha = (sum((alpha.lm$residuals)^2) - n1*(alpha.coef[2]^2)*gamma1[length(gamma1)])/n
  alpha = c(alpha.coef,s2alpha)
  # Outcome model
  Xt0 = as.matrix(data.main)
  Xt0[,Aname] = c(model.matrix(ME2.formula,data.main) %*% gamma2[-length(gamma2)])
  Xt0=as.data.frame(Xt0)
  beta = coxph(Outcome.formula,data=Xt0,ties="breslow")$coefficients
  # now we calculate the sandwich variance
  param = c(alpha,beta,gamma)
  param_sep = c(0,cumsum(c(length(alpha),length(beta),length(gamma1),length(gamma2))))
  # joint estimating equation for alpha and beta
  # type 1, bread; type 2, meat
  ee_alpha_beta = function(param, type=1) {
    alpha = param[(param_sep[1]+1):param_sep[2]]
    beta = param[(param_sep[2]+1):param_sep[3]]
    gamma1 = param[(param_sep[3]+1):param_sep[4]]
    gamma2 = param[(param_sep[4]+1):param_sep[5]]
    # ee for mediator model
    Xm[1:n1,Aname] = c(model.matrix(ME1.formula,data.main) %*% gamma1[-length(gamma1)])
    M = data[,all.vars(Mediator.formula)[1]]
    alpha.coef = alpha[-length(alpha)]
    s2alpha= alpha[length(alpha)]
    resi=M - c(Xm %*% alpha.coef)
    out1=Xm * resi
    offset.term = c(rep((alpha.coef[2]^2)*gamma1[length(gamma1)],n1),rep(0,n2))
    out2= s2alpha+offset.term - resi^2
    out = cbind(out1,out2)
    # ee for outcome model
    Xt0 = as.matrix(data.main)
    Xt0[,Aname] = c(model.matrix(ME2.formula,data.main) %*% gamma2[-length(gamma2)])
    Xt0=as.data.frame(Xt0)
    Xt = model.matrix(Outcome.formula,Xt0)[,-1]
    Time = data.main[,Tname]
    D = data.main[,Dname]
    unique_t<-unique(Time)
    index_t1 = num_t <- c()
    for (i in (1:length(unique_t))) {
      index=which(Time==unique_t[i])
      num_t[i] = length(index)
      index_t1[i] = index[1]
    }
    exp_LP = exp(c(Xt %*% beta))
    E_no = rev(cumsum(rev(exp_LP)))
    E_deno = apply(apply(apply(Xt * exp_LP,2,rev),2,cumsum),2,rev)
    E = E_deno/E_no
    # This is for the ties based on Breslow's approximation
    E = matrix(unlist(rep(as.data.frame(t(E[index_t1,])),num_t)),ncol=length(beta),byrow=T)
    out3=rbind((Xt - E) * (D %*% t(rep(1,dim(Xt)[2]))),matrix(0,ncol=length(beta),nrow=n2))
    out = cbind(out,out3)
    if (type==1) return(apply(out,2,sum))
    if (type==2) return(out)
  }
  Uab = ee_alpha_beta(param=param,type=2)
  Iabg = numDeriv::jacobian(ee_alpha_beta, x=param, method="simple", type=1)/n
  Iabg_ab = Iabg[,(param_sep[1]+1):param_sep[3]]
  Iabg_g  = Iabg[,(param_sep[3]+1):param_sep[5]]
  Meat = crossprod(Uab + t(Iabg_g %*% solve(Ig) %*% t(Ug)))
  cov_ab = (solve(Iabg_ab) %*% (Meat/n) %*% t(solve(Iabg_ab)))/n
  return(list(Point=c(alpha,beta),COV=cov_ab))
}



rrc.f=function(Outcome.formula,Mediator.formula,ME1.formula,ME2.formula,data,Time_sep,
                Iname,Tname,Dname,Aname,Asname,Mname,Wname) {
  data.main = data[which(data[,Iname]==0),]
  data.main = data.main[order(data.main[,Tname]),]
  data.va   = data[which(data[,Iname]==1),]
  data.va = data.va[order(data.va[,Tname]),]
  data=rbind(data.main,data.va)
  n=dim(data)[1]
  n1 = dim(data.main)[1]
  n2=n-n1
  #######################################
  # split the datasets based on Time_sep
  #######################################
  index_t_main = index_t_va = c()
  n_rrc = length(Time_sep)
  for (i in (1:n_rrc)) {
    index_t_main[i] = which(data.main[,Tname]>=Time_sep[i])[1]
    index_t_va[i] = which(data.va[,Tname]>=Time_sep[i])[1]
  }
  #index_t_main[i+1] = n1;index_t_va[i+1] = n2;
  ##########################
  # measurement error models
  ##########################
  ## measurement error model 1
  gamma1.lm = lm(ME1.formula,data=data.va)
  gamma1.coef = gamma1.lm$coefficient
  s2gamma1 = mean((gamma1.lm$residuals)^2)
  ## measurement error model 2
  gamma2 = c()
  for (j in (1:n_rrc)) {
    gamma2.now = lm(ME2.formula,data=data.va[index_t_va[j]:n2,])$coefficients
    gamma2 = c(gamma2,gamma2.now)
  }
  # design matrices
  Xe1 = model.matrix(ME1.formula,data.va)
  Xe2 = model.matrix(ME2.formula,data.va)
  Ae = data.va[,all.vars(ME1.formula)[1]]
  ### all parameters in ME models
  gamma = c(gamma1.coef,s2gamma1,gamma2)
  gamma_sep = c(0,cumsum(c(length(gamma1.coef)+1,length(gamma2))))
  gamma1 = gamma[(gamma_sep[1]+1):gamma_sep[2]]
  gamma2 = gamma[(gamma_sep[2]+1):gamma_sep[3]]
  n_gamma2 = length(gamma2)/n_rrc
  # estimating equation for gamma
  # type 1, bread; type 2, meat
  ee_gamma = function(gamma,type=1) {
    gamma1 = gamma[(gamma_sep[1]+1):gamma_sep[2]]
    gamma2 = gamma[(gamma_sep[2]+1):gamma_sep[3]]
    gamma1.coef = gamma1[-length(gamma1)]
    s2gamma1= gamma1[length(gamma1)]
    # ee for gamma1
    resi=Ae - c(Xe1 %*% gamma1.coef)
    out1=Xe1 * resi
    out2= s2gamma1 - resi^2
    out3 = cbind(out1,out2)
    out4=matrix(0,ncol=dim(out3)[2],nrow=n)
    out4[(n1+1):n,] = out3
    # ee for gamma2
    for (j in (1:n_rrc)) {
      index.now = index_t_va[j]:n2
      gamma2.now = gamma2[((j-1)*n_gamma2+1):(n_gamma2*j)]
      resi=rep(0,n2)
      resi[index.now] = c(Ae - c(Xe2 %*% gamma2.now))[index.now]
      if (j == 1) {
        out5 = Xe2 * resi
      } else {
        out5 = cbind(out5,Xe2 * resi)
      }
    }
    out6=matrix(0,ncol=dim(out5)[2],nrow=n)
    out6[(n1+1):n,] = out5
    out=cbind(out4,out6)
    if (type==1) return(apply(out,2,sum))
    if (type==2) return(out)
  }
  Ug = ee_gamma(gamma=gamma,type=2)
  Ig = numDeriv::jacobian(ee_gamma, x=gamma, method="simple", type=1)
  # mediator model
  Xm = model.matrix(Mediator.formula,data)
  Xm[1:n1,Aname] = c(model.matrix(ME1.formula,data.main) %*% gamma1[-length(gamma1)])
  M = data[,all.vars(Mediator.formula)[1]]
  alpha.lm = lm(M~Xm-1)
  alpha.coef = alpha.lm$coefficients
  s2alpha = (sum((alpha.lm$residuals)^2) - n1*(alpha.coef[2]^2)*gamma1[length(gamma1)])/n
  alpha = c(alpha.coef,s2alpha)
  # Outcome model
  Xt0 = as.matrix(data.main)
  index_t_main=c(index_t_main,n1)
  for (j in (1:n_rrc)) {
    index.now = index_t_main[j]:index_t_main[j+1]
    gamma2.now = gamma2[((j-1)*n_gamma2+1):(n_gamma2*j)]
    Xt0[index.now,Aname] = c(model.matrix(ME2.formula,data.main) %*% gamma2.now)[index.now]
  }
  Xt0=as.data.frame(Xt0)
  beta = coxph(Outcome.formula,data=Xt0,ties="breslow")$coefficients
  # now we calculate the sandwich variance
  param = c(alpha,beta,gamma)
  param_sep = c(0,cumsum(c(length(alpha),length(beta),length(gamma1),length(gamma2))))
  # joint estimating equation for alpha and beta
  # type 1, bread; type 2, meat
  ee_alpha_beta = function(param, type=1) {
    alpha = param[(param_sep[1]+1):param_sep[2]]
    beta = param[(param_sep[2]+1):param_sep[3]]
    gamma1 = param[(param_sep[3]+1):param_sep[4]]
    gamma2 = param[(param_sep[4]+1):param_sep[5]]
    # ee for mediator model
    Xm[1:n1,Aname] = c(model.matrix(ME1.formula,data.main) %*% gamma1[-length(gamma1)])
    M = data[,all.vars(Mediator.formula)[1]]
    alpha.coef = alpha[-length(alpha)]
    s2alpha= alpha[length(alpha)]
    resi=M - c(Xm %*% alpha.coef)
    out1=Xm * resi
    offset.term = c(rep((alpha.coef[2]^2)*gamma1[length(gamma1)],n1),rep(0,n2))
    out2= s2alpha+offset.term - resi^2
    out = cbind(out1,out2)
    # ee for outcome model
    Xt0 = as.matrix(data.main)
    for (j in (1:n_rrc)) {
      index.now = index_t_main[j]:index_t_main[j+1]
      gamma2.now = gamma2[((j-1)*n_gamma2+1):(n_gamma2*j)]
      Xt0[index.now,Aname] = c(model.matrix(ME2.formula,data.main) %*% gamma2.now)[index.now]
    }
    Xt0=as.data.frame(Xt0)
    Xt = model.matrix(Outcome.formula,Xt0)[,-1]
    Time = data.main[,Tname]
    D = data.main[,Dname]
    unique_t<-unique(Time)
    index_t1 = num_t <- c()
    for (i in (1:length(unique_t))) {
      index=which(Time==unique_t[i])
      num_t[i] = length(index)
      index_t1[i] = index[1]
    }
    exp_LP = exp(c(Xt %*% beta))
    E_no = rev(cumsum(rev(exp_LP)))
    E_deno = apply(apply(apply(Xt * exp_LP,2,rev),2,cumsum),2,rev)
    E = E_deno/E_no
    # This is for the ties based on Breslow's approximation
    E = matrix(unlist(rep(as.data.frame(t(E[index_t1,])),num_t)),ncol=length(beta),byrow=T)
    out3=rbind((Xt - E) * (D %*% t(rep(1,dim(Xt)[2]))),matrix(0,ncol=length(beta),nrow=n2))
    out = cbind(out,out3)
    if (type==1) return(apply(out,2,sum))
    if (type==2) return(out)
  }
  Uab = ee_alpha_beta(param=param,type=2)
  Iabg = numDeriv::jacobian(ee_alpha_beta, x=param, method="simple", type=1)/n
  Iabg_ab = Iabg[,(param_sep[1]+1):param_sep[3]]
  Iabg_g  = Iabg[,(param_sep[3]+1):param_sep[5]]
  Meat = crossprod(Uab + t(Iabg_g %*% solve(Ig) %*% t(Ug)))
  cov_ab = (solve(Iabg_ab) %*% (Meat/n) %*% t(solve(Iabg_ab)))/n
  return(list(Point=c(alpha,beta),COV=cov_ab))
}

get_estimated_effect = function(res,Mname,Aname,Wname,a0=0,a1=1,cw=0) {
  #res=gold.e;a0=0;a1=1;cw=0
  n_alpha = length(Wname)+3
  n_beta  = length(Wname)+3
  param = res$Point
  param_cov = res$COV
  names(param)=c(paste("alpha",0:(n_alpha-2),sep=""),"s2alpha",paste("beta",1:n_beta,sep=""))
  colnames(param_cov)=rownames(param_cov) = names(param)
  name_alpha_w = names(param)[3:(n_alpha-1)]
  name_beta_w   = names(param)[(n_alpha+4):(n_alpha+n_beta)]
  
  alpha0=param[1];alpha1=param[2];s2alpha = param[n_alpha]
  beta1 = param[n_alpha+1];beta2 = param[n_alpha+2];beta3 = param[n_alpha+3]
  for (i in (1:length(Wname))) {
    assign(name_alpha_w[i],param[2+i])
    assign(name_beta_w[i],param[n_alpha+3+i])
  }
  
  # calculate the mediation measures
  .A = "(beta2*alpha1+beta3*alpha1*a1)*(a1-a0)"
  .B = paste("beta1*(a1-a0)+beta3*(alpha0+alpha1*a0+beta2*s2alpha)*(a1-a0)+0.5*(beta3^2)*s2alpha*(a1^2-a0^2)",
             paste("+beta3*(",paste(paste(name_alpha_w,"*",cw),collapse = "+"),")*(a1-a0)"))
  # NIE
  NIE_fun = function() {output = .A;return(output)}
  variable=names(param)
  NIE_D=deriv(parse(text=NIE_fun()),variable)
  NIE_D = eval(NIE_D)
  NIE_p = NIE_D[1]
  lambda= t(attr(NIE_D,"gradient"))
  V_NIE = as.vector(t(lambda) %*% param_cov %*% lambda)
  
  # NDE
  NDE_fun = function() {output = .B;return(output)}
  variable=names(param)
  NDE_D=deriv(parse(text=NDE_fun()),variable)
  NDE_D = eval(NDE_D)
  NDE_p = NDE_D[1]
  lambda= t(attr(NDE_D,"gradient"))
  V_NDE = as.vector(t(lambda) %*% param_cov %*% lambda)
  
  # TE
  TE_fun = function() {output = paste(.A,"+",.B,sep="");return(output)}
  variable=names(param)
  TE_D=deriv(parse(text=TE_fun()),variable)
  TE_D = eval(TE_D)
  TE_p = TE_D[1]
  lambda= t(attr(TE_D,"gradient"))
  V_TE = as.vector(t(lambda) %*% param_cov %*% lambda)
  
  # MP
  MP_fun = function() {output = paste("(",.A,")/(",.A,"+",.B,")",sep="");return(output)}
  variable=names(param)
  MP_D=deriv(parse(text=MP_fun()),variable)
  MP_D = eval(MP_D)
  MP_p = MP_D[1]
  lambda= t(attr(MP_D,"gradient"))
  V_MP = as.vector(t(lambda) %*% param_cov %*% lambda)
  list(Point=c(NIE=NIE_p,NDE=NDE_p,TE=TE_p,MP=MP_p),
       SE=c(NIE=sqrt(V_NIE),NDE=sqrt(V_NDE),TE=sqrt(V_TE),MP=sqrt(V_MP)))
}


combined.f=function(Outcome.formula,Mediator.formula,ME1.formula,ME2.formula,data,Time_sep,
                    Iname,Tname,Dname,Aname,Asname,Mname,Wname) {
  
  naive.e = gold.f(Outcome.formula,Mediator.formula,data,Iname,Tname,Dname,is.gold=0)
  gold.e =gold.f(Outcome.formula,Mediator.formula,data,Iname,Tname,Dname,is.gold=1)
  orc1.e=orc1.f(Outcome.formula,Mediator.formula,ME1.formula,data,Iname,Tname,Dname,Aname,Asname,Mname,Wname)
  orc2.e=orc2.f(Outcome.formula,Mediator.formula,ME1.formula,ME2.formula,data,Iname,Tname,Dname,Aname,Asname,Mname,Wname)
  rrc.e=rrc.f(Outcome.formula,Mediator.formula,ME1.formula,ME2.formula,
              data,Time_sep,Iname,Tname,Dname,Aname,Asname,Mname,Wname)
  
  naive.res = get_estimated_effect(res=naive.e,Mname,Aname,Wname,a0=0,a1=1,cw=0)
  gold.res = get_estimated_effect(res=gold.e,Mname,Aname,Wname,a0=0,a1=1,cw=0)
  orc1.res = get_estimated_effect(res=orc1.e,Mname,Aname,Wname,a0=0,a1=1,cw=0)
  orc2.res = get_estimated_effect(res=orc2.e,Mname,Aname,Wname,a0=0,a1=1,cw=0)
  rrc.res = get_estimated_effect(res=rrc.e,Mname,Aname,Wname,a0=0,a1=1,cw=0)
  
  Naive.res = rbind(naive.res$Point,naive.res$Point-1.96*naive.res$SE,naive.res$Point+1.96*naive.res$SE)
  Gold.res  = rbind(gold.res$Point,gold.res$Point-1.96*gold.res$SE,gold.res$Point+1.96*gold.res$SE)
  ORC1.res = rbind(orc1.res$Point,orc1.res$Point-1.96*orc1.res$SE,orc1.res$Point+1.96*orc1.res$SE)
  ORC2.res = rbind(orc2.res$Point,orc2.res$Point-1.96*orc2.res$SE,orc2.res$Point+1.96*orc2.res$SE)
  RRC.res = rbind(rrc.res$Point,rrc.res$Point-1.96*rrc.res$SE,rrc.res$Point+1.96*rrc.res$SE)
  
  rownames(Naive.res)=rownames(Gold.res)=rownames(ORC1.res)=rownames(ORC2.res)=rownames(RRC.res)=c("point","CI.lower","CI.upper")
  colnames(Naive.res)=colnames(Gold.res)=colnames(ORC1.res)=colnames(ORC2.res)=colnames(RRC.res)=c("NIE","NDE","TE","MP")
  list(Naive=Naive.res,Gold=Gold.res,ORC1=ORC1.res,ORC2=ORC2.res,RRC=RRC.res)
}

simulation.f=function(n=10000,n2=500,replication=1000,rho=0.3) {
  #n=10000;n2=500;replication=5;rho=0.3
  output1=output2=output3=matrix(0,ncol=20,nrow=replication)
  
  Iname="I"
  Tname="Time"
  Dname="Event"
  Aname="A"
  Asname = "As"
  Mname = "M"
  Wname="W"
  Time_sep=seq(0,40,by=10)
  Outcome.formula= as.formula(Surv(Time,Event) ~ A + M + I(A*M) + W)
  Mediator.formula = as.formula(M ~ A + W)
  ME1.formula = as.formula(A ~ As + W)
  ME2.formula = as.formula(A ~ As + M + W)

  for (i in (1:replication)) {
    data.raw = gen_data(n=n,n2=n2,rho=rho,corma=0.4,beta123=c(log(1.5),log(1.5),log(1.2)))
    data = data.raw$data
    true_effect = get_true_effect.f(data=data.raw,a0=0, a1=1, cw=0)
    allres =  combined.f(Outcome.formula,Mediator.formula,ME1.formula,ME2.formula,data,Time_sep,
                         Iname,Tname,Dname,Aname,Asname,Mname,Wname)
    # point estimator
    output1[i,]=c(allres$Naive[1,],allres$Gold[1,],allres$ORC1[1,],allres$ORC2[1,],allres$RRC[1,])
    # CI.lower
    output2[i,]=c(allres$Naive[2,],allres$Gold[2,],allres$ORC1[2,],allres$ORC2[2,],allres$RRC[2,])
    output3[i,]=c(allres$Naive[3,],allres$Gold[3,],allres$ORC1[3,],allres$ORC2[3,],allres$RRC[3,])
  }
  # Bias SE and CR

  BIAS0=(apply(output1,2,median,na.rm=T)-rep(true_effect,5))/rep(true_effect,5)*100
  SE0 = apply(output1,2,sd,na.rm=T)
  par.mat=matrix(rep(true_effect,5*replication),ncol=20,nrow=replication,byrow=T)
  CR0 = apply(((par.mat-output2)>0)*((par.mat-output3)<0),2,mean,na.rm=T)*100
  res0=round(rbind(BIAS0,SE0,CR0),3)
  
  BIAS1=matrix(res0[1,],ncol=5,nrow=4,byrow=F)
  SE1=matrix(res0[2,],ncol=5,nrow=4,byrow=F)
  CR1=matrix(res0[3,],ncol=5,nrow=4,byrow=F)
  colnames(BIAS1)=colnames(SE1)=colnames(CR1)=c("Naive","Gold","ORC1","ORC2","RRC")
  rownames(BIAS1)=rownames(SE1)=rownames(CR1)=c("NIE","NDE","TE","MP")
  list(BIAS=BIAS1,SE=SE1,CR=CR1)
}



aid <- paste(Sys.getenv("PBS_ARRAYID"),
             Sys.getenv("SLURM_ARRAY_TASK_ID"),
             sep = ""
)
aid = as.numeric(aid)[1]

set.seed(aid+62622)

Sys.time()
res1=simulation.f(n=10000,n2=250,replication=100,rho=0.25)
print(res1)
Sys.time()

Sys.time()
res2=simulation.f(n=10000,n2=250,replication=100,rho=0.5)
print(res2)
Sys.time()

Sys.time()
res3=simulation.f(n=10000,n2=250,replication=100,rho=0.75)
print(res3)
Sys.time()

