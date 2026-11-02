#### Functions for Different Latency Weight Curve #### 

library(npreg)

##### Constant Latency Weight Model #####

lag1<-function(x,t1,t2,w){
  y<-c()
  for(i in 1:length(x)){
    if(x[i]<t1||x[i]>t2){
      y_add=0
    }
    else{
      y_add=w
    }
    y<-c(y,y_add)
  }
  return(y)
}

##### Linear Latency Weight Model #####

lag2<-function(x,t1,t2,tp,w1,w2){
  y<-c()
  for(i in 1:length(x)){
    if(x[i]<t1||x[i]>t2){
      y_add<-0
    }
    else if(x[i]>=t1&&x[i]<=tp){
      y_add<-(1-w1)/(tp-t1)*x[i]+(w1*tp-t1)/(tp-t1)
    }
    else if(x[i]>tp&&x[i]<=t2){
      y_add<-(w2-1)/(t2-tp)*x[i]+(t2-w2*tp)/(t2-tp)
    }
    y<-c(y,y_add)
  }
  return(y)
}

##### Polynomial Latency Weight Model Up to Quadratic Terms #####

lag3<-function(x,a0,a1,a2){
  y<-c()
  for(i in 1:length(x)){
    y_add<-a0+a1*x[i]+a2*x[i]^2
    y<-c(y,y_add)
  }
  return(y)
}

##### Exponential Latency Weight Model Up to Cubic Terms With No Intercept #####

lag4<-function(x,l1,l2,l3){
  y<-c()
  for(i in 1:length(x)){
    y_add=exp(2*(l1*x[i]+l2*x[i]^2+l3*x[i]^3))
    y<-c(y,y_add)
  }
  return(y)
}

##### Exponential Latency Weight Model Up to Cubic Terms With Intercept #####

lag4_1<-function(x,l0,l1,l2,l3){
  y<-c()
  for(i in 1:length(x)){
    y_add=exp(2*(l0+l1*x[i]+l2*x[i]^2+l3*x[i]^3))
    y<-c(y,y_add)
  }
  return(y)
}

##### Exponential Latency Weight Model Up to Quadratic Terms With Intercept #####

lag4_2<-function(x,l0,l1,l2){
  y<-c()
  for(i in 1:length(x)){
    y_add=exp(2*(l0+l1*x[i]+l2*x[i]^2))
    y<-c(y,y_add)
  }
  return(y)
}

##### Should be Spline Model but does not need self-written function #####

lag5<-function(){
  
}

##### Two-Stage Constant Latency Weight Model #####

lag6<-function(x1,w1,w2,t1,t2,t3){
y<-c()
for(i in 1:t1){
  y<-c(y,0)
}
for(i in (t1+1):t2){
  y<-c(y,1)
}
for(i in (t2+1):t3){
  y<-c(y,w1)
}
for(i in (t3+1):56){
  y<-c(y,w2)
}
return(y)
}


x1<-seq(0,55,by=1)

#### Visualization Trial ####
### peak at 8-15, almost 0 at 30, no lag Cox corrected multivariate=4.69


plot(x1,rep(0,length(x1)),type='n',
     xlab="Latency Time",ylab="Weighted Effect Coefficient",main="Visualization of Various Latency Curves",
     xlim=c(0,55),ylim=c(0,3))
lines(x1,log(4.69)*rep(1,length(x1)),type='l',col='black')
lines(x1,log(5.25)*lag1(x1,t1=2,t2=51,w=1),type='l',col='red')
lines(x1,log(10.51)*lag2(x1,t1=0,t2=53,tp=12,w1=0.3,w2=0),type='l',col='green')
#lines(x1,pmax(lag3(x1,a0=2.9,a1=-0.035,a2=11),0),type='l',col='blue')
lines(x1,lag4_2(x1,l0=0.31,l1=0.041,l2=-0.0021),type='l',col='orange')
mod.ss<-ss(x1,lag4_2(x1,l0=0.38,l1=0.0350,l2=-0.0019),nknots=6)
plot(mod.ss,add=T,col="purple",ci=F)

lines(x1,2.74*lag6(x1,w1=0.48,w2=0.14,t1=3,t2=17,t3=28),type='l',col='blue')
legend("topright",
       legend=c("No Latency","Latency Model 1","Latency Model 2","Latency Model 4-2", "Latency Spline Model","Latency Model 6"),
       col=c("black","red","green","orange","purple",'blue'),lty=1,cex=0.7)


