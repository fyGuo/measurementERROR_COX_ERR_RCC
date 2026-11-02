library(haven)

valid_data<-read_sas("~/OneDrive - Harvard University/JUL2015-FASRC/data/valid_2009_short.sas7bdat")

valid_data$VALUE100<-valid_data$VALUE/100
valid_data$EVALUE100<-valid_data$EVALUE/100
valid_data$yy1<-ifelse(valid_data$YEAR<70,1,0)
valid_data$yy2<-ifelse(valid_data$YEAR>=70&valid_data$YEAR<80,1,0)
valid_data$yy3<-ifelse(valid_data$YEAR>=80&valid_data$YEAR<90,1,0)
valid_data$yy4<-ifelse(valid_data$YEAR>=90&valid_data$YEAR<100,1,0)
valid_data$yy5<-ifelse(valid_data$YEAR>=100,1,0)

# Some preliminary analysis
length(unique(valid_data$ID))
# 862
hist(valid_data$VALUE100)
# peak at 0, decreasing slope --> probably not normal, unless normal around 0 or <0
hist(valid_data$EVALUE100)
# double peak at 0 and 0.05, down-up-down, decreasing slope --> not normal

sum(valid_data$VALUE100==0)
# 336
sum(valid_data$EVALUE100==0)
# 66
sum(valid_data$VALUE100==0&valid_data$EVALUE100==0)
# 9
View(valid_data[which(valid_data$VALUE100==0&valid_data$EVALUE100==0),])

# Separate into training and testing group
# Fit model on training group
# Assess model performance on testing group
# Report MSE, R-squared
# Draw residual plot

library(caret)
# Additive model
ctrl.add<-trainControl(method="cv",number=5)
model.add<-train(VALUE~EVALUE+SMOKE+NATIVE+HISPANIC+COMP5+COMP6+yy2, data=valid_data, method="lm",trControl=ctrl.add)
print(model.add)
#Linear Regression 
#2833 samples
#7 predictor
#No pre-processing
#Resampling: Cross-Validated (5 fold) 
#Summary of sample sizes: 2267, 2267, 2266, 2266, 2266 
#Resampling results:
#  RMSE     Rsquared   MAE     
#1.87208  0.2456883  1.338731
#Tuning parameter 'intercept' was held constant at a value of TRUE

emcov_data<-valid_data[,c("SMOKE","NATIVE","HISPANIC","COMP5","COMP6","yy2","EVALUE_LOGPOS")]
VALUE_LOGPOS_PREDICT<-predict(model.mul$finalModel, data=emcov_data)
VALUE_PREDICT<-exp(VALUE_LOGPOS_PREDICT)-min(valid_data$VALUE[valid_data$VALUE>0])/2
RMSE.mul<-mean((valid_data$VALUE-VALUE_PREDICT)^2)^(0.5)
# 2.104287

em.lm1<-lm()

# transform to log

# Multiplicative model
valid_data$VALUE_LOGPOS<-log(valid_data$VALUE+min(valid_data$VALUE[valid_data$VALUE>0])/2)
valid_data$EVALUE_LOGPOS<-log(valid_data$EVALUE+min(valid_data$EVALUE[valid_data$EVALUE>0])/2)
ctrl.mul<-trainControl(method="cv",number=5)
model.mul<-train(VALUE_LOGPOS~EVALUE_LOGPOS+SMOKE+NATIVE+HISPANIC+COMP5+COMP6+yy2, data=valid_data, method="lm",trControl=ctrl.mul)
print(model.mul)
#Linear Regression 
#2833 samples
#7 predictor
#No pre-processing
#Resampling: Cross-Validated (5 fold) 
#Summary of sample sizes: 2266, 2266, 2266, 2267, 2267 
#Resampling results:
#  RMSE      Rsquared   MAE     
#1.356156  0.1016481  1.085785
#Tuning parameter 'intercept' was held constant at a value of TRUE

emcov_data<-valid_data[,c("SMOKE","NATIVE","HISPANIC","COMP5","COMP6","yy2","EVALUE_LOGPOS")]
VALUE_LOGPOS_PREDICT<-predict(model.mul$finalModel, data=emcov_data)
VALUE_PREDICT<-exp(VALUE_LOGPOS_PREDICT)-min(valid_data$VALUE[valid_data$VALUE>0])/2
RMSE.mul<-mean((valid_data$VALUE-VALUE_PREDICT)^2)^(0.5)
# 2.104287

# Multiplicative model for positive exposure range
valid_data_positive<-valid_data[which(valid_data$VALUE>0&valid_data$EVALUE>0),]
valid_data_positive$VALUE_LOGPOS<-log(valid_data_positive$VALUE)
valid_data_positive$EVALUE_LOGPOS<-log(valid_data_positive$EVALUE)
# 2440 rows
ctrl.mul.pos<-trainControl(method="cv",number=5)
model.mul.pos<-train(VALUE_LOGPOS~EVALUE_LOGPOS+SMOKE+NATIVE+HISPANIC+COMP5+COMP6+yy2, data=valid_data_positive, method="lm",trControl=ctrl.mul.pos)
print(model.mul.pos)
#Linear Regression 
#2440 samples
#7 predictor
#No pre-processing
#Resampling: Cross-Validated (5 fold) 
#Summary of sample sizes: 1952, 1953, 1951, 1951, 1953 
#Resampling results:
#  RMSE       Rsquared   MAE      
#0.9891032  0.1498305  0.7800505
#Tuning parameter 'intercept' was held constant at a value of TRUE

plot(model.mul.pos$finalModel)

emcov_data.pos<-valid_data_positive[,c("SMOKE","NATIVE","HISPANIC","COMP5","COMP6","yy2","EVALUE_LOGPOS")]
VALUE_LOGPOS_PREDICT<-predict(model.mul.pos$finalModel, data=emcov_data.pos)
VALUE_PREDICT<-exp(VALUE_LOGPOS_PREDICT)
RMSE.mul.pos<-mean((valid_data_positive$VALUE-VALUE_PREDICT)^2)^(0.5)
# 2.025913

plot(valid_data_positive$EVALUE, valid_data_positive$VALUE-VALUE_PREDICT, type="p")
plot(valid_data_positive$EVALUE_LOGPOS, valid_data_positive$VALUE_LOGPOS-VALUE_LOGPOS_PREDICT, type="p")


# Multiplicative model for positive, with squared term
valid_data_positive<-valid_data[which(valid_data$VALUE>0&valid_data$EVALUE>0),]
valid_data_positive$VALUE_LOGPOS<-log(valid_data_positive$VALUE)
valid_data_positive$EVALUE_LOGPOS<-log(valid_data_positive$EVALUE)
# 2440 rows
ctrl.mul.pos<-trainControl(method="cv",number=5)
model.mul.pos.sq<-train(VALUE_LOGPOS~EVALUE_LOGPOS+I(EVALUE_LOGPOS^2)+SMOKE+NATIVE+HISPANIC+COMP5+COMP6+yy2, data=valid_data_positive, method="lm",trControl=ctrl.mul.pos)
print(model.mul.pos.sq)

plot(model.mul.pos.sq$finalModel)

emcov_data.pos<-valid_data_positive[,c("SMOKE","NATIVE","HISPANIC","COMP5","COMP6","yy2","EVALUE_LOGPOS")]
VALUE_LOGPOS_PREDICT<-predict(model.mul.pos$finalModel, data=emcov_data.pos)
VALUE_PREDICT<-exp(VALUE_LOGPOS_PREDICT)
RMSE.mul.pos<-mean((valid_data_positive$VALUE-VALUE_PREDICT)^2)^(0.5)



# additive ME
model.add.nocov<-lm(VALUE~EVALUE, data=valid_data)
# R-squared: 0.1116
valid_data$VALUE.add.nocov<-fitted(model.add.nocov)
model.add.cov<-lm(VALUE~EVALUE+SMOKE+NATIVE+HISPANIC+COMP5+COMP6+yy2, data=valid_data)
# R-squared: 0.2511, NATIVE not significant at 0.05
valid_data$VALUE.add.cov<-fitted(model.add.cov)

# multiplicative ME
valid_data$logVALUE.pos<-log(valid_data$VALUE+min(valid_data$VALUE[valid_data$VALUE>0])/2)
valid_data$logEVALUE.pos<-log(valid_data$EVALUE+min(valid_data$EVALUE[valid_data$EVALUE>0])/2)
model.mul.nocov<-lm(logVALUE.pos~logEVALUE.pos,data=valid_data)
# R-squared: 0.03764
valid_data$VALUE.mul.nocov<-exp(fitted(model.mul.nocov))-min(valid_data$VALUE[valid_data$VALUE>0])/2
model.mul.cov<-lm(logVALUE.pos~logEVALUE.pos+SMOKE+NATIVE+HISPANIC+COMP5+COMP6+yy2, data=valid_data)
# R-squared: 0.1058
# Will log-transformation decrease R-squared? 
valid_data$VALUE.mul.cov<-exp(fitted(model.mul.cov))-min(valid_data$VALUE[valid_data$VALUE>0])/2

# plot of fitted VS true
valid_data.fitted<-data.frame(EVALUE=valid_data$EVALUE,
                              VALUE=valid_data$VALUE,
                              VALUE.add.nocov=valid_data$VALUE.add.nocov,
                              VALUE.add.cov=valid_data$VALUE.add.cov,
                              VALUE.mul.nocov=valid_data$VALUE.mul.nocov,
                              VALUE.mul.cov=valid_data$VALUE.mul.cov)
plot.fittedVSevalue<-valid_data.fitted %>% tidyr::gather("id", "fitted", 2:6) %>% 
  ggplot(., aes(EVALUE, fitted))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE, color="blue")+
  ggtitle("True and Fitted Value VS Error-prone Value")+
  facet_wrap(~id)

plot.fittedVSevalue

#plot(model.add.nocov)
plot(model.add.cov)
#plot(model.mul.nocov)
plot(model.mul.cov)

RMSEtotal<-c(0)
for(j in 3:6){
  RMSEtotal<-c(RMSEtotal,sqrt(mean((valid_data.fitted[,j]-valid_data.fitted$VALUE)^2)))
}
RMSEtotal<-RMSEtotal[-c(1)]
RMSEtotal
#add.nocov    add.cov  mul.nocov    mul.cov
#2.030557 1.864325 2.306218 2.104287
sqrt(mean(model.mul.cov$residuals^2))
# 1.352035
sqrt(mean(model.mul.nocov$residuals^2))
# 1.402653




# fit by calendar years
valid_data_2<-valid_data
valid_data_2$logVALUE.pos<-log(valid_data_2$VALUE+min(valid_data_2$VALUE[valid_data_2$VALUE>0])/2)
valid_data_2$logEVALUE.pos<-log(valid_data_2$EVALUE+min(valid_data_2$EVALUE[valid_data_2$EVALUE>0])/2)
valid_data_2$VALUE.add.nocov<-0
valid_data_2$VALUE.add.cov<-0
valid_data_2$VALUE.mul.nocov<-0
valid_data_2$VALUE.mul.cov<-0

RMSE_2<-c(0,0,0,0)
Rsquared_2<-c(0,0,0,0)
AdjRsquared_2<-c(0,0,0,0)
for(i in 67:76){
  rowselect<-which(valid_data_2$YEAR==i)
  if(length(rowselect)==0){
    cat("skip year ", i, " because sample size=0\n")
    next
  }
  else{
    cat("year", i, " sample size=", length(rowselect),"\n")
  }
  valid_data_2.select<-valid_data_2[rowselect,]
  
  model.add.nocov<-lm(VALUE~EVALUE, data=valid_data_2.select)
  valid_data_2$VALUE.add.nocov[rowselect]<-fitted(model.add.nocov)
  model.add.cov<-lm(VALUE~EVALUE+SMOKE+NATIVE+HISPANIC+COMP5+COMP6+yy2, data=valid_data_2.select)
  valid_data_2$VALUE.add.cov[rowselect]<-fitted(model.add.cov)
  
  model.mul.nocov<-lm(logVALUE.pos~logEVALUE.pos,data=valid_data_2.select)
  valid_data_2$VALUE.mul.nocov[rowselect]<-exp(fitted(model.mul.nocov))-min(valid_data_2$VALUE[valid_data_2$VALUE>0])/2
  model.mul.cov<-lm(logVALUE.pos~logEVALUE.pos+SMOKE+NATIVE+HISPANIC+COMP5+COMP6+yy2, data=valid_data_2.select)
  valid_data_2$VALUE.mul.cov[rowselect]<-exp(fitted(model.mul.cov))-min(valid_data_2$VALUE[valid_data_2$VALUE>0])/2
  
  RMSE_2<-rbind(RMSE_2, c(sqrt(mean(model.add.nocov$residuals^2)),
                          sqrt(mean(model.add.cov$residuals^2)),
                          sqrt(mean(model.mul.nocov$residuals^2)),
                          sqrt(mean(model.mul.cov$residuals^2))))
  Rsquared_2<-rbind(Rsquared_2, c(summary(model.add.nocov)$r.squared,
                                  summary(model.add.cov)$r.squared,
                                  summary(model.mul.nocov)$r.squared,
                                  summary(model.mul.cov)$r.squared))
  AdjRsquared_2<-rbind(Rsquared_2, c(summary(model.add.nocov)$adj.r.squared,
                                     summary(model.add.cov)$adj.r.squared,
                                     summary(model.mul.nocov)$adj.r.squared,
                                     summary(model.mul.cov)$adj.r.squared))
}

valid_data_2.fitted<-data.frame(EVALUE=valid_data_2$EVALUE,
                              VALUE=valid_data_2$VALUE,
                              VALUE.add.nocov=valid_data_2$VALUE.add.nocov,
                              VALUE.add.cov=valid_data_2$VALUE.add.cov,
                              VALUE.mul.nocov=valid_data_2$VALUE.mul.nocov,
                              VALUE.mul.cov=valid_data_2$VALUE.mul.cov)


library(ggplot2)
plot.fittedVSevalue_2<-valid_data_2.fitted %>% tidyr::gather("id", "fitted", 2:6) %>% 
  ggplot(., aes(EVALUE, fitted))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE, color="blue")+
  ggtitle("True and Fitted Value in Calendar Year Model VS Error-prone Value")+
  facet_wrap(~id)

plot.fittedVSevalue_2

plot.fittedVStrue_2<-valid_data_2.fitted %>% tidyr::gather("id", "fitted", 3:6) %>% 
  ggplot(., aes(VALUE, fitted))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE, color="blue")+
  ggtitle("Fitted Value in Calendar Year Model VS True Value")+
  facet_wrap(~id)

plot.fittedVStrue_2

RMSE_2<-RMSE_2[-c(1),]
Rsquared_2<-Rsquared_2[-c(1),]
AdjRsquared_2<-AdjRsquared_2[-c(1),]
colnames(RMSE_2)<-c("add.nocov","add.cov","mul.nocov","mul.cov")
colnames(Rsquared_2)<-c("add.nocov","add.cov","mul.nocov","mul.cov")
colnames(AdjRsquared_2)<-c("add.nocov","add.cov","mul.nocov","mul.cov")
colMeans(Rsquared_2)
#add.nocov    add.cov  mul.nocov    mul.cov 
#0.09939920 0.23878039 0.05728835 0.15325118
colMeans(AdjRsquared_2)
#add.nocov    add.cov  mul.nocov    mul.cov 
#0.08945928 0.20116932 0.05155951 0.11481654

RMSEtotal_2<-c(0)
for(j in 3:6){
  RMSEtotal_2<-c(RMSEtotal_2,sqrt(mean((valid_data_2.fitted[,j]-valid_data_2.fitted$VALUE)^2)))
}
RMSEtotal_2<-RMSEtotal_2[-c(1)]
RMSEtotal_2
#add.nocov    add.cov  mul.nocov    mul.cov
#1.876166 1.640864 2.153154 1.862894



# fit by calendar years

valid_data_3<-valid_data
valid_data_3$logVALUE.pos<-log(valid_data_3$VALUE+min(valid_data_3$VALUE[valid_data_3$VALUE>0])/2)
valid_data_3$logEVALUE.pos<-log(valid_data_3$EVALUE+min(valid_data_3$EVALUE[valid_data_3$EVALUE>0])/2)
valid_data_3$VALUE.add.nocov<-0
valid_data_3$VALUE.add.cov<-0
valid_data_3$VALUE.mul.nocov<-0
valid_data_3$VALUE.mul.cov<-0

library(haven)
final_data<-read_sas("~/OneDrive - Harvard University/JUL2015-FASRC/data/final_data.sas7bdat")
final_data$YR_LAST<-as.numeric(final_data$year)-1900
final_data_1<-final_data[which((is.na(final_data$EXPSTRT)==FALSE&final_data$YR_TV>=final_data$EXPSTRT&final_data$YR_TV<=final_data$YR_LAST)|(is.na(final_data$EXPSTRT)==TRUE&final_data$CUMEXP!=0)),]
DT<-unique(final_data_1$AGE_TV[which(final_data_1$CASE==1)])

RMSE_3<-c(0,0,0,0)
Rsquared_3<-c(0,0,0,0)
AdjRsquared_3<-c(0,0,0,0)
for(i in 1:length(DT)){
  rowselect<-which(valid_data_3$AGE_TV==DT[i])
  if(length(rowselect)==0){
    cat("skip risk set ",i," by AGE_TV ", DT[i], " because sample size=0\n")
    next
  }
  else{
    cat("risk set", i, " by AGE_TV ", DT[i]," sample size=", length(rowselect),"\n")
  }
  valid_data_3.select<-valid_data_3[rowselect,]
  
  model.add.nocov<-lm(VALUE~EVALUE, data=valid_data_3.select)
  valid_data_3$VALUE.add.nocov[rowselect]<-fitted(model.add.nocov)
  model.add.cov<-lm(VALUE~EVALUE+SMOKE+NATIVE+HISPANIC+COMP5+COMP6+yy2, data=valid_data_3.select)
  valid_data_3$VALUE.add.cov[rowselect]<-fitted(model.add.cov)
  
  model.mul.nocov<-lm(logVALUE.pos~logEVALUE.pos,data=valid_data_3.select)
  valid_data_3$VALUE.mul.nocov[rowselect]<-exp(fitted(model.mul.nocov))-min(valid_data_3$VALUE[valid_data_3$VALUE>0])/2
  model.mul.cov<-lm(logVALUE.pos~logEVALUE.pos+SMOKE+NATIVE+HISPANIC+COMP5+COMP6+yy2, data=valid_data_3.select)
  valid_data_3$VALUE.mul.cov[rowselect]<-exp(fitted(model.mul.cov))-min(valid_data_3$VALUE[valid_data_3$VALUE>0])/2
  
  RMSE_3<-rbind(RMSE_3, c(sqrt(mean(model.add.nocov$residuals^2)),
                          sqrt(mean(model.add.cov$residuals^2)),
                          sqrt(mean(model.mul.nocov$residuals^2)),
                          sqrt(mean(model.mul.cov$residuals^2))))
  Rsquared_3<-rbind(Rsquared_3, c(summary(model.add.nocov)$r.squared,
                                  summary(model.add.cov)$r.squared,
                                  summary(model.mul.nocov)$r.squared,
                                  summary(model.mul.cov)$r.squared))
  AdjRsquared_3<-rbind(Rsquared_3, c(summary(model.add.nocov)$adj.r.squared,
                                     summary(model.add.cov)$adj.r.squared,
                                     summary(model.mul.nocov)$adj.r.squared,
                                     summary(model.mul.cov)$adj.r.squared))
}
# risk set 1  by AGE_TV  58  sample size= 12 
# skip risk set  2  by AGE_TV  71  because sample size=0
# risk set 3  by AGE_TV  57  sample size= 16 
# risk set 4  by AGE_TV  63  sample size= 7 
# skip risk set  5  by AGE_TV  73  because sample size=0
# skip risk set  6  by AGE_TV  77  because sample size=0
# skip risk set  7  by AGE_TV  70  because sample size=0
# skip risk set  8  by AGE_TV  75  because sample size=0
# risk set 9  by AGE_TV  66  sample size= 1 
# skip risk set  10  by AGE_TV  67  because sample size=0
# risk set 11  by AGE_TV  62  sample size= 5 
# skip risk set  12  by AGE_TV  83  because sample size=0
# risk set 13  by AGE_TV  56  sample size= 22 
# skip risk set  14  by AGE_TV  68  because sample size=0
# risk set 15  by AGE_TV  65  sample size= 4 
# skip risk set  16  by AGE_TV  76  because sample size=0
# risk set 17  by AGE_TV  55  sample size= 23 
# risk set 18  by AGE_TV  61  sample size= 7 
# skip risk set  19  by AGE_TV  69  because sample size=0
# skip risk set  20  by AGE_TV  72  because sample size=0
# skip risk set  21  by AGE_TV  84  because sample size=0
# skip risk set  22  by AGE_TV  80  because sample size=0
# risk set 23  by AGE_TV  41  sample size= 77 
# risk set 24  by AGE_TV  50  sample size= 38 
# risk set 25  by AGE_TV  43  sample size= 77 
# risk set 26  by AGE_TV  59  sample size= 9 
# risk set 27  by AGE_TV  54  sample size= 26 
# risk set 28  by AGE_TV  64  sample size= 4 
# skip risk set  29  by AGE_TV  78  because sample size=0
# skip risk set  30  by AGE_TV  79  because sample size=0
# risk set 31  by AGE_TV  53  sample size= 28 
# skip risk set  32  by AGE_TV  86  because sample size=0
# risk set 33  by AGE_TV  46  sample size= 55 
# risk set 34  by AGE_TV  60  sample size= 10 
# skip risk set  35  by AGE_TV  74  because sample size=0
# risk set 36  by AGE_TV  47  sample size= 48 
# skip risk set  37  by AGE_TV  88  because sample size=0
# risk set 38  by AGE_TV  44  sample size= 67 
# risk set 39  by AGE_TV  48  sample size= 50

valid_data_3.fitted<-data.frame(EVALUE=valid_data_3$EVALUE,
                                VALUE=valid_data_3$VALUE,
                                VALUE.add.nocov=valid_data_3$VALUE.add.nocov,
                                VALUE.add.cov=valid_data_3$VALUE.add.cov,
                                VALUE.mul.nocov=valid_data_3$VALUE.mul.nocov,
                                VALUE.mul.cov=valid_data_3$VALUE.mul.cov)


library(ggplot2)
plot.fittedVSevalue_3<-valid_data_3.fitted %>% tidyr::gather("id", "fitted", 2:6) %>% 
  ggplot(., aes(EVALUE, fitted))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE, color="blue")+
  ggtitle("True and Fitted Value in Risk Set Model VS Error-prone Value")+
  facet_wrap(~id)

plot.fittedVSevalue_3

plot.fittedVStrue_3<-valid_data_3.fitted %>% tidyr::gather("id", "fitted", 3:6) %>% 
  ggplot(., aes(VALUE, fitted))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE, color="blue")+
  ggtitle("Fitted Value in Risk Set Model VS True Value")+
  facet_wrap(~id)

plot.fittedVStrue_3

RMSE_3<-RMSE_3[-c(1),]
Rsquared_3<-Rsquared_3[-c(1),]
AdjRsquared_3<-AdjRsquared_3[-c(1),]
colnames(RMSE_3)<-c("add.nocov","add.cov","mul.nocov","mul.cov")
colnames(Rsquared_3)<-c("add.nocov","add.cov","mul.nocov","mul.cov")
colnames(AdjRsquared_3)<-c("add.nocov","add.cov","mul.nocov","mul.cov")
colMeans(Rsquared_3)
#add.nocov   add.cov mul.nocov   mul.cov 
#0.2780830 0.5957379 0.2193154 0.5047457
colMeans(AdjRsquared_3)
#add.nocov   add.cov mul.nocov   mul.cov 
#0.2699346 0.5731589 0.2100205 0.4818471

RMSEtotal_3<-c(0)
for(j in 3:6){
  RMSEtotal_3<-c(RMSEtotal_3,sqrt(mean((valid_data_3.fitted[,j]-valid_data_3.fitted$VALUE)^2)))
}
RMSEtotal_3<-RMSEtotal_3[-c(1)]
RMSEtotal_3
#add.nocov    add.cov  mul.nocov    mul.cov
#2.695278 2.657610 2.751629 2.689666


