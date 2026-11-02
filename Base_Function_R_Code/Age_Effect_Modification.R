## Effect modification

setwd("~/OneDrive - Harvard University/JUL2015-FASRC/data")

library(haven)
final_data<-read_sas("final_data.sas7bdat")
# 88843 rows, 40 columns

summary(final_data)
# 438 NAs for EXPSTRT, EXPEND, AGESTRT
# 1 NA for CEXP
# 42783 NAs for LASTE

final_data$YR_LAST<-as.numeric(final_data$year)-1900

final_data_1<-final_data[which((is.na(final_data$EXPSTRT)==FALSE&final_data$YR_TV>=final_data$EXPSTRT&final_data$YR_TV<=final_data$YR_LAST)|(is.na(final_data$EXPSTRT)==TRUE&final_data$CUMEXP!=0)),]
# 88405 rows, 41 columns
# manual double check
# View(final_data[which(is.na(final_data$EXPSTRT)==TRUE),])

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
# final_data_1: 88405*51

final_data_1$AGE_OUT<-final_data_1$AGE_TV+final_data_1$time

library(survival)
# these two COXPH models output match with draft results
model.coxph.naive<-coxph(Surv(time,CASE)~(cumexp100)+strata(AGE_TV),data=final_data_1)
summary(model.coxph.naive)
model.coxph.adjusted<-coxph(Surv(time,CASE)~(cumexp100+yy2+yy3+yy4+yy5+SMOKE+native+hispanic+COMP5+COMP6)+strata(AGE_TV),data=final_data_1)
summary(model.coxph.adjusted)

model.coxph.modified<-coxph(Surv(time,CASE)~(cumexp100*AGE_TV+AGE_TV+cumexp100+yy2+yy3+yy4+yy5+SMOKE+native+hispanic+COMP5+COMP6)+strata(AGE_TV),data=final_data_1)
summary(model.coxph.modified)
coef(model.coxph.modified)

anova(model.coxph.modified,model.coxph.adjusted)

x2<-seq(19,106,by=1)
plot(x2,rep(0,length(x2)),type='n',
     xlab="Attained Age",ylab="RR-scale Effect Coefficient",main="Effect Modification of Age in UMiner Data",
     xlim=c(19,106),ylim=c(0,3))
lines(x2,exp(-0.005015325*x2+0.612936726),type='l',col='blue')
lines(x2,(-0.006820258*x2+1.176266),type='l',col='red')
legend("topright",
       legend=c("Cox","ERR"),
       col=c("blue","red"),lty=1,cex=0.7)


