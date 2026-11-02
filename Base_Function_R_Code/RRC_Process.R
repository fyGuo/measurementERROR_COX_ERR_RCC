#### Load Uminer Data ####

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

valid_data<-read_sas("~/OneDrive - Harvard University/JUL2015-FASRC/data/valid_2009_short.sas7bdat")
valid_data$VALUE100<-valid_data$VALUE/100
valid_data$EVALUE100<-valid_data$EVALUE/100
valid_data$yy1<-ifelse(valid_data$YEAR<70,1,0)
valid_data$yy2<-ifelse(valid_data$YEAR>=70&valid_data$YEAR<80,1,0)
valid_data$yy3<-ifelse(valid_data$YEAR>=80&valid_data$YEAR<90,1,0)
valid_data$yy4<-ifelse(valid_data$YEAR>=90&valid_data$YEAR<100,1,0)
valid_data$yy5<-ifelse(valid_data$YEAR>=100,1,0)


#### Create Risk Sets by Unique Death AGE_TV ####

dtime_year<-sort(unique(final_data_1$AGE_TV[which(final_data_1$CASE==1)]))
# 39
# minimum death time: 41, maximum death time: 88
dtime_month<-sort(unique(final_data_1$AGE_OUT[which(final_data_1$CASE==1)]))
# 113

mv_exp<-read_sas("~/OneDrive - Harvard University/JUL2015-FASRC/analysis/mv_exp.sas7bdat")
mv_exp$YREXP100<-mv_exp$YREXP/100
mv_exp_temp<-data.frame(EVALUE100=mv_exp$YREXP100,SMOKE=mv_exp$SMOKE,NATIVE=mv_exp$native,HISPANIC=mv_exp$hispanic,COMP5=mv_exp$COMP5,COMP6=mv_exp$COMP6,yy2=mv_exp$yy2,yy3=mv_exp$yy3,yy4=mv_exp$yy4,yy5=mv_exp$yy5)

mv_exp$YREXP100_em<-NA
mv_exp$YREXP100_em[which(mv_exp$surf!=1)]<-mv_exp$YREXP100[which(mv_exp$surf!=1)]

summary(valid_data_1$AGE_TV)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#41.00   43.00   47.00   47.87   51.00   66.00 
# age distribution in rows which are used to build the prediction model
summary(final_data_1$AGE_TV[which(final_data_1$surf==1)])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#11.00   27.00   33.00   33.57   39.00   70.00 
# age distribution in rows which need to be predicted

# no rows in valid_data with AGE_TV>88
valid_data_1<-valid_data[which(valid_data$AGE_TV>=41&valid_data$AGE_TV<=88),]
RS_v<-matrix(data=NA,nrow=nrow(valid_data_1),ncol=length(dtime_year))
# 839 rows
RS_m<-matrix(data=NA,nrow=nrow(final_data_1),ncol=length(dtime_year))
RS_predict_em<-matrix(data=NA,nrow=nrow(final_data_1),ncol=length(dtime_year))
for(k in 1:length(dtime_year)){
  if(k==length(dtime_year)){
    RS_v[,k]<-as.numeric(valid_data_1$AGE_TV>=dtime_year[k])
    RS_m[,k]<-as.numeric(final_data_1$AGE_TV>=dtime_year[k])
    #em<-lm(VALUE100~EVALUE100+SMOKE+NATIVE+HISPANIC+COMP5+COMP6+yy2+yy3+yy4+yy5, data=valid_data_1[which(RS_v[,k]==1),])
    #RS_predict_em[,k]<-predict(em,mv_exp_temp)
    #mv_exp$YREXP100_em[which(mv_exp$surf==1&RS_m[,k]==1)]<-RS_predict_em[which(mv_exp$surf==1&RS_m[,k]==1),k]
  }
  else{
    RS_v[,k]<-as.numeric(valid_data_1$AGE_TV>=dtime_year[k]&valid_data_1$AGE_TV<dtime_year[k+1])
    RS_m[,k]<-as.numeric(final_data_1$AGE_TV>=dtime_year[k]&final_data_1$AGE_TV<dtime_year[k+1])
    #em<-lm(VALUE100~EVALUE100+SMOKE+NATIVE+HISPANIC+COMP5+COMP6+yy2+yy3+yy4+yy5, data=valid_data_1[which(RS_v[,k]==1),])
    #RS_predict_em[,k]<-predict(em,mv_exp_temp)
    #mv_exp$YREXP100_em[which(mv_exp$surf==1&RS_m[,k]==1)]<-RS_predict_em[which(mv_exp$surf==1&RS_m[,k]==1),k]
  }
}

# colSums(RS_v)
#[1] 156  77 125  55  48  98 106  28  26  23  22  16  12   9  10   7   5   7   4   4   1   0   0   0   0   0
#[27]   0   0   0   0   0   0   0   0   0   0   0   0   0
# 0 entries in RS_v[,k] when k>=22, 1 entry in k=21
# combine risk set from k>=20, 5 entries in total
# colSums(RS_m)
# 30542 entries with AGE_TV<41, set RS_m[,1]=1


#### Adjust for Error in individual Risk Sets ####

dtime_year<-sort(unique(final_data_1$AGE_TV[which(final_data_1$CASE==1)]))[1:20]
# 20
RS_v<-matrix(data=NA,nrow=nrow(valid_data_1),ncol=length(dtime_year))
RS_m<-matrix(data=NA,nrow=nrow(final_data_1),ncol=length(dtime_year))
RS_predict_em<-matrix(data=NA,nrow=nrow(final_data_1),ncol=length(dtime_year))
for(k in 1:length(dtime_year)){
  if(k==length(dtime_year)){
    RS_v[,k]<-as.numeric(valid_data_1$AGE_TV>=dtime_year[k])
    RS_m[,k]<-as.numeric(final_data_1$AGE_TV>=dtime_year[k])
    em<-lm(VALUE100~EVALUE100+SMOKE+NATIVE+HISPANIC+COMP5+COMP6+yy2+yy3+yy4+yy5, data=valid_data_1[which(RS_v[,k]==1),])
    RS_predict_em[,k]<-predict(em,mv_exp_temp)
    mv_exp$YREXP100_em[which(mv_exp$surf==1&RS_m[,k]==1)]<-RS_predict_em[which(mv_exp$surf==1&RS_m[,k]==1),k]
  }
  else if(k==1){
    RS_v[,k]<-as.numeric(valid_data_1$AGE_TV>=dtime_year[k]&valid_data_1$AGE_TV<dtime_year[k+1])
    RS_m[,k]<-as.numeric(final_data_1$AGE_TV<dtime_year[k+1])
    em<-lm(VALUE100~EVALUE100+SMOKE+NATIVE+HISPANIC+COMP5+COMP6+yy2+yy3+yy4+yy5, data=valid_data_1[which(RS_v[,k]==1),])
    RS_predict_em[,k]<-predict(em,mv_exp_temp)
    mv_exp$YREXP100_em[which(mv_exp$surf==1&RS_m[,k]==1)]<-RS_predict_em[which(mv_exp$surf==1&RS_m[,k]==1),k]
  }
  else{
    RS_v[,k]<-as.numeric(valid_data_1$AGE_TV>=dtime_year[k]&valid_data_1$AGE_TV<dtime_year[k+1])
    RS_m[,k]<-as.numeric(final_data_1$AGE_TV>=dtime_year[k]&final_data_1$AGE_TV<dtime_year[k+1])
    em<-lm(VALUE100~EVALUE100+SMOKE+NATIVE+HISPANIC+COMP5+COMP6+yy2+yy3+yy4+yy5, data=valid_data_1[which(RS_v[,k]==1),])
    RS_predict_em[,k]<-predict(em,mv_exp_temp)
    mv_exp$YREXP100_em[which(mv_exp$surf==1&RS_m[,k]==1)]<-RS_predict_em[which(mv_exp$surf==1&RS_m[,k]==1),k]
  }
}

mv_exp$yrpredexpxm<-ave(mv_exp$predexpxm, factor(mv_exp$id), FUN=function(x) c(NA,diff(x)))
mv_exp$yrpredexpxm[which(is.na(mv_exp$yrpredexpxm)==TRUE)]<-mv_exp$predexpxm[which(is.na(mv_exp$yrpredexpxm)==TRUE)]
summary(mv_exp$yrpredexpxm)

mv_exp$CUMEXP100_em<-ave(mv_exp$YREXP100_em,factor(mv_exp$id),FUN=cumsum)
mv_exp$age_out<-mv_exp$age_tv+mv_exp$time


#### Fit Cox/ERR after RRC Adjustment ####

library(survival)
model.coxph.rrc.adjusted.b<-coxph(Surv(time,CASE)~(CUMEXP100_em+SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5)+strata(age_tv),data=mv_exp,ties="breslow")

model.coxph.adjusted.b<-coxph(Surv(time,CASE)~(cumexp100+SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5)+strata(age_tv),data=mv_exp,ties="breslow")
round(model.coxph.rrc.adjusted.b$coefficients,digits=3)
#CUMEXP100_em        SMOKE       native     hispanic        COMP5        COMP6          yy2          yy3 
#1.543        1.620       -0.822       -1.053        0.409        0.673        0.072        0.088 
#yy4          yy5 
#-0.247       -0.683 
round(model.coxph.adjusted.b$coefficients,digits=3)
#cumexp100     SMOKE    native  hispanic     COMP5     COMP6       yy2       yy3       yy4       yy5 
#0.297     1.671    -0.861    -0.953     0.302     0.485     0.140     0.105    -0.265    -0.709


model.err.rrc.adjusted<-ERR(data=mv_exp,outcome=c("CASE"),
                            covariates=c("SMOKE","native","hispanic","COMP5","COMP6","yy2","yy3","yy4","yy5"),
                            exposure=c("CUMEXP100_em"),times=c("age_tv","age_out"),strata=c("age_tv"),
                            start=c(as.numeric(model.coxph.adjusted.b$coefficients[2:10]),1))

model.err.adjusted<-ERR(data=mv_exp,outcome=c("CASE"),
                        covariates=c("SMOKE","native","hispanic","COMP5","COMP6","yy2","yy3","yy4","yy5"),
                        exposure=c("cumexp100"),times=c("age_tv","age_out"),strata=c("age_tv"),
                        start=c(as.numeric(model.coxph.adjusted.b$coefficients[2:10]),1))

round(model.err.adjusted$estimate,digits=3)

round(model.err.rrc.adjusted$estimate,digits=3)
#1.636 -0.853 -1.043  0.501  0.643  0.123  0.094 -0.250 -0.700  3.016

