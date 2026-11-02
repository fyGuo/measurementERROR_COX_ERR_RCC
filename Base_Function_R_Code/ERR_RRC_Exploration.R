#### Load and Explore Uminer Data with Cox/ERR and RRC ####


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

# data check
summary(final_data_1$YR_TV[which(final_data_1$CASE==1)])
# all cases happen at or after 1967, when we have true yearly exposure
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#67.00   85.00   95.00   93.11  102.00  111.00
which(final_data_1$YR_TV==67&final_data_1$CASE==1)
#[1] 26851
View(final_data_1[26851,])

summary(final_data_1$YR_TV[which(final_data_1$surf==1)])
# all exposure values with error happen at or before 1967, when we have error-prone yearly exposure
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#57.00   61.00   63.00   62.83   65.00   67.00 
sum(final_data_1$surf)
#[1] 5242
length(which(final_data_1$YR_TV>=57&final_data_1$YR_TV<=67))
#[1] 7693
# not every exposure between 1957 and 1967 is error-prone, surf==1 is the correct way to identify it
# Conclusion: cumulative exposure at death time = cumulative corrected exposure from entry to 1967 (surf==1) + cumulative true exposure from 1967 to current year (surf==0)


valid_data<-read_sas("~/OneDrive - Harvard University/JUL2015-FASRC/data/valid_2009_short.sas7bdat")
#Two risk sets with risk set width s>=180
#300 and 481 person year in each, 126 and 156 unique ids in each
#length(which(valid_data$AGE_TV>=41&valid_data$AGE_TV<=44))
#[1] 300
#length(which(valid_data$AGE_TV>=46&valid_data$AGE_TV<=88))
#[1] 481
#length(unique(valid_data$ID[which(valid_data$AGE_TV>=41&valid_data$AGE_TV<=44)]))
#[1] 126
#length(unique(valid_data$ID[which(valid_data$AGE_TV>=46&valid_data$AGE_TV<=88)]))
#[1] 146

valid_data$VALUE100<-valid_data$VALUE/100
valid_data$EVALUE100<-valid_data$EVALUE/100
valid_data$yy1<-ifelse(valid_data$YEAR<70,1,0)
valid_data$yy2<-ifelse(valid_data$YEAR>=70&valid_data$YEAR<80,1,0)
valid_data$yy3<-ifelse(valid_data$YEAR>=80&valid_data$YEAR<90,1,0)
valid_data$yy4<-ifelse(valid_data$YEAR>=90&valid_data$YEAR<100,1,0)
valid_data$yy5<-ifelse(valid_data$YEAR>=100,1,0)
#summary(valid_data$YEAR)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#67.00   69.00   70.00   70.52   72.00   76.00
# only yy1 and yy2 will have non-zero values

RS1<-valid_data[which(valid_data$AGE_TV>=41&valid_data$AGE_TV<=44),]
RS2<-valid_data[which(valid_data$AGE_TV>=46&valid_data$AGE_TV<=88),]
em1<-lm(VALUE100~EVALUE100+SMOKE+NATIVE+HISPANIC+COMP5+COMP6, data=RS1)
em2<-lm(VALUE100~EVALUE100+SMOKE+NATIVE+HISPANIC+COMP5+COMP6, data=RS2)
#mean(em1$residuals^2)
#[1] 0.0002968713
#mean(em2$residuals^2)
#[1] 0.0002733588
em1<-lm(VALUE~EVALUE+SMOKE+NATIVE+HISPANIC+COMP5+COMP6+yy2, data=RS1)
em2<-lm(VALUE~EVALUE+SMOKE+NATIVE+HISPANIC+COMP5+COMP6+yy2, data=RS2)
#yy3,yy4,yy5 will have no coefficients due to sigularity
#adding yy2 will improve mean squared error
#mean(em1$residuals^2)
#[1] 0.0002816023
#mean(em2$residuals^2)
#[1] 0.0002624194

RS1$logVALUE<-log(RS1$VALUE+1)
RS1$logEVALUE<-log(RS1$EVALUE+1)
RS2$logVALUE<-log(RS2$VALUE+1)
RS2$logEVALUE<-log(RS2$EVALUE+1)

#log(Value+1)=log(EVALUE+1)+Z+epsilon
#Value=[(Evalue+1)*exp(Z+epsilon)]*exp(sigma_epsilon^2/2)-1

em1.log<-lm(logVALUE~logEVALUE+SMOKE+NATIVE+HISPANIC+COMP5+COMP6+yy2, data=RS1)
em2.log<-lm(logVALUE~logEVALUE+SMOKE+NATIVE+HISPANIC+COMP5+COMP6+yy2, data=RS2)
mv_exp$logYREXP<-log(mv_exp$YREXP+1)
mv_exp_temp.log<-data.frame(logEVALUE=mv_exp$logYREXP,SMOKE=mv_exp$SMOKE,NATIVE=mv_exp$native,HISPANIC=mv_exp$hispanic,COMP5=mv_exp$COMP5,COMP6=mv_exp$COMP6,yy2=mv_exp$yy2,yy3=mv_exp$yy3,yy4=mv_exp$yy4,yy5=mv_exp$yy5)

mv_exp$YREXP100_em1.log1<-exp(predict(em1.log,mv_exp_temp.log))/100
mv_exp$YREXP100_em1.log2<-(exp(predict(em1.log,mv_exp_temp.log))*exp(summary(em1.log)$sigma^2/2)-1)/100
mv_exp$YREXP100_em2.log1<-exp(predict(em2.log,mv_exp_temp.log))/100
mv_exp$YREXP100_em2.log2<-(exp(predict(em2.log,mv_exp_temp.log))*exp(summary(em2.log)$sigma^2/2)-1)/100

sum(mv_exp$YREXP100_em1.log1[which(mv_exp$surf==1&mv_exp$age_tv<=45)])
sum(mv_exp$YREXP100_em1.log2[which(mv_exp$surf==1&mv_exp$age_tv<=45)])
sum(mv_exp$yrpredexpxm[which(mv_exp$surf==1&mv_exp$age_tv<=45)])
sum(mv_exp$YREXP100_em1[which(mv_exp$surf==1&mv_exp$age_tv<=45)])

sum(mv_exp$YREXP100_em2.log1[which(mv_exp$surf==1&mv_exp$age_tv>45)])
sum(mv_exp$YREXP100_em2.log2[which(mv_exp$surf==1&mv_exp$age_tv>45)])
sum(mv_exp$yrpredexpxm[which(mv_exp$surf==1&mv_exp$age_tv>45)])
sum(mv_exp$YREXP100_em2[which(mv_exp$surf==1&mv_exp$age_tv>45)])


mv_exp$YREXP100_em.log<-NA
mv_exp$YREXP100_em.log[which(mv_exp$surf!=1)]<-mv_exp$YREXP100[which(mv_exp$surf!=1)]
mv_exp$YREXP100_em.log[which(mv_exp$surf==1&mv_exp$age_tv<=45)]<-mv_exp$YREXP100_em1.log1[which(mv_exp$surf==1&mv_exp$age_tv<=45)]
mv_exp$YREXP100_em.log[which(mv_exp$surf==1&mv_exp$age_tv>45)]<-mv_exp$YREXP100_em2.log1[which(mv_exp$surf==1&mv_exp$age_tv>45)]

mv_exp$CUMEXP100_em.log<-ave(mv_exp$YREXP100_em.log,factor(mv_exp$id),FUN=cumsum)
mv_exp$CUMEXP100_em.log<-abs(mv_exp$CUMEXP100_em.log)

model.coxph.rrc.log.adjusted.b<-coxph(Surv(time,CASE)~(CUMEXP100_em.log+SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5)+strata(age_tv),data=mv_exp,ties="breslow")
model.coxph.rrc.adjusted.b<-coxph(Surv(time,CASE)~(CUMEXP100_em+SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5)+strata(age_tv),data=mv_exp,ties="breslow")
model.coxph.adjusted.b<-coxph(Surv(time,CASE)~(cumexp100+SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5)+strata(age_tv),data=mv_exp,ties="breslow")
round(exp(as.numeric(model.coxph.adjusted.b$coefficients)),digits=3)
#[1] 1.346 5.317 0.423 0.385 1.352 1.624 1.150 1.110 0.767 0.492
round(exp(as.numeric(model.coxph.rrc.adjusted.b$coefficients)),digits=3)
#[1] 5.347 5.218 0.424 0.346 1.692 2.112 1.127 1.190 0.874 0.563
round(exp(as.numeric(model.coxph.rrc.log.adjusted.b$coefficients)),digits=3)
#[1] 8.642 5.002 0.420 0.329 1.956 2.319 1.038 0.978 0.686 0.427

model.coxph.rrc.log.naive.b<-coxph(Surv(time,CASE)~(CUMEXP100_em.log)+strata(age_tv),data=mv_exp,ties="breslow")
model.coxph.rrc.naive.b<-coxph(Surv(time,CASE)~(CUMEXP100_em)+strata(age_tv),data=mv_exp,ties="breslow")
model.coxph.naive.b<-coxph(Surv(time,CASE)~(cumexp100)+strata(age_tv),data=mv_exp,ties="breslow")
round(exp(as.numeric(model.coxph.naive.b$coefficients)),digits=3)
#[1] 1.461
round(exp(as.numeric(model.coxph.rrc.naive.b$coefficients)),digits=3)
#[1] 7.737
round(exp(as.numeric(model.coxph.rrc.log.naive.b$coefficients)),digits=3)
#[1] 10.856


model.err.adjusted<-ERR_breslow(data=mv_exp,outcome=c("CASE"),
                                covariates=c("SMOKE","native","hispanic","COMP5","COMP6","yy2","yy3","yy4","yy5"),
                                exposure=c("cumexp100"),times=c("age_tv","age_out"),strata=c("age_tv"),
                                start=c(as.numeric(model.coxph.rrc.adjusted.b$coefficients[2:10]),3))
model.err.rrc.adjusted<-ERR_breslow(data=mv_exp,outcome=c("CASE"),
                                    covariates=c("SMOKE","native","hispanic","COMP5","COMP6","yy2","yy3","yy4","yy5"),
                                    exposure=c("CUMEXP100_em"),times=c("age_tv","age_out"),strata=c("age_tv"),
                                    start=c(as.numeric(model.coxph.rrc.adjusted.b$coefficients[2:10]),3))
model.err.rrc.log.adjusted<-ERR_breslow(data=mv_exp,outcome=c("CASE"),
                                    covariates=c("SMOKE","native","hispanic","COMP5","COMP6","yy2","yy3","yy4","yy5"),
                                    exposure=c("CUMEXP100_em.log"),times=c("age_tv","age_out"),strata=c("age_tv"),
                                    start=c(as.numeric(model.coxph.rrc.adjusted.b$coefficients[2:10]),3))
round(model.err.adjusted$estimate, 3)
#[1]  1.673 -0.808 -0.945  0.258  0.399  0.184  0.169 -0.164 -0.578  0.750
round(model.err.rrc.adjusted$estimate, 3)
#[1]  1.655 -0.892 -1.050  0.577  0.692  0.164  0.157 -0.177 -0.635  3.105
round(model.err.rrc.log.adjusted$estimate, 3)
#[1]  1.633 -0.926 -1.073  0.665  0.731  0.101  0.029 -0.340 -0.830  3.421

model.err.unadjusted<-ERR_breslow(data=mv_exp,outcome=c("CASE"),
                                      covariates=NULL,
                                      exposure=c("cumexp100"),times=c("age_tv","age_out"),strata=c("age_tv"),
                                      start=c(3))
model.err.rrc.unadjusted<-ERR_breslow(data=mv_exp,outcome=c("CASE"),
                                      covariates=NULL,
                                      exposure=c("CUMEXP100_em"),times=c("age_tv","age_out"),strata=c("age_tv"),
                                      start=c(3))
model.err.rrc.log.unadjusted<-ERR_breslow(data=mv_exp,outcome=c("CASE"),
                                      covariates=NULL,
                                      exposure=c("CUMEXP100_em.log"),times=c("age_tv","age_out"),strata=c("age_tv"),
                                      start=c(3))
round(model.err.unadjusted$estimate, 3)
#1.609
round(model.err.rrc.unadjusted$estimate, 3)
#[1] 5.522
round(model.err.rrc.log.unadjusted$estimate, 3)
#[1] 5.459


#### Explore: Fit GEE Error Model ####

library(geepack)
library(gee)
mv_exp$wave<-ave(mv_exp$YR_TV, factor(mv_exp$id), FUN=function(x) c(NA,diff(x)))
mv_exp$wave[which(is.na(mv_exp$wave)==TRUE)]<-mv_exp$predexpxm[which(is.na(mv_exp$yrpredexpxm)==TRUE)]

gee1<-geeglm(formula=VALUE100~EVALUE100+SMOKE+NATIVE+HISPANIC+COMP5+COMP6+yy2,
          family=gaussian,data=RS1,id=ID,waves=YEAR,corstr="independence")

gee1.ar1<-geeglm(formula=VALUE100~EVALUE100+SMOKE+NATIVE+HISPANIC+COMP5+COMP6+yy2,
             family=gaussian,data=RS1,id=ID,waves=YEAR,corstr="ar1")


mv_exp<-read_sas("~/OneDrive - Harvard University/JUL2015-FASRC/analysis/mv_exp.sas7bdat")
mv_exp$YREXP100<-mv_exp$YREXP/100
mv_exp_temp<-data.frame(EVALUE100=mv_exp$YREXP100,SMOKE=mv_exp$SMOKE,NATIVE=mv_exp$native,HISPANIC=mv_exp$hispanic,COMP5=mv_exp$COMP5,COMP6=mv_exp$COMP6,yy2=mv_exp$yy2,yy3=mv_exp$yy3,yy4=mv_exp$yy4,yy5=mv_exp$yy5)
mv_exp$YREXP100_em1<-predict(em1,mv_exp_temp)
mv_exp$YREXP100_em2<-predict(em2,mv_exp_temp)

mv_exp$yrpredexpxm<-ave(mv_exp$predexpxm, factor(mv_exp$id), FUN=function(x) c(NA,diff(x)))
mv_exp$yrpredexpxm[which(is.na(mv_exp$yrpredexpxm)==TRUE)]<-mv_exp$predexpxm[which(is.na(mv_exp$yrpredexpxm)==TRUE)]

# Use prediction value only for 1163 subjects, 5242 rows, between 1957-1967, indexed by surf==1
View(mv_exp[which(mv_exp$surf==1),c("id","age_tv","yrpredexpxm","YREXP100","YREXP100_em1","YREXP100_em2")])
sum(mv_exp$yrpredexpxm[which(mv_exp$surf==1&mv_exp$age_tv<=45)])
sum(mv_exp$YREXP100_em1[which(mv_exp$surf==1&mv_exp$age_tv<=45)])

View(mv_exp[which(mv_exp$surf==1&mv_exp$age_tv>45),c("id","age_tv","yrpredexpxm","YREXP100","YREXP100_em1","YREXP100_em2")])
sum(mv_exp$yrpredexpxm[which(mv_exp$surf==1&mv_exp$age_tv>45)])
sum(mv_exp$YREXP100_em2[which(mv_exp$surf==1&mv_exp$age_tv>45)])

mv_exp$YREXP100_em<-NA
mv_exp$YREXP100_em[which(mv_exp$surf!=1)]<-mv_exp$YREXP100[which(mv_exp$surf!=1)]
mv_exp$YREXP100_em[which(mv_exp$surf==1&mv_exp$age_tv<=45)]<-mv_exp$YREXP100_em1[which(mv_exp$surf==1&mv_exp$age_tv<=45)]
mv_exp$YREXP100_em[which(mv_exp$surf==1&mv_exp$age_tv>45)]<-mv_exp$YREXP100_em2[which(mv_exp$surf==1&mv_exp$age_tv>45)]

mv_exp$CUMEXP100_em<-ave(mv_exp$YREXP100_em,factor(mv_exp$id),FUN=cumsum)
mv_exp$CUMEXP100_em<-abs(mv_exp$CUMEXP100_em)
mv_exp$age_out<-mv_exp$age_tv+mv_exp$time

library(survival)
model.coxph.rrc.adjusted.b<-coxph(Surv(time,CASE)~(CUMEXP100_em+SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5)+strata(age_tv),data=mv_exp,ties="breslow")
model.coxph.adjusted.b<-coxph(Surv(time,CASE)~(cumexp100+SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5)+strata(age_tv),data=mv_exp,ties="breslow")
round(as.numeric(model.coxph.adjusted.b$coefficients),digits=3)
round(as.numeric(model.coxph.rrc.adjusted.b$coefficients),digits=3)

model.coxph.rrc.naive.b<-coxph(Surv(time,CASE)~(CUMEXP100_em)+strata(age_tv),data=mv_exp,ties="breslow")
model.coxph.naive.b<-coxph(Surv(time,CASE)~(cumexp100)+strata(age_tv),data=mv_exp,ties="breslow")
round(as.numeric(model.coxph.naive.b$coefficients),digits=3)
round(as.numeric(model.coxph.rrc.naive.b$coefficients),digits=3)

model.err.rrc.adjusted<-ERR_breslow(data=mv_exp,outcome=c("CASE"),
                            covariates=c("SMOKE","native","hispanic","COMP5","COMP6","yy2","yy3","yy4","yy5"),
                            exposure=c("CUMEXP100_em"),times=c("age_tv","age_out"),strata=c("age_tv"),
                            start=c(as.numeric(model.coxph.rrc.adjusted.b$coefficients[2:10]),3))
model.err.rrc.adjusted$estimate
model.err.rrc.unadjusted<-ERR_breslow(data=mv_exp,outcome=c("CASE"),
                            covariates=NULL,
                            exposure=c("CUMEXP100_em"),times=c("age_tv","age_out"),strata=c("age_tv"),
                            start=c(3))
model.err.rrc.unadjusted$estimate

model.coxph.adjusted.b<-coxph(Surv(time,CASE)~(cumexp100+SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5)+strata(AGE_TV),data=final_data_1,ties="breslow")
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


summary(model.err.rrc.adjusted)
model.err.rrc.adjusted$estimate
# when em1 and em2 are fit without yy2-yy5
#[1]  1.6564452 -0.8962619 -1.0488900  0.5642577  0.6873010  0.1756176  0.1809412 -0.1475585 -0.6007430
#[10]  3.0020554
# when em1 and em2 are fit with yy2-yy5
#[1]  1.6553962 -0.8920066 -1.0498169  0.5766766  0.6918478  0.1638341  0.1572038 -0.1769214 -0.6345456
#[10]  3.1057013

model.coxph.rrc.adjusted.b$coefficients
# when em1 and em2 are fit without yy2-yy5
#CUMEXP100_em        SMOKE       native     hispanic        COMP5        COMP6          yy2          yy3 
#1.55591239   1.65805750  -0.87135359  -1.05181524   0.50988980   0.72819328   0.13668896   0.20661693 
#yy4          yy5 
#-0.09893083  -0.53948784 
# when em1 and em2 are fit with yy2-yy5
#CUMEXP100_em        SMOKE       native     hispanic        COMP5        COMP6          yy2          yy3 
#1.6765346    1.6521924   -0.8571185   -1.0602918    0.5259610    0.7474950    0.1196825    0.1735995 
#yy4          yy5 
#-0.1342674   -0.5745794 


#[1] "id"           "age_tv"       "predexpxm"    "surtotxm"     "BIRTH_YY"     "LUNG_CA"      "can"         
#[8] "EXPSTRT"      "EXPEND"       "SMOKE"        "ETHNIC"       "age_new"      "year"         "mnth"        
#[15] "AGE_STRT"     "FYR"          "COMP5"        "COMP6"        "COMP12"       "_NAME_"       "CUMEXP"      
#[22] "YR_TV"        "YREXP"        "comp_m"       "yrn"          "AGE"          "CEXP"         "CASE"        
#[29] "y1"           "y2"           "y3"           "y4"           "y5"           "y6"           "y7"          
#[36] "LASTE"        "type"         "type1"        "type2"        "type3"        "time"         "etype"       
#[43] "yr_last0"     "yr_last"      "surf"         "white"        "hispanic"     "native"       "cumexp100"   
#[50] "cumexp500"    "yy1"          "yy2"          "yy3"          "yy4"          "yy5"          "casenew"     
#[57] "i"            "YREXP100"     "YREXP100_em1" "YREXP100_em2" "yrpredexpxm"  "YREXP100_em"  "CUMEXP100_em"
#[64] "age_out" 





#### Explore: Covariance Matrix after RRC Hand Calculation Turned Function ####

# Calculate covariance matrix for beta_RRC, incorporating covariance matrix for em1,em2 estimation

# main study sample size
n1<-length(unique(mv_exp$id))
# validation study sample size, consisting of unique subjects in RS1 combined with RS2
n2<-length(union(unique(RS1$ID),unique(RS2$ID)))
# number of covariates and exposures
p<-1
q<-9
# beta_RRC from ERR model performed on RRC corrected exposures
beta_rrc<-as.numeric(model.err.rrc.adjusted$estimate)
# alpha from measurement error model
# applies to age_tv<=45, and >45
alpha1<-as.numeric(em1$coefficients)
alpha2<-as.numeric(em2$coefficients)

# I_beta
case<-which(mv_exp$CASE==1)
# covariate matrix
Z<-as.matrix(mv_exp[,c("SMOKE","native","hispanic","COMP5","COMP6","yy2","yy3","yy4","yy5")])
Xe<-mv_exp$cumexp100
X<-mv_exp$CUMEXP100_em
# calculate risk function
risk<-exp(Z%*%beta_rrc[1:q])*(1+X*beta_rrc[(q+1):(q+p)])
adjX<-X/(1+X*beta_rrc[(q+1):(q+p)])
design<-cbind(Z,adjX)

for(i in 1:length(case)){
  case_i<-case[i]
  # create risk set indicator variable Ym = I(Tj<=Ti)
  Risk_set<-which(mv_exp$age_tv==mv_exp$age_tv[case_i]&mv_exp$age_out>=mv_exp$age_out[case_i])
  Ym<-as.numeric(mv_exp$id %in% Risk_set)
  # calculate I_beta term
  Ii<--(sum(risk[Risk_set])*t(design[Risk_set,]*risk[Risk_set])%*%design[Risk_set,]-as.matrix(colSums(design[Risk_set,]*risk[Risk_set]))%*%colSums(design[Risk_set,]*risk[Risk_set]))/(sum(risk[Risk_set])^2)
  Ii[p+q,p+q]<--adjX[case_i]^2+(sum(adjX[Risk_set]*risk[Risk_set])^2)/(sum(risk[Risk_set])^2)
  I<-I+Ii
  # calculate H_beta term
  # U_beta
  Ubeta<-design[case_i,]-as.matrix(colSums(design[Risk_set,]*risk[Risk_set]))/sum(risk[Risk_set])
  term<-0
  for(j in 1:length(case)){
    case_j<-case[j]
    Risk_set_j<-which(mv_exp$age_tv==mv_exp$age_tv[case_j]&mv_exp$age_out>=mv_exp$age_out[case_j])
    if(mv_exp$age_tv[case_i]==mv_exp$age_tv[case_j]&mv_exp$age_out[case_i]>=mv_exp$age_out[case_j]){
      term<-term+risk[case_j]/sum(risk[Risk_set_j])*(design[case_j,]-as.matrix(colSums(design[Risk_set_j,]*risk[Risk_set_j]))/sum(risk[Risk_set_j]))
    }
  }
  Ubeta_tilde<-Ubeta-term
  Hbeta<-t(Ubeta_tilde)%*%Ubeta_tilde
  # calculate H_beta adjustment term
  # V_alpha
  
  # U_star
}


# covariance matrix after RCC hand calculation

# basic parameters
# number of exposure variables
p<-1
# number of risk factor variables
q<-9
# main study sample size
n1<-length(unique(mv_exp$id))
# validation study sample size, consisting of unique subjects in RS1 combined with RS2
n2<-length(union(unique(RS1$ID),unique(RS2$ID)))

# estimated ERR model coefficients after RRC
beta_rrc<-as.numeric(model.err.rrc.adjusted$estimate)
# alpha from measurement error model
# applies to age_tv<=45, and >45
alpha1<-as.numeric(em1$coefficients)
alpha2<-as.numeric(em2$coefficients)
# due to singularity, set alpha1/2[9,10,11]=0
alpha1[c(9,10,11)]<-c(0,0,0)
alpha2[c(9,10,11)]<-c(0,0,0)

# covariates matrix
Z<-as.matrix(mv_exp[,c("SMOKE","native","hispanic","COMP5","COMP6","yy2","yy3","yy4","yy5")])
Xe<-mv_exp$cumexp100
X<-mv_exp$CUMEXP100_em
# calculate risk function
risk<-exp(Z%*%beta_rrc[1:q])*(1+X*beta_rrc[(q+1):(q+p)])
adjX<-X/(1+X*beta_rrc[(q+1):(q+p)])
design<-cbind(adjX,Z)
design_em<-cbind(rep(1,nrow(mv_exp)),mv_exp$YREXP100_em,Z)

# for 132 cases, calculate S0,S1, only dependent on death time
# row numbers for 132 case
case_index<-which(mv_exp$CASE==1)
S0<-c()
S1<-matrix(NA,nrow=p+q,ncol=length(case_index))
A<-matrix(NA,nrow=1+p+q,ncol=length(case_index))
for(i in 1:length(case_index)){
  # row number of case i
  case_i<-case_index[i]
  # create risk set of case i (row numbers)
  Risk_set_i<-which(mv_exp$age_tv==mv_exp$age_tv[case_i]&mv_exp$age_out>=mv_exp$age_out[case_i])
  # calculate S0 = sum of risk in Risk_set_i
  S0<-c(S0,sum(risk[Risk_set_i]))
  # calculate S1 = sum of (adjX,Z)^T * risk in Risk_set_i
  S1[,i]<-colSums(design[Risk_set_i,]*risk[Risk_set_i])
  # calculate A = sum of (1,YREXP100_em,Z) * surf from entry till 1967 for case i
  A[,i]<-colSums(design_em[which(mv_exp$id==mv_exp$id[case_i]),]*mv_exp$surf[which(mv_exp$id==mv_exp$id[case_i])])
}

I<-matrix(0,nrow=p+q,ncol=p+q)
Hbeta<-matrix(0,nrow=p+q,ncol=p+q)
for(i in 1:length(case_index)){
  # row number of case i
  case_i<-case_index[i]
  # create risk set of case i (row numbers)
  Risk_set_i<-which(mv_exp$age_tv==mv_exp$age_tv[case_i]&mv_exp$age_out>=mv_exp$age_out[case_i])
  
  # Calculate I_beta_inverse --> Calculate U_dbeta
  Ii<--(S0[i]*t(design[Risk_set_i,]*risk[Risk_set_i])%*%design[Risk_set_i,]-as.matrix(colSums(design[Risk_set_i,]*risk[Risk_set_i]))%*%colSums(design[Risk_set_i,]*risk[Risk_set_i]))/(S0[i]^2)
  Ii[1,1]<--adjX[case_i]^2+(sum(adjX[Risk_set_i]*risk[Risk_set_i])^2)/(S0[i]^2)
  I<-I+Ii
  
  # Calculate H_beta_alpha
  # 1. Calculate H_beta
  # (1) Calculate U_beta
  Ubetai<-design[case_i,]-S1[,i]/S0[i]
  # (2) Calculate U_beta_adj
  # find case index set j, such that Ti>=Tj, i is in Risk_set_j
  # Do we still need to stratify by age_tv? Now, yes
  Ubetai_adj<-0
  for(j in 1:length(case_index)){
    case_j<-case_index[j]
    Risk_set_j<-which(mv_exp$age_tv==mv_exp$age_tv[case_j]&mv_exp$age_out>=mv_exp$age_out[case_j])
    if(mv_exp$age_tv[case_i]==mv_exp$age_tv[case_j]&mv_exp$age_out[case_i]>=mv_exp$age_out[case_j]){
      Ubetai_adj<-Ubetai_adj+risk[case_i]/(n1*S0[j])*(design[case_j,]-S1[,j]/S0[j])
    }
  }
  # (3) Calculate U_beta_tilde=U_beta-U_beta_adj, H_beta=U_beta_tilde%*%t(U_beta_tilde)
  Ubetai_tilde<-Ubetai-Ubetai_adj
  Hbetai<-as.matrix(Ubetai_tilde)%*%t(as.matrix(Ubetai_tilde))
  Hbeta<-Hbeta+Hbetai
  
  # 2. Calculate H_alpha
  # (1) Calculate Ubeta_star = U_beta_dalpha
  
  # (2) Calculate V_alpha: separate process, see below
  
}

design_v1<-cbind(rep(1,nrow(RS1)),RS1$EVALUE100,RS1$SMOKE,RS1$NATIVE,RS1$HISPANIC,RS1$COMP5,RS1$COMP6,RS1$yy2)
design_v2<-cbind(rep(1,nrow(RS2)),RS2$EVALUE100,RS2$SMOKE,RS2$NATIVE,RS2$HISPANIC,RS2$COMP5,RS2$COMP6,RS2$yy2)
Ualpha_v1<-rep(0,1+p+q-3)
Ualpha2_v1<-matrix(0,nrow=1+p+q-3,ncol=1+p+q-3)
for(i in 1:nrow(RS1)){
  Ualphai<-design_v1[i,]*(RS1$VALUE100[i]-sum(alpha1[1:8]*design_v1[i,]))
  Ualpha_v1<-Ualpha_v1+Ualphai
  Ualpha2i<-as.matrix(design_v1[i,])%*%t(as.matrix(design_v1[i,]))
  Ualpha2_v1<-Ualpha2_v1+Ualpha2i
}
Ualpha_v2<-rep(0,1+p+q-3)
Ualpha2_v2<-matrix(0,nrow=1+p+q-3,ncol=1+p+q-3)
for(i in 1:nrow(RS2)){
  Ualphai<-design_v2[i,]*(RS2$VALUE100[i]-sum(alpha2[1:8]*design_v2[i,]))
  Ualpha_v2<-Ualpha_v2+Ualphai
  Ualpha2i<-as.matrix(design_v2[i,])%*%t(as.matrix(design_v2[i,]))
  Ualpha2_v2<-Ualpha2_v2+Ualpha2i
}
Ualpha<-c(Ualpha_v1,Ualpha_v2)
Ualpha_sq<-as.matrix(Ualpha)%*%t(as.matrix(Ualpha))
Ualpha2<-rbind(cbind(Ualpha2_v1,matrix(0,nrow=nrow(Ualpha2_v1),ncol=ncol(Ualpha2_v2))),cbind(matrix(0,nrow=nrow(Ualpha2_v2),ncol=ncol(Ualpha2_v1)),Ualpha2_v2))
Valpha<-solve(Ualpha2)%*%Ualpha_sq%*%solve(Ualpha2)


RS1$rsi<-0
RS2$rsi<-1
RS<-rbind(RS1,RS2)
em<-lm(VALUE100~(EVALUE100+SMOKE+NATIVE+HISPANIC+COMP5+COMP6)*rsi, data=RS)



#### Exploration: Manual RRC ####

# Manual RRC
valid_data$VALUE100<-valid_data$VALUE/100
valid_data$EVALUE100<-valid_data$EVALUE/100
valid_data$yy1<-ifelse(valid_data$YEAR<70,1,0)
valid_data$yy2<-ifelse(valid_data$YEAR>=70&valid_data$YEAR<80,1,0)
valid_data$yy3<-ifelse(valid_data$YEAR>=80&valid_data$YEAR<90,1,0)
valid_data$yy4<-ifelse(valid_data$YEAR>=90&valid_data$YEAR<100,1,0)
valid_data$yy5<-ifelse(valid_data$YEAR>=100,1,0)
valid_data$YREXP<-valid_data$yrexp
valid_data$hispanic<-valid_data$HISPANIC
valid_data$native<-valid_data$NATIVE

final_data_1$PREDYREXP<-rep(0,nrow(final_data_1))

death_times<-unique(final_data_1$AGE_TV[which(final_data_1$CASE==1)])
for(i in 1:length(death_times)){
  d_time<-death_times[i]
  RS_m<-which(final_data_1$AGE_TV==d_time&final_data_1$AGE_OUT>=d_time)
  RS_v<-valid_data[which(as.matrix(valid_data$AGE_TV)==d_time),]
  
  emodel<-lm(VALUE~YREXP+SMOKE+native+hispanic+COMP5+COMP6+yy2, data=RS_v)
  design<-final_data_1[RS_m,c("YREXP","SMOKE","hispanic","native","COMP5","COMP6","yy2")]
  design$PREDYREXP<-predict(emodel,design)
  
  final_data_1$PREDYREXP[RS_m]<-design$PREDYREXP
}

for(i in 1:length(death_times)){
  d_time<-death_times[i]
  RS_v<-valid_data[which(as.matrix(valid_data$AGE_TV)==d_time),]
  print(nrow(RS_v))
}
# many 0 risk sets


mv_exp$YREXP100_em<-NA
mv_exp$YREXP100_em[which(mv_exp$surf!=1)]<-mv_exp$YREXP100[which(mv_exp$surf!=1)]
mv_exp$YREXP100_em[which(mv_exp$surf==1&mv_exp$age_tv<=45)]<-mv_exp$YREXP100_em1[which(mv_exp$surf==1&mv_exp$age_tv<=45)]
mv_exp$YREXP100_em[which(mv_exp$surf==1&mv_exp$age_tv>45)]<-mv_exp$YREXP100_em2[which(mv_exp$surf==1&mv_exp$age_tv>45)]

mv_exp$CUMEXP100_em<-ave(mv_exp$YREXP100_em,factor(mv_exp$id),FUN=cumsum)
mv_exp$age_out<-mv_exp$age_tv+mv_exp$time

final_data_1$YREXP_RRC[which(final_data_1$surf!=1)]<-final_data_1$YREXP[which(final_data_1$surf!=1)]
final_data_1$YREXP_RRC[which(final_data_1$surf==1)]<-final_data_1$PREDYREXP[which(final_data_1$surf==1)]
final_data_1$CUMEXP_RRC[which(final_data_1$surf==1)]<-ave(final_data_1$YREXP_RRC,factor(final_data_1$ID),FUN=cumsum)
final_data_1$CUMEXP_RRC100<-final_data_1$CUMEXP_RRC/100

model.coxph.rrc.naive.b<-coxph(Surv(time,CASE)~(CUMEXP_RRC100+SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5)+strata(AGE_TV),data=final_data_1,ties="breslow")
model.coxph.rrc.adjusted.b<-coxph(Surv(time,CASE)~(CUMEXP_RRC100+SMOKE+native+hispanic+COMP5+COMP6+yy2+yy3+yy4+yy5)+strata(AGE_TV),data=final_data_1,ties="breslow")
