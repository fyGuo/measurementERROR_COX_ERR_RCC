/*===================================================================================
   Macro RRC - measurement error correction for time-varying exposures
   for Cox Proportional Hazard model.

  This macro can be used for
    
  1: time-varying cumulative average exposure
  2: b~a-month moving average exposure
  3: time-varying cumulative total exposure
  4: time-varying point exposure
  5: time-independent exposure

  Residual method is not recommended for the analysis anymore. Please use residual=0 if possible.  
    
===================================================================================*/


%macro rrc
(
id=,                /* Variable name of ID, required */
main=,              /* Name of main dataset, required */
validation=,        /* Name of validation dataset, required */
surrogate=,         /* exposure with error in main dataset, required */
true=,              /* exposure without error in main and validation dataset, required */
confounder= ,       /* covariates adjusted in the model, optional */
extravariableinval=,   /* extra variables included in the measurement error model, NOT in the primary regression model*/  
extravariableinmain=, /* extra variables included in the primary regression model, NOT in the measurement error model */    
case= ,             /* variable that indicates whether an event occurs or not, required */
time= ,             /* variable for time-to-event outcome, required */
missing=,           /* missing indicator, new added */
period=,            /* period indicator for exposure, new added */
groupnum=5,         /* the minimum number of subjects in each risk set, optional,  default=5 */
increments=1,       /* the unit of expsoure, optional, default=1 */
type=1,             /* type of data,
                       1: time-varying cumulative average exposure
                       2: b~a-month moving average exposure
                       3: time-varying cumulative total exposure
                       4: time-varying point exposure
                       5: time-independent exposure */
month=12,           /* for 12-month average exposure with month1=1, could be any other integer */
month1=1,           /* for b~a-month moving average exposure. could be any other integer and month1<month */   
emodel=1,           /* if model=1, it assumes general linear measurement error model;
                       if model=2, it assumes additive measurement error model.
                       only required for type=1. */
icc=,               /* the correlation of repeated measurements of true exposure,
                       only required for type=1, emodel=2 */
dampfactor=,        /* the damp factor of the covariance matrix for true expsoure, only required for type=1, emodel=2*/
timestart=,         /* the start of follow-up time, for time-independent exposure, only required for type=5*/
timeend=,           /* the end of follow-up time, for time-independent exposure, only required for type=5*/
residual=0,         /* if residual=1,do residual analysis for dimension reduction; otherwise, residual=0,*/
initialvalue=,     /* the initial value for all the exposure and covariates */
xupper=,          /* the upper bound for all the exposure and covariates */
xlower=,         /* the lower bound for all the exposure and covariates */
filename=RRCOutput.txt,   /*filename for the output file*/
csvfilename=none,       /* filename for CSV output */
path=/tmp,                 /* path to directory to put data files, you can use . to output everything in your own directory */
    
/* The following option is only for Type=2 now */
    
truepredict=, /* Name of the variable for predicted true value */
surrogatefunc=,  /* Name of the variable for the function (such as 12 month) of surroagte exposure for output */
outdata=  /* Name of the output data for the main study with both the surrogate and the predicted true exposure */    
);

/* Check some arguments */
%let nex=%sysfunc(countw(&true,' '));
%let nexs=%sysfunc(countw(&surrogate,' '));
%let ncovm=%sysfunc(countw(&confounder,' '));
%let ncove=%sysfunc(countw(&extravariableinval,' '));
%let ncovem=%sysfunc(countw(&extravariableinmain,' '));

%if &type ne 2 %then %do;
    %let ncov=%sysfunc(countw(&confounder,' '));
%end;

%if &nex^=&nexs %then %goto out1;
%if &nex=0 or &nexs=0 %then %goto out4;
%let in_main=%sysfunc(countw(&main,' '));
%if &in_main=0 %then %goto out2;
%let in_vali=%sysfunc(countw(&validation,' '));
%if &in_vali=0 %then %goto out3;

/*indicator for residual analysis*/
%let isres=0;

data mainy0;
%if &type=5 %then %do;
    set &main;
    do &time=&timestart to &timeend;
        if &timestart<&timeend then do;
            if &time<&timeend then casenew=0;
            if &time=&timeend then casenew=&case;
            end;
        else do;
            casenew=&case;
            end;
        output;
    end;
%end;
%else %do;
 set &main;
 casenew=&case;
%end;
run;

/*sort the data*/
proc sort data=mainy0;
    by &id &time;
   run;

proc sort data=&validation;
    by &id &time;
   run;

data validationy0;
    set &validation;
    nindex=_N_;
run;


    

/* Not recommending residual analysis for now */

%if &residual=1 %then %do;

    %put ****** The residual analysis is not recommended for this analysis ********;
    %put ****** Please use residual=0 **********;
    endsas;

%end;


/* this else means if residual = 0 */
%else %do;

 %if &type=2 %then %do;

        data mainy0;
            set mainy0;
            array sur{&nex} &surrogate;

            do i=1 to &nex;
                if sur{i}=. then sur{i}=0;
                end;
        run;

        %if &ncovm>0 %then %do;
            data mainy0;
                set mainy0;
                array cov{&ncovm} &confounder;        

                do i=1 to &ncovm;
                    if cov{i}=. then delete;
                    end;
            run;
            
            %end;
        
        %if &ncove>0 %then %do;
            data mainy0;
                set mainy0;
                array ecov{&ncove} &extravariableinval;
                
                do i=1 to &ncove;
                    if ecov{i}=. then delete;
                    end;
            run;
            
            %end;

        %if &ncovem>0 %then %do;
            data mainy0;
                set mainy0;
                array ecov{&ncovem} &extravariableinmain;
                
                do i=1 to &ncovem;
                    if ecov{i}=. then delete;
                    end;
            run;
            
            %end;

 %end;

/*corresponding to other types*/

 %else %do;

        data mainy0;
            set mainy0;
            array sur{&nex} &surrogate;
            do i=1 to &nex;
                if sur{i}=. then delete;
                end;
        run;

    %if &ncovm>0 %then %do;
        data mainy0;
            set mainy0;
            array cov{&ncovm} &confounder;

            do i=1 to &ncovm;
                if cov{i}=. then delete;
                end;
            
            run;
            
     %end;

 
     %if &ncove>0 %then %do;
         data mainy0;
             set mainy0;
             array ecov{&ncove} &extravariableinval;
                
             do i=1 to &ncove;
                if ecov{i}=. then delete;
                end;
         run;
            
      %end;

      %if &ncovem>0 %then %do;
          data mainy0;
              set mainy0;
              array ecov{&ncovem} &extravariableinmain;
                
              do i=1 to &ncovem;
                 if ecov{i}=. then delete;
                 end;
          run;
            
       %end;
 
        
 %end;


/*convert sas dataset into .dat file for Fortran to read data*/
%let fname=main;
%if &type=2 %then %do;

  data mainy;
  set mainy0;
  file "&path./&fname..&sysjobid..dat";
  put &id &time casenew &surrogate &extravariableinval &confounder &extravariableinmain &period &missing;
  run;
  
%end;
%else %if &type=1 and &emodel=2 %then %do;

  data mainy;
  set mainy0;
  file "&path./&fname..&sysjobid..dat";
  put &id &time casenew &surrogate &confounder;
  run;

%end;
%else %do;

  data mainy;
  set mainy0;
  file "&path./&fname..&sysjobid..dat";
  put &id &time casenew &surrogate &extravariableinval &confounder &extravariableinmain;
  run;
  
%end;

        

%let fname=validation;
%if &emodel=1 %then %do;

   data validationy;
   set validationy0;
   file "&path./&fname..&sysjobid..dat";
   put &id &time &true &surrogate &extravariableinval &confounder;
   run;
   
%end;

%else %do;

   data validationy;
   set validationy0;
   file "&path./&fname..&sysjobid..dat";
   put &id &time &true &surrogate &confounder;
   run;

%end;


/*corresponding to the end of the residual analysis part*/
%end;

title "Main study";
proc print data=mainy0(obs=100);
    var &id &time casenew &surrogate &extravariableinval &extravariableinmain &period &missing;
run;
title;

title "Validation study";
proc print data=validationy0;
    var &id &time &true &surrogate &extravariableinval;
run;
title;





proc sql;
select count(&id) into: mainobs
from mainy;
quit;

proc sql;
create table mainid as
select unique(&id),count(&time) as idobs
from mainy group by &id;
quit;

proc sql;
select count(&id) into: nmain
from mainid;
quit;

proc sql;
select max(idobs) into: nage_main
from mainid;
quit;

data mainoutcome;
set mainy;
where casenew=1;
run;

proc sort data=mainoutcome nodupkey;
by &time;
run;

proc sql;
select count(&time) into: nageobs
from mainoutcome;
quit;

proc sql;
select count(&id) into: valobs
from validationy;
quit;

proc sql;
create table valid as
select unique(&id),count(&time) as validobs
from validationy group by &id;
quit;

proc sql;
select count(&id) into: nval
from valid;
quit;

proc sql;
select max(validobs) into: nage_val
from valid;
quit;

    

/*This part is particularly for the correlation of repeated measurements.*/

%let fname=covpara;
%if &type=1 and &emodel=2 %then %do;
data covpara;
rho=&icc;
theta=&dampfactor;
run;
data _null_;
set covpara;
file "&path./&fname..&sysjobid..dat";
put rho theta;
run;
%end;


/* This part is to output the initial values for the exposure and covariates */
%let fname=xinitial;
data initial;
do i=1 to (&nex+&ncovm+&ncovem);
inival=scan("&initialvalue",i," ");
if inival=. then inival=0;
output;
end;
run;

data _null_;
set initial;
file "&path./&fname..&sysjobid..dat";
put inival;
run;

/* This part is to output the upper bound values for the exposure and covariates */
%let fname=xupbound;
data xup;
do i=1 to (&nex+&ncovm+&ncovem);
upval=scan("&xupper",i," ");
if upval=. then upval=5;
output;
end;
run;

data _null_;
set xup;
file "&path./&fname..&sysjobid..dat";
put upval;
run;


/* This part is to output the lower bound values for the exposure and covariates */
%let fname=xlowbound;
data xlow;
do i=1 to (&nex+&ncovm+&ncovem);
lowval=scan("&xlower",i," ");
if lowval=. then lowval=-5;
output;
end;
run;

data _null_;
set xlow;
file "&path./&fname..&sysjobid..dat";
put lowval;
run;

    

/*This part is to output the units for the exposures*/
%let fname=exunit;
data increment;
do i=1 to &nex;
unit=scan("&increments",i, " ");
output;
end;
run;
data _null_;
set increment;
file "&path./&fname..&sysjobid..dat";
put unit;
run;


/*This part is to output variable names*/
%let fname=expnames;
data expnames;
do i=1 to &nex;
truename=scan("&surrogate",i);
output;
end;
run;
data _null_;
set expnames;
file "&path./&fname..&sysjobid..dat";
put truename;
run;


%let fname=covnames;
data covnames;
do i=1 to &ncovm;
covname=scan("&confounder",i);
output;
end;

do i=1 to &ncovem;
covname=scan("&extravariableinmain",i);
output;
end;
run;

data _null_;
set covnames;
file "&path./&fname..&sysjobid..dat";
put covname;
run;

    

/* Output run parameters */
%let fname=runparams;
data runparams;
p1 = scan("&nmain",1);
p2 = scan("&nval",1);
p3 = scan("&nageobs",1);
p4 = scan("&nage_val",1);
p5 = scan("&nage_main",1);
p6 = scan("&nex",1);
%if &type=1 and &emodel=2 %then %do;
p7 = scan("&ncov",1);
%end;
p8 = scan("&isres",1);
%if &type=1 and &emodel=1 %then
  %do;
    p9 = scan("&type",1);
  %end;
%else  %if &type=1 and &emodel=2 %then
  %do;
    p9 = scan("&emodel",1);
  %end;
%else
  %do;
    p9 = scan("&type",1);
  %end;
p10 = scan("&mainobs",1);
p11 = scan("&valobs",1);
p12 = scan("&groupnum",1);
%if &type=2 %then %do;
p13 = scan("&month",1);
p14 = scan("&month1",1);
%end;
%if &emodel=1 %then %do;
p15 = scan("&ncovm",1);
p16 = scan("&ncove",1);
p17 = scan("&ncovem",1);
%end;
output;
run;
data _null_;
set runparams;
file "&path./&fname..&sysjobid..dat";
put p1;
put p2;
put p3;
put p4;
put p5;
put p6;
%if &type=1 and &emodel=2 %then %do;
put p7;
%end;
put p8;
put p9;
put p10;
put p11;
put p12;
%if &type=2 %then %do;
put p13;
put p14;
%end;
%if &emodel=1 %then %do;
put p15;
put p16;
put p17;    
%end;
run;

/*Call Unix command to modify RRC Fortran program based on arguments from the users*/
/*
                       1: time-varying cumulative average exposure
                       2: b~a-month moving average exposure
                       3: time-varying cumulative total exposure
                       4: time-varying point exposure
                       5: time-independent exposure
*/

%let mydir=/udd/n2bch/SAS/macros/RRC/code;

/*%let mydir=/udd/resjp/RRC-1.3-ALL/;*/
    

%if &type=1 and &emodel=1 %then
  %do;
    %let myprog=&mydir.rrc_cum;
    %put &myprog &path &sysjobid;
    x "&myprog &path &sysjobid";
  %end;
%else  %if &type=1 and &emodel=2 %then
  %do;
    %let myprog=&mydir.rrc_timevarying_method_GENEnew;
    x "&myprog &path &sysjobid";
  %end;
%else  %if &type=2 %then
  %do;
    %let myprog=&mydir.rrc_abmon;
    x "&myprog &path &sysjobid";
  %end;
%else  %if &type=3 %then
  %do;
    %let myprog=&mydir.rrc_cum;
    x "&myprog &path &sysjobid";
  %end;
%else
  %do;
    %let myprog=&mydir.rrc_new;
    x "&myprog &path &sysjobid";
  %end;


%if &type=2 and "&outdata" ne "" %then %do;

/* merge the xsurrogate data with the xhat data */
%let fname=xhat;
data _main_xhat;
    infile "&path./&fname..&sysjobid..dat";
    input &id &time &truepredict;
run;

proc sort data=_main_xhat;
    by &id &time;
run;


%let fname=xsurrogate;
data _main_xsurr;
    infile "&path./&fname..&sysjobid..dat";
    input &id &time &surrogatefunc;
run;

proc sort data=_main_xsurr;
    by &id &time;
run;
      

data &outdata;
    merge _main_xhat _main_xsurr;
    by &id &time;
run;


data &outdata;
    merge &outdata(in=small) mainy0;
    by &id &time;
    if small;
run;    
    
data "&outdata";
    set &outdata;
run;

%end;

%if &type=2 %then %do;

%let fname=xhat;
x "rm &path./&fname..&sysjobid..dat";

%let fname=xsurrogate;
x "rm &path./&fname..&sysjobid..dat";  

%end;

%let fname=main;
x "rm &path./&fname..&sysjobid..dat";
%let fname=validation;
x "rm &path./&fname..&sysjobid..dat";
%let fname=exunit;
x "rm &path./&fname..&sysjobid..dat";
%let fname=expnames;
x "rm &path./&fname..&sysjobid..dat";

%let fname=covnames;
x "rm &path./&fname..&sysjobid..dat";

%let fname=runparams;
x "rm &path./&fname..&sysjobid..dat";
%let fname=covpara;
%if &type=1 and &emodel=2 %then %do;
x "rm &path./&fname..&sysjobid..dat";
%end;


%let fname=xinitial;
x "rm &path./&fname..&sysjobid..dat";

%let fname=xupbound;
x "rm &path./&fname..&sysjobid..dat";

%let fname=xlowbound;                 
x "rm &path./&fname..&sysjobid..dat";


%let out=Output;
x "mv &path./&out..&sysjobid..txt &filename";
x "mv &path./&out..&sysjobid..csv &csvfilename";
x "rm -f none";

%goto out;

%out1: %err_msg(1);
%out2: %err_msg(2);
%out3: %err_msg(3);
%out4: %err_msg(4);
%out:;
quit;
run;
%mend;

%macro err_msg(err_num);
   %put **********************  ERROR IN SAS MACRO RRC:    **********************;
   %put ******                                                                  ******;

   %if &err_num = 1 %then %do;
     %put ******    The number of true predictors and the number of surrogate predictors should be equal.    ******;
         %end;
   %if &err_num = 2 %then %do;
     %put ******    The argument main= is required. You need to specify the name of main dataset.    ******;
         %end;
   %if &err_num = 3 %then %do;
     %put ******    The argument validation= is required. You need to specify the name of validation dataset.    ******;
         %end;
   %if &err_num = 4 %then %do;
     %put ******    The argument true= and surrogate= are required. You need to have at least one true and one surrogare in your datas
ets.    ******;
   %end;

%mend;

