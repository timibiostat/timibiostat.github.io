/* 
#############################################################################
** AUTHOR:  Michael Palazzolo
** DATE:    July 7, 2023
** Reference: Klien et al. SAS and R functions to Compute Pseudo-values for Censored Data Regression. CMPB. 2008
** ----------------------------------------------------------------
** PURPOSE: This macro derives pseudo-values for each observation based on Kaplan-Meier estimates for the survival curve.
Pseudo-values are based on the leave-one-out KM estimator.
This macro derives pseudo-values for a specific endpoint.
The output dataset will contain the variable 'pseudo' which contains the pseudo-values.
This macros run time is dependent on the size of the input dataset.
** ----------------------------------------------------------------
##############################################################################
*/


*** ===============================;
*** Description of MACRO variables ;
*** ===============================;
*** INDATA = SAS dataset for pseudo values;
*** ID = Unique subject ID;
*** TIME = Time variable in survival analysis;
*** CENSOR = Censoring indicator variable;
*** TIMEPOINT = Timepoint to calculate the pseudo values at (can me multiple timepoints);
*** HOWMANY = Number of observations in the dataset;
*** OUTDATA = Output dataset;
*** ================================;




*** #############################################################;
*** MACRO starts here;
*** #############################################################;
%macro pseudoD (indata=,id=,time=,censor=,timepoint=,howmany=,outdata=);
/***************************************************************************************/
/************************ OVERALL KM ESTIMATE (THETA) **********************************/
/***************************************************************************************/
proc lifetest data=&indata. noprint plots=none timelist=&timepoint. reduceout outsurv=sall;
time &time.*&censor.(0);
run;
/*proc print data=sall; run;*/

data sall;
set sall;
theta = survival;
keep theta;
run;

data sout;
set &indata.;
run;
*=======================================================================================;

/***************************************************************************************/
/************ KM ESTIMATE BASED ON LEAVE-ONE-OUT ESTIMATOR (THETAMINI) *****************/
/***************************************************************************************/
%do ip=1 %to &howmany.;
proc lifetest data=&indata. noprint plots=none timelist=&timepoint. reduceout outsurv=salli;
time &time.*&censor.(0);
where &id. ne &ip.;
run;

data salli2;
set salli;
thetamini = survival;
&id. = &ip.;
keep &id. thetamini;
run;
/*proc print data=salli2; run;*/

data souti;
merge salli2 sall;
run;
/*proc print data=souti; run;*/

proc sort data=sout; by &id.; run;
proc sort data=souti; by &id.; run;
data sout;
merge sout souti;
by &id.;
run;
/*proc print data=sout; run;*/
%end;
*=======================================================================================;

/***************************************************************************************/
/*************************** PSEUDO-VALUE CALCULATION **********************************/
/***************************************************************************************/
data sfinal;
set sout;
pseudo = 1 - (&howmany. * theta - (&howmany.-1) * thetamini);
run;

data &outdata.;
	set sfinal;
run;
*=======================================================================================;
%mend pseudoD;
