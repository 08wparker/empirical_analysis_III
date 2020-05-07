log using "~/Desktop/UChicago/Empirical_Analysis/EAIII/partA/ps/ps4/ps4.log", replace

*****************************
******* Problem Set 4 *******
*****************************
**Futing Chen, Hongfan Chen, Will Parker **


clear all
set more off
cd "~/Desktop/UChicago/Empirical_Analysis/EAIII/partA/ps/ps4/data"
use final5.dta


replace avgmath= avgmath-100 if avgmath>100
replace avgmath=. if mathsize==0


*** 1 ***
reg avgmath classize, robust cluster(schlcode)              
reg avgmath classize tipuach c_size, robust cluster(schlcode)

/* Interpretation:
The estimated coefficient of classize drops from .317 (p=.000) to .017 (p=.686).
This result suggests that (1) the percentage of disadvantaged students and 
enrollment are important explanatory variables and (2) omitting them will lead 
to bias in the estimated class size effect. For example, because disadvantaged
students are overrepresented in small classes, the OLS estimate without 
controling for enrollment will be upwardly biased and suggests that students in
larger classes have higher scores.*/


*** 2 ***
* limit the sample to schools with enrollment between 20-60 students
keep if c_size>=20 & c_size<=60
* Generate a (predicted) large class dummy: largeclass=I{enrollment<=40}
gen largeclass1 = c_size <= 40 
* OLS with a linear trend in enrollment
gen diff_40=c_size-40
reg avgmath i.largeclass##c.diff_40 tipuach, robust cluster(schlcode)
/*largeclass: -4.3 (se=1.592, p=.007)*/


*** 3 ***
* Use Local Linear Regression to point estimate large class effect
gen x=0 in 1
lpoly avgmath diff_40 if diff_40 <= 0, deg(1) gen(L) at(x) nograph 
lpoly avgmath diff_40 if diff_40 >= 0, deg(1) gen(R) at(x) nograph
gen rd_largeclass = L - R
list rd_largeclass in 1  // -1.647234 

* Use a nonparametric bootstrap to estimate the standard error 
program rd, rclass
	tempvar x L R 
	gen `x' = 0 in 1 
	lpoly avgmath diff_40 if diff_40 <= 0, deg(1) gen(`L') at(`x') nograph 
	lpoly avgmath diff_40 if diff_40 >= 0, deg(1) gen(`R') at(`x') nograph
	return scalar rd = `L' - `R' in 1
end

bootstrap rd=r(rd), rep(100) seed(123) nowarn: rd  // bootstrap se= 2.534

* Compare these results to the estimates you obtained with OLS.
/* The estimate given by the nonparametric approach is -1.647 (se=2.534, p=.516).
This is much smaller in absolute magnitude than the OLS estimate 
-4.3 (se=1.592, p=.007).*/


*** 4 ***
* Estimate the effect of class size on math scores using fuzzy RDD. 
ivregress 2sls avgmath (classize c.classize#c.diff_40=largeclass ///
	i.largeclass#c.diff_40) diff_40 tipuach, robust cluster(schlcode)
/*-.394(se=.191, p=.039)*/


*** 5 ***
* Use -rdrobust- to estimate the effect of class size on math scores
* Sharp RD
rdrobust avgmath c_size, c(40) all  //1.806(se=5.055, p=.721)
* Fuzzy RD
rdrobust avgmath c_size, c(40) fuzzy(classize) all //-.369(se=.634, p=.56)



*** 6 ***
clear 
use final5.dta
gen pclassize = c_size/(int((c_size - 1)/40) + 1)

label var pclassize "Predicted class size"
bysort c_size: egen avgclassize=mean(classize)
label var avgclassize "Average class size"
tw (line avgclassize c_size) (line pclassize c_size)


*** 7 ***
ivregress 2sls avgmath (classize=pclassize) c_size, robust cluster(schlcode)


*** 8 ***
ivregress 2sls avgmath (classize=pclassize) c_size tipuach, robust cluster(schlcode)
est store e1
/* After controlling for the percentage of disadvantaged students, the estimated
coefficient on class size changes from -.075 (se=.11, p= .499) to
-.223 (se=.1, p=.025). It increases in absolute magnitude and becomes 
statistically significant, suggesting that it is better to condition the RDD on 
the percent of disadvantaged students.*/


*** 9 ***
hist c_size,  width(5) kdensity xline(40 80 120 160 200, lcol(red)) 
/* There is a jump at 40, suggesting that there may be sorting into 
schools that split students into smaller classes.*/


*** 10 ***
gen bin = autocode(c_size , 12, 0 , 240)
bysort bin: egen binclassize=mean(classize)
bysort bin: egen binmath=mean(avgmath)
label var binclassize "Average class size"
label var binmath "Average math score"
label var bin "Enrollment"

tw (conn binclassize bin) (conn binmath bin, yaxis(2)), xtick(0(20)240) xlabel(0(20)240) 



*** 11 ***
gen c_sizesq=c_size^2
reg avgmath c_size c_sizesq
predict mathq
reg classize c_size c_sizesq
predict classizeq
tw (lfit classize c_size) (line classizeq c_size) (scatter binclassize bin), ///
	xtick(0(20)240) xlabel(0(20)240) ///
	xtitle("Enrollment") ytitle("Class size") ///
	legend(lab(1 "Linear") lab(2 "Quadratic")) name(plot3, replace)
tw (lfit avgmath c_size) (line mathq c_size) (scatter binmath bin), ///
	xtick(0(20)240) xlabel(0(20)240) ///
	xtitle("Enrollment") ytitle("Math scores") ///
	legend(lab(1 "Linear") lab(2 "Quadratic")) name(plot4, replace)
graph combine plot3 plot4
/*The polynomial approximation captures the non-linearities well*/



*** 12 ***
* (1) check sensitivity to bandwidths
gen byte discont1= (c_size>=36 & c_size<=45) | (c_size>=76 & c_size<=85) | ///
	(c_size>=116 & c_size<=125) | (c_size>=156 & c_size<=165) | ///
	(c_size>=196 & c_size<=205) | (c_size>=236 & c_size<=245) 
gen byte discont2= (c_size>=38 & c_size<=43) | (c_size>=78 & c_size<=83) | ///
	(c_size>=118 & c_size<=123) | (c_size>=158 & c_size<=163) | ///
	(c_size>=198 & c_size<=203) | (c_size>=238 & c_size<=243) 
ivregress 2sls avgmath (classize=pclassize) c_size tipuach if discont1==1, ///
	robust cluster(schlcode)  //+-5 discontinuity sample
est store d1
ivregress 2sls avgmath (classize=pclassize) c_size tipuach if discont2==1, ///
	robust cluster(schlcode)  //+-3
est store d2
esttab d1 d2 using Check1.rtf, b(3) se(3) star(* 0.05 ** 0.01 *** 0.001) nogaps replace 

* (2) check sensitivity to control for enrollment
gen c_sizecube=c_size^3
ivregress 2sls avgmath (classize=pclassize) c_size tipuach c_sizesq, robust cluster(schlcode)  
est store e2
ivregress 2sls avgmath (classize=pclassize) c_size tipuach c_sizesq c_sizecube, robust cluster(schlcode)  
est store e3
esttab e1 e2 e3 using Check2.rtf, b(3) se(3) star(* 0.05 ** 0.01 *** 0.001) nogaps replace 


*** 13 ***
ivregress 2sls tipuach (classize=pclassize) c_size c_sizesq, robust cluster(schlcode)  


log close
