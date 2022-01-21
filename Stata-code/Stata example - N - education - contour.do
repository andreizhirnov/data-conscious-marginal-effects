
*** A contourplot for the marginal effects of education for the analysis presented in
*** Nagler, Jonathan. 1991. "The Effect of Registration Laws and Education on U.S. Voter Turnout." American Political Science Review 85(4): 1393–1405.
*** The dataset was made public by William D. Berry, Jacqueline H. R. DeMeritt, and Justin Esarey as part of the replication materials for their article entitled:
*** "Testing for Interaction in Binary Logit and Probit Models: Is a Product Term Essential?" and can be downloaded from: https://jdemeritt.weebly.com/uploads/2/2/7/7/22771764/bde.zip (named "scobit.dta")

clear all

* Specify the function that returns differences in the predicted values of the dependent variable for two matrices with the values of covariates (x and x_new).
* If x_new includes an additional increment added to the main explanatory variable, it can be used for the first-difference method.
mata
real matrix me(coef, x, x_new) {
dydx=normal(coef*x_new')-normal(coef*x')								/* Replace normal() with the appropriate function as needed */
return(dydx)
}
end

* Specify the function that returns a vector of marginal effects by row; this function uses me() internally
mata
real matrix me_byrow(coef, X, Z) {
dydx=me(coef, X, Z)
means=mean(dydx)'
ra=mm_quantile(dydx, 1, (0.025 \ 0.975))'						        /* Confidence level can be changed here */
return((means,ra))
}
end

* Load the data and estimate the model
use "scobit.dta",clear
drop if newvote==-1
probit newvote closing neweduc educ2 cloeduc cloeduc2 age age2 south gov

keep if e(sample)
matrix beta=e(b)[.,e(depvar) + ":"]
matrix vcov=e(V)[e(depvar) + ":",e(depvar) + ":"]

preserve

* Simulate the coefficients
drawnorm coef1-coef`=colsof(beta)', n(10000) means(beta) cov(vcov) clear
putmata coef=(*), replace
restore

* Creating the necessary datasets
egen mclo=mean(closing)
egen medu=median(neweduc)
egen age1=mean(age)
egen south1=mode(south)
egen gov1=mode(gov)
collapse (mean) age=age1 south=south1 gov=gov1 mclo medu (count) obs=age1, by(closing neweduc)

gen age2=age^2
gen educ2=neweduc^2
gen cloeduc=closing*neweduc
gen cloeduc2=closing*neweduc^2

* Push the data to mata
putmata X=(closing neweduc educ2 cloeduc cloeduc2 age age2 south gov 1), replace
preserve
replace neweduc = neweduc+1
replace educ2 = neweduc^2
replace cloeduc = closing*neweduc
replace cloeduc2=closing*neweduc^2
putmata X1=(closing neweduc educ2 cloeduc cloeduc2 age age2 south gov 1), replace
restore

mata: mebr=me_byrow(coef, X, X1)
getmata (me_est lb ub)=mebr

gen significant=(lb>0 & ub>0)|(lb<0 & ub<0)

* Add additional observations to anchor marker sizes
gen counter=_n
qui sum counter
loc extra1=`=r(max)'+1
loc extra2=`=r(max)'+2
loc extra3=`=r(max)'+3
loc extra4=`=r(max)'+4
set obs `extra4'

qui sum obs
replace significant=1 in `extra1'/`extra2'
replace significant=0 in `extra3'/`extra4'
replace obs=`=r(min)' in `extra1'/`extra4'
replace obs=`=r(max)' in `extra2'/`extra3'

* Break the ME values into steps
qui sum me_est,detail
loc locut=`r(min)' +(`r(max)'-`r(min)')*1/4
loc medcut=`r(min)'+(`r(max)'-`r(min)')*2/4
loc hicut=`r(min)' +(`r(max)'-`r(min)')*3/4
loc minest =`r(min)'
loc maxest =`r(max)'

local colr = "white*.5 yellow*.5 orange*.5 red*.5" /* Color ramp from less intense to more intense colors */

/* Note that the replication do file uses additional graphical parameters, which leads to different axis and legend labels from this minimal example. */
twoway (contour me_est neweduc closing if me_est!=., ccuts(`locut' `medcut' `hicut') ccolors(`colr')) ///
(scatter neweduc closing [fw=obs] if significant==0, msymbol(oh) mlcolor(black%95) mlwidth(vthin) msize(*.25)) /// 
(scatter neweduc closing [fw=obs] if significant==1, msymbol(o) mfcolor(black%95) mlwidth(none) msize(*.25)), ///
xtitle(Closing Date) ytitle(Education) ztitle("") zlabel(`minest' `locut' `medcut' `hicut' `maxest') ///
legend(off)  clegend(title(`"Effect Size"', size(medsmall) pos(12) justification(right)) width(5) height(25))