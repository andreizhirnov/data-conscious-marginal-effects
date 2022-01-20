
*** A plot for the DAME and MEM estimates for the analysis presented in
*** Nagler, Jonathan. 1991. "The Effect of Registration Laws and Education on U.S. Voter Turnout." American Political Science Review 85(4): 1393–1405.
*** The dataset was made public by William D. Berry, Jacqueline H. R. DeMeritt, and Justin Esarey as part of the replication materials for their article entitled:
*** "Testing for Interaction in Binary Logit and Probit Models: Is a Product Term Essential?" and can be downloaded from: https://jdemeritt.weebly.com/uploads/2/2/7/7/22771764/bde.zip (named "scobit.dta")

clear all

* Specify the function that returns differences in the predicted values of the dependent variable for two matrices with the values of covariates (x and x_new).
* If x_new includes an additional increment added to the main explanatory variable, it can be used for the first-difference method.
mata
real matrix me(coef, x, x_new) {
dydx=normal(coef*x_new')-normal(coef*x')								/* Replace normal() with the appropriate link function as needed */
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

* Specify the function that returns a vector of weighted average marginal effects; this function uses me() internally
mata
real matrix me_wt(coef, X, Z, group_id, weight) {
dydx=me(coef, X, Z)
groups=uniqrows(group_id)
wtm=J(cols(dydx), rows(groups), .)
obs=J(rows(groups),1,.)
for (i=1; i<=rows(groups); i++) {
	wtmc=(group_id:==groups[i]):*weight
	obs[i]=sum(wtmc)
	wtm[.,i]=wtmc/sum(wtmc)
	}
dydxw=dydx*wtm
means=mean(dydxw)'
ra=mm_quantile(dydxw, 1, (0.025 \ 0.975))'						        /* Confidence level can be changed here */
return((groups,obs,means,ra))
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
collapse (count) wt=newvote, by(closing neweduc age south gov)
gen age2=age^2
gen educ2=neweduc^2
gen cloeduc=closing*neweduc
gen cloeduc2=closing*neweduc^2

* Push the data to mata
putmata wt=wt group_id=neweduc X=(closing neweduc educ2 cloeduc cloeduc2 age age2 south gov 1), replace
preserve
replace closing = closing+1
replace cloeduc = closing*neweduc
replace cloeduc2 = closing*educ2
putmata X1=(closing neweduc educ2 cloeduc cloeduc2 age age2 south gov 1), replace
restore

mata: dame=me_wt(coef, X, X1, group_id, wt)

** Marginal effects at means
qui sum neweduc
loc mn=r(min)
loc mx=r(max)

collapse (mean) age closing neweduc (median) south gov [fw=wt]
expand 21
replace neweduc=`mn' + (_n-1)*(`mx'-`mn')/20
gen educ2=neweduc^2
gen cloeduc=closing*neweduc
gen cloeduc2=closing*educ2
gen age2=age^2

putmata X=(closing neweduc educ2 cloeduc cloeduc2 age age2 south gov 1), replace
preserve
replace closing = closing+1
replace cloeduc = closing*neweduc
replace cloeduc2 = closing*educ2
putmata X1=(closing neweduc educ2 cloeduc cloeduc2 age age2 south gov 1), replace
restore
mata: mem=me_byrow(coef, X, X1)
getmata (mem lbm ubm)=mem
getmata (midpoint obs dame_est lb ub)=dame, force

* Plot the DAME and MEM estimates
/* Note that the replication do file uses additional graphical parameters, which leads to different axis and legend labels from this minimal example. */
twoway (line mem neweduc, lpattern(solid)) ///
(rline lbm ubm neweduc, lpattern(dash)) ///
(rspike lb ub midpoint) ///
(scatter dame_est midpoint [fw=obs], msymbol(o) msize(*.25)), /// 
yline(0, lcolor(red)) ytitle("ME of Closing Date") xtitle("Education") legend(off)
local tn="n-clo-dame"
graph export `tn'.png
graph export `tn'.svg
graph export `tn'.pdf