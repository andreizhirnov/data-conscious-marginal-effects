

cd "C:/Users/az310/Dropbox/ME paper X/replication Stata vignette/data"

*** A contourplot with the marginal effects of the registration closing date for an analysis presented in
*** Nagler, Jonathan. 1991. "The Effect of Registration Laws and Education on U.S. Voter Turnout." _American Political Science Review_ 85(4): 1393–1405.
*** The dataset was made public by William D. Berry, Jacqueline H. R. DeMeritt, and Justin Esarey in the replication materials to their article 
*** "Testing for Interaction in Binary Logit and Probit Models: Is a Product Term Essential?" and
*** can be downloaded from https://jdemeritt.weebly.com/uploads/2/2/7/7/22771764/bde.zip (scobit.dta)

clear all

* specify the function that returns a vector of marginal effects for a given matrix and a function
mata
real matrix me(coef, x, x_new) {
dydx=normal(coef*x_new')-normal(coef*x')								/* Replace normal() with the appropriate function as needed */
return(dydx)
}
end
* specify the function that returns a vector of marginal effects by row; this function uss me() internally
mata
real matrix me_byrow(coef, X, Z) {
dydx=me(coef, X, Z)
means=mean(dydx)'
ra=mm_quantile(dydx, 1, (0.025 \ 0.975))'						        /* Confidence level can be changed here */
return((means,ra))
}
end

** load the data and estimate the model
use N.dta,clear
drop if newvote==-1
probit newvote closing neweduc educ2 cloeduc cloeduc2 age age2 south gov

keep if e(sample)
matrix beta=e(b)[.,e(depvar) + ":"]
matrix vcov=e(V)[e(depvar) + ":",e(depvar) + ":"]

preserve
* simulate coefficients
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

* push data to mata
putmata X=(closing neweduc educ2 cloeduc cloeduc2 age age2 south gov 1), replace
preserve
replace closing = closing+1
replace cloeduc = closing*neweduc
replace cloeduc2 = closing*educ2
putmata X1=(closing neweduc educ2 cloeduc cloeduc2 age age2 south gov 1), replace
restore

mata: mebr=me_byrow(coef, X, X1)
getmata (me_est lb ub)=mebr

gen significant=(lb>0 & ub>0)|(lb<0 & ub<0)

** add extra observations to ancor the scatter sizes
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

** break the ME values into steps
qui sum me_est,detail
loc locut=`r(min)' +(`r(max)'-`r(min)')*1/4
loc medcut=`r(min)'+(`r(max)'-`r(min)')*2/4
loc hicut=`r(min)' +(`r(max)'-`r(min)')*3/4
loc minest =`r(min)'
loc maxest =`r(max)'

local colr = `""red*.5" "orange*.5" "yellow*.5" "white*.5""' /* color ramp: from more intense to less intense colors */

twoway (contour me_est closing neweduc if me_est!=., ccuts(`locut' `medcut' `hicut') ccolors(`colr')) ///
(scatter closing neweduc [fw=obs] if significant==0, msymbol(oh) mlcolor(black%95) mlwidth(vthin) msize(*.25)) /// 
(scatter closing neweduc [fw=obs] if significant==1, msymbol(o) mfcolor(black%95) mlwidth(none) msize(*.25)), ///
xtitle(Education) ytitle(Closing Date) ztitle("") zlabel(`minest' `locut' `medcut' `hicut' `maxest') ///
legend(off)  clegend(title(`"Effect Size"', size(medsmall) pos(12) justification(right)) width(5) height(25)) 

