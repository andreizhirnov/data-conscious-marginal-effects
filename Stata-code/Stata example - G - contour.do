
*** A contourplot for the marginal effects of polarization for the analysis presented in
*** Golder, Sona. 2006. The Logic of Pre-Electoral Coalition Formation. Columbus: Ohio State University Press.
*** The dataset is available from Matt Golder's website at: http://mattgolder.com/files/interactions/interaction3.zip" (named "interaction3.dta")

clear all

* Specify the function that returns the partial derivative of the predicted values of the dependent variable; 
* it uses a matrix of covariate values needed to compute the linear prediction of the model (x) and a matrix of covariate values needed to compute the linear component of the first derivative of the predicted value of the dependent variable (z).
mata
real matrix me(coef, x, z) {
int_coef_names = ("polarization","polarization_threshold")				/* The coefficients used to compute the derivative of the linear prediction */
nam = st_matrixcolstripe("beta")
k = J(cols(int_coef_names),1,.)
for (j=1; j<=cols(int_coef_names); j++) {
 k[j] = selectindex(nam[.,2]:==int_coef_names[j])
}
dydx = (coef[.,k]*z'):*normalden(coef*x') 								/* Replace normalden() with the derivative of the inverse link function as needed */
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
use "interaction3.dta",clear

xtprobit pec polarization threshold polarization_threshold seatshare seatshare_2 incompatibility asymmetry asym_seat, re i(ident) 
keep if e(sample)
matrix beta=e(b)[.,e(depvar) + ":"]
matrix vcov=e(V)[e(depvar) + ":",e(depvar) + ":"]

preserve

* Simulate the coefficients
drawnorm coef1-coef`=colsof(beta)', n(10000) means(beta) cov(vcov) clear
putmata coef=(*), replace
restore

* Create the necessary datasets
foreach var of varlist seatshare incompatibility asymmetry {
qui sum `var'
replace `var'=`=r(mean)'
}

collapse (mean) seatshare incompatibility asymmetry (count) obs=seatshare, by(polarization threshold)

gen polarization_threshold=polarization*threshold
gen seatshare_2=seatshare^2
gen asym_seat=seatshare*asymmetry

putmata wt=obs X=(polarization threshold polarization_threshold seatshare seatshare_2 incompatibility asymmetry asym_seat 1) Z=(1 threshold), replace
mata: mebr=me_byrow(coef, X, Z)
getmata (me_est lb ub)=mebr

gen significant=(lb>0 & ub>0)|(lb<0 & ub<0)

* Add additional observations to anchor the marker sizes
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

local colr = "white*.5 yellow*.5 orange*.5 red*.5" /* Color ramp goes from less intense colors to more intense colors */

/* Note that the replication do file uses additional graphical parameters, which leads to different axis and legend labels from this minimal example. */
twoway (contour me_est polarization threshold if me_est!=., ccuts(`locut' `medcut' `hicut') ccolors(`colr')) ///
(scatter polarization threshold [fw=obs] if significant==0, msymbol(oh) mlcolor(black%95) mlwidth(vthin) msize(*.25)) /// 
(scatter polarization threshold [fw=obs] if significant==1, msymbol(o) mfcolor(black%95) mlwidth(none) msize(*.25)), ///
xtitle(Threshold) ytitle(Polarization) ztitle("") zlabel(`minest' `locut' `medcut' `hicut' `maxest') ///
legend(off)  clegend(title("Effect Size", size(medsmall) pos(12) justification(right)) width(5) height(25))
