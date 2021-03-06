

*** A contourplot for the marginal effects of the proximity of the next election for the analysis presented in
*** Arceneaux, Kevin, Martin Johnson, Rene Lindstädt, and Ryan J. Vander Wielen. 2016. "The Influence of News Media on Political Elites: Investigating Strategic Responsiveness in Congress." American Journal of Political Science 60(1): 5–29.
*** The dataset is part of the published replication materials and can be downloaded from: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/27597" (named "FoxNews_Master.dta")

clear all

* Specify the function that returns differences in the predicted values of the dependent variable for two matrices with the values of covariates (x and x_new).
* If x_new includes an additional increment added to the main explanatory variable, it can also be used for the first-difference method.
mata
real matrix me(coef, x, x_new) {
dydx=logistic(coef*x_new') - logistic(coef*x')			/* Replace logistic() with the appropriate link function as needed */
return(dydx)
}
end

* Specify the function that returns a vector of marginal effects by row; this function uses me() internally
mata
real matrix me_byrow(coef, X, Z) {
dydx=me(coef, X, Z)
means=mean(dydx)'
ra=mm_quantile(dydx, 1, (0.025 \ 0.975))'				/* Confidence level can be changed here */
return((means,ra))
}
end

* Load the data and estimate the model
use "FoxNews_Master.dta",clear
gen dvprop=dv/100
gen daysdv=daystoelection*dvprop
gen days2dv=daystoelection2*dvprop
gen days3dv=daystoelection3*dvprop 
logit PartyVote daystoelection daystoelection2 daystoelection3 dvprop daysdv days2dv days3dv Retirement seniorit qualchal_lag qualchal spendgap_lag spendgap distpart_lag RegPass Susp OtherPass Amend ProPart if PresencePartyUnity==1 & Republican==1 & FoxNews==1, cluster(dist2)

keep if e(sample)
matrix beta=e(b)[.,e(depvar) + ":"]
matrix vcov=e(V)[e(depvar) + ":",e(depvar) + ":"]

preserve

* Simulate the coefficients
drawnorm coef1-coef`=colsof(beta)', n(10000) means(beta) cov(vcov) clear
putmata coef=(*), replace
restore

* Create the necessary datasets
foreach x of varlist qualchal qualchal_lag Retirement {
qui sum `x'
replace `x'= (r(mean)>0.5)
}

foreach x of varlist seniorit spendgap_lag spendgap distpart_lag {
qui sum `x'
replace `x'=r(mean)
}

* Find the modal vote type. The user can also set Amend=1 and all other vote type dummies to 0.
local dummies Amend OtherPass ProPart RegPass Susp
egen baseline = rowmax(`dummies')
replace baseline = 1-baseline
tabstat `dummies' baseline, save
mata props = st_matrix("r(StatTotal)")
mata st_local("modal", st_matrixcolstripe("r(StatTotal)")[selectindex(props :== max(props))[1,1],2])
foreach v in `dummies' {
	replace `v'=0
}
replace `modal'=1 if "`modal'"! = "baseline"

* Round the values of the constitutive terms to reduce overplotting
replace dvprop=round(dvprop,0.02) 
replace daystoelection=round(daystoelection,10)
collapse (mean) Amend OtherPass ProPart qualchal qualchal_lag RegPass Retirement Susp seniorit spendgap_lag spendgap distpart_lag (count) obs=Amend, by(daystoelection dvprop) 
gen daystoelection2=daystoelection^2
gen daystoelection3=daystoelection^3
gen daysdv=daystoelection*dvprop
gen days2dv=daystoelection2*dvprop
gen days3dv=daystoelection3*dvprop

putmata wt=obs X=(daystoelection daystoelection2 daystoelection3 dvprop daysdv days2dv days3dv Retirement seniorit qualchal_lag qualchal spendgap_lag spendgap distpart_lag RegPass Susp OtherPass Amend ProPart 1), replace
preserve 
replace daystoelection = daystoelection+1
replace daystoelection2=daystoelection^2 
replace daystoelection3=daystoelection^3
replace daysdv=daystoelection*dvprop
replace days2dv=daystoelection2*dvprop
replace days3dv=daystoelection3*dvprop
putmata X1=(daystoelection daystoelection2 daystoelection3 dvprop daysdv days2dv days3dv Retirement seniorit qualchal_lag qualchal spendgap_lag spendgap distpart_lag RegPass Susp OtherPass Amend ProPart 1), replace
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
loc locut=`r(min)'+(`r(max)'-`r(min)')*1/5 
loc lmedcut=`r(min)'+(`r(max)'-`r(min)')*2/5
loc hmedcut=`r(min)'+(`r(max)'-`r(min)')*3/5
loc hicut=`r(min)'+(`r(max)'-`r(min)')*4/5
loc minest =`r(min)'
loc maxest =`r(max)'

loc colr= "navy*.5 ltblue*.5 white*.5 orange*.5 red*.5" 	/* Color ramp has intense colors at both ends */

/* Note that the replication do file uses additional graphical parameters, which leads to different axis and legend labels from this minimal example. */
twoway (contour me_est daystoelection dvprop if me_est!=., ccuts(`locut' `lmedcut' `hmedcut' `hicut') ccolors(`colr')) ///
(scatter daystoelection dvprop [fw=obs] if significant==0, msymbol(oh) mlcolor(black%95) mlwidth(vthin) msize(*.25)) /// 
(scatter daystoelection dvprop [fw=obs] if significant==1, msymbol(o) mfcolor(black%95) mlwidth(none) msize(*.25)), ///
xtitle(Democratic Vote Share) ytitle(Days to Election) ztitle("") zlabel(`minest' `locut' `lmedcut' `hmedcut' `hicut' `maxest') ///
legend(off)  clegend(title(`"Effect Size"', size(medsmall) pos(12) justification(right)) width(5) height(25)) 