
*** A plot for the DAME and MEM estimates for the analysis presented in
*** Arceneaux, Kevin, Martin Johnson, Rene Lindstädt, and Ryan J. Vander Wielen. 2016. "The Influence of News Media on Political Elites: Investigating Strategic Responsiveness in Congress." American Journal of Political Science 60(1): 5–29.
*** The dataset is part of the published replication materials and can be downloaded from: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/27597" (named "FoxNews_Master.dta")

clear all

* Specify the function that returns differences in the predicted values of the dependent variable for two matrices with the values of covariates (x and x_new).
* If x_new includes an additional increment added to the main explanatory variable, it can be used in the first-difference method as well.
mata
real matrix me(coef, x, x_new) {
dydx=logistic(coef*x_new') - logistic(coef*x')							/* Replace logistic() with the appropriate link function as needed */
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
qui sum dvprop
loc mn=r(min)
loc mx=r(max)
xtile group_id = dvprop, nq(10)
egen midpoint=median(dvprop),by(group_id)

collapse (count) wt=PartyVote, by(midpoint daystoelection dvprop Amend OtherPass ProPart qualchal qualchal_lag RegPass Retirement Susp seniorit spendgap_lag spendgap distpart_lag ) 
gen daystoelection2=daystoelection^2
gen daystoelection3=daystoelection^3
gen daysdv=daystoelection*dvprop
gen days2dv=daystoelection2*dvprop
gen days3dv=daystoelection3*dvprop

putmata wt=wt group_id=midpoint X=(daystoelection daystoelection2 daystoelection3 dvprop daysdv days2dv days3dv Retirement seniorit qualchal_lag qualchal spendgap_lag spendgap distpart_lag RegPass Susp OtherPass Amend ProPart 1), replace
preserve 
replace daystoelection = daystoelection+1
replace daystoelection2=daystoelection^2 
replace daystoelection3=daystoelection^3
replace daysdv=daystoelection*dvprop
replace days2dv=daystoelection2*dvprop
replace days3dv=daystoelection3*dvprop
putmata X1=(daystoelection daystoelection2 daystoelection3 dvprop daysdv days2dv days3dv Retirement seniorit qualchal_lag qualchal spendgap_lag spendgap distpart_lag RegPass Susp OtherPass Amend ProPart 1), replace
restore

mata: dame=me_wt(coef, X, X1, group_id, wt)

** Marginal effects at means
* Find the modal vote type. The user can also set Amend=1 and all other vote type dummies to 0.
local dummies Amend OtherPass ProPart RegPass Susp
egen baseline = rowmax(`dummies')
replace baseline = 1-baseline
tabstat `dummies' baseline [fw=wt], save
mata props = st_matrix("r(StatTotal)")
mata st_local("modal", st_matrixcolstripe("r(StatTotal)")[selectindex(props :== max(props))[1,1],2])

* Set all other variables to their means or modes
collapse (mean) qualchal qualchal_lag Retirement daystoelection seniorit (median) spendgap_lag spendgap distpart_lag [fw=wt]
foreach v in `dummies' {
	gen `v'=0
}
replace `modal'=1 if "`modal'"! = "baseline"

expand 21
gen dvprop=`mn' + (_n-1)*(`mx'-`mn')/20
gen daystoelection2=daystoelection^2
gen daystoelection3=daystoelection^3
gen daysdv=daystoelection*dvprop
gen days2dv=daystoelection2*dvprop
gen days3dv=daystoelection3*dvprop

putmata X=(daystoelection daystoelection2 daystoelection3 dvprop daysdv days2dv days3dv Retirement seniorit qualchal_lag qualchal spendgap_lag spendgap distpart_lag RegPass Susp OtherPass Amend ProPart 1), replace
preserve 
replace daystoelection = daystoelection+1
replace daystoelection2=daystoelection^2 
replace daystoelection3=daystoelection^3
replace daysdv=daystoelection*dvprop
replace days2dv=daystoelection2*dvprop
replace days3dv=daystoelection3*dvprop
putmata X1=(daystoelection daystoelection2 daystoelection3 dvprop daysdv days2dv days3dv Retirement seniorit qualchal_lag qualchal spendgap_lag spendgap distpart_lag RegPass Susp OtherPass Amend ProPart 1), replace
restore
 
mata: mem=me_byrow(coef, X, X1)
getmata (mem lbm ubm)=mem
getmata (midpoint obs dame_est lb ub)=dame, force

* Plot the DAME and MEM estimates
/* Note that the replication do file uses additional graphical parameters, which leads to different axis and legend labels from this minimal example. */
twoway (line mem dvprop, lpattern(solid)) ///
(rline lbm ubm dvprop, lpattern(dash)) ///
(rspike lb ub midpoint) ///
(scatter dame_est midpoint [fw=obs], msymbol(o) msize(*.25)), /// 
yline(0, lcolor(red)) ytitle("ME of Days to Election") xtitle("Democratic Vote Share") legend(off) ylab(-.0001(.00008).00022)