
*** A plot for the DAME and MEM estimates for the analysis presented in
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

gen wt=1
qui sum threshold
loc mn=r(min)
loc mx=r(max)
xtile group_id = threshold, nq(10)
egen midpoint=median(threshold),by(group_id)
putmata wt=wt group_id=midpoint X=(polarization threshold polarization_threshold seatshare seatshare_2 incompatibility asymmetry asym_seat 1) Z=(1 threshold), replace
mata: dame=me_wt(coef, X, Z, group_id, wt)

** Marginal effects at means
collapse (mean) polarization seatshare incompatibility asymmetry [fw=wt]
expand 21
gen threshold=`mn' + (_n-1)*(`mx'-`mn')/20
gen polarization_threshold=polarization*threshold
gen seatshare_2=seatshare^2
gen asym_seat=seatshare*asymmetry

putmata X=(polarization threshold polarization_threshold seatshare seatshare_2 incompatibility asymmetry asym_seat 1) Z=(1 threshold), replace
mata: mem=me_byrow(coef, X, Z)
getmata (mem lbm ubm)=mem
mata: group = dame[.,1]
getmata (midpoint obs dame_est lb ub)=dame, force

* Plot the DAME and MEM estimates
/* Note that the replication do file uses additional graphical parameters, which leads to different axis and legend labels from this minimal example. */
twoway (line mem threshold, lpattern(solid)) ///
(rline lbm ubm threshold, lpattern(dash)) ///
(rspike lb ub midpoint) ///
(scatter dame_est midpoint [fw=obs], msymbol(o) msize(*.25)), /// 
yline(0, lcolor(red)) ytitle("ME of polarization") xtitle("Effective Electoral Threshold") legend(off) ylab(-.0006(.0021).0078)

