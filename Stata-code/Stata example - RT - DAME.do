
*** A plot for the DAME and MEM estimates for the analysis presented in
*** Robertson, Graeme B., and Emmanuel Teitelbaum. 2011. "Foreign Direct Investment, Regime Type, and Labor Protest in Developing Countries." American Journal of Political Science 55(3): 665–677.
*** The dataset is part of the replication materials and can be downloaded from Emmanuel Teitelbaum's website: https://home.gwu.edu/~ejt/pages/Data_files/Robertson%20Teitelbaum%202011.dta"

clear all

* Specify the function that returns the partial derivative of the predicted values of the dependent variable; 
* it uses a matrix of covariate values needed to compute the linear prediction of the model (x) and a matrix of covariate values needed to compute the linear component of the first derivative of the predicted value of the dependent variable (z).
mata 
real matrix me(coef, X, Z) {
int_coef_names = ("l_l_flows","l_demflows") 	/* The coefficients used to compute the derivative of the linear prediction */
nam=st_matrixcolstripe("beta")
k=J(cols(int_coef_names),1,.)
for (j=1; j<=cols(int_coef_names); j++) {
 k[j]=selectindex(nam[.,2]:==int_coef_names[j])
}
dydx=(coef[.,k]*Z'):*exp(coef*X') 				/* Replace exp() with the derivative of the inverse link function as needed */
return(dydx)
}
end

* Specify the function that returns a vector of marginal effects by row; this function uses me() internally
mata
real matrix me_byrow(coef, X, Z) {
dydx=me(coef, X, Z)
means=mean(dydx)'
ra=mm_quantile(dydx, 1, (0.025 \ 0.975))'		 /* Confidence level can be changed here */
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
ra=mm_quantile(dydxw, 1, (0.025 \ 0.975))'		/* Confidence level can be changed here */
return((groups,obs,means,ra))
}
end

* Load the data and estimate the model
use "Robertson Teitelbaum 2011.dta",clear

tsset country year 
gen l_l_flows=L.l_flows
gen l_polity2=L.polity2
gen l_dispute=L.dispute
gen l_demflows=l_l_flows*l_polity2
 
xtnbreg dispute l_l_flows l_polity2 l_demflows l_dispute open_penn l_gdp_pc_penn gdp_grth inflation_1 urban xratchg l_pop time, re 
keep if e(sample)
matrix beta=e(b)[.,e(depvar) + ":"]
matrix vcov=e(V)[e(depvar) + ":",e(depvar) + ":"]

preserve

* Simulate the coefficients
drawnorm coef1-coef`=colsof(beta)', n(10000) means(beta) cov(vcov) clear
putmata coef=(*), replace
restore

gen wt=1
qui sum l_polity2
loc mn=r(min)
loc mx=r(max)
xtile group_id = l_polity2, nq(4)
egen midpoint=median(l_polity2),by(group_id)
putmata wt=wt group_id=midpoint X=(l_l_flows l_polity2 l_demflows l_dispute open_penn l_gdp_pc_penn gdp_grth inflation_1 urban xratchg l_pop time 1) Z=(1 l_polity2) , replace
mata: dame=me_wt(coef, X, Z, group_id, wt)

** Marginal effects at means
collapse (mean) l_dispute open_penn l_gdp_pc_penn gdp_grth inflation_1 urban xratchg l_pop time l_l_flows [fw=wt]
expand 21
gen l_polity2=`mn' + (_n-1)*(`mx'-`mn')/20
gen l_demflows = l_l_flows*l_polity2

putmata X=(l_l_flows l_polity2 l_demflows l_dispute open_penn l_gdp_pc_penn gdp_grth inflation_1 urban xratchg l_pop time 1) Z=(1 l_polity2) , replace
mata: mem=me_byrow(coef, X, Z)
getmata (mem lbm ubm)=mem
mata: group = dame[.,1]
getmata (midpoint obs dame_est lb ub)=dame, force

* Plot the DAME and MEM estimates
/* Note that the replication do file uses additional graphical parameters, which leads to different axis and legend labels from this minimal example. */
twoway (line mem l_polity2, lpattern(solid)) ///
(rline lbm ubm l_polity2, lpattern(dash)) ///
(rspike lb ub midpoint) ///
(scatter dame_est midpoint [fw=obs], msymbol(o) msize(*.25)), /// 
yline(0, lcolor(red)) ytitle("ME of ln(FDI flows)") xtitle("Polity 2") legend(off)