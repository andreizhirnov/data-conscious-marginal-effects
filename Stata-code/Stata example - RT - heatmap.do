
*** A heatmap for the marginal effects of logged FDI flows for the analysis presented in
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
ra=mm_quantile(dydx, 1, (0.025 \ 0.975))'		/* Confidence level can be changed here */
return((means,ra))
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

* Create the necessary datasets
global varlist "l_dispute open_penn l_gdp_pc_penn gdp_grth inflation_1 urban xratchg l_pop time"
foreach y of global varlist {
qui sum `y'
replace `y'=r(mean)
}

egen l_flows1=cut(l_l_flows), at(-9.22, -7.22, -5.22, -3.22, -1.22, 1.22, 3.22, 5.22, 7.22, 9.22)
egen flows_mean=mean(l_l_flows), by(l_flows1)
replace l_l_flows=flows_mean
collapse (mean) l_dispute open_penn l_gdp_pc_penn gdp_grth inflation_1 urban xratchg l_pop time (count) obs=dispute, by(l_l_flows l_polity2) 
gen l_demflows = l_l_flows*l_polity2

putmata X=(l_l_flows l_polity2 l_demflows l_dispute open_penn l_gdp_pc_penn gdp_grth inflation_1 urban xratchg l_pop time 1) Z=(1 l_polity2) , replace
mata: mebr=me_byrow(coef, X, Z)
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

* Color ramp from less intense to more intense colors
loc scolr="yellow*.25"
loc ecolr="red*.95"

* Variable grid
qui sum l_l_flows [fw=obs] 
loc x1max: disp %9.1f r(max)
loc x1min: disp %9.1f r(min)
loc s1=round((`x1max'-`x1min')/4, 0.1)				/* Number of ticks on the y axis can be changed here or changing the "ylab" option below */		
qui sum l_polity2 [fw=obs]
loc x2max: disp %9.4f r(max)
loc x2min: disp %9.4f r(min)
loc s2=round((`x2max'-`x2min')/4, 0.0001)		

/* Note that the replication do file uses additional graphical parameters, which leads to different axis and legend labels from this minimal example. */
twoway (contour me_est l_l_flows l_polity2 if me_est!=., levels(100) crule(linear) scolor(`scolr') ecolor(`ecolr') zlab(#5, labsize(medsmall))) ///
(scatter l_l_flows l_polity2 [fw=obs] if significant==0, msymbol(oh) mlcolor(black%95) mlwidth(vthin) msize(*.25)) ///
(scatter l_l_flows l_polity2 [fw=obs] if significant==1, msymbol(o) mfcolor(black%95) mlwidth(none) msize(*.25)), ///
xsca(alt) ysca(alt) xtitle("") ytitle("") ztitle("") ///
ylab(`x1min'(`s1')`x1max', grid gmax labsize(medsmall)) xlab(`x2min'(`s2')`x2max', labsize(medsmall) grid gmax) /// 
legend(off) clegend(title("Effect Size", size(medsmall) pos(12) justification(right)) ring(0) width(5) height(25)) nodraw name(yx, replace)

twoway histogram l_polity2 [fw=obs], frac ysca(alt reverse) xtitle("Polity 2", size(medsmall)) ytitle("") ///
xlab(`x2min'(`s2')`x2max') ylab(#4, nogrid labsize(medsmall)) ///
fysize(20) fcolor(black%95) lwidth(vthin) lcolor(white%25) nodraw name(hy, replace) 

twoway histogram l_l_flows  [fw=obs], frac xsca(alt reverse) horiz ytitle("ln(FDI flows)", size(medsmall)) xtitle("") ///
ylab(`x1min'(`s1')`x1max', grid gmax labsize(medsmall)) xlab(#4, nogrid labsize(medsmall))  ///
fxsize(20) fcolor(black%95) lwidth(vthin) lcolor(white%25) nodraw name(hx, replace)

gr combine hx yx hy, hole(3) imargin(zero) scale(1.1) xsize(5.5) ysize(5.5)
