
*** A heatmap for the marginal effects of polarization for the analysis presented in
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

collapse (count) obs=pec, by(polarization threshold seatshare incompatibility asymmetry)

gen polarization_threshold=polarization*threshold
gen seatshare_2=seatshare^2
gen asym_seat=seatshare*asymmetry

putmata wt=obs X=(polarization threshold polarization_threshold seatshare seatshare_2 incompatibility asymmetry asym_seat 1) Z=(1 threshold), replace
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
qui sum polarization [fw=obs] 
loc x1max: disp %9.1f r(max)
loc x1min: disp %9.1f r(min)
loc s1=round((`x1max'-`x1min')/4, 0.1)				/* Number of ticks on the y axis can be changed here or changing the "ylab" option below */		
qui sum threshold [fw=obs]
loc x2max: disp %9.4f r(max)
loc x2min: disp %9.4f r(min)
loc s2=round((`x2max'-`x2min')/4, 0.0001)		

/* Note that the replication do file uses additional graphical parameters, which leads to different axis and legend labels from this minimal example. */
twoway (contour me_est polarization threshold if me_est!=., levels(100) crule(linear) scolor(`scolr') ecolor(`ecolr') zlab(#5, labsize(medsmall))) ///
(scatter polarization threshold [fw=obs] if significant==0, msymbol(oh) mlcolor(black%95) mlwidth(vthin) msize(*.25)) ///
(scatter polarization threshold [fw=obs] if significant==1, msymbol(o) mfcolor(black%95) mlwidth(none) msize(*.25)), ///
xsca(alt) ysca(alt) xtitle("") ytitle("") ztitle("") ///
ylab(`x1min'(`s1')`x1max', grid gmax labsize(medsmall)) xlab(`x2min'(`s2')`x2max', labsize(medsmall) grid gmax) /// 
legend(off) clegend(title("Effect Size", size(medsmall) pos(12) justification(right)) ring(0) width(5) height(25)) nodraw name(yx, replace)

twoway histogram threshold [fw=obs], frac ysca(alt reverse) xtitle("Effective Electoral Threshold", size(medsmall)) ytitle("") ///
xlab(`x2min'(`s2')`x2max') ylab(#4, nogrid labsize(medsmall)) ///
fysize(20) fcolor(black%95) lwidth(vthin) lcolor(white%25) nodraw name(hy, replace) 

twoway histogram polarization [fw=obs], frac xsca(alt reverse) horiz ytitle("Polarization", size(medsmall)) xtitle("") ///
ylab(`x1min'(`s1')`x1max', grid gmax labsize(medsmall)) xlab(#4, nogrid labsize(medsmall))  ///
fxsize(20) fcolor(black%95) lwidth(vthin) lcolor(white%25) nodraw name(hx, replace)

gr combine hx yx hy, hole(3) imargin(zero) scale(1.1) xsize(5.5) ysize(5.5)	