
*** A contourplot for the marginal effects of registration closing date for the analysis presented in
*** Nagler, Jonathan. 1991. "The Effect of Registration Laws and Education on U.S. Voter Turnout." American Political Science Review 85(4): 1393–1405.
*** The dataset was made public by William D. Berry, Jacqueline H. R. DeMeritt, and Justin Esarey as part of the replication materials for their article entitled:
*** "Testing for Interaction in Binary Logit and Probit Models: Is a Product Term Essential?" and can be downloaded from: https://jdemeritt.weebly.com/uploads/2/2/7/7/22771764/bde.zip (named "scobit.dta")

clear all

* Load the data and estimate the model
use "scobit.dta",clear
drop if newvote==-1
probit newvote c.closing##c.neweduc##c.neweduc c.age##c.age ib(freq).south ib(freq).gov

keep if e(sample)

** points for the background
foreach v of varlist closing neweduc {
	qui sum `v'
	local `v'_s "(`r(min)'(`=(r(max)-r(min))/15')`r(max)')"
}
margins, dydx(closing) at(closing=`closing_s' neweduc=`neweduc_s' (mean) _continuous (base) _factor) saving(temp_bg,replace)

** points for the foreground
foreach v of varlist closing neweduc {
	qui sum `v' 
	egen `v'_r = cut(`v'), at(`=r(min)-0.5'(`=(1+r(max)-r(min))/15')`=r(max)+0.5')
}
egen group_id=group(closing_r neweduc_r)
margins, dydx(closing) over(group_id) at((omean) _continuous (base) _factor (mean) closing (median) neweduc) saving(temp_me,replace)

collapse (count) obs=closing (mean) closing neweduc, by(group_id)

rename group_id _by1
merge 1:1 _by1 using temp_me, nogenerate keepusing(_margin _pvalue)
gen significant=(_pvalue>=0.975)|(_pvalue<=0.025)
append using temp_bg, keep(_at? _margin)
rename _margin me_est
foreach v of varlist closing neweduc {
	foreach a of varlist _at? {
		local l: variable label `a'
		if "`l'"=="`v'" {
			replace `v'=`a' if mi(`v')
			break
		}
	}
}

* Additional observations to anchor the marker sizes
gen counter=_n
qui sum counter
loc coreobs=r(max)
set obs `=`coreobs'+4'

qui sum obs
replace significant=1 in `=`coreobs'+1'/`=`coreobs'+2'
replace significant=0 in `=`coreobs'+3'/`=`coreobs'+4'
replace obs=r(min) in `=`coreobs'+1'/`=`coreobs'+4'
replace obs=r(max) in `=`coreobs'+2'/`=`coreobs'+3'

* Break the ME values into steps
qui sum me_est 
matrix mimx = (r(min), r(max))*(4,3,2,1,0\0,1,2,3,4)/4
local ccuts = "`=mimx[1,2]' `=mimx[1,3]' `=mimx[1,4]'"
local zlabs = "`=mimx[1,1]' `ccuts' `=mimx[1,5]'"

local colr = "red*.5 orange*.5 yellow*.5 white*.5" /* Color ramp from more intense to less intense colors */

/* Note that the replication do file uses additional graphical parameters, which leads to different axis and legend labels from this minimal example. */
twoway (contour me_est closing neweduc if !mi(me_est), ccuts(`ccuts') ccolors(`colr')) ///
(scatter closing neweduc [fw=obs] if significant==0, msymbol(oh) mlcolor(black%95) mlwidth(vthin) msize(*.25)) /// 
(scatter closing neweduc [fw=obs] if significant==1, msymbol(o) mfcolor(black%95) mlwidth(none) msize(*.25)), ///
xtitle(Education) ytitle(Closing Date) ztitle("") zlabel(`zlabs') ///
legend(off)  clegend(title("Effect Size", size(medsmall) pos(12) justification(right)) width(5) height(25)) 