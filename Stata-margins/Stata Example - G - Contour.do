
*** A contourplot for the marginal effects of polarization for the analysis presented in
*** Golder, Sona. 2006. The Logic of Pre-Electoral Coalition Formation. Columbus: Ohio State University Press.
*** The dataset is available from Matt Golder's website at: http://mattgolder.com/files/interactions/interaction3.zip" (named "interaction3.dta")

clear all

* Load the data and estimate the model
use "interaction3.dta",clear

xtprobit pec c.polarization##c.threshold c.seatshare##c.seatshare incompatibility c.asymmetry##c.seatshare, re i(ident) 
keep if e(sample)

** aim to create a 15 x 15 grid
foreach v of varlist polarization threshold {
	qui sum `v'
	local `v'_s "`r(min)'(`=(r(max)-r(min))/15')`r(max)'" 
}
margins, dydx(polarization) at(polarization=(`polarization_s') threshold=(`threshold_s') (mean) _all) saving(temp_bg,replace)

foreach v of varlist polarization threshold {
	qui sum `v' 
	egen `v'_r = cut(`v'), at(`=r(min)-0.5'(`=(1+r(max)-r(min))/15')`=r(max)+0.5')
}
egen group_id=group(polarization_r threshold_r)
margins, dydx(polarization) over(group_id) at((omean) _all (mean) polarization threshold) saving(temp_me,replace)

collapse (count) obs=pec (mean) polarization threshold, by(group_id)

rename group_id _by1
merge 1:1 _by1 using temp_me, nogenerate keepusing(_margin _pvalue)
gen significant=(_pvalue>=0.975)|(_pvalue<=0.025)
append using temp_bg, keep(_at? _margin)
rename _margin me_est
foreach v of varlist polarization threshold {
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

local colr = "white*.5 yellow*.5 orange*.5 red*.5" /* Color ramp goes from less intense colors to more intense colors */

/* Note that the replication do file uses additional graphical parameters, which leads to different axis and legend labels from this minimal example. */
twoway (contour me_est polarization threshold if me_est!=., ccuts(`ccuts') ccolors(`colr')) ///
(scatter polarization threshold [fw=obs] if significant==0, msymbol(oh) mlcolor(black%95) mlwidth(vthin) msize(*.25)) /// 
(scatter polarization threshold [fw=obs] if significant==1, msymbol(o) mfcolor(black%95) mlwidth(none) msize(*.25)), ///
xtitle(Threshold) ytitle(Polarization) ztitle("") zlabel(`zlabs') ///
legend(off)  clegend(title("Effect Size", size(medsmall) pos(12) justification(right)) width(5) height(25))
