

*** A contourplot for the marginal effects of the proximity of the next election for the analysis presented in
*** Arceneaux, Kevin, Martin Johnson, Rene Lindstädt, and Ryan J. Vander Wielen. 2016. "The Influence of News Media on Political Elites: Investigating Strategic Responsiveness in Congress." American Journal of Political Science 60(1): 5–29.
*** The dataset is part of the published replication materials and can be downloaded from: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/27597" (named "FoxNews_Master.dta")

clear all

* Load the data and estimate the model
use "FoxNews_Master.dta",clear
gen dvprop=dv/100
gen vtype_char = "OtherProc" if !mi(RegPass)
foreach v of varlist RegPass Amend OtherPass ProPart Susp {
	replace vtype_char="`v'" if `v'
}
encode vtype_char,gen(votetype)

logit PartyVote c.daystoelection##c.daystoelection##c.daystoelection##c.dvprop seniorit spendgap_lag spendgap distpart_lag ib(freq).votetype ib(freq).Retirement ib(freq).qualchal_lag ib(freq).qualchal if PresencePartyUnity==1 & Republican==1 & FoxNews==1, cluster(dist2)

keep if e(sample)

** points for the background
foreach v of varlist daystoelection dvprop {
	qui sum `v'
	local `v'_s "(`r(min)'(`=(r(max)-r(min))/15')`r(max)')"
}
margins, dydx(daystoelection) at(daystoelection=`daystoelection_s' dvprop=`dvprop_s' (mean) _continuous (base) _factor) saving(temp_bg,replace)

** points for the foreground
foreach v of varlist daystoelection dvprop {
	qui sum `v' 
	egen `v'_r = cut(`v'), at(`=r(min)-0.5'(`=(1+r(max)-r(min))/15')`=r(max)+0.5')
}
egen group_id=group(daystoelection_r dvprop_r)
margins, dydx(daystoelection) over(group_id) at((omean) _continuous (base) _factor (mean) daystoelection dvprop) saving(temp_me,replace)

collapse (count) obs=daystoelection (mean) daystoelection dvprop, by(group_id)

rename group_id _by1
merge 1:1 _by1 using temp_me, nogenerate keepusing(_margin _pvalue)
gen significant=(_pvalue>=0.975)|(_pvalue<=0.025)
append using temp_bg, keep(_at? _margin)
rename _margin me_est
foreach v of varlist daystoelection dvprop {
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
matrix mimx = (r(min), r(max))*(5,4,3,2,1,0\0,1,2,3,4,5)/5
local ccuts = "`=mimx[1,2]' `=mimx[1,3]' `=mimx[1,4]' `=mimx[1,5]'"
local zlabs = "`=mimx[1,1]' `ccuts' `=mimx[1,6]'"

loc colr= "navy*.5 ltblue*.5 white*.5 orange*.5 red*.5"
/* Color ramp has intense colors at both ends */

/* Note that the replication do file uses additional graphical parameters, which leads to different axis and legend labels from this minimal example. */
twoway (contour me_est daystoelection dvprop if me_est!=., ccuts(`ccuts') ccolors(`colr')) ///
(scatter daystoelection dvprop [fw=obs] if significant==0, msymbol(oh) mlcolor(black%95) mlwidth(vthin) msize(*.25)) /// 
(scatter daystoelection dvprop [fw=obs] if significant==1, msymbol(o) mfcolor(black%95) mlwidth(none) msize(*.25)), ///
xtitle(Democratic Vote Share) ytitle(Days to Election) ztitle("") zlabel(`zlabs') ///
legend(off)  clegend(title(`"Effect Size"', size(medsmall) pos(12) justification(right)) width(5) height(25))  
