
*** A heatmap for the marginal effects of logged FDI flows for the analysis presented in
*** Robertson, Graeme B., and Emmanuel Teitelbaum. 2011. "Foreign Direct Investment, Regime Type, and Labor Protest in Developing Countries." American Journal of Political Science 55(3): 665–677.
*** The dataset is part of the replication materials and can be downloaded from Emmanuel Teitelbaum's website: https://home.gwu.edu/~ejt/pages/Data_files/Robertson%20Teitelbaum%202011.dta"

clear all

* Load the data and estimate the model
use "Robertson Teitelbaum 2011.dta",clear

tsset country year 
gen l_l_flows=L.l_flows
gen l_polity2=L.polity2
gen l_dispute=L.dispute
gen l_demflows=l_l_flows*l_polity2
 
xtnbreg dispute c.l_l_flows##c.l_polity2 l_dispute open_penn l_gdp_pc_penn gdp_grth inflation_1 urban xratchg l_pop time, re
keep if e(sample)

*** points for the background
foreach v of varlist l_l_flows l_polity2 {
	qui sum `v'
	local `v'_s "(`r(min)'(`=(r(max)-r(min))/15')`r(max)')"
}
margins, dydx(l_l_flows) at(l_l_flows=`l_l_flows_s' l_polity2=`l_polity2_s' (mean) _all) saving(temp_bg,replace)

*** points for the foreground
foreach v of varlist l_l_flows l_polity2 {
	qui sum `v' 
	egen `v'_r = cut(`v'), at(`=r(min)-0.5'(`=(1+r(max)-r(min))/10')`=r(max)+0.5')
}
egen group_id=group(l_l_flows_r l_polity2_r)
margins, dydx(l_l_flows) over(group_id) at((omean) _all (mean) l_l_flows l_polity2) saving(temp_me,replace)

collapse (count) obs=l_l_flows (mean) l_l_flows l_polity2, by(group_id)

rename group_id _by1
merge 1:1 _by1 using temp_me, nogenerate keepusing(_margin _pvalue)
gen significant=(_pvalue>=0.975)|(_pvalue<=0.025)
append using temp_bg, keep(_at? _margin)
rename _margin me_est
foreach v of varlist l_l_flows l_polity2 {
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

* Color ramp from less intense to more intense colors
loc scolr="yellow*.25"
loc ecolr="red*.95"

* Labs: numlist specification for 5 equally spaced values for each variable
local nra = 5-1
foreach v of varlist l_l_flows l_polity2 {
qui sum `v'
local s = (r(max)-r(min))/`nra'
local r = 10^(floor(ln(`s')/ln(10))-1)
local s = round(`s',`r') 
local from = round(r(min),`r')
local to = `from' + `nra'*`s'
local `v'_r = "`from'(`s')`to'"
} 

/* Note that the replication do file uses additional graphical parameters, which leads to different axis and legend labels from this minimal example. */
twoway (contour me_est l_l_flows l_polity2 if me_est!=., levels(100) crule(linear) scolor(`scolr') ecolor(`ecolr') zlab(#5, labsize(medsmall))) ///
(scatter l_l_flows l_polity2 [fw=obs] if significant==0, msymbol(oh) mlcolor(black%95) mlwidth(vthin) msize(*.25)) ///
(scatter l_l_flows l_polity2 [fw=obs] if significant==1, msymbol(o) mfcolor(black%95) mlwidth(none) msize(*.25)), ///
xsca(alt) ysca(alt) xtitle("") ytitle("") ztitle("") ///
ylab(`l_l_flows_r', grid gmax labsize(medsmall)) xlab(`l_polity2_r', labsize(medsmall) grid gmax) /// 
legend(off) clegend(title("Effect Size", size(medsmall) pos(12) justification(right)) ring(0) width(5) height(25)) nodraw name(yx, replace)

twoway histogram l_polity2 [fw=obs], frac ysca(alt reverse) xtitle("Polity 2", size(medsmall)) ytitle("") ///
xlab(`l_polity2_r') ylab(#4, nogrid labsize(medsmall)) ///
fysize(20) fcolor(black%95) lwidth(vthin) lcolor(white%25) nodraw name(hy, replace) 

twoway histogram l_l_flows  [fw=obs], frac xsca(alt reverse) horiz ytitle("ln(FDI flows)", size(medsmall)) xtitle("") ///
ylab(`l_l_flows_r', grid gmax labsize(medsmall)) xlab(#4, nogrid labsize(medsmall))  ///
fxsize(20) fcolor(black%95) lwidth(vthin) lcolor(white%25) nodraw name(hx, replace)

gr combine hx yx hy, hole(3) imargin(zero) scale(1.1) xsize(5.5) ysize(5.5)
