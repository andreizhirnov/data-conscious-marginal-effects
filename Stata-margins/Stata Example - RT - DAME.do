
*** A plot for the DAME and MEM estimates for the analysis presented in
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

** DAME
xtile group_id = l_polity2, nq(4)
margins, dydx(l_l_flows) over(group_id) saving(temp_dame, replace) 
** marginal effects at means
qui sum l_polity2
loc cuts="`=r(min)'(1)`=r(max)'"
margins, dydx(l_l_flows) at(l_polity2=(`cuts')) atmeans saving(temp_mem, replace)
** combine information
collapse (median) l_polity2 (count) obs=l_polity2, by(group_id)
rename group_id _by1
merge 1:1 _by1 using temp_dame, nogenerate keepusing(_margin _ci_lb _ci_ub)
rename (_margin _ci_lb _ci_ub)(dame lb ub)
append using temp_mem, keep(_at? _margin _ci_lb _ci_ub) 
foreach v of varlist l_polity2 {
	foreach a of varlist _at? {
		local l: variable label `a'
		if "`l'"=="`v'" {
			replace `v'=`a' if mi(`v')
			break
		}
	}
}
rename (_margin _ci_lb _ci_ub)(mem lbm ubm)

* Plot the DAME and MEM estimates
/* Note that the replication do file uses additional graphical parameters, which leads to different axis and legend labels from this minimal example. */
twoway (line mem l_polity2, lpattern(solid)) ///
(rline lbm ubm l_polity2, lpattern(dash)) ///
(rspike lb ub l_polity2) ///
(scatter dame l_polity2 [fw=obs], msymbol(o) msize(*.25)), /// 
yline(0, lcolor(red)) ytitle("ME of ln(FDI flows)") xtitle("Polity 2") legend(off)