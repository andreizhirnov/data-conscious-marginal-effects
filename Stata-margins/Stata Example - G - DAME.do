
*** A plot for the DAME and MEM estimates for the analysis presented in
*** Golder, Sona. 2006. The Logic of Pre-Electoral Coalition Formation. Columbus: Ohio State University Press.
*** The dataset is available from Matt Golder's website at: http://mattgolder.com/files/interactions/interaction3.zip" (named "interaction3.dta")

clear all

* Load the data and estimate the model
use "interaction3.dta",clear

xtprobit pec c.polarization##c.threshold c.seatshare##c.seatshare incompatibility c.asymmetry##c.seatshare, re i(ident) 
keep if e(sample)

** DAME
xtile group_id = threshold, nq(10)
margins, dydx(polarization) over(group_id) saving(temp_dame, replace) 
** marginal effects at means
qui sum threshold
loc cuts="`=r(min)'(`=(r(max)-r(min))/20')`=r(max)'"
margins, dydx(polarization) at(threshold=(`cuts') (mean) _all) saving("temp_mem", replace)

** combine information
collapse (median) threshold (count) obs=threshold, by(group_id)
rename group_id _by1
merge 1:1 _by1 using temp_dame, nogenerate keepusing(_margin _ci_lb _ci_ub)
rename (_margin _ci_lb _ci_ub)(dame lb ub)
append using temp_mem, keep(_at? _margin _ci_lb _ci_ub)
foreach v of varlist threshold {
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
twoway (line mem threshold, lpattern(solid)) ///
(rline lbm ubm threshold, lpattern(dash)) ///
(rspike lb ub threshold) ///
(scatter dame threshold [fw=obs], msymbol(o) msize(*.25)), /// 
yline(0, lcolor(red)) ytitle("ME of polarization") xtitle("Effective Electoral Threshold") legend(off)

