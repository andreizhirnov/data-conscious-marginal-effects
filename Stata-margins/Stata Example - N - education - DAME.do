
*** A plot for the DAME and MEM estimates for the analysis presented in
*** Nagler, Jonathan. 1991. "The Effect of Registration Laws and Education on U.S. Voter Turnout." American Political Science Review 85(4): 1393–1405.
*** The dataset was made public by William D. Berry, Jacqueline H. R. DeMeritt, and Justin Esarey as part of the replication materials for their article entitled:
*** "Testing for Interaction in Binary Logit and Probit Models: Is a Product Term Essential?" and can be downloaded from: https://jdemeritt.weebly.com/uploads/2/2/7/7/22771764/bde.zip (named "scobit.dta")

clear all

* Load the data and estimate the model
use "scobit.dta",clear
drop if newvote==-1
probit newvote c.closing##c.neweduc##c.neweduc c.age##c.age ib(freq).south ib(freq).gov

keep if e(sample)

** DAME
xtile group_id = closing, nq(10)
margins, dydx(neweduc) over(group_id) saving(temp_dame, replace) 
** marginal effects at means
qui sum closing
loc cuts="`=r(min)'(`=(r(max)-r(min))/20')`=r(max)'"
margins, dydx(neweduc) at(closing=(`cuts') (mean) _continuous (base) _factor) saving(temp_mem, replace)
** combine information
collapse (median) closing (count) obs=closing, by(group_id)
rename group_id _by1
merge 1:1 _by1 using temp_dame, nogenerate keepusing(_margin _ci_lb _ci_ub)
rename (_margin _ci_lb _ci_ub)(dame lb ub)
append using temp_mem, keep(_at? _margin _ci_lb _ci_ub)
foreach v of varlist closing {
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
twoway (line mem closing, lpattern(solid)) ///
(rline lbm ubm closing, lpattern(dash)) ///
(rspike lb ub closing) ///
(scatter dame closing [fw=obs], msymbol(o) msize(*.25)), /// 
yline(0, lcolor(red)) ytitle("ME of Education") xtitle("Closing Date") legend(off) 