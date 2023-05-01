
*** A plot for the DAME and MEM estimates for the analysis presented in
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

** DAME
xtile group_id = dvprop, nq(10)
margins, dydx(daystoelection) over(group_id) vce(unconditional) saving(temp_dame, replace)
** marginal effects at means
qui sum dvprop
loc cuts="`=r(min)'(`=(r(max)-r(min))/20')`=r(max)'"
margins, dydx(daystoelection) at(dvprop=(`cuts') (mean) _continuous (base) _factor) vce(unconditional) saving(temp_mem, replace)
** combine information
collapse (median) dvprop (count) obs=dvprop, by(group_id)
rename group_id _by1
merge 1:1 _by1 using temp_dame, nogenerate keepusing(_margin _ci_lb _ci_ub)
rename (_margin _ci_lb _ci_ub)(dame lb ub)
append using temp_mem, keep(_at? _margin _ci_lb _ci_ub)
foreach v of varlist dvprop {
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
twoway (line mem dvprop, lpattern(solid)) ///
(rline lbm ubm dvprop, lpattern(dash)) ///
(rspike lb ub dvprop) ///
(scatter dame dvprop [fw=obs], msymbol(o) msize(*.25)), /// 
yline(0, lcolor(red)) ytitle("ME of Days to Election") xtitle("Democratic Vote Share") legend(off)