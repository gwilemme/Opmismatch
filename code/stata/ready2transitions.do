/***** Transition rates *****/

/**** core program ****/
capture program drop build_rates
program define build_rates
	args i

	keep emp0 emp1 occ0 occ1 wag0 wag1 sameemployer
	replace occ1 = occ0 if (emp0==1) & (emp1==0) /* correction for consistent E to U transitions (diagonal transition matrix) */
	replace occ1 = occ0 if (emp0==0) & (emp1==0) /* correction for consistent U to U transitions (diagonal transition matrix) */
	drop if (occ0==-1) | (occ1==-1)
	
	gen state0 = ""
	replace state0 = "U" + string(occ0) if (emp0==0)
	replace state0 = "E" + string(occ0) if (emp0==1)
	gen state1 = ""
	replace state1 = "U" + string(occ1) if (emp1==0)
	replace state1 = "Enew" + string(occ1) if (emp1==1) 
	replace state1 = "Esame" + string(occ1) if (emp1==1) & (sameemployer==1) & (occ1==occ0) & (emp0==1)

	contract state0 state1, zero /** count the number of observations per combination state0-state1 **/
	bysort state0: egen totfreq = sum(_freq)
	gen rate`i' = _freq / totfreq /** transition rate **/
	rename _freq N`i' /** number of observations **/
	drop totfreq

end


/*** initialize ***/
use ".....\data\ready\cpsfeb17.dta", clear
drop emp0 occ0 wag0

/* merge with previous month */
merge 1:1 id using ".....\data\ready\cpsjan17.dta", keepusing(emp0 occ0 wag0)
keep if _merge == 3
drop _merge 

build_rates 1

tempfile mycopy
save `mycopy'


/*** loop ***/
local mth1 "feb17 mar17 apr17 may17 jun17 jul17 aug17 sep17 oct17 nov17 dec17 jan18 feb18 mar18 apr18 may18 jun18 jul18 aug18 sep18 oct18 nov18 dec18"
local mth0 "jan17 feb17 mar17 apr17 may17 jun17 jul17 aug17 sep17 oct17 nov17 dec17 jan18 feb18 mar18 apr18 may18 jun18 jul18 aug18 sep18 oct18 nov18"

quietly forval i=2/23 {
	local m1 `:word `i' of `mth1''
	local m0 `:word `i' of `mth0''
	
	use ".....\data\ready\cps`m1'.dta", clear
	drop emp0 occ0 wag0
	merge 1:1 id using ".....\data\ready\cps`m0'.dta", keepusing(emp0 occ0 wag0)
	keep if _merge == 3
	drop _merge 

	build_rates `i'

	merge 1:1 state0 state1 using `mycopy', nogenerate
	saveold `mycopy', replace
}


gen rate = 0
gen Nrate = 0
gen N = 0


quietly forval i=1/23 {
	replace rate = rate + rate`i' if !missing(rate`i')
	replace Nrate = Nrate + 1 if !missing(rate`i')
	replace N = N + N`i'  if !missing(N`i')

}
replace rate = rate/Nrate
keep state0 state1 rate N


outsheet using "C:\Users\LaWile\Dropbox\Economics\These\paper1\data\moments\rates.csv", comma replace


