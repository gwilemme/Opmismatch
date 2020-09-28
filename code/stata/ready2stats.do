/***** Employment stocks and wages *****/

/**** core program ****/
capture program drop build_E
program define build_E
	args i

	keep emp1 occ1 wag1
	gen N`i' = 1 /* number of observations */
	gen Nw`i'  = !(missing(wag1)) /* number of observed wages */
	collapse (sum) wag1 N`i' Nw`i', by(emp1 occ1)
	gen w`i' = wag1 / Nw`i'
	drop wag1

	egen tot = sum(N`i')
	gen share`i' = N`i'/tot /*share of workers by occupation */
	drop tot
end

/********** do not consider jan17, start with feb17 ****/ 
/*** initialize ***/
use ".....\data\ready\cpsfeb17.dta", clear
build_E "1"
	
tempfile mycopy
save `mycopy'


/*** loop ***/
local mths "feb17 mar17 apr17 may17 jun17 jul17 aug17 sep17 oct17 nov17 dec17 jan18 feb18 mar18 apr18 may18 jun18 jul18 aug18 sep18 oct18 nov18 dec18"
quietly forval i=2/23 {
	local mth `:word `i' of `mths''

use ".....\data\ready\cps`mth'.dta", clear
build_E `i'

merge 1:1 emp1 occ1 using `mycopy', nogenerate
saveold `mycopy', replace
}

gen share = 0
gen w = 0
gen N = 0
gen Nw = 0

quietly forval i=1/23 {
	replace share = share + share`i'
	replace w = w + w`i'
	replace N = N + N`i'
	replace Nw = Nw + Nw`i'
}
replace share = share/23
replace w = w/23
keep emp1 occ1 share w N Nw

outsheet using ".....\data\moments\stats.csv", comma replace


