*==============================================================*
* 00_setup.do
* Shared preamble for all Stata replication scripts.
* Do NOT run this file directly — it is called by each script.
*==============================================================*

*---------------------------------------------------*
* 0) LOAD DATA
*---------------------------------------------------*
use "data/20yr_threshold_12pairs.dta", clear

* A handful of women (6 as of the last build) are exact full-row duplicates,
* almost certainly from a boundary-edge case in the grid/ethnic-homeland
* spatial join (a GPS point landing exactly on the line between two adjacent
* polygons can get matched to both). Drop the repeat copies.
duplicates drop

*---------------------------------------------------*
* 1) SAMPLE RESTRICTIONS + IDS
*---------------------------------------------------*
drop if v012 < 18 & v511 == .

* Fully-exposed cohorts only. fully_exposed marks women whose entire at-risk
* window -- from the control country's minimum legal age to the treated
* country's -- fell inside the years the legal gap was in force. The pair
* files apply only an overlap rule, which lets in women who spent part of
* their teenage years under identical law on both sides; those dilute the
* treatment. The flag is built once in build/02_pool_pairs.R so that this
* file and the four R scripts share one definition.
capture confirm variable fully_exposed
if _rc {
    di as error "fully_exposed not found -- rebuild the analysis dataset with"
    di as error "build/02_pool_pairs.R, which constructs the flag."
    exit 111
}
keep if fully_exposed == 1

capture drop bp eth_id multiple_unions not_working domestic_work paid_cash ///
    white_collar white_collar_cond educ_any educ_years_pos

* v021 (DHS PSU number) is only unique WITHIN a single country-survey round;
* it restarts from 1 every round, so raw v021 collides across different
* countries and different rounds of the same country. Build a cluster ID
* that's actually unique to one physical PSU before clustering SEs on it.
capture drop cluster_id
egen cluster_id = group(country_name survey_year v021)

* A single DHS cluster can fall within the bandwidth of more than one
* border pair (e.g. eastern DRC is near Rwanda, Burundi, Uganda, and
* Zambia all at once), so the same women are duplicated across pair
* files. Keep only each cluster's comparison against its NEAREST border
* (smallest |rd|); drop its appearances in the other, more distant pairs.
capture drop abs_rd min_abs_rd
gen abs_rd = abs(rd)
egen min_abs_rd = min(abs_rd), by(cluster_id)
drop if abs_rd > min_abs_rd
drop abs_rd min_abs_rd

encode country_pair, gen(bp)
encode eth_fe, gen(eth_id)

drop if missing(child_marriage, treated, rd, bp, survey_year, v021)

*---------------------------------------------------*
* 2) RD SHAPE + FE + CONTROLS
*---------------------------------------------------*
global RDBP "i.bp c.rd#i.bp c.rd#i.bp#i.treated"
global FE   "i.survey_year nightlights_mean population_density_mean agri_suitability_mean malaria_ecology_mean capdist water_wetland_pct elevation_mean traveltime_mean rainfall_mean temperature_mean dist_coast_km i.eth_id"
global X    "c.v012 i.v025"

* Region FE (Southern/Western/Eastern/Central Africa dummies) is excluded:
* strictly redundant with border-pair FE (i.bp), since each of the 12 pairs
* sits entirely within a single region. Every other pixel-level geographic
* control from Anderson (2018)'s spec is kept.

*---------------------------------------------------*
* 3) CLEAN KEY VARIABLES
*---------------------------------------------------*
replace v133 = . if inlist(v133, 97, 98, 99)
replace v106 = . if inlist(v106, 9)
replace v525 = . if inlist(v525, 0, 96, 97, 98, 99)
replace v741 = . if v741 == 9
replace v717 = . if inlist(v717, 10, 96, 98, 99)

*---------------------------------------------------*
* 4) CONSTRUCT OUTCOMES
*---------------------------------------------------*

* Extensive margin: any education
gen educ_any = .
replace educ_any = 1 if v133 > 0 & !missing(v133)
replace educ_any = 0 if v133 == 0
label var educ_any "Any education"

* Intensive margin: years conditional on having any schooling
gen educ_years_pos = .
replace educ_years_pos = v133 if v133 > 0 & !missing(v133)
label var educ_years_pos "Education years | schooling > 0"

* multiple_unions: ever been in more than one union
* v503 = 1 (one union only)   -> 0
* v503 = 2+ (more than one)   -> 1
* v503 = 9 (don't know)       -> .
* v503 = . (never in union)   -> 0
gen multiple_unions = 0
replace multiple_unions = 1 if v503 >= 2 & !missing(v503)
replace multiple_unions = . if v503 == 9
label var multiple_unions "More than one union"

gen not_working = .
replace not_working = 1 if v714 == 0
replace not_working = 0 if v714 == 1
label var not_working "Not working"

* Domestic work, conditional on currently working and on occupation being known.
* Women who are working but whose occupation is missing stay missing: coding them
* 0 would assert they are not domestic workers when v717 does not say so. The
* earlier version did exactly that, because in Stata missing != 6 evaluates true.
gen domestic_work = .
replace domestic_work = 1 if v714 == 1 & v717 == 6
replace domestic_work = 0 if v714 == 1 & !missing(v717) & v717 != 6
label var domestic_work "Works in domestic jobs"

gen paid_cash = .
replace paid_cash = 1 if v714 == 1 & inlist(v741,1,2)
replace paid_cash = 0 if v714 == 1 & inlist(v741,0,3)
label var paid_cash "Currently working and paid in cash (cash only or cash and in-kind)"

* white_collar extensive margin: anyone whose recorded occupation is
* professional, technical or managerial (v717 in 1,2) counts as white collar,
* whether or not she is working right now -- v717 is still reported for women
* who worked recently but are not currently employed. Order matters: the
* inlist() line comes LAST so it overrides the v714==0 zero.
* Working women whose occupation is missing stay missing.
gen white_collar = .
replace white_collar = 0 if v714 == 0
replace white_collar = 0 if v714 == 1 & !missing(v717) & !inlist(v717,1,2)
replace white_collar = 1 if inlist(v717,1,2)
label var white_collar "White-collar employment"

gen white_collar_cond = .
replace white_collar_cond = 1 if v714 == 1 & inlist(v717,1,2)
replace white_collar_cond = 0 if v714 == 1 & !missing(v717) & !inlist(v717,1,2)
label var white_collar_cond "White-collar among workers"

