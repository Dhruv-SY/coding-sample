****************************************************
* 12_greg_ethnicity.do
* Robustness: replace the Murdock (1959) ethnic-homeland fixed
* effects with fixed effects from the GREG dataset (Weidmann,
* Rod & Cederman 2010; digitized Soviet Atlas Narodov Mira),
* following Anderson (2018)'s GREG robustness check. The two
* maps come from independent cartographic traditions, so
* agreement across them shows the results do not hinge on
* Murdock's particular boundary choices.
* Requires data/greg_crosswalk.dta (built by 03_greg_crosswalk.R).
* Reduced form only, 50/100/200km bandwidths.
****************************************************

do "code/00_setup.do"

* attach GREG group at each cluster coordinate
merge m:1 gps_lat gps_lon using "data/greg_crosswalk.dta", keep(master match) nogenerate
quietly count if missing(greg_fe)
di "Unmatched rows (should be 0): " r(N)
drop if missing(greg_fe)

encode greg_fe, gen(greg_id)

* swap Murdock ethnicity FE for GREG FE in the control set
global FE "i.survey_year nightlights_mean population_density_mean agri_suitability_mean malaria_ecology_mean capdist water_wetland_pct elevation_mean traveltime_mean rainfall_mean temperature_mean dist_coast_km i.greg_id"

local bandwidths "50 100 200"
local outcomes  "child_marriage v133 educ_any educ_years_pos v106 v212 v201 v525 multiple_unions not_working domestic_work paid_cash white_collar white_collar_cond"

tempname memhold
tempfile results_bw
postfile `memhold' str40 outcome ///
    double rf_b50 rf_se50 rf_p50 rf_N50 ///
    double rf_b100 rf_se100 rf_p100 rf_N100 ///
    double rf_b200 rf_se200 rf_p200 rf_N200 ///
    using `results_bw', replace

foreach y of local outcomes {

    di "========================================"
    di "Outcome: `y'"
    di "========================================"

    local rf_b50  .
    local rf_se50 .
    local rf_p50  .
    local rf_N50  .
    local rf_b100  .
    local rf_se100 .
    local rf_p100  .
    local rf_N100  .
    local rf_b200  .
    local rf_se200 .
    local rf_p200  .
    local rf_N200  .

    foreach h of local bandwidths {

        preserve
            keep if abs(rd) <= `h'

            if "`y'" == "white_collar_cond" {
                keep if v714 == 1
                keep if !missing(white_collar_cond)
            }
            else if "`y'" == "educ_years_pos" {
                keep if v133 > 0 & !missing(v133)
                keep if !missing(educ_years_pos)
            }
            else {
                keep if !missing(`y')
            }

            count
            if r(N) > 0 {
                quietly reg `y' treated $RDBP $X $FE, vce(cluster cluster_id)
                local rf_b  = _b[treated]
                local rf_se = _se[treated]
                local rf_t  = _b[treated] / _se[treated]
                local rf_p  = 2*ttail(e(df_r), abs(`rf_t'))
                local rf_N  = e(N)
                if `h' == 50 {
                    local rf_b50  = `rf_b'
                    local rf_se50 = `rf_se'
                    local rf_p50  = `rf_p'
                    local rf_N50  = `rf_N'
                }
                else if `h' == 100 {
                    local rf_b100  = `rf_b'
                    local rf_se100 = `rf_se'
                    local rf_p100  = `rf_p'
                    local rf_N100  = `rf_N'
                }
                else if `h' == 200 {
                    local rf_b200  = `rf_b'
                    local rf_se200 = `rf_se'
                    local rf_p200  = `rf_p'
                    local rf_N200  = `rf_N'
                }
            }
        restore
    }

    post `memhold' ("`y'") ///
        (`rf_b50') (`rf_se50') (`rf_p50') (`rf_N50') ///
        (`rf_b100') (`rf_se100') (`rf_p100') (`rf_N100') ///
        (`rf_b200') (`rf_se200') (`rf_p200') (`rf_N200')
}

postclose `memhold'

use `results_bw', clear

replace outcome = "Child marriage (first stage)"   if outcome == "child_marriage"
replace outcome = "Education (years)"              if outcome == "v133"
replace outcome = "Any education"                  if outcome == "educ_any"
replace outcome = "Education (years | schooling>0)" if outcome == "educ_years_pos"
replace outcome = "Education level"                if outcome == "v106"
replace outcome = "Age at first birth"             if outcome == "v212"
replace outcome = "Children ever born"             if outcome == "v201"
replace outcome = "Age at first sex"               if outcome == "v525"
replace outcome = "Multiple unions"                if outcome == "multiple_unions"
replace outcome = "Not working"                    if outcome == "not_working"
replace outcome = "Domestic work"                  if outcome == "domestic_work"
replace outcome = "Paid in cash"                   if outcome == "paid_cash"
replace outcome = "White-collar employment"        if outcome == "white_collar"
replace outcome = "White-collar among workers"     if outcome == "white_collar_cond"

gen rf_star50 = cond(rf_p50<.01,"***",cond(rf_p50<.05,"**",cond(rf_p50<.1,"*","")))
gen rf_star100 = cond(rf_p100<.01,"***",cond(rf_p100<.05,"**",cond(rf_p100<.1,"*","")))
gen rf_star200 = cond(rf_p200<.01,"***",cond(rf_p200<.05,"**",cond(rf_p200<.1,"*","")))

gen rf_coef50  = string(rf_b50,"%9.3f") + rf_star50
gen rf_coef100 = string(rf_b100,"%9.3f") + rf_star100
gen rf_coef200 = string(rf_b200,"%9.3f") + rf_star200

gen rf_se50_s  = "(" + string(rf_se50,"%9.3f") + ")"
gen rf_se100_s = "(" + string(rf_se100,"%9.3f") + ")"
gen rf_se200_s = "(" + string(rf_se200,"%9.3f") + ")"

gen rf_N50_s   = string(rf_N50,"%9.0f")
gen rf_N100_s  = string(rf_N100,"%9.0f")
gen rf_N200_s  = string(rf_N200,"%9.0f")

export delimited using "results/results_greg_reducedform.csv", replace

file open tex_rf using "results/results_greg_reducedform.tex", write replace

file write tex_rf "\begin{table}[htbp]" _n
file write tex_rf "\centering" _n
file write tex_rf "\caption{Reduced-Form Effect of Treatment: GREG Ethnic Fixed Effects (Atlas Narodov Mira)}" _n
file write tex_rf "\label{tab:rf_bw_greg}" _n
file write tex_rf "\begin{tabular}{lccc}" _n
file write tex_rf "\hline\hline" _n
file write tex_rf " & \multicolumn{3}{c}{Bandwidth (km)} \\" _n
file write tex_rf "\cline{2-4}" _n
file write tex_rf "Outcome & 50 & 100 & 200 \\" _n
file write tex_rf "\hline" _n

count
local N = r(N)

forvalues i = 1/`N' {
    local out    = outcome[`i']
    local rfc50  = rf_coef50[`i']
    local rfc100 = rf_coef100[`i']
    local rfc200 = rf_coef200[`i']
    local rfs50  = rf_se50_s[`i']
    local rfs100 = rf_se100_s[`i']
    local rfs200 = rf_se200_s[`i']
    local rfn50  = rf_N50_s[`i']
    local rfn100 = rf_N100_s[`i']
    local rfn200 = rf_N200_s[`i']

    file write tex_rf "`out'" " & " "`rfc50'" " & " "`rfc100'" " & " "`rfc200'" " \\" _n
    file write tex_rf " " " & " "`rfs50'" " & " "`rfs100'" " & " "`rfs200'" " \\" _n
    file write tex_rf "Obs." " & " "`rfn50'" " & " "`rfn100'" " & " "`rfn200'" " \\" _n
    file write tex_rf "\addlinespace" _n
}

file write tex_rf "\hline" _n
file write tex_rf "\multicolumn{4}{p{0.9\linewidth}}{\footnotesize Notes: Identical to the main reduced-form specification except that the Murdock (1959) ethnic-homeland fixed effects are replaced with fixed effects from the Geo-Referencing of Ethnic Groups dataset (GREG; Weidmann, Rod and Cederman 2010), the digitized Soviet Atlas Narodov Mira, following the analogous robustness check in Anderson (2018). Each survey cluster is assigned the primary GREG group at its coordinate (46 groups in the estimation sample, versus 102 Murdock homelands). Standard errors clustered at the DHS-cluster level. Stars: * \$p<0.10\$, ** \$p<0.05\$, *** \$p<0.01\$.}" _n
file write tex_rf "\end{tabular}" _n
file write tex_rf "\end{table}" _n
file close tex_rf
