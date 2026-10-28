****************************************************
* 10_robustness_leave_one_pair_out.do
* Reduced-form-only robustness: drop each of the 12 border
* pairs one at a time and re-estimate the reduced-form effect
* of treatment on every main outcome (100km bandwidth).
* Generalizes 05_robustness_no_namibia.do / 06_robustness_no_niger.do
* to all 12 pairs systematically.
****************************************************

do "code/00_setup.do"

local bandwidths "100"
local outcomes  "child_marriage v133 educ_any educ_years_pos v106 v212 v201 v525 multiple_unions not_working domestic_work paid_cash white_collar white_collar_cond"

* short labels for column headers
local pair_BenBurk    "Benin-Burkina Faso"
local pair_BenNig     "Benin-Niger"
local pair_BurTan     "Burundi-Tanzania"
local pair_CivBurk    "Cote d'Ivoire-Burkina Faso"
local pair_CivMali    "Cote d'Ivoire-Mali"
local pair_GhaBurk    "Ghana-Burkina Faso"
local pair_NamAng     "Namibia-Angola"
local pair_NamSaf     "Namibia-South Africa"
local pair_NamZam     "Namibia-Zambia"
local pair_NamZim     "Namibia-Zimbabwe"
local pair_RwaTan     "Rwanda-Tanzania"
local pair_UgaTan     "Uganda-Tanzania"

local pair_codes "BenBurk BenNig BurTan CivBurk CivMali GhaBurk NamAng NamSaf NamZam NamZim RwaTan UgaTan"

*---------------------------------------------------*
* STORE REDUCED-FORM COEFFICIENTS: baseline (all pairs) + each exclusion,
* at each bandwidth
*---------------------------------------------------*
tempname memhold
tempfile results_loo
postfile `memhold' str40 outcome str40 excluded int bw ///
    double rf_b rf_se rf_p rf_N ///
    using `results_loo', replace

foreach h of local bandwidths {
    foreach y of local outcomes {

        * ---- baseline: all 12 pairs, no exclusion ----
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
                local b  = _b[treated]
                local se = _se[treated]
                local t  = _b[treated] / _se[treated]
                local p  = 2*ttail(e(df_r), abs(`t'))
                local N  = e(N)
                post `memhold' ("`y'") ("All Pairs") (`h') (`b') (`se') (`p') (`N')
            }
        restore

        * ---- drop each pair in turn ----
        foreach code of local pair_codes {
            local pname "`pair_`code''"

            preserve
                keep if abs(rd) <= `h'
                drop if country_pair == "`pname'"

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
                    local b  = _b[treated]
                    local se = _se[treated]
                    local t  = _b[treated] / _se[treated]
                    local p  = 2*ttail(e(df_r), abs(`t'))
                    local N  = e(N)
                    post `memhold' ("`y'") ("No `pname'") (`h') (`b') (`se') (`p') (`N')
                }
            restore
        }
    }
}

postclose `memhold'

*---------------------------------------------------*
* BETTER OUTCOME LABELS + EXPORT WIDE CSV (outcome x excluded-pair x bw)
*---------------------------------------------------*
use `results_loo', clear

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

gen rf_star = cond(rf_p<.01,"***",cond(rf_p<.05,"**",cond(rf_p<.1,"*","")))
gen rf_coef = string(rf_b,"%9.3f") + rf_star

export delimited outcome excluded bw rf_b rf_se rf_p rf_N using "results/leave_one_pair_out_reducedform.csv", replace

*---------------------------------------------------*
* LATEX TABLES: one per bandwidth, outcome (rows) x excluded pair (columns)
*---------------------------------------------------*
levelsof outcome, local(outs)
local col_order `" "All Pairs" "No Benin-Burkina Faso" "No Benin-Niger" "No Burundi-Tanzania" "No Cote d'Ivoire-Burkina Faso" "No Cote d'Ivoire-Mali" "No Ghana-Burkina Faso" "No Namibia-Angola" "No Namibia-South Africa" "No Namibia-Zambia" "No Namibia-Zimbabwe" "No Rwanda-Tanzania" "No Uganda-Tanzania" "'

foreach h of local bandwidths {
    file open tex using "results/leave_one_pair_out_reducedform.tex", write replace
    file write tex "\begin{table}[htbp]" _n
    file write tex "\centering" _n
    file write tex "\tiny" _n
    file write tex "\caption{Leave-One-Pair-Out Robustness: Reduced-Form Effect of Treatment (`h'km Bandwidth)}" _n
    file write tex "\label{tab:leave_one_pair_out_bw`h'}" _n
    file write tex "\begin{tabular}{l*{13}{c}}" _n
    file write tex "\hline\hline" _n
    file write tex "Outcome & All Pairs & No BEN-BFA & No BEN-NER & No BDI-TZA & No CIV-BFA & No CIV-MLI & No GHA-BFA & No NAM-AGO & No NAM-ZAF & No NAM-ZMB & No NAM-ZWE & No RWA-TZA & No UGA-TZA \\" _n
    file write tex "\hline" _n

    foreach o of local outs {
        local row `"`o'"'
        foreach c of local col_order {
            count if outcome == `"`o'"' & excluded == `"`c'"' & bw == `h'
            if r(N) > 0 {
                quietly levelsof rf_coef if outcome == `"`o'"' & excluded == `"`c'"' & bw == `h', local(cellval) clean
                local row `"`row' & `cellval'"'
            }
            else {
                local row `"`row' & --"'
            }
        }
        file write tex `"`row' \\"' _n
    }

    file write tex "\hline" _n
    file write tex "\multicolumn{14}{p{0.95\linewidth}}{\footnotesize Notes: Each cell reports the reduced-form coefficient on treated from a separate regression at the `h'km bandwidth, with the column's named border pair excluded from the sample (``All Pairs'' is the unrestricted baseline). All regressions include border-pair fixed effects, border-pair-specific slopes and kinks, geographic controls, ethnicity fixed effects, survey-year fixed effects, and individual controls, with standard errors clustered at the DHS-cluster level. Full coefficients, standard errors, and observation counts are in leave\_one\_pair\_out\_reducedform.csv. Stars: * \$p<0.10\$, ** \$p<0.05\$, *** \$p<0.01\$.}" _n
    file write tex "\end{tabular}" _n
    file write tex "\end{table}" _n
    file close tex
}
