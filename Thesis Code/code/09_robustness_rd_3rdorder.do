****************************************************
* 09_robustness_rd_3rdorder.do
* Robustness: 3rd-order RD polynomial in distance to border
* (matches Anderson (2018) robustness Table A7, column 3)
****************************************************

do "code/00_setup.do"

* Add quadratic and cubic terms in rd, interacted with border-pair FE
* and border-pair x treated, on top of the linear terms already in $RDBP
global RDBP "i.bp c.rd#i.bp c.rd#i.bp#i.treated c.rd#c.rd#i.bp c.rd#c.rd#i.bp#i.treated c.rd#c.rd#c.rd#i.bp c.rd#c.rd#c.rd#i.bp#i.treated"

*---------------------------------------------------*
* 5) BANDWIDTHS + OUTCOMES
*---------------------------------------------------*
local bandwidths "50 100 200"
local outcomes  "v133 educ_any educ_years_pos v106 v212 v201 v525 multiple_unions not_working domestic_work paid_cash white_collar white_collar_cond"

*---------------------------------------------------*
* 6) STORE REDUCED FORM + OBSERVATIONS
*---------------------------------------------------*
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
            if r(N) == 0 {
                restore
                continue
            }

            *-----------------------------*
            * Reduced form
            *-----------------------------*
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
        restore
    }

    post `memhold' ("`y'") ///
        (`rf_b50') (`rf_se50') (`rf_p50') (`rf_N50') ///
        (`rf_b100') (`rf_se100') (`rf_p100') (`rf_N100') ///
        (`rf_b200') (`rf_se200') (`rf_p200') (`rf_N200')
}

postclose `memhold'

*---------------------------------------------------*
* 7) OPEN RESULTS
*---------------------------------------------------*
use `results_bw', clear

*---------------------------------------------------*
* 8) BETTER OUTCOME LABELS
*---------------------------------------------------*
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

*---------------------------------------------------*
* 9) ADD STARS
*---------------------------------------------------*
gen rf_star50 = cond(rf_p50<.01,"***",cond(rf_p50<.05,"**",cond(rf_p50<.1,"*","")))
gen rf_star100 = cond(rf_p100<.01,"***",cond(rf_p100<.05,"**",cond(rf_p100<.1,"*","")))
gen rf_star200 = cond(rf_p200<.01,"***",cond(rf_p200<.05,"**",cond(rf_p200<.1,"*","")))

*---------------------------------------------------*
* 10) FORMAT STRINGS
*---------------------------------------------------*
gen rf_coef50  = string(rf_b50,"%9.3f") + rf_star50
gen rf_coef100 = string(rf_b100,"%9.3f") + rf_star100
gen rf_coef200 = string(rf_b200,"%9.3f") + rf_star200

gen rf_se50_s  = "(" + string(rf_se50,"%9.3f") + ")"
gen rf_se100_s = "(" + string(rf_se100,"%9.3f") + ")"
gen rf_se200_s = "(" + string(rf_se200,"%9.3f") + ")"

gen rf_N50_s   = string(rf_N50,"%9.0f")
gen rf_N100_s  = string(rf_N100,"%9.0f")
gen rf_N200_s  = string(rf_N200,"%9.0f")

*---------------------------------------------------*
* EXPORT LATEX TABLE: REDUCED FORM
*---------------------------------------------------*
file open tex_rf using "results/results_rd3rdorder_reducedform.tex", write replace

file write tex_rf "\begin{table}[htbp]" _n
file write tex_rf "\centering" _n
file write tex_rf "\caption{Reduced-Form Effect of Treatment Across Bandwidths (3rd-Order RD Polynomial)}" _n
file write tex_rf "\label{tab:rf_bw_rd3}" _n
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
file write tex_rf "\multicolumn{4}{p{0.9\linewidth}}{\footnotesize Notes: For each outcome and bandwidth, the table reports the reduced-form coefficient on treated with clustered standard errors in parentheses and the number of observations used in the regression. All regressions include border-pair fixed effects, a 3rd-order border-pair-specific RD polynomial in distance to border interacted with treatment, geographic controls, ethnicity fixed effects, survey-year fixed effects, and individual controls. Stars indicate significance: * \$p<0.10\$, ** \$p<0.05\$, *** \$p<0.01\$.}" _n
file write tex_rf "\end{tabular}" _n
file write tex_rf "\end{table}" _n
file close tex_rf
