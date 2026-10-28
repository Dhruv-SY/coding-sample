****************************************************
* 02_first_stage.do
* First-stage estimates across outcome-specific samples
****************************************************

do "code/00_setup.do"

*---------------------------------------------------*
* 5) BANDWIDTHS + OUTCOME SAMPLES
*---------------------------------------------------*
local bandwidths "50 100 200"
local outcomes  "v133 educ_any educ_years_pos v106 v212 v201 v525 multiple_unions not_working domestic_work paid_cash white_collar white_collar_cond"

*---------------------------------------------------*
* 6) STORE FIRST-STAGE RESULTS ONLY
*---------------------------------------------------*
tempname memhold
tempfile results_firststage_bw

postfile `memhold' str40 outcome ///
    double fs_b50 fs_se50 fs_p50 fs_N50 ///
    double fs_b100 fs_se100 fs_p100 fs_N100 ///
    double fs_b200 fs_se200 fs_p200 fs_N200 ///
    using `results_firststage_bw', replace

foreach y of local outcomes {

    di "========================================"
    di "First-stage outcome sample: `y'"
    di "========================================"

    local fs_b50  .
    local fs_se50 .
    local fs_p50  .
    local fs_N50  .

    local fs_b100  .
    local fs_se100 .
    local fs_p100  .
    local fs_N100  .

    local fs_b200  .
    local fs_se200 .
    local fs_p200  .
    local fs_N200  .

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
            * First stage
            *-----------------------------*
            quietly reg child_marriage treated $RDBP $X $FE, vce(cluster cluster_id)

            local fs_b  = _b[treated]
            local fs_se = _se[treated]
            local fs_t  = _b[treated] / _se[treated]
            local fs_p  = 2*ttail(e(df_r), abs(`fs_t'))
            local fs_N  = e(N)

            if `h' == 50 {
                local fs_b50  = `fs_b'
                local fs_se50 = `fs_se'
                local fs_p50  = `fs_p'
                local fs_N50  = `fs_N'
            }
            else if `h' == 100 {
                local fs_b100  = `fs_b'
                local fs_se100 = `fs_se'
                local fs_p100  = `fs_p'
                local fs_N100  = `fs_N'
            }
            else if `h' == 200 {
                local fs_b200  = `fs_b'
                local fs_se200 = `fs_se'
                local fs_p200  = `fs_p'
                local fs_N200  = `fs_N'
            }

        restore
    }

    post `memhold' ("`y'") ///
        (`fs_b50') (`fs_se50') (`fs_p50') (`fs_N50') ///
        (`fs_b100') (`fs_se100') (`fs_p100') (`fs_N100') ///
        (`fs_b200') (`fs_se200') (`fs_p200') (`fs_N200')
}

postclose `memhold'

*---------------------------------------------------*
* 7) OPEN RESULTS
*---------------------------------------------------*
use `results_firststage_bw', clear

*---------------------------------------------------*
* 8) BETTER OUTCOME SAMPLE LABELS
*---------------------------------------------------*
replace outcome = "Education (years)"                 if outcome == "v133"
replace outcome = "Any education"                     if outcome == "educ_any"
replace outcome = "Education (years | schooling>0)"   if outcome == "educ_years_pos"
replace outcome = "Education level"                   if outcome == "v106"
replace outcome = "Age at first birth"                if outcome == "v212"
replace outcome = "Children ever born"                if outcome == "v201"
replace outcome = "Age at first sex"                  if outcome == "v525"
replace outcome = "Multiple unions"                   if outcome == "multiple_unions"
replace outcome = "Not working"                       if outcome == "not_working"
replace outcome = "Domestic work"                     if outcome == "domestic_work"
replace outcome = "Paid in cash"                      if outcome == "paid_cash"
replace outcome = "White-collar employment"           if outcome == "white_collar"
replace outcome = "White-collar among workers"        if outcome == "white_collar_cond"

*---------------------------------------------------*
* 9) ADD STARS
*---------------------------------------------------*
gen fs_star50  = cond(fs_p50<.01,"***",cond(fs_p50<.05,"**",cond(fs_p50<.1,"*","")))
gen fs_star100 = cond(fs_p100<.01,"***",cond(fs_p100<.05,"**",cond(fs_p100<.1,"*","")))
gen fs_star200 = cond(fs_p200<.01,"***",cond(fs_p200<.05,"**",cond(fs_p200<.1,"*","")))

*---------------------------------------------------*
* 10) FORMAT STRINGS
*---------------------------------------------------*
gen fs_coef50  = string(fs_b50,"%9.3f") + fs_star50
gen fs_coef100 = string(fs_b100,"%9.3f") + fs_star100
gen fs_coef200 = string(fs_b200,"%9.3f") + fs_star200

gen fs_se50_s  = "(" + string(fs_se50,"%9.3f") + ")"
gen fs_se100_s = "(" + string(fs_se100,"%9.3f") + ")"
gen fs_se200_s = "(" + string(fs_se200,"%9.3f") + ")"

gen fs_N50_s   = string(fs_N50,"%9.0f")
gen fs_N100_s  = string(fs_N100,"%9.0f")
gen fs_N200_s  = string(fs_N200,"%9.0f")

*---------------------------------------------------*
* 11) EXPORT LATEX TABLE: FIRST STAGE ONLY
*---------------------------------------------------*
file open tex_fs using "results/results_firststage.tex", write replace

file write tex_fs "\begin{table}[htbp]" _n
file write tex_fs "\centering" _n
file write tex_fs "\caption{First-Stage Effect of Treatment on Child Marriage Across Bandwidths}" _n
file write tex_fs "\label{tab:firststage_bw}" _n
file write tex_fs "\begin{tabular}{lccc}" _n
file write tex_fs "\hline\hline" _n
file write tex_fs " & \multicolumn{3}{c}{Bandwidth (km)} \\" _n
file write tex_fs "\cline{2-4}" _n
file write tex_fs "Outcome sample & 50 & 100 & 200 \\" _n
file write tex_fs "\hline" _n

count
local N = r(N)

forvalues i = 1/`N' {
    local out    = outcome[`i']
    local fsc50  = fs_coef50[`i']
    local fsc100 = fs_coef100[`i']
    local fsc200 = fs_coef200[`i']
    local fss50  = fs_se50_s[`i']
    local fss100 = fs_se100_s[`i']
    local fss200 = fs_se200_s[`i']
    local fsn50  = fs_N50_s[`i']
    local fsn100 = fs_N100_s[`i']
    local fsn200 = fs_N200_s[`i']

    file write tex_fs "`out'" " & " "`fsc50'" " & " "`fsc100'" " & " "`fsc200'" " \\" _n
    file write tex_fs " " " & " "`fss50'" " & " "`fss100'" " & " "`fss200'" " \\" _n
    file write tex_fs "Obs." " & " "`fsn50'" " & " "`fsn100'" " & " "`fsn200'" " \\" _n
    file write tex_fs "\addlinespace" _n
}

file write tex_fs "\hline" _n
file write tex_fs "\multicolumn{4}{p{0.9\linewidth}}{\footnotesize Notes: For each outcome sample and bandwidth, the table reports the first-stage coefficient on treatment from a regression of child marriage on treatment. Standard errors clustered at the DHS cluster level are in parentheses. The first stage is estimated on the same outcome-specific sample used in the corresponding IV regression. All regressions include border-pair fixed effects, border-pair-specific slopes and kinks, geographic controls, ethnicity fixed effects, survey-year fixed effects, and individual controls. Stars indicate significance: * \$p<0.10\$, ** \$p<0.05\$, *** \$p<0.01\$.}" _n
file write tex_fs "\end{tabular}" _n
file write tex_fs "\end{table}" _n
file close tex_fs
