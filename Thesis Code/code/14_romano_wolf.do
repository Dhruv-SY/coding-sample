****************************************************
* 14_romano_wolf.do
* Romano-Wolf (2005, 2016) stepdown multiple-hypothesis
* correction for the reduced-form estimates at the 100km
* bandwidth, adjusted within three pre-specified outcome
* families (matching the thesis appendix structure):
*   - Education: years, level
*   - Fertility/timing: age at first birth, children ever born,
*     age at first sex, multiple unions
*   - Labour market: not working, domestic work, white-collar,
*     white-collar among workers
* The procedure controls the family-wise error rate under
* arbitrary dependence across outcomes within a family, via
* cluster bootstrap (resampling cluster_id).
****************************************************

do "code/00_setup.do"

keep if abs(rd) <= 100

* rwolf2 requires a common estimation command per equation; each
* outcome keeps its own estimation sample via the if conditions.
capture drop wcc_sample eyp_sample
gen byte wcc_sample = v714 == 1 & !missing(white_collar_cond)
gen byte eyp_sample = v133 > 0 & !missing(v133) & !missing(educ_years_pos)

tempname memhold
tempfile results_rw
postfile `memhold' str40 outcome str20 family double p_model p_resample p_rw using `results_rw', replace

*---------------------------------------------------*
* Family 1: education
*---------------------------------------------------*
rwolf2 (reg v133 treated $RDBP $X $FE, vce(cluster cluster_id)) ///
       (reg v106 treated $RDBP $X $FE, vce(cluster cluster_id)), ///
    indepvars(treated, treated) cluster(cluster_id) seed(42) reps(1000) ///
    usevalid nodots

matrix RW = e(RW)
local eqs "v133 v106"
local i = 1
foreach y of local eqs {
    post `memhold' ("`y'") ("Education") (RW[`i',1]) (RW[`i',2]) (RW[`i',3])
    local i = `i' + 1
}

*---------------------------------------------------*
* Family 2: fertility / timing
*---------------------------------------------------*
rwolf2 (reg v212 treated $RDBP $X $FE, vce(cluster cluster_id)) ///
       (reg v201 treated $RDBP $X $FE, vce(cluster cluster_id)) ///
       (reg v525 treated $RDBP $X $FE, vce(cluster cluster_id)) ///
       (reg multiple_unions treated $RDBP $X $FE, vce(cluster cluster_id)), ///
    indepvars(treated, treated, treated, treated) cluster(cluster_id) seed(42) reps(1000) ///
    usevalid nodots

matrix RW = e(RW)
local eqs "v212 v201 v525 multiple_unions"
local i = 1
foreach y of local eqs {
    post `memhold' ("`y'") ("Fertility") (RW[`i',1]) (RW[`i',2]) (RW[`i',3])
    local i = `i' + 1
}

*---------------------------------------------------*
* Family 3: labour market
*---------------------------------------------------*
rwolf2 (reg not_working treated $RDBP $X $FE, vce(cluster cluster_id)) ///
       (reg domestic_work treated $RDBP $X $FE, vce(cluster cluster_id)) ///
       (reg white_collar treated $RDBP $X $FE, vce(cluster cluster_id)) ///
       (reg white_collar_cond treated $RDBP $X $FE if wcc_sample, vce(cluster cluster_id)), ///
    indepvars(treated, treated, treated, treated) cluster(cluster_id) seed(42) reps(1000) ///
    usevalid nodots

matrix RW = e(RW)
local eqs "not_working domestic_work white_collar white_collar_cond"
local i = 1
foreach y of local eqs {
    post `memhold' ("`y'") ("Labour market") (RW[`i',1]) (RW[`i',2]) (RW[`i',3])
    local i = `i' + 1
}

postclose `memhold'

*---------------------------------------------------*
* FORMAT + EXPORT
*---------------------------------------------------*
use `results_rw', clear

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

gen star = cond(p_rw < .01, "***", cond(p_rw < .05, "**", cond(p_rw < .1, "*", "")))
gen pm_s  = string(p_model, "%9.3f")
gen pr_s  = string(p_resample, "%9.3f")
gen prw_s = string(p_rw, "%9.3f") + star

export delimited using "results/romano_wolf.csv", replace

file open tex using "results/romano_wolf.tex", write replace
file write tex "\begin{table}[htbp]" _n
file write tex "\centering" _n
file write tex "\caption{Romano-Wolf Multiple-Hypothesis Correction (100km Bandwidth)}" _n
file write tex "\label{tab:romano_wolf}" _n
file write tex "\begin{tabular}{llccc}" _n
file write tex "\hline\hline" _n
file write tex "Family & Outcome & Model p-value & Resample p-value & Romano-Wolf p-value \\" _n
file write tex "\hline" _n

count
local N = r(N)
forvalues i = 1/`N' {
    local fam = family[`i']
    local out = outcome[`i']
    local pm  = pm_s[`i']
    local pr  = pr_s[`i']
    local prw = prw_s[`i']
    file write tex "`fam'" " & " "`out'" " & " "`pm'" " & " "`pr'" " & " "`prw'" " \\" _n
}

file write tex "\hline" _n
file write tex "\multicolumn{5}{p{0.9\linewidth}}{\footnotesize Notes: Romano-Wolf (2005, 2016) stepdown adjusted p-values for the reduced-form coefficient on treated at the 100km bandwidth, controlling the family-wise error rate within three pre-specified outcome families (education, fertility/timing, labour market) under arbitrary within-family dependence. Computed via cluster bootstrap resampling DHS clusters (1,000 replications, rwolf2). Stars are based on the Romano-Wolf p-value: * \$p<0.10\$, ** \$p<0.05\$, *** \$p<0.01\$.}" _n
file write tex "\end{tabular}" _n
file write tex "\end{table}" _n
file close tex
