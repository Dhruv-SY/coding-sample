****************************************************
* 11_lee_bounds.do
* Lee (2009) bounds for CONDITIONAL outcomes whose sample
* composition the treatment itself shifts.
*
* Applied to the two timing outcomes, which are only defined for
* a treatment-affected subpopulation:
*   - Age at first birth: conditional on ever having given birth,
*     and treatment reduces fertility
*   - Age at first sex: conditional on ever having had sex
* For these, the conditional treatment-control comparison mixes
* the causal effect with a composition change (e.g. marginal
* girls pushed INTO school by the law). Lee bounds trim the
* group with the higher selection rate to equalize selection,
* bounding the effect for always-selected women.
*
* NOTE: Lee bounds are defined for a simple binary-treatment
* comparison; they do not accommodate the RD terms or fixed
* effects. Bounds below are computed on the raw treated/control
* contrast within the 100km bandwidth and should be read as
* bounding/supporting evidence for the conditional outcomes,
* not as replacements for the main RD estimates.
****************************************************

do "code/00_setup.do"

keep if abs(rd) <= 100

* Selection is encoded by missingness of the outcome in the full
* 100km sample: each variable below is missing exactly for women
* outside its conditioning group (plus item nonresponse).

tempname memhold
tempfile results_lee
postfile `memhold' str40 outcome double sel_c sel_t lower upper N using `results_lee', replace

* outcome list: conditional outcomes only
local leevars "v212 v525"

foreach y of local leevars {

    di "========================================"
    di "Lee bounds: `y'"
    di "========================================"

    * selection rates by arm (share of women with the outcome defined)
    quietly count if treated == 0
    local n0 = r(N)
    quietly count if treated == 0 & !missing(`y')
    local s0 = r(N) / `n0'
    quietly count if treated == 1
    local n1 = r(N)
    quietly count if treated == 1 & !missing(`y')
    local s1 = r(N) / `n1'

    capture noisily leebounds `y' treated
    if _rc == 0 {
        matrix b = e(b)
        local lo = b[1,1]
        local hi = b[1,2]
        local n  = e(N)
        post `memhold' ("`y'") (`s0') (`s1') (`lo') (`hi') (`n')
    }
    else {
        di as error "leebounds failed for `y', rc=" _rc
        post `memhold' ("`y'") (`s0') (`s1') (.) (.) (.)
    }
}

postclose `memhold'

use `results_lee', clear

replace outcome = "Age at first birth"               if outcome == "v212"
replace outcome = "Age at first sex"                 if outcome == "v525"

gen sel_c_s = string(sel_c, "%9.3f")
gen sel_t_s = string(sel_t, "%9.3f")
gen lo_s    = string(lower, "%9.3f")
gen hi_s    = string(upper, "%9.3f")
gen N_s     = string(N, "%9.0f")

export delimited using "results/lee_bounds.csv", replace

file open tex using "results/lee_bounds.tex", write replace
file write tex "\begin{table}[htbp]" _n
file write tex "\centering" _n
file write tex "\caption{Lee (2009) Bounds for Conditional Outcomes (100km Bandwidth)}" _n
file write tex "\label{tab:lee_bounds}" _n
file write tex "\begin{tabular}{lccccc}" _n
file write tex "\hline\hline" _n
file write tex "Outcome & Selection (control) & Selection (treated) & Lower bound & Upper bound & N \\" _n
file write tex "\hline" _n

count
local N = r(N)
forvalues i = 1/`N' {
    local out = outcome[`i']
    local sc  = sel_c_s[`i']
    local st  = sel_t_s[`i']
    local lo  = lo_s[`i']
    local hi  = hi_s[`i']
    local nn  = N_s[`i']
    file write tex "`out'" " & " "`sc'" " & " "`st'" " & " "`lo'" " & " "`hi'" " & " "`nn'" " \\" _n
}

file write tex "\hline" _n
file write tex "\multicolumn{6}{p{0.95\linewidth}}{\footnotesize Notes: Lee (2009) bounds for the timing outcomes, which are defined only on a treatment-affected subpopulation (ever gave birth; ever had sex). Selection columns report the share of women in each arm for whom the outcome is defined; the group with higher selection is trimmed from above/below to equalize selection rates, bounding the average treatment effect for always-selected women. Bounds are computed on the raw treated/control contrast within the 100km bandwidth (Lee bounds do not accommodate the RD terms or fixed effects) and complement, rather than replace, the main RD estimates. Computed with the leebounds Stata package (Tauchmann 2014).}" _n
file write tex "\end{tabular}" _n
file write tex "\end{table}" _n
file close tex
