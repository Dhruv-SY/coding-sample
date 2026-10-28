****************************************************
* 15_summary_stats.do
* Descriptive statistics by country (mean and median) on the FULL
* country DHS samples -- analogous to thesis Tables 16 and 17,
* which report "summary statistics for the full country samples,"
* not the border-bandwidth-restricted RD estimation sample.
****************************************************

clear all
set more off

*---------------------------------------------------*
* 0) LOAD + APPEND ALL 12 PRE-SPATIAL-MERGE COUNTRY-PAIR FILES
*    (full national DHS samples, before any border/bandwidth
*    restriction). Countries appearing in multiple pairs (e.g.
*    Tanzania, Namibia, Burkina Faso) are duplicated across files
*    and de-duplicated below.
*---------------------------------------------------*
tempfile combined

use "../Pre Spatial Merge/Cote d'Ivoire-Mali.dta", clear
save `combined', replace

use "../Pre Spatial Merge/Rwanda-Tanzania.dta", clear
append using `combined'
save `combined', replace

use "../Pre Spatial Merge/Burundi-Tanzania.dta", clear
append using `combined'
save `combined', replace

use "../Pre Spatial Merge/Cote d'Ivoire-Burkina Faso.dta", clear
append using `combined'
save `combined', replace

use "../Pre Spatial Merge/Uganda-Tanzania.dta", clear
append using `combined'
save `combined', replace

use "../Pre Spatial Merge/Namibia-South Africa.dta", clear
append using `combined'
save `combined', replace

use "../Pre Spatial Merge/Namibia-Angola.dta", clear
append using `combined'
save `combined', replace

use "../Pre Spatial Merge/Namibia-Zambia.dta", clear
append using `combined'
save `combined', replace

use "../Pre Spatial Merge/Namibia-Zimbabwe.dta", clear
append using `combined'
save `combined', replace

use "../Pre Spatial Merge/Ghana-Burkina Faso.dta", clear
append using `combined'
save `combined', replace

use "../Pre Spatial Merge/Benin-Niger.dta", clear
append using `combined'
save `combined', replace

use "../Pre Spatial Merge/Benin-Burkina Faso.dta", clear
append using `combined'
save `combined', replace

use `combined', clear

* countries appearing in multiple pair files are exact full-row
* duplicates of the same DHS respondent; keep one copy each
duplicates drop

*---------------------------------------------------*
* 1) SAME SAMPLE RESTRICTION AND VARIABLE CONSTRUCTION AS 00_setup.do
*---------------------------------------------------*
drop if v012 < 18 & v511 == .

replace v133 = . if inlist(v133, 97, 98, 99)
replace v525 = . if inlist(v525, 0, 96, 97, 98, 99)
replace v741 = . if v741 == 9
replace v717 = . if inlist(v717, 10, 96, 98, 99)

gen not_working = .
replace not_working = 1 if v714 == 0
replace not_working = 0 if v714 == 1

gen child_marriage = 0
replace child_marriage = 1 if v511 < 18 & !missing(v511)

gen married  = !missing(v511)
gen employed = 1 - not_working
label var married  "Ever married/in union"
label var employed "Employed"
label var child_marriage "Married before 18"

*---------------------------------------------------*
* 2) MEAN TABLE (Table 16 analogue)
*---------------------------------------------------*
preserve
    collapse (mean) v012 married child_marriage v212 v133 employed ///
             (count) obs=v012, by(country_name)

    gen treat_grp = .
    replace treat_grp = 1 if inlist(country_name, ///
        "Cote d'Ivoire", "Rwanda", "Burundi", "Uganda", "Namibia", "Ghana", "Benin")
    replace treat_grp = 0 if inlist(country_name, ///
        "Mali", "Tanzania", "Burkina Faso", "South Africa", "Angola", ///
        "Zambia", "Zimbabwe", "Niger")

    tempfile means
    save `means'

    collapse (mean) v012 married child_marriage v212 v133 employed [aw=obs], by(treat_grp)
    gen country_name = cond(treat_grp==1, "Treatment avg.", "Control avg.")
    drop treat_grp
    tempfile avgs
    save `avgs'

    use `means', clear
    drop treat_grp
    append using `avgs'

    file open tex using "results/summary_stats_mean.tex", write replace
    file write tex "\begin{table}[htbp]" _n
    file write tex "\centering" _n
    file write tex "\caption{Mean Summary Statistics by Country (Full Country Samples)}" _n
    file write tex "\label{tab:sumstats_mean}" _n
    file write tex "\begin{tabular}{lccccccc}" _n
    file write tex "\hline\hline" _n
    file write tex "Country & Obs. & Age & Married & Before 18 & First birth & Educ. (yrs) & Employed \\" _n
    file write tex "\hline" _n

    count
    local N = r(N)
    forvalues i = 1/`N' {
        local cty  = country_name[`i']
        local o    = obs[`i']
        local age  = v012[`i']
        local mar  = married[`i']
        local b18  = child_marriage[`i']
        local fb   = v212[`i']
        local ed   = v133[`i']
        local emp  = employed[`i']

        if "`cty'" == "Treatment avg." | "`cty'" == "Control avg." {
            file write tex "`cty' & & " %9.2f (`age') " & " %9.2f (`mar') " & " %9.2f (`b18') ///
                " & " %9.2f (`fb') " & " %9.2f (`ed') " & " %9.2f (`emp') " \\" _n
        }
        else {
            file write tex "`cty' & " %9.0fc (`o') " & " %9.2f (`age') " & " %9.2f (`mar') ///
                " & " %9.2f (`b18') " & " %9.2f (`fb') " & " %9.2f (`ed') " & " %9.2f (`emp') " \\" _n
        }
    }

    file write tex "\hline" _n
    file write tex "\multicolumn{8}{p{0.95\linewidth}}{\footnotesize Notes: Mean values by country, full national DHS samples (all women 15--49 in the survey rounds used, not restricted to the border-bandwidth RD estimation sample). Treatment countries have minimum legal marriage age \$\geq\$18; control countries have minimum legal marriage age \$<\$18. Child marriage (Before 18) is defined as first marriage or cohabitation before age 18 and is defined for all women (never-married adult women are coded 0). Employed is one minus the not-working indicator. Age is respondent age at interview (15--49).}" _n
    file write tex "\end{tabular}" _n
    file write tex "\end{table}" _n
    file close tex

    export delimited using "results/summary_stats_mean.csv", replace
restore

*---------------------------------------------------*
* 3) MEDIAN TABLE (Table 17 analogue)
*---------------------------------------------------*
preserve
    collapse (median) v012 married child_marriage v212 v133 employed ///
             (count) obs=v012, by(country_name)

    file open tex using "results/summary_stats_median.tex", write replace
    file write tex "\begin{table}[htbp]" _n
    file write tex "\centering" _n
    file write tex "\caption{Median Summary Statistics by Country (Full Country Samples)}" _n
    file write tex "\label{tab:sumstats_median}" _n
    file write tex "\begin{tabular}{lccccccc}" _n
    file write tex "\hline\hline" _n
    file write tex "Country & Obs. & Age & Married & Before 18 & First birth & Educ. (yrs) & Employed \\" _n
    file write tex "\hline" _n

    count
    local N = r(N)
    forvalues i = 1/`N' {
        local cty  = country_name[`i']
        local o    = obs[`i']
        local age  = v012[`i']
        local mar  = married[`i']
        local b18  = child_marriage[`i']
        local fb   = v212[`i']
        local ed   = v133[`i']
        local emp  = employed[`i']

        file write tex "`cty' & " %9.0fc (`o') " & " %9.0f (`age') " & " %9.0f (`mar') ///
            " & " %9.0f (`b18') " & " %9.0f (`fb') " & " %9.0f (`ed') " & " %9.0f (`emp') " \\" _n
    }

    file write tex "\hline" _n
    file write tex "\multicolumn{8}{p{0.95\linewidth}}{\footnotesize Notes: Median values by country, full national DHS samples (all women 15--49 in the survey rounds used, not restricted to the border-bandwidth RD estimation sample). Age refers to respondent age at interview (15--49). Child marriage (Before 18) is defined as marriage before age 18 and is defined for all women (never-married adult women are coded 0). Medians are computed using available observations for each variable.}" _n
    file write tex "\end{tabular}" _n
    file write tex "\end{table}" _n
    file close tex

    export delimited using "results/summary_stats_median.csv", replace
restore
