****************************************************
* 13_twoway_cluster.do
* Reduced-form results with the updated pixel controls
* (rainfall, temperature, distance to coast),
* comparing our standard DHS-cluster SEs against two-way
* ethnicity x country clustering (matching Anderson (2018),
* via ivreg2's Cameron-Gelbach-Miller implementation).
*
* CAVEAT: with only 15 countries and ~100 ethnic homelands in
* our 12-pair sample, many ethnicities are nested within a
* single country. ivreg2 flags the two-way cluster covariance
* matrix as not of full rank ("number of clusters insufficient")
* -- point estimates are unaffected, but the two-way SEs should
* be read with real caution, not as a drop-in replacement for
* the standard clustering.
****************************************************

do "code/00_setup.do"

encode country_name, gen(country_id)

local bandwidths "50 100 200"
local outcomes  "v133 educ_any educ_years_pos v106 v212 v201 v525 multiple_unions not_working domestic_work paid_cash white_collar white_collar_cond"

tempname memhold
tempfile results_bw
postfile `memhold' str40 outcome ///
    double b50 se_std50 se_2way50 N50 ///
    double b100 se_std100 se_2way100 N100 ///
    double b200 se_std200 se_2way200 N200 ///
    using `results_bw', replace

foreach y of local outcomes {

    local b50 .
    local se_std50 .
    local se_2way50 .
    local N50 .
    local b100 .
    local se_std100 .
    local se_2way100 .
    local N100 .
    local b200 .
    local se_std200 .
    local se_2way200 .
    local N200 .

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
                local b  = _b[treated]
                local se_std = _se[treated]
                local n  = e(N)

                quietly ivreg2 `y' treated $RDBP $X $FE, cluster(eth_id country_id)
                local se_2way = _se[treated]

                if `h' == 50 {
                    local b50 = `b'
                    local se_std50 = `se_std'
                    local se_2way50 = `se_2way'
                    local N50 = `n'
                }
                else if `h' == 100 {
                    local b100 = `b'
                    local se_std100 = `se_std'
                    local se_2way100 = `se_2way'
                    local N100 = `n'
                }
                else if `h' == 200 {
                    local b200 = `b'
                    local se_std200 = `se_std'
                    local se_2way200 = `se_2way'
                    local N200 = `n'
                }
            }
        restore
    }

    post `memhold' ("`y'") (`b50') (`se_std50') (`se_2way50') (`N50') (`b100') (`se_std100') (`se_2way100') (`N100') (`b200') (`se_std200') (`se_2way200') (`N200')
}

postclose `memhold'

use `results_bw', clear

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

export delimited using "results/twoway_cluster_comparison.csv", replace

file open tex using "results/twoway_cluster_comparison.tex", write replace
file write tex "\begin{table}[htbp]" _n
file write tex "\centering" _n
file write tex "\caption{Reduced-Form Effect of Treatment: Standard vs. Two-Way (Ethnicity x Country) Clustering}" _n
file write tex "\label{tab:twoway_cluster}" _n
file write tex "\begin{tabular}{lcccccc}" _n
file write tex "\hline\hline" _n
file write tex " & \multicolumn{2}{c}{50km} & \multicolumn{2}{c}{100km} & \multicolumn{2}{c}{200km} \\" _n
file write tex "Outcome & Coef. & SE(std) / SE(2way) & Coef. & SE(std) / SE(2way) & Coef. & SE(std) / SE(2way) \\" _n
file write tex "\hline" _n

count
local N = r(N)
forvalues i = 1/`N' {
    local out = outcome[`i']
    local b50 = b50[`i']
    local ses50 = se_std50[`i']
    local se2_50 = se_2way50[`i']
    local b100 = b100[`i']
    local ses100 = se_std100[`i']
    local se2_100 = se_2way100[`i']
    local b200 = b200[`i']
    local ses200 = se_std200[`i']
    local se2_200 = se_2way200[`i']

    local line = "`out'" + " & " + string(`b50',"%9.3f") + " & " + string(`ses50',"%9.3f") + " / " + string(`se2_50',"%9.3f") + ///
        " & " + string(`b100',"%9.3f") + " & " + string(`ses100',"%9.3f") + " / " + string(`se2_100',"%9.3f") + ///
        " & " + string(`b200',"%9.3f") + " & " + string(`ses200',"%9.3f") + " / " + string(`se2_200',"%9.3f") + " \\"
    file write tex `"`line'"' _n
}

file write tex "\hline" _n
file write tex "\multicolumn{7}{p{0.95\linewidth}}{\footnotesize Notes: For each outcome and bandwidth, the table reports the reduced-form coefficient on treated, the standard error clustered at the DHS-cluster level (SE std), and the standard error from two-way clustering at the ethnicity x country level via ivreg2 (SE 2way), matching Anderson (2018). CAUTION: with only 15 countries and roughly 100 ethnic homelands in this sample, many ethnicities are nested within a single country, and ivreg2 flags the two-way cluster covariance matrix as not of full rank. The two-way SEs should be read as indicative, not a reliable replacement for the standard clustering. All regressions include border-pair fixed effects, border-pair-specific slopes and kinks, the full geographic/pixel control set (including rainfall, temperature, and distance to coast), ethnicity fixed effects, survey-year fixed effects, and individual controls.}" _n
file write tex "\end{tabular}" _n
file write tex "\end{table}" _n
file close tex
