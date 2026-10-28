****************************************************
* 16_balance_test.do
* Balance table for RD covariates
* Formatted like Anderson (2017) "Legal Origins and
* Female HIV" Appendix Table A5: a single column of
* "coefficient on treatment" per covariate, rather than
* separate control/treated means.
****************************************************

clear all
set more off

capture mkdir "results"

* Use the same data load, sample restrictions, cluster_id construction, and
* nearest-border dedup as every other script, instead of a separate hardcoded
* copy of this logic (which had drifted: stale dataset path, and missing the
* cluster_id variable entirely despite the regressions below clustering on it).
do "code/00_setup.do"

*---------------------------------------------------*
* 2) RD SPECIFICATION
*---------------------------------------------------*

global RDBP "i.bp c.rd#i.bp c.rd#i.bp#i.treated"

* Balance sample: 100 km bandwidth
capture drop balance_sample
gen byte balance_sample = abs(rd) <= 100

*---------------------------------------------------*
* 3) CLEAN / CONSTRUCT BALANCE VARIABLES
*---------------------------------------------------*

capture drop no_educ

* Wife no education (mirrors Anderson Table A5's "Wife no education")
gen byte no_educ = .
replace no_educ = 1 if v133 == 0
replace no_educ = 0 if v133 > 0 & !missing(v133)
label var no_educ "Wife no education"

* Labels (population_density_mean is already people/km^2, extracted directly
* from the GHSL raster at each woman's own coordinate -- no pixel-area
* conversion needed, unlike the old PRIO-CSV-based construction)
label var v012                   "Wife age"
label var population_density_mean "Population density"
label var nightlights_mean        "Luminosity"
label var water_wetland_pct       "Water area"
label var elevation_mean          "Elevation"
label var agri_suitability_mean   "Soil suitability"
label var malaria_ecology_mean    "Malaria index"
label var capdist                 "Distance to capital"
label var traveltime_mean         "Travel time to market"
label var rainfall_mean           "Rainfall"
label var temperature_mean        "Temperature"
label var dist_coast_km           "Distance to coast"
label var diamsec_s               "Diamond deposit in pixel"
label var dist_colonial_railway_km   "Distance to colonial railway"
label var dist_protestant_mission_km "Distance to Protestant mission"
label var dist_catholic_mission_km   "Distance to Catholic mission"
label var dist_trade_route_km        "Distance to trade route"
label var gdppc_2000              "GDP per capita (2000, country)"
label var ruggedness              "Terrain ruggedness (country)"
label var dist_coast_ctry         "Distance to coast (country-level)"

*---------------------------------------------------*
* 4) VARIABLES TO TEST (Table-A5 style: one row per covariate)
*---------------------------------------------------*

local balvars ///
    v012 ///
    no_educ ///
    population_density_mean ///
    nightlights_mean ///
    water_wetland_pct ///
    elevation_mean ///
    agri_suitability_mean ///
    malaria_ecology_mean ///
    capdist ///
    traveltime_mean ///
    rainfall_mean ///
    temperature_mean ///
    dist_coast_km ///
    diamsec_s ///
    dist_colonial_railway_km ///
    dist_protestant_mission_km ///
    dist_catholic_mission_km ///
    dist_trade_route_km ///
    gdppc_2000 ///
    ruggedness ///
    dist_coast_ctry

*---------------------------------------------------*
* 5) STORE BALANCE RESULTS
*---------------------------------------------------*

tempfile balance_results
tempname balhold

postfile `balhold' str60 varname double b se p N using `balance_results', replace

foreach y of local balvars {

    preserve

        keep if balance_sample == 1
        keep if !missing(`y', treated, rd, bp, cluster_id)

        capture quietly reg `y' treated $RDBP, vce(cluster cluster_id)

        if _rc == 0 {
            local b  = _b[treated]
            local se = _se[treated]
            local t  = `b' / `se'
            local p  = 2 * ttail(e(df_r), abs(`t'))
            local N  = e(N)

            post `balhold' ("`y'") (`b') (`se') (`p') (`N')
        }
        else {
            di as error "Balance regression failed for `y', rc=" _rc
        }

    restore
}

postclose `balhold'

*---------------------------------------------------*
* 6) STORE JOINT TESTS (fixed-effects balance)
*---------------------------------------------------*

tempfile joint_results
tempname jhold

postfile `jhold' str60 testname double p N using `joint_results', replace

preserve

    keep if balance_sample == 1

    capture quietly reg treated i.survey_year $RDBP, vce(cluster cluster_id)
    if _rc == 0 {
        capture quietly testparm i.survey_year
        if _rc == 0 {
            post `jhold' ("Survey-year fixed effects") (r(p)) (e(N))
        }
    }

    capture quietly reg treated i.eth_id $RDBP, vce(cluster cluster_id)
    if _rc == 0 {
        capture quietly testparm i.eth_id
        if _rc == 0 {
            post `jhold' ("Ethnicity fixed effects") (r(p)) (e(N))
        }
    }

restore

postclose `jhold'

*---------------------------------------------------*
* 7) FORMAT
*---------------------------------------------------*

use `balance_results', clear

replace varname = "Wife age"               if varname == "v012"
replace varname = "Wife no education"      if varname == "no_educ"
replace varname = "Population density"     if varname == "population_density_mean"
replace varname = "Luminosity"             if varname == "nightlights_mean"
replace varname = "Water area"             if varname == "water_wetland_pct"
replace varname = "Elevation"              if varname == "elevation_mean"
replace varname = "Soil suitability"       if varname == "agri_suitability_mean"
replace varname = "Malaria index"          if varname == "malaria_ecology_mean"
replace varname = "Distance to capital"    if varname == "capdist"
replace varname = "Travel time to market"  if varname == "traveltime_mean"
replace varname = "Rainfall"               if varname == "rainfall_mean"
replace varname = "Temperature"            if varname == "temperature_mean"
replace varname = "Distance to coast"      if varname == "dist_coast_km"
replace varname = "Diamond deposit in pixel" if varname == "diamsec_s"
replace varname = "Distance to colonial railway"   if varname == "dist_colonial_railway_km"
replace varname = "Distance to Protestant mission" if varname == "dist_protestant_mission_km"
replace varname = "Distance to Catholic mission"   if varname == "dist_catholic_mission_km"
replace varname = "Distance to trade route"        if varname == "dist_trade_route_km"
replace varname = "GDP per capita (country-level)"     if varname == "gdppc_2000"
replace varname = "Terrain ruggedness (country-level)" if varname == "ruggedness"
replace varname = "Distance to coast (country-level)"  if varname == "dist_coast_ctry"

gen star = cond(p < .01, "***", cond(p < .05, "**", cond(p < .1, "*", "")))
gen coef_s = string(b, "%9.3f") + star
gen se_s   = "(" + string(se, "%9.3f") + ")"
gen N_s    = string(N, "%9.0f")

tempfile panelA
save `panelA', replace

use `joint_results', clear
gen p_s = string(p, "%9.3f")
gen N_s = string(N, "%9.0f")

tempfile panelB
save `panelB', replace

*---------------------------------------------------*
* 8) EXPORT LATEX BALANCE TABLE (Anderson Table A5 style)
*---------------------------------------------------*

use `panelA', clear

file open tex using "results/balance_table.tex", write replace

file write tex "\begin{table}[htbp]" _n
file write tex "\centering" _n
file write tex "\caption{Estimations of Control Variables (RD Balance Test)}" _n
file write tex "\label{tab:balance}" _n
file write tex "\begin{tabular}{lc}" _n
file write tex "\hline\hline" _n
file write tex "Variable & Coefficient on Treated \\" _n
file write tex "\hline" _n

count
local N = r(N)

forvalues i = 1/`N' {
    local v  = varname[`i']
    local cf = coef_s[`i']
    local se = se_s[`i']

    file write tex "`v'" " & " "`cf'" " \\" _n
    file write tex " " " & " "`se'" " \\" _n
}

file write tex "\hline" _n

use `panelB', clear
count
local N2 = r(N)

forvalues i = 1/`N2' {
    local v  = testname[`i']
    local pp = p_s[`i']
    local nn = N_s[`i']

    file write tex "`v'" " & p-value: `pp' (N=`nn') \\" _n
}

file write tex "\hline" _n
file write tex "\multicolumn{2}{p{0.75\linewidth}}{\footnotesize Notes: Each row reports the coefficient (standard error) on treatment from a separate regression of the covariate on treatment, including border-pair fixed effects, border-pair-specific RD slopes, and treated-side kinks in distance to the border. Standard errors clustered at the DHS-cluster level (country x survey-round x PSU). Estimates use observations within 100 km of the border cutoff. GDP per capita, terrain ruggedness, and this country-level distance to coast (distinct from the pixel-level distance-to-coast row above, which is measured at each respondent's own coordinate) are COUNTRY-level variables: within a border pair they take only two values, so their rows are within-pair country differences whose clustered standard errors are understated; read them as descriptive. Stars: * \$p<0.10\$, ** \$p<0.05\$, *** \$p<0.01\$.}" _n
file write tex "\end{tabular}" _n
file write tex "\end{table}" _n

file close tex
