*==============================================================*
* MASTER DO-FILE
* Child Marriage, Legal Regimes, and Life Outcomes
* Dhruv Yadav and Pedro Silva
*
* Runs every Stata script in the analysis in sequence.
* Set the working directory to the root of this replication
* package (the folder containing master.do) before running.
*
* NOTE: code/12_greg_ethnicity.do needs data/greg_crosswalk.dta. It ships
* with the package; rebuild it with build/03_greg_crosswalk.R if needed.
*==============================================================*

*--------------------------------------------------------------*
* SET WORKING DIRECTORY
* Edit this line to point to your local copy of the package
*--------------------------------------------------------------*
* cd "path/to/Replication Package"

*--------------------------------------------------------------*
* CONFIRM DATA EXISTS
*--------------------------------------------------------------*
capture confirm file "data/20yr_threshold_12pairs.dta"
if _rc {
    di as error "Analysis dataset not found."
    di as error "Place 20yr_threshold_12pairs.dta in data/, or rebuild it"
    di as error "by running build/01_build_pairs.R then build/02_pool_pairs.R."
    exit 601
}

*--------------------------------------------------------------*
* RUN SCRIPTS
* Each script calls code/00_setup.do to load the data, construct
* outcomes and set the specification globals, then applies its
* own filters and writes LaTeX tables to results/.
*--------------------------------------------------------------*

* --- MAIN RESULTS ---
* Reduced-form effect of treatment, 13 outcomes x 3 bandwidths
do "code/01_main_reducedform.do"

* First-stage: effect of treatment on child marriage
do "code/02_first_stage.do"

* Heterogeneity by urban/rural residence
do "code/03_heterogeneity_urban_rural.do"

* Exclude migrants who moved in the same period as their marriage
do "code/04_robustness_migration.do"

* --- ROBUSTNESS: SAMPLE AND SPECIFICATION ---
* Exclude the 4 Namibia border pairs
do "code/05_robustness_no_namibia.do"

* Exclude the Benin-Niger border pair
do "code/06_robustness_no_niger.do"

* Drop pixel-level geographic controls and ethnicity FE
do "code/07_robustness_no_pixel_controls.do"

* 2nd-order RD polynomial (matches Anderson (2018) main spec)
do "code/08_robustness_rd_2ndorder.do"

* 3rd-order RD polynomial (matches Anderson (2018) Table A7 col. 3)
do "code/09_robustness_rd_3rdorder.do"

* Leave-one-pair-out, reduced form only (100km)
do "code/10_robustness_leave_one_pair_out.do"


* Lee (2009) bounds for conditional outcomes
do "code/11_lee_bounds.do"

* GREG (Atlas Narodov Mira) ethnic FE instead of Murdock
* (requires data/greg_crosswalk.dta -- see build/03_greg_crosswalk.R)
do "code/12_greg_ethnicity.do"

* --- ROBUSTNESS: INFERENCE WITH FEW CLUSTERS ---
* Two-way ethnicity x country clustering (Cameron-Gelbach-Miller,
* via ivreg2). See the caveat at the top of the script: with 15
* countries and many ethnicities nested within one country, the
* two-way covariance matrix is not of full rank.
do "code/13_twoway_cluster.do"


* Romano-Wolf multiple-hypothesis correction across three outcome families.
* Slow: a 1,000-replication cluster bootstrap per family.
do "code/14_romano_wolf.do"

* --- DESCRIPTIVES AND FIGURES ---
* Summary statistics on the full country DHS samples
* (reads ../Pre Spatial Merge/ -- see build/README.md)
do "code/15_summary_stats.do"

* RD covariate balance table
do "code/16_balance_test.do"

* RD figures, all 13 outcomes plus child marriage
do "code/17_graphs.do"

di ""
di "============================================================"
di " All Stata scripts completed. LaTeX tables in results/"
di "============================================================"
