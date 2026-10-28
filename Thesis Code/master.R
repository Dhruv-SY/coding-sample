# ==============================================================================
# MASTER R SCRIPT
# Child Marriage, Legal Regimes, and Life Outcomes
# Dhruv Yadav and Pedro Silva
#
# Runs every R script in the analysis. Set the working directory to the
# root of this replication package (the folder containing master.R)
# before running.
#
# Order does not matter relative to master.do: the one cross-language
# dependency (data/greg_crosswalk.dta, needed by code/12_greg_ethnicity.do)
# is now produced by build/03_greg_crosswalk.R, not by this file.
# ==============================================================================

# ------------------------------------------------------------------------------
# SET WORKING DIRECTORY
# Edit this line to point to your local copy of the package
# ------------------------------------------------------------------------------
# setwd("path/to/Replication Package")

# ------------------------------------------------------------------------------
# CONFIRM DATA EXISTS
# ------------------------------------------------------------------------------
data_file <- "data/20yr_threshold_12pairs.dta"
if (!file.exists(data_file)) {
  stop(
    "Analysis dataset not found: ", data_file, "\n",
    "Place 20yr_threshold_12pairs.dta in data/, or rebuild it by running\n",
    "build/01_build_pairs.R then build/02_pool_pairs.R."
  )
}

# ------------------------------------------------------------------------------
# RUN SCRIPTS
# ------------------------------------------------------------------------------

# Conley spatial HAC SEs at 50, 100 and 200 km cutoffs
source("code/18_conley_se_robustness.R")

# Conley 50 km table and Lehner data-driven spatial cutoff
source("code/19_lehner_conley.R")



# Sample map, Anderson (2018) Figure 1 style.
# Reads ../Robustness Data/Murdock_shapefile_v2_2024/ and downloads
# GADM boundaries into ../.gadm_cache on first run.
source("code/20_sample_map.R")

cat("\n")
cat("============================================================\n")
cat(" R scripts completed. Outputs in results/\n")
cat("============================================================\n")
