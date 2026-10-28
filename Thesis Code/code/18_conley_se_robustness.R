# ==============================================================================
# 18_conley_se_robustness.R
#
# Conley / Spatial-HAC Standard Error Robustness Check
# Reduced-Form Border RD Estimates
#
# Reduced-form specification:
#   Y = alpha + rho * treated
#       + border-pair-specific RD slopes
#       + border-pair-specific treated-side kinks
#       + age + urban/rural
#       + geographic controls
#       + border-pair FE + ethnicity FE + survey-year FE
#
# The coefficient of interest is rho: the effect of residence on the strict-law
# side of the border.
# ==============================================================================

rm(list = ls())

# ------------------------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------------------------

packages <- c(
  "haven",
  "dplyr",
  "stringr",
  "purrr",
  "tibble",
  "readr",
  "fixest"
)

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(haven)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(readr)
library(fixest)

# ------------------------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------------------------

DATA_PATH <- "data/20yr_threshold_12pairs.dta"

RESULTS_DIR <- "results"

CSV_OUT <- file.path(RESULTS_DIR, "conley_reduced_form_results.csv")
TEX_OUT <- file.path(RESULTS_DIR, "conley_reduced_form_results.tex")

# ------------------------------------------------------------------------------
# 2. Variable names
# ------------------------------------------------------------------------------

RD_VAR      <- "rd"
TREAT_VAR   <- "treated"
CLUSTER_VAR <- "v021"

# Change these manually if your GPS variables have different names.
LAT_VAR <- "lat"
LON_VAR <- "lon"

BANDWIDTHS <- c(50, 100, 200)
CONLEY_CUTOFFS <- c(50, 100, 200)

# ------------------------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------------------------

first_existing <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

se_from_vcov <- function(vc, coef_name) {
  if (is.null(vc) || !coef_name %in% rownames(vc)) return(NA_real_)
  v <- vc[coef_name, coef_name]
  if (is.na(v) || v < 0) return(NA_real_)
  sqrt(v)
}

pval_z <- function(b, se) {
  if (anyNA(c(b, se)) || se <= 0) return(NA_real_)
  2 * pnorm(-abs(b / se))
}

sig_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.01 ~ "***",
    p < 0.05 ~ "**",
    p < 0.10 ~ "*",
    TRUE ~ ""
  )
}

fmt <- function(x, d = 3) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = d))
}

fmt_n <- function(x) {
  formatC(x, format = "d", big.mark = ",")
}

clean_label <- function(x) {
  x %>%
    str_replace_all("&", "\\\\&") %>%
    str_replace_all("%", "\\\\%") %>%
    str_replace_all("\\|", "$\\\\mid$")
}

get_conley_vcov <- function(model, cutoff_km, lat_var, lon_var) {
  tryCatch(
    fixest::vcov_conley(
      model,
      lat = lat_var,
      lon = lon_var,
      cutoff = cutoff_km,
      distance = "spherical"
    ),
    error = function(e) {
      warning("Conley vcov failed at ", cutoff_km, " km: ", conditionMessage(e))
      nm <- names(coef(model))
      matrix(NA_real_, length(nm), length(nm), dimnames = list(nm, nm))
    }
  )
}

make_rf_formula <- function(outcome_var, controls, fe_str) {
  control_str <- if (length(controls) > 0) {
    paste(controls, collapse = " + ")
  } else {
    "1"
  }
  
  # Correct translation of Stata:
  #   treated
  #   c.rd#i.bp
  #   c.rd#i.bp#i.treated
  #
  # In R/fixest:
  #   treated + rd:bp + rd:bp:treated
  rhs <- paste(
    TREAT_VAR,
    paste0(RD_VAR, ":bp"),
    paste0(RD_VAR, ":bp:", TREAT_VAR),
    control_str,
    sep = " + "
  )
  
  as.formula(
    paste0(outcome_var, " ~ ", rhs, " | ", fe_str),
    env = parent.frame()
  )
}

# ------------------------------------------------------------------------------
# 4. Load data
# ------------------------------------------------------------------------------

df <- read_dta(DATA_PATH) %>%
  as.data.frame()

# Auto-detect latitude and longitude if needed
if (!LAT_VAR %in% names(df)) {
  LAT_VAR <- first_existing(df, c("lat", "latitude", "LAT", "LATNUM", "y", "gps_lat"))
}

if (!LON_VAR %in% names(df)) {
  LON_VAR <- first_existing(df, c("lon", "longitude", "LONG", "LONGNUM", "x", "gps_lon"))
}

if (is.na(LAT_VAR) || is.na(LON_VAR)) {
  stop("Could not find latitude/longitude variables. Set LAT_VAR and LON_VAR manually.")
}

message("Using latitude variable: ", LAT_VAR)
message("Using longitude variable: ", LON_VAR)

# ------------------------------------------------------------------------------
# 5. Sample restrictions and IDs
# ------------------------------------------------------------------------------

df <- df %>%
  mutate(
    bp = as.factor(country_pair),
    eth_id = as.factor(eth_fe),
    survey_year = as.factor(survey_year),
    v025 = as.factor(v025)
  ) %>%
  filter(!(v012 < 18 & is.na(v511))) %>%
  # Fully-exposed cohorts only: the whole at-risk window fell inside the
  # years the legal gap was in force. Flag built in build/02_pool_pairs.R
  # so Stata and R share one definition.
  filter(fully_exposed == 1) %>%
  filter(
    !is.na(child_marriage),
    !is.na(treated),
    !is.na(rd),
    !is.na(bp),
    !is.na(survey_year),
    !is.na(v021)
  )

# v021 (DHS PSU number) is only unique WITHIN a single country-survey round;
# it restarts from 1 every round, so raw v021 collides across different
# countries and different rounds of the same country. Build a cluster ID
# that's actually unique to one physical PSU before clustering SEs on it.
df <- df %>%
  mutate(cluster_id = as.factor(paste(country_name, survey_year, v021, sep = "_")))
CLUSTER_VAR <- "cluster_id"

# Make treated numeric 0/1 if imported as labelled/factor-like
df <- df %>%
  mutate(
    treated = as.numeric(treated),
    rd = as.numeric(rd)
  )

# A single DHS cluster can fall within the bandwidth of more than one
# border pair (e.g. eastern DRC is near Rwanda, Burundi, Uganda, and
# Zambia all at once), so the same women are duplicated across pair
# files. Keep only each cluster's comparison against its NEAREST border
# (smallest |rd|); drop its appearances in the other, more distant pairs.
df <- df %>%
  group_by(cluster_id) %>%
  filter(abs(rd) == min(abs(rd))) %>%
  ungroup()

# ------------------------------------------------------------------------------
# 6. Clean variables and construct outcomes
# ------------------------------------------------------------------------------

df <- df %>%
  mutate(
    v133 = ifelse(v133 %in% c(97, 98, 99), NA, v133),
    v106 = ifelse(v106 %in% c(9), NA, v106),
    v525 = ifelse(v525 %in% c(0, 96, 97, 98, 99), NA, v525),
    v741 = ifelse(v741 %in% c(9), NA, v741),
    v717 = ifelse(v717 %in% c(10, 96, 98, 99), NA, v717),
    
    educ_any = case_when(
      !is.na(v133) & v133 > 0 ~ 1,
      !is.na(v133) & v133 == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    
    educ_years_pos = case_when(
      !is.na(v133) & v133 > 0 ~ v133,
      TRUE ~ NA_real_
    ),
    
    age_first_sex = v525,
    
    multiple_unions = case_when(
      !is.na(v503) & v503 > 1 & v503 < 9 ~ 1,
      !is.na(v503) & v503 == 1 ~ 0,
      is.na(v503) & is.na(v511) ~ 0,
      !is.na(v503) & v503 == 9 ~ NA_real_,
      TRUE ~ NA_real_
    ),
    
    not_working = case_when(
      v714 == 0 ~ 1,
      v714 == 1 ~ 0,
      TRUE ~ NA_real_
    ),
    
    domestic_work = case_when(
      v714 == 1 & v717 == 6 ~ 1,
      v714 == 1 & !is.na(v717) & v717 != 6 ~ 0,
      TRUE ~ NA_real_
    ),
    
    paid_cash = case_when(
      v714 == 1 & v741 %in% c(1, 2) ~ 1,
      v714 == 1 & v741 %in% c(0, 3) ~ 0,
      TRUE ~ NA_real_
    ),
    
    white_collar = case_when(
      v717 %in% c(1, 2) ~ 1,
      v714 == 1 & !is.na(v717) & !(v717 %in% c(1, 2)) ~ 0,
      v714 == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    
    white_collar_cond = case_when(
      v714 == 1 & v717 %in% c(1, 2) ~ 1,
      v714 == 1 & !is.na(v717) & !(v717 %in% c(1, 2)) ~ 0,
      TRUE ~ NA_real_
    )
  )

# ------------------------------------------------------------------------------
# 7. Controls and fixed effects
# ------------------------------------------------------------------------------

possible_geo_controls <- c(
  "nightlights_mean",
  "population_density_mean",
  "agri_suitability_mean",
  "malaria_ecology_mean",
  "capdist",
  "water_wetland_pct",
  "elevation_mean",
  "traveltime_mean",
  "rainfall_mean",
  "temperature_mean",
  "dist_coast_km"
)

geo_controls <- possible_geo_controls[possible_geo_controls %in% names(df)]

controls <- c(
  "v012",
  "v025",
  geo_controls
)

FE_STR <- "bp + eth_id + survey_year"

message("Controls used: ", paste(controls, collapse = ", "))
message("Fixed effects: ", FE_STR)

# ------------------------------------------------------------------------------
# 8. Outcomes
# ------------------------------------------------------------------------------

outcomes <- tribble(
  ~tag,                    ~var,                 ~label,
  "educ_years",             "v133",               "Education (years)",
  "educ_any",               "educ_any",           "Any education",
  "educ_years_pos",         "educ_years_pos",     "Education years | schooling > 0",
  "educ_level",             "v106",               "Education level",
  "age_first_birth",        "v212",               "Age at first birth",
  "children_ever_born",     "v201",               "Children ever born",
  "age_first_sex",          "age_first_sex",      "Age at first sex",
  "multiple_unions",        "multiple_unions",    "Multiple unions",
  "not_working",            "not_working",        "Not working",
  "domestic_work",          "domestic_work",      "Domestic work",
  "paid_cash",              "paid_cash",          "Paid in cash",
  "white_collar",           "white_collar",       "White-collar employment",
  "white_collar_cond",      "white_collar_cond",  "White-collar employment among workers"
)

# Benchmarks from the current Stata reduced-form table
# (results/results_withedumargins_reducedform.tex, 2022-nightlights build).
stata_rf_benchmarks <- tribble(
  ~tag,                    ~bw, ~stata_coef,

  "educ_years",               50,  0.962,
  "educ_years",              100,  1.205,
  "educ_years",              200,  1.298,

  "educ_any",                 50,  0.122,
  "educ_any",                100,  0.161,
  "educ_any",                200,  0.158,

  "educ_level",               50,  0.282,
  "educ_level",              100,  0.334,
  "educ_level",              200,  0.334,

  "age_first_birth",          50,  0.857,
  "age_first_birth",         100,  0.931,
  "age_first_birth",         200,  0.755,

  "children_ever_born",       50,  -0.266,
  "children_ever_born",      100,  -0.330,
  "children_ever_born",      200,  -0.280,

  "age_first_sex",            50,  0.527,
  "age_first_sex",           100,  0.633,
  "age_first_sex",           200,  0.480,


  "multiple_unions",          50,  0.010,
  "multiple_unions",         100,  -0.008,
  "multiple_unions",         200,  -0.011,

  "not_working",              50,  0.103,
  "not_working",             100,  0.113,
  "not_working",             200,  0.149,

  "domestic_work",            50,  0.005,
  "domestic_work",           100,  -0.016,
  "domestic_work",           200,  -0.008,

  "paid_cash",                50,  0.045,
  "paid_cash",               100,  0.057,
  "paid_cash",               200,  0.067,

  "white_collar",             50,  0.026,
  "white_collar",            100,  0.029,
  "white_collar",            200,  0.026,
  "educ_years_pos",           50,  0.282,
  "educ_years_pos",          100,  0.444,
  "educ_years_pos",          200,  0.647,

  "white_collar_cond",        50,  0.033,
  "white_collar_cond",       100,  0.044,
  "white_collar_cond",       200,  0.044
)

# ------------------------------------------------------------------------------
# 9. Estimate reduced-form models with Conley SEs
# ------------------------------------------------------------------------------

results <- list()

for (bw in BANDWIDTHS) {
  
  message("\n==================================================")
  message("Bandwidth: ", bw, " km")
  message("==================================================")
  
  df_bw <- df %>%
    filter(abs(.data[[RD_VAR]]) <= bw)
  
  for (j in seq_len(nrow(outcomes))) {
    
    tag   <- outcomes$tag[j]
    y_var <- outcomes$var[j]
    label <- outcomes$label[j]
    
    message("Estimating: ", label)
    
    if (!y_var %in% names(df_bw)) {
      warning("Outcome variable missing: ", y_var)
      next
    }
    
    vars_needed <- c(
      y_var,
      TREAT_VAR,
      RD_VAR,
      CLUSTER_VAR,
      LAT_VAR,
      LON_VAR,
      "bp",
      "eth_id",
      "survey_year",
      controls
    )
    
    vars_needed <- unique(vars_needed)
    
    df_est <- df_bw %>%
      filter(if_all(all_of(vars_needed), ~ !is.na(.x)))
    
    N <- nrow(df_est)
    
    if (N == 0) {
      warning("No complete observations for ", label, " at ", bw, " km.")
      next
    }
    
    fml <- make_rf_formula(y_var, controls, FE_STR)
    
    model <- tryCatch(
      feols(
        fml,
        data = df_est,
        cluster = as.formula(paste0("~", CLUSTER_VAR)),
        warn = FALSE,
        notes = FALSE
      ),
      error = function(e) {
        warning("feols failed for ", label, " at ", bw, " km: ", conditionMessage(e))
        NULL
      }
    )
    
    if (is.null(model)) next
    
    coef_name <- TREAT_VAR
    
    if (!coef_name %in% names(coef(model))) {
      warning("Coefficient on treated not found for ", label, " at ", bw, " km.")
      next
    }
    
    rf_b <- unname(coef(model)[coef_name])
    
    vc_dhs <- tryCatch(
      vcov(model, cluster = as.formula(paste0("~", CLUSTER_VAR))),
      error = function(e) NULL
    )
    
    se_dhs <- se_from_vcov(vc_dhs, coef_name)
    p_dhs  <- pval_z(rf_b, se_dhs)
    
    vc_bp <- tryCatch(
      vcov(model, cluster = ~bp),
      error = function(e) NULL
    )
    
    se_bp <- se_from_vcov(vc_bp, coef_name)
    p_bp  <- pval_z(rf_b, se_bp)
    
    conley_ses <- setNames(
      rep(NA_real_, length(CONLEY_CUTOFFS)),
      paste0("se_conley_", CONLEY_CUTOFFS, "km")
    )
    
    for (cut in CONLEY_CUTOFFS) {
      vc_con <- get_conley_vcov(model, cut, LAT_VAR, LON_VAR)
      conley_ses[paste0("se_conley_", cut, "km")] <-
        se_from_vcov(vc_con, coef_name)
    }
    
    se_c100 <- conley_ses["se_conley_100km"]
    p_c100  <- pval_z(rf_b, se_c100)
    
    row <- tibble(
      outcome_tag = tag,
      outcome = label,
      outcome_var = y_var,
      bandwidth_km = bw,
      N = N,
      rf_coef = rf_b,
      se_dhs_cluster = se_dhs,
      p_dhs_cluster = p_dhs,
      se_bp_cluster = se_bp,
      p_bp_cluster = p_bp,
      se_conley_50km = conley_ses["se_conley_50km"],
      se_conley_100km = conley_ses["se_conley_100km"],
      se_conley_200km = conley_ses["se_conley_200km"],
      p_conley_100km = p_c100,
      stars_conley_100 = sig_stars(p_c100)
    )
    
    results[[length(results) + 1]] <- row
  }
}

results_df <- bind_rows(results)

# ------------------------------------------------------------------------------
# 10. Benchmark comparison with Stata reduced-form table
# ------------------------------------------------------------------------------

results_df <- results_df %>%
  left_join(stata_rf_benchmarks, by = c("outcome_tag" = "tag", "bandwidth_km" = "bw")) %>%
  mutate(
    diff_vs_stata = rf_coef - stata_coef
  )

message("\n==================================================")
message("Benchmark comparison: R reduced-form vs Stata reduced-form")
message("==================================================")

print(
  results_df %>%
    select(outcome, bandwidth_km, rf_coef, stata_coef, diff_vs_stata) %>%
    arrange(bandwidth_km, outcome)
)

# Stop if R does not match Stata reasonably well.
# You can loosen/tighten this threshold.
max_abs_diff <- max(abs(results_df$diff_vs_stata), na.rm = TRUE)

message("Maximum absolute difference vs Stata benchmark: ", round(max_abs_diff, 6))

if (max_abs_diff > 0.01) {
  warning(
    "R coefficients differ from Stata by more than 0.01. ",
    "Do not use Conley results until the specification/sample match is fixed."
  )
}

# ------------------------------------------------------------------------------
# 11. Save CSV
# ------------------------------------------------------------------------------

write_csv(results_df, CSV_OUT)
message("Saved CSV to: ", CSV_OUT)

# ------------------------------------------------------------------------------
# 12. Export LaTeX
# ------------------------------------------------------------------------------

con <- file(TEX_OUT, open = "wt")

writeLines("% Auto-generated by 18_conley_se_robustness.R", con)
writeLines("% Reduced-form estimates with Conley spatial standard errors", con)
writeLines("", con)

for (bw in BANDWIDTHS) {
  
  bw_data <- results_df %>%
    filter(bandwidth_km == bw) %>%
    mutate(
      outcome_latex = clean_label(outcome),
      coef_latex = paste0(fmt(rf_coef), stars_conley_100)
    )
  
  writeLines(c(
    "\\begin{table}[htbp]",
    "\\centering",
    paste0("\\caption{Reduced-Form Estimates with Conley Spatial Standard Errors (Bandwidth ", bw, " km)}"),
    paste0("\\label{tab:conley_rf_bw", bw, "}"),
    "\\small",
    "\\begin{tabular}{l r r r r r r r}",
    "\\hline\\hline",
    "Outcome & $N$ & RF coef. & DHS & BP & Conley 50 & Conley 100 & Conley 200 \\\\",
    "\\hline"
  ), con)
  
  for (i in seq_len(nrow(bw_data))) {
    r <- bw_data[i, ]
    
    line <- paste(
      r$outcome_latex,
      fmt_n(r$N),
      r$coef_latex,
      fmt(r$se_dhs_cluster),
      fmt(r$se_bp_cluster),
      fmt(r$se_conley_50km),
      fmt(r$se_conley_100km),
      fmt(r$se_conley_200km),
      sep = " & "
    )
    
    writeLines(paste0(line, " \\\\"), con)
  }
  
  writeLines(c(
    "\\hline",
    paste0(
      "\\multicolumn{8}{p{0.98\\linewidth}}{\\footnotesize ",
      "\\textit{Notes:} The table reports reduced-form estimates at the ",
      bw,
      " km bandwidth. The coefficient is the estimated effect of residence on the strict-law side of the border. ",
      "All regressions include border-pair fixed effects, border-pair-specific slopes and treated-side kinks in distance to the border, ",
      "ethnicity fixed effects, survey-year fixed effects, respondent age, urban/rural status, and geographic controls. ",
      "DHS denotes standard errors clustered at the DHS cluster level. BP denotes standard errors clustered at the border-pair level. ",
      "Conley 50, 100, and 200 denote Conley spatial-HAC standard errors using the corresponding distance cutoff in kilometres. ",
      "Stars are based on Conley 100 km standard errors: * $p<0.10$, ** $p<0.05$, *** $p<0.01$.}"
    ),
    "\\end{tabular}",
    "\\end{table}",
    ""
  ), con)
}

close(con)

message("Saved LaTeX to: ", TEX_OUT)

# ------------------------------------------------------------------------------
# 13. Done
# ------------------------------------------------------------------------------

message("Done.")
