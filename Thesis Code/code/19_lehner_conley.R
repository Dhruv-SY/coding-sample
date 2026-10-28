# ==============================================================================
# 19_lehner_conley.R
#
# Creates two appendix tables:
#
# Table A: Conley 50 km cutoff
#   - Columns: RD bandwidths 50, 100, 200 km
#   - SEs: always Conley 50 km
#
# Table B: Lehner-style data-driven cutoff
#   - Columns: RD bandwidths 50, 100, 200 km
#   - SEs: Conley SEs using residual-covariogram selected cutoff
#
# Reduced-form specification:
#   Y = alpha + rho * treated
#       + border-pair-specific RD slopes
#       + border-pair-specific treated-side kinks
#       + respondent age + urban/rural
#       + geographic controls
#       + border-pair FE + ethnicity FE + survey-year FE
# ==============================================================================

rm(list = ls())

# ------------------------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------------------------

packages <- c(
  "haven",
  "dplyr",
  "stringr",
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
library(tibble)
library(readr)
library(fixest)

# ------------------------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------------------------

DATA_PATH <- "data/20yr_threshold_12pairs.dta"

RESULTS_DIR <- "results"

CSV_OUT <- file.path(RESULTS_DIR, "conley_50_and_lehner_results.csv")
COV_OUT <- file.path(RESULTS_DIR, "conley_lehner_covariogram_diagnostics.csv")

TEX_OUT_50 <- file.path(RESULTS_DIR, "conley_50_table.tex")
TEX_OUT_LEHNER <- file.path(RESULTS_DIR, "conley_lehner_table.tex")

# ------------------------------------------------------------------------------
# 2. Variable names and settings
# ------------------------------------------------------------------------------

RD_VAR      <- "rd"
TREAT_VAR   <- "treated"
CLUSTER_VAR <- "v021"

LAT_VAR <- "gps_lat"
LON_VAR <- "gps_lon"

BANDWIDTHS <- c(50, 100, 200)

# Fixed Conley cutoff table
FIXED_CONLEY_CUTOFF <- 50

# Lehner-style residual covariogram selector settings
MAX_CUTOFF_KM <- 300
DIST_BIN_KM   <- 10
MAX_PAIRS     <- 300000

set.seed(12345)

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

clean_label <- function(x) {
  x %>%
    str_replace_all("&", "\\\\&") %>%
    str_replace_all("%", "\\\\%") %>%
    str_replace_all("\\|", "$\\\\mid$")
}

haversine_km <- function(lat1, lon1, lat2, lon2) {
  r <- 6371.0088
  to_rad <- pi / 180
  
  lat1 <- lat1 * to_rad
  lon1 <- lon1 * to_rad
  lat2 <- lat2 * to_rad
  lon2 <- lon2 * to_rad
  
  dlat <- lat2 - lat1
  dlon <- lon2 - lon1
  
  a <- sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  
  r * c
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
# 4. Lehner-style residual covariogram cutoff selector
# ------------------------------------------------------------------------------

select_covariogram_cutoff <- function(df_est,
                                      residuals,
                                      lat_var,
                                      lon_var,
                                      cluster_var,
                                      max_cutoff_km = 300,
                                      bin_width_km = 10,
                                      max_pairs = 300000) {
  
  clus <- df_est %>%
    mutate(.resid = residuals) %>%
    group_by(.data[[cluster_var]]) %>%
    summarise(
      lat = mean(.data[[lat_var]], na.rm = TRUE),
      lon = mean(.data[[lon_var]], na.rm = TRUE),
      resid = mean(.resid, na.rm = TRUE),
      n_obs = dplyr::n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(lat), !is.na(lon), !is.na(resid))
  
  n <- nrow(clus)
  
  if (n < 5) {
    return(list(
      cutoff_km = 50,
      covariogram = tibble(),
      status = "too_few_clusters_default_50"
    ))
  }
  
  clus$resid_c <- clus$resid - mean(clus$resid, na.rm = TRUE)
  
  total_pairs <- n * (n - 1) / 2
  
  if (total_pairs <= max_pairs) {
    pair_idx <- utils::combn(n, 2)
    i <- pair_idx[1, ]
    j <- pair_idx[2, ]
  } else {
    i <- sample.int(n, max_pairs, replace = TRUE)
    j <- sample.int(n, max_pairs, replace = TRUE)
    
    keep <- i != j
    i <- i[keep]
    j <- j[keep]
  }
  
  d <- haversine_km(
    clus$lat[i], clus$lon[i],
    clus$lat[j], clus$lon[j]
  )
  
  cov_prod <- clus$resid_c[i] * clus$resid_c[j]
  
  pair_df <- tibble(
    dist_km = d,
    cov_prod = cov_prod
  ) %>%
    filter(
      !is.na(dist_km),
      !is.na(cov_prod),
      dist_km <= max_cutoff_km
    ) %>%
    mutate(
      bin = floor(dist_km / bin_width_km) * bin_width_km,
      bin_mid = bin + bin_width_km / 2
    )
  
  if (nrow(pair_df) == 0) {
    return(list(
      cutoff_km = 50,
      covariogram = tibble(),
      status = "no_pairs_default_50"
    ))
  }
  
  covariogram <- pair_df %>%
    group_by(bin, bin_mid) %>%
    summarise(
      cov_mean = mean(cov_prod, na.rm = TRUE),
      n_pairs = dplyr::n(),
      .groups = "drop"
    ) %>%
    arrange(bin_mid) %>%
    mutate(
      cov_smooth = (
        dplyr::lag(cov_mean, default = first(cov_mean)) +
          cov_mean +
          dplyr::lead(cov_mean, default = last(cov_mean))
      ) / 3
    )
  
  positive_seen <- FALSE
  selected <- NA_real_
  
  for (r in seq_len(nrow(covariogram))) {
    cm <- covariogram$cov_smooth[r]
    
    if (!is.na(cm) && cm > 0) {
      positive_seen <- TRUE
    }
    
    if (positive_seen && !is.na(cm) && cm <= 0) {
      selected <- covariogram$bin_mid[r]
      break
    }
  }
  
  if (is.na(selected)) {
    selected <- min(max_cutoff_km, max(covariogram$bin_mid, na.rm = TRUE))
    status <- "no_zero_crossing_use_max_positive_range"
  } else {
    status <- "selected_first_zero_crossing"
  }
  
  selected <- max(selected, bin_width_km)
  
  list(
    cutoff_km = selected,
    covariogram = covariogram,
    status = status
  )
}

# ------------------------------------------------------------------------------
# 5. Load data
# ------------------------------------------------------------------------------

df <- read_dta(DATA_PATH) %>%
  as.data.frame()

if (!LAT_VAR %in% names(df)) {
  LAT_VAR <- first_existing(df, c("gps_lat", "lat", "latitude", "LAT", "LATNUM", "ycoord", "y", "cluster_lat"))
}

if (!LON_VAR %in% names(df)) {
  LON_VAR <- first_existing(df, c("gps_lon", "lon", "long", "longitude", "LONG", "LONGNUM", "xcoord", "x", "cluster_lon"))
}

if (is.na(LAT_VAR) || is.na(LON_VAR)) {
  cat("\nNo obvious latitude/longitude variables found. Possible coordinate variables:\n")
  print(names(df)[str_detect(tolower(names(df)), "lat|lon|long|coord|gps|xcoord|ycoord")])
  stop("Set LAT_VAR and LON_VAR manually.")
}

message("Using latitude variable: ", LAT_VAR)
message("Using longitude variable: ", LON_VAR)

# ------------------------------------------------------------------------------
# 6. Sample restrictions and IDs
# ------------------------------------------------------------------------------

df <- df %>%
  mutate(
    bp = as.factor(country_pair),
    eth_id = as.factor(eth_fe),
    survey_year = as.factor(survey_year),
    v025 = as.factor(v025),
    treated = as.numeric(treated),
    rd = as.numeric(rd)
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
# 7. Clean variables and construct outcomes
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
# 8. Controls and fixed effects
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
# 9. Outcomes
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

# ------------------------------------------------------------------------------
# 10. Estimate reduced-form models
# ------------------------------------------------------------------------------

results <- list()
covariogram_outputs <- list()

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
    
    # Fixed Conley 50 SE
    vc_50 <- get_conley_vcov(model, FIXED_CONLEY_CUTOFF, LAT_VAR, LON_VAR)
    se_50 <- se_from_vcov(vc_50, coef_name)
    p_50  <- pval_z(rf_b, se_50)
    
    # Lehner-style selected cutoff
    cutoff_obj <- select_covariogram_cutoff(
      df_est = df_est,
      residuals = resid(model),
      lat_var = LAT_VAR,
      lon_var = LON_VAR,
      cluster_var = CLUSTER_VAR,
      max_cutoff_km = MAX_CUTOFF_KM,
      bin_width_km = DIST_BIN_KM,
      max_pairs = MAX_PAIRS
    )
    
    selected_cutoff <- cutoff_obj$cutoff_km
    message("Selected cutoff: ", round(selected_cutoff, 1), " km")
    
    if (nrow(cutoff_obj$covariogram) > 0) {
      covariogram_outputs[[length(covariogram_outputs) + 1]] <-
        cutoff_obj$covariogram %>%
        mutate(
          outcome_tag = tag,
          outcome = label,
          bandwidth_km = bw,
          selected_cutoff_km = selected_cutoff,
          selector_status = cutoff_obj$status
        )
    }
    
    vc_lehner <- get_conley_vcov(model, selected_cutoff, LAT_VAR, LON_VAR)
    se_lehner <- se_from_vcov(vc_lehner, coef_name)
    p_lehner  <- pval_z(rf_b, se_lehner)
    
    row <- tibble(
      outcome_tag = tag,
      outcome = label,
      outcome_var = y_var,
      bandwidth_km = bw,
      N = N,
      rf_coef = rf_b,
      
      se_conley_50 = se_50,
      p_conley_50 = p_50,
      stars_conley_50 = sig_stars(p_50),
      
      selected_cutoff_km = selected_cutoff,
      se_lehner = se_lehner,
      p_lehner = p_lehner,
      stars_lehner = sig_stars(p_lehner),
      selector_status = cutoff_obj$status
    )
    
    results[[length(results) + 1]] <- row
  }
}

results_df <- bind_rows(results)
covariogram_df <- bind_rows(covariogram_outputs)

# ------------------------------------------------------------------------------
# 11. Save CSVs
# ------------------------------------------------------------------------------

write_csv(results_df, CSV_OUT)
message("Saved results CSV to: ", CSV_OUT)

if (nrow(covariogram_df) > 0) {
  write_csv(covariogram_df, COV_OUT)
  message("Saved covariogram diagnostics to: ", COV_OUT)
}

# ------------------------------------------------------------------------------
# 12. LaTeX export helper
# ------------------------------------------------------------------------------

write_three_bandwidth_table <- function(results_df,
                                        tex_out,
                                        caption,
                                        label,
                                        coef_star_col,
                                        se_col,
                                        note_extra) {
  
  con <- file(tex_out, open = "wt")
  
  writeLines("% Auto-generated by 19_lehner_conley.R", con)
  writeLines("", con)
  
  writeLines(c(
    "\\begin{table}[H]",
    "\\centering",
    "\\small",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\begin{tabular}{l>{\\centering\\arraybackslash}p{2.2cm}>{\\centering\\arraybackslash}p{2.2cm}>{\\centering\\arraybackslash}p{2.2cm}}",
    "\\toprule",
    " & \\multicolumn{3}{c}{Bandwidth (km)} \\\\",
    "\\cmidrule(lr){2-4}",
    "Outcome & 50 & 100 & 200 \\\\",
    "\\midrule"
  ), con)
  
  get_coef <- function(tag, bw) {
    r <- results_df %>% filter(outcome_tag == tag, bandwidth_km == bw)
    if (nrow(r) == 0) return("")
    paste0("$", fmt(r$rf_coef[1]), r[[coef_star_col]][1], "$")
  }
  
  get_se <- function(tag, bw) {
    r <- results_df %>% filter(outcome_tag == tag, bandwidth_km == bw)
    if (nrow(r) == 0) return("")
    paste0("(", fmt(r[[se_col]][1]), ")")
  }
  
  write_panel <- function(panel_title, tags) {
    writeLines(paste0("\\textit{", panel_title, "} \\\\[2pt]"), con)
    
    for (tag in tags) {
      row_data <- results_df %>% filter(outcome_tag == tag)
      if (nrow(row_data) == 0) next
      
      lab <- clean_label(row_data$outcome[1])
      
      coef_line <- paste(
        lab,
        get_coef(tag, 50),
        get_coef(tag, 100),
        get_coef(tag, 200),
        sep = " & "
      )
      
      se_line <- paste(
        "",
        get_se(tag, 50),
        get_se(tag, 100),
        get_se(tag, 200),
        sep = " & "
      )
      
      writeLines(paste0(coef_line, " \\\\"), con)
      writeLines(paste0(se_line, " \\\\[4pt]"), con)
    }
    
    writeLines("\\\\[-2pt]", con)
  }
  
  write_panel(
    "Human capital",
    c("educ_years", "educ_any", "educ_years_pos", "educ_level")
  )
  
  write_panel(
    "Fertility",
    c("age_first_birth", "children_ever_born", "living_children")
  )
  
  write_panel(
    "Sexual behaviour",
    c("age_first_sex")
  )
  
  write_panel(
    "Household and labour market",
    c(
      "multiple_unions",
      "not_working",
      "domestic_work",
      "paid_cash",
      "white_collar",
      "white_collar_cond"
    )
  )
  
  writeLines(c(
    "\\bottomrule",
    paste0(
      "\\multicolumn{4}{p{0.92\\textwidth}}{\\footnotesize ",
      "\\textit{Notes:} Reduced-form estimates. ",
      "The coefficient is the estimated effect of residence on the strict-law side of the border. ",
      "Columns report estimates using 50, 100, and 200 km bandwidths around the border cutoff. ",
      note_extra,
      " All regressions include border-pair fixed effects, border-pair-specific slopes and treated-side kinks in distance to the border, ",
      "pixel-level geographic controls, ethnicity fixed effects, survey-year fixed effects, respondent age, and urban/rural residence. ",
      "$^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$.}"
    ),
    "\\end{tabular}",
    "\\end{table}",
    ""
  ), con)
  
  close(con)
}

# ------------------------------------------------------------------------------
# 13. Export Table A: fixed Conley 50
# ------------------------------------------------------------------------------

write_three_bandwidth_table(
  results_df = results_df,
  tex_out = TEX_OUT_50,
  caption = "Reduced-Form Estimates with Conley Spatial Standard Errors: 50 km Cutoff",
  label = "tab:conley_50",
  coef_star_col = "stars_conley_50",
  se_col = "se_conley_50",
  note_extra = "Standard errors in parentheses are Conley spatial-HAC standard errors computed using a fixed 50 km spatial cutoff in all columns."
)

message("Saved Conley 50 table to: ", TEX_OUT_50)

# ------------------------------------------------------------------------------
# 14. Export Table B: Lehner-style data-driven cutoff
# ------------------------------------------------------------------------------

write_three_bandwidth_table(
  results_df = results_df,
  tex_out = TEX_OUT_LEHNER,
  caption = "Reduced-Form Estimates with Data-Driven Conley Spatial Standard Errors",
  label = "tab:conley_lehner",
  coef_star_col = "stars_lehner",
  se_col = "se_lehner",
  note_extra = "Standard errors in parentheses are Conley spatial-HAC standard errors using a data-driven spatial cutoff selected from the empirical residual covariogram. The selected cutoff is the first distance bin at which the smoothed residual covariance becomes non-positive, with residuals aggregated to the DHS-cluster level."
)

message("Saved Lehner-style table to: ", TEX_OUT_LEHNER)

# ------------------------------------------------------------------------------
# 15. Print selected cutoffs
# ------------------------------------------------------------------------------

message("\nSelected Lehner-style cutoffs:")
print(
  results_df %>%
    select(outcome, bandwidth_km, selected_cutoff_km, selector_status) %>%
    arrange(bandwidth_km, outcome)
)

message("Done.")