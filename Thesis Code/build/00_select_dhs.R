# =====================================================
# 00_select_dhs.R
# STEP 0 of the data build: raw DHS country files -> one file per border pair.
#
#   ../Women/<country>/<round>/   ->   ../Pre Spatial Merge/<pair>.dta
#
# For every country in every border pair this script:
#   1. reads the DHS individual recode (IR) .DTA for each survey round,
#   2. keeps the analysis variables (KEEP_VARS below, including the full
#      domestic violence module), filling any the round does not carry
#      with missing,
#   3. tags country_name and survey_year,
#   4. merges cluster GPS from the round's GE shapefile on v001 = DHSCLUST,
#   5. stacks the two countries of the pair and applies the cohort
#      restriction (see below),
#   6. writes ../Pre Spatial Merge/<pair>.dta.
#
# Run from the package root (the folder containing master.do), not from
# build/. This script is NOT called by master.do or master.R.
#
# ---------------------------------------------------------------
# PROVENANCE
# ---------------------------------------------------------------
# The original Step 0 was done interactively and no script survived. This
# file was reconstructed from its outputs: the variable list, the round
# coverage and the cohort cut were all recovered by inspecting the
# existing ../Pre Spatial Merge/*.dta files and matching them back to the
# raw DHS rounds. See build/README.md for what was verified and how.
# =====================================================

suppressMessages({
  library(haven)
  library(dplyr)
  library(sf)
})

WOMEN_DIR <- "../Women"
OUT_DIR   <- "../Pre Spatial Merge"

# Writing to OUT_DIR overwrites the pair files the rest of the build reads.
# This script reproduces them exactly (see build/README.md, "Verification"),
# but overwriting the canonical inputs should still be deliberate: set
# OVERWRITE to TRUE, or point OUT_DIR somewhere else and diff the result.
OVERWRITE <- FALSE

if (!OVERWRITE && dir.exists(OUT_DIR) &&
    length(list.files(OUT_DIR, pattern = "\\.dta$")) > 0) {
  stop("\n", OUT_DIR, " already holds pair files.\n",
       "Set OVERWRITE <- TRUE to replace them, or change OUT_DIR to write ",
       "elsewhere and compare.\n")
}

# -----------------------------------------------------
# 1) BORDER PAIRS
#
# gap_start is the first year the two countries' female minimum legal
# marriage ages differed, from border_pairs_15yr_threshold.xlsx. These
# are the 12 pairs whose gap held for at least 20 years.
# -----------------------------------------------------
# gap_end is the last year of the gap; ctrl_age is the CONTROL country's
# female minimum legal marriage age with parental consent. Both feed the
# cohort restriction below.
PAIRS <- tibble::tribble(
  ~pair,                        ~country_a,       ~country_b,     ~gap_start, ~gap_end, ~ctrl_age,
  "Cote d'Ivoire-Mali",         "Cote d'Ivoire",  "Mali",          1995,       2023,     15,
  "Rwanda-Tanzania",            "Rwanda",         "Tanzania",      1995,       2023,     15,
  "Burundi-Tanzania",           "Burundi",        "Tanzania",      1995,       2023,     15,
  "Cote d'Ivoire-Burkina Faso", "Cote d'Ivoire",  "Burkina Faso",  1995,       2023,     17,
  "Uganda-Tanzania",            "Uganda",         "Tanzania",      1995,       2023,     15,
  "Namibia-South Africa",       "Namibia",        "South Africa",  1996,       2023,     15,
  "Namibia-Angola",             "Namibia",        "Angola",        1996,       2023,     15,
  "Namibia-Zambia",             "Namibia",        "Zambia",        1996,       2021,     16,
  "Namibia-Zimbabwe",           "Namibia",        "Zimbabwe",      1996,       2021,     16,
  "Ghana-Burkina Faso",         "Ghana",          "Burkina Faso",  1999,       2023,     17,
  "Benin-Niger",                "Benin",          "Niger",         2004,       2023,     15,
  "Benin-Burkina Faso",         "Benin",          "Burkina Faso",  2004,       2023,     17
)

# country_name as written into the output files -> folder under ../Women/
COUNTRY_DIR <- c(
  "Cote d'Ivoire" = "Ivory Coast", "Mali" = "Mali", "Rwanda" = "Rwanda",
  "Tanzania" = "Tanzania", "Burundi" = "Burundi", "Burkina Faso" = "Burkina Faso",
  "Uganda" = "Uganda", "Namibia" = "Namibia", "South Africa" = "South Africa",
  "Angola" = "Angola", "Zambia" = "Zambia", "Zimbabwe" = "Zimbabwe",
  "Ghana" = "Ghana", "Benin" = "Benin", "Niger" = "Niger"
)

# -----------------------------------------------------
# 2) VARIABLES KEPT
# Rounds that do not carry a variable get it as missing, so the stacked
# pair file has the same columns throughout.
# -----------------------------------------------------
KEEP_VARS <- c(
  # domestic violence module
  "d101a","d101b","d101c","d101d","d101e","d101f",
  "d105a","d105b","d105c","d105d","d105e","d105f",
  "d105g","d105h","d105i","d105j","d105k",
  "d106","d107","d108","d113","d121","d126","d130a",
  # identifiers, weights, dates
  "v001","v005","v008","v011","v012","v021",
  # residence, region, background
  "v025","v103","v104","v105","v106","v107","v130","v133","v135","v136",
  "v149","v190",
  # fertility, marriage, sexual activity
  "v201","v212","v218","v501","v503","v508","v511","v525",
  # employment
  "v704","v714","v716","v717","v731","v732","v739","v741",
  # decision making and attitudes
  "v743a","v743b","v743d","v743f",
  "v744a","v744b","v744c","v744d","v744e",
  "v745a","v745b",
  "v761","v763a","v781","v822",
  # occasionally absent
  "v542","v861"
)

cmc <- function(year, month = 1) (year - 1900) * 12 + month

# -----------------------------------------------------
# 3) LOCATE AND READ ONE SURVEY ROUND
# -----------------------------------------------------

# Round folders are named like "BJ_2011-12_DHS_07042026_1213_252852".
# survey_year is the first year in that name.
round_year <- function(round_dir) {
  as.integer(sub("^[A-Z]{2}_(\\d{4}).*$", "\\1", basename(round_dir)))
}

read_round <- function(country, round_dir) {
  ir <- list.files(round_dir, pattern = "IR.*FL\\.DTA$",
                   recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  ge <- list.files(round_dir, pattern = "GE.*FL\\.shp$",
                   recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(ir) == 0) { warning("no IR file in ", round_dir); return(NULL) }
  if (length(ge) == 0) { warning("no GE shapefile in ", round_dir); return(NULL) }

  d <- read_dta(ir[1])
  names(d) <- tolower(names(d))

  # keep what exists; create what does not, so every round has the same columns
  present <- intersect(KEEP_VARS, names(d))
  d <- d[, present, drop = FALSE]
  for (v in setdiff(KEEP_VARS, present)) d[[v]] <- NA_real_
  d <- d[, KEEP_VARS, drop = FALSE]

  # strip DHS value labels: rounds label the same variable differently and
  # bind_rows() would otherwise fail or silently keep one round's labels
  d[] <- lapply(d, function(x) { attributes(x) <- NULL; x })

  d$country_name <- country
  d$survey_year  <- round_year(round_dir)

  # cluster GPS
  g <- suppressWarnings(st_read(ge[1], quiet = TRUE))
  g <- as.data.frame(g)[, c("DHSCLUST", "LATNUM", "LONGNUM", "URBAN_RURA")]
  names(g) <- c("dhsclust", "gps_lat", "gps_lon", "urban_rura")
  g <- g[!duplicated(g$dhsclust), ]
  # NOTE: DHS writes 0,0 for clusters with no GPS fix. Those are kept here,
  # coordinates and all, because the original build kept them. They are
  # removed in Step 1 (01_build_pairs.R): a point at 0,0 sits in the Gulf
  # of Guinea, falls inside no country polygon, and drops out there.

  d <- left_join(d, g, by = c("v001" = "dhsclust"))
  d$dhsclust <- d$v001
  d <- d[!is.na(d$gps_lat) & !is.na(d$gps_lon), ]

  d
}

read_country <- function(country) {
  dir <- file.path(WOMEN_DIR, COUNTRY_DIR[[country]])
  if (!dir.exists(dir)) stop("No folder for ", country, " at ", dir)
  rounds <- list.dirs(dir, recursive = FALSE)
  out <- lapply(rounds, function(r) {
    x <- read_round(country, r)
    if (!is.null(x)) cat(sprintf("     %s %d: %d women\n",
                                 country, round_year(r), nrow(x)))
    x
  })
  bind_rows(out)
}

# -----------------------------------------------------
# 4) BUILD EACH PAIR
#
# COHORT RESTRICTION: keep a woman only if her marriage-risk window --
# from the control country's minimum legal age (ctrl_age) up to 18 --
# overlaps the years the legal gap was actually in force.
#
#   LOWER  v011 >  cmc(gap_start, 1)  - 19 * 12
#          she had not yet turned 19 when the gap opened, so she was still
#          at or below the top of the risk window.
#
#   UPPER  v011 <= cmc(gap_end, 12) - ctrl_age * 12
#          she reached ctrl_age on or before the gap closed, so she entered
#          the risk window while the laws still differed.
#
# Women outside that window were never exposed to the differential regime.
# Rounds too early to contain any eligible woman (e.g. Burkina Faso 1998
# for the pairs whose gap starts in 2004) empty out here rather than being
# excluded by hand. The upper bound binds only where a survey runs past the
# end of the gap -- Mali 2023, Zambia 2024, Burkina Faso 2021.
# -----------------------------------------------------
dir.create(OUT_DIR, showWarnings = FALSE)

for (i in seq_len(nrow(PAIRS))) {
  p <- PAIRS[i, ]
  cat(sprintf("=== %s (gap from %d) ===\n", p$pair, p$gap_start))

  d <- bind_rows(read_country(p$country_a), read_country(p$country_b))

  lo <- cmc(p$gap_start, 1) - 19 * 12
  hi <- cmc(p$gap_end, 12) - p$ctrl_age * 12
  before <- nrow(d)
  d <- d[!is.na(d$v011) & d$v011 > lo & d$v011 <= hi, ]
  cat(sprintf("   cohort cut (%d < v011 <= %d): %d -> %d women\n",
              lo, hi, before, nrow(d)))

  # Drop women under 18 who have never married (v511 missing). Their
  # marriage outcome is still censored -- they may yet marry before 18 --
  # so they carry no information on the outcome. code/00_setup.do applies
  # the identical restriction, so this is belt-and-braces, but it is what
  # the original Step 0 did and it keeps the pair files consistent.
  before <- nrow(d)
  d <- d[!(d$v012 < 18 & is.na(d$v511)), ]
  cat(sprintf("   drop under-18 never-married: %d -> %d women\n",
              before, nrow(d)))

  out <- file.path(OUT_DIR, paste0(p$pair, ".dta"))
  write_dta(d, out)
  cat(sprintf("   wrote %s (%d rows, %d cols)\n\n", out, nrow(d), ncol(d)))
}

cat("Done.\n")
