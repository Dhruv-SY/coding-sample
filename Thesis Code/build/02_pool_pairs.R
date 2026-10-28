# =====================================================
# 02_pool_pairs.R
# Pool the 12 per-pair files written by 01_build_pairs.R into the
# single analysis dataset used by every script in code/.
#
#   ../Pixel Data/<pair>_pixel.dta   x12   ->   data/20yr_threshold_12pairs.dta
#
# The 12 pairs are those whose female minimum-legal-marriage-age gap
# was in force for at least 20 years (the "20yr threshold" in the
# filename). Pairs at the 15yr and 10yr thresholds were explored but
# are not the estimation sample.
#
# This is a plain row append: no deduplication, no bandwidth filter,
# no variable selection. All of that happens downstream in
# code/00_setup.do (duplicate drop, nearest-border rule, sample
# restrictions). 01_build_pairs.R has already applied the 200km
# export bandwidth to each pair file.
#
# Three pairs (Benin-Niger, Namibia-South Africa, Namibia-Zimbabwe)
# were built from DHS rounds that do not carry v861 and v542. Those
# two columns append as missing for those pairs. Neither is used in
# any regression in code/, so this is harmless -- but do not add
# either to a specification without checking coverage first.
# =====================================================

suppressMessages({
  library(haven)
  library(dplyr)
})

PIXEL_DIR <- "../Pixel Data"
OUT_PATH  <- "data/20yr_threshold_12pairs.dta"

# Pair metadata. gap_start/gap_end are the years the two countries' female
# minimum legal marriage ages differed; ctrl_age and treat_age are those ages
# on the control and treated side. From border_pairs_15yr_threshold.xlsx.
PAIR_INFO <- tibble::tribble(
  ~pair,                        ~gap_start, ~gap_end, ~ctrl_age, ~treat_age,
  "Cote d'Ivoire-Mali",          1995,       2023,     15,        18,
  "Rwanda-Tanzania",             1995,       2023,     15,        21,
  "Burundi-Tanzania",            1995,       2023,     15,        18,
  "Cote d'Ivoire-Burkina Faso",  1995,       2023,     17,        18,
  "Uganda-Tanzania",             1995,       2023,     15,        18,
  "Namibia-South Africa",        1996,       2023,     15,        18,
  "Namibia-Angola",              1996,       2023,     15,        18,
  "Namibia-Zambia",              1996,       2021,     16,        18,
  "Namibia-Zimbabwe",            1996,       2021,     16,        18,
  "Ghana-Burkina Faso",          1999,       2023,     17,        18,
  "Benin-Niger",                 2004,       2023,     15,        18,
  "Benin-Burkina Faso",          2004,       2023,     17,        18
)

PAIRS <- c(
  "Cote d'Ivoire-Mali", "Rwanda-Tanzania", "Burundi-Tanzania",
  "Cote d'Ivoire-Burkina Faso", "Uganda-Tanzania", "Namibia-South Africa",
  "Namibia-Angola", "Namibia-Zambia", "Namibia-Zimbabwe",
  "Ghana-Burkina Faso", "Benin-Niger", "Benin-Burkina Faso"
)

pieces <- vector("list", length(PAIRS))

for (i in seq_along(PAIRS)) {
  f <- file.path(PIXEL_DIR, paste0(PAIRS[i], "_pixel.dta"))
  if (!file.exists(f)) {
    stop("Missing pair file: ", f, "\nRun build/01_build_pairs.R first.")
  }
  d <- read_dta(f)
  cat(sprintf("  %-28s %7d rows, %3d cols\n", PAIRS[i], nrow(d), ncol(d)))
  pieces[[i]] <- d
}

pooled <- bind_rows(pieces)

# -----------------------------------------------------
# FULLY-EXPOSED FLAG
#
# A woman's at-risk window runs from the control country's minimum legal age
# to the treated country's (e.g. 15 to 18; 15 to 21 for Rwanda-Tanzania) --
# the ages over which the two laws actually differ. fully_exposed == 1 means
# that whole window fell inside the years the gap was in force, so she was
# subject to the differential regime for all of it.
#
# The pair files themselves apply only an OVERLAP rule (the window need merely
# touch the gap), which admits partially-exposed cohorts -- a woman who was 17
# when the gap opened spent ages 15 and 16 under identical laws on both sides.
# Those women dilute the treatment. The flag is computed here, once, so the
# Stata and R analysis scripts share one definition rather than each
# re-implementing it.
#
# It is only a flag: no rows are dropped here. code/00_setup.do and the R
# scripts subset on it, so the unrestricted sample stays available.
# -----------------------------------------------------
cmc <- function(y, m = 1) (y - 1900) * 12 + m

pooled <- pooled %>%
  left_join(PAIR_INFO, by = c("country_pair" = "pair")) %>%
  mutate(
    .w0 = v011 + ctrl_age  * 12,          # enters the at-risk window
    .w1 = v011 + treat_age * 12,          # leaves it
    .g0 = cmc(gap_start, 1),              # gap opens
    .g1 = cmc(gap_end, 12),               # gap closes
    fully_exposed = as.integer(.w0 >= .g0 & .w1 <= .g1)
  ) %>%
  select(-.w0, -.w1, -.g0, -.g1)

cat(sprintf("\nfully_exposed: %d of %d rows (%.1f%%)\n",
            sum(pooled$fully_exposed, na.rm = TRUE), nrow(pooled),
            100 * mean(pooled$fully_exposed, na.rm = TRUE)))

cat(sprintf("\nPooled: %d rows, %d cols across %d pairs\n",
            nrow(pooled), ncol(pooled), length(PAIRS)))

dir.create("data", showWarnings = FALSE)
write_dta(pooled, OUT_PATH)
cat("Written:", OUT_PATH, "\n")
