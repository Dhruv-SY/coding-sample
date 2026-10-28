# Thesis Code

## Child Marriage, Legal Regimes, and Life Outcomes

**Dhruv Yadav — Pedro Silva**

---

## If you are reading this as a code sample

The analysis cannot be run without licensed DHS microdata, so three places
show the work most directly:

- **`build/`** — the data-construction chain, raw DHS survey files through to
  the estimation dataset. `00_select_dhs.R` is a reconstruction: the original
  step was done interactively and no script survived, so it was reverse-
  engineered from its outputs and verified to reproduce all twelve
  intermediate files exactly, row for row and column for column.
  `build/README.md` documents every recovered rule and what was checked.
- **`code/00_setup.do`** — the shared preamble every Stata script calls. It
  handles the non-obvious sample construction: DHS cluster identifiers that
  collide across survey rounds, clusters falling within range of more than one
  border, and the cohort restriction that defines exposure to the legal
  discontinuity.
- **Known issues**, at the bottom of this file — what is broken, what is
  disabled and why, and which parts of the build are reconstructed rather
  than original.

Specifications are duplicated across Stata and R and agree to five decimal
places; `code/18_conley_se_robustness.R` carries that cross-check as an
assertion that halts the run if the two implementations diverge.

---

## Paper

Analysis code and data construction code for the working paper *Child
Marriage, Legal Regimes, and Life Outcomes*. The paper estimates the effect
of stricter minimum-legal-marriage-age regimes on women's education,
fertility and labour market outcomes, using variation in the legal age
across international borders in Sub-Saharan Africa.

---

## Empirical strategy

A **border regression discontinuity design**. Women living near the same
international border are compared across 12 border pairs where the female
minimum legal marriage age differs between the two countries, and where
that gap was in force for at least 20 years.

Estimates are **reduced form**: the coefficient on `treated`, an indicator
for living on the stricter-law side of the border. This is an
intention-to-treat effect of the legal regime, not a 2SLS estimate of the
effect of child marriage itself. Child marriage is reported as a first
stage (`code/02_first_stage.do`) and as an outcome, but is not used to
instrument anything.

The estimation sample is restricted to **fully-exposed cohorts**: women whose
entire at-risk window (from the control country's minimum legal age up to the
treated country's) fell inside the years the two laws actually differed.
Partially-exposed cohorts -- who passed part of their teenage years under
identical law on both sides -- would dilute the treatment, so they are excluded
via the `fully_exposed` flag built in `build/02_pool_pairs.R`.

**12 border pairs:** AGO-NAM, BEN-BFA, BEN-NER, BFA-CIV, BFA-GHA, BDI-TZA,
CIV-MLI, NAM-ZAF, NAM-ZMB, NAM-ZWE, RWA-TZA, TZA-UGA

**Core specification (Stata):**

```stata
global RDBP "i.bp c.rd#i.bp c.rd#i.bp#i.treated"
global FE   "i.survey_year nightlights_mean population_density_mean ///
             agri_suitability_mean malaria_ecology_mean capdist ///
             water_wetland_pct elevation_mean traveltime_mean ///
             rainfall_mean temperature_mean dist_coast_km i.eth_id"
global X    "c.v012 i.v025"

reg Y treated $RDBP $X $FE, vce(cluster cluster_id)
```

Controls in every specification:

- Border-pair fixed effects, with pair-specific RD slopes and treated-side kinks
- Ethnicity fixed effects (Murdock 1959 homelands)
- Survey-year fixed effects
- Pixel-level geographic controls: nightlights (VIIRS 2022), population
  density, agricultural suitability, malaria ecology, distance to capital,
  wetland share, elevation, travel time, rainfall, temperature, distance
  to coast (PRIO-GRID and FAO)
- Individual controls: age, urban/rural

Region fixed effects are deliberately excluded — each of the 12 pairs sits
entirely within one region, so they are collinear with the pair fixed effects.

Standard errors are clustered on `cluster_id`, a PSU identifier built as
`group(country_name survey_year v021)`. Raw DHS `v021` restarts at 1 every
survey round and so collides across countries and rounds; clustering on it
directly would be wrong.

### A note on inference

With only 12 pairs across 15 countries, the number of clusters at the level
treatment actually varies is small, and PSU-clustered standard errors are
optimistic. Several conservative alternatives are included — two-way
ethnicity x country clustering (`13`), Conley spatial HAC at 50/100/200 km
and a Lehner data-driven cutoff (`18`, `19`), and a Romano-Wolf stepdown
correction across three outcome families (`14`). These give materially wider
intervals than the main tables. Any write-up should report them alongside
the headline results rather than in place of them.

---

## Requirements

**Stata** 15 or later, plus three SSC packages:

```stata
ssc install ivreg2      // 13_twoway_cluster.do   (two-way clustered SEs)
ssc install leebounds   // 11_lee_bounds.do       (Lee 2009 bounds)
ssc install rwolf2      // 14_romano_wolf.do      (Romano-Wolf stepdown)
```

**R** 4.0 or later:

```r
install.packages(c("haven", "dplyr", "fixest", "readr", "stringr",
                   "purrr", "tibble", "sf", "geodata", "ggplot2"))
```

`fixest` 0.10.0 or later is required for formula-based Conley spatial
standard errors.

---

## Quickstart

Set the working directory to this folder and run the two master scripts, in
either order — they share no dependencies.

Scripts are numbered in the order they run: `master.do` executes `01`-`17`,
`master.R` executes `18`-`20`. `00_setup.do` is the shared preamble every
Stata script calls, not a step of its own.

```r
setwd("path/to/Thesis Code")
source("master.R")
```

```stata
cd "path/to/Thesis Code"
do master.do
```

All outputs are written to `results/`.

The analysis dataset is **not** distributed: it is derived from DHS
microdata, which cannot be redistributed. A replicator needs their own
DHS account and must rebuild it from the raw survey files — see
`build/README.md` — those scripts are not called by either master script.

---

## Folder structure

```text
Thesis Code/
|-- README.md
|-- master.do              all Stata analysis
|-- master.R               all R analysis
|-- build/                 data construction (merge code, not run by masters)
|   |-- README.md
|   |-- 00_select_dhs.R    raw DHS         -> Pre Spatial Merge (vars, cohorts, GPS)
|   |-- 01_build_pairs.R   Pre Spatial Merge -> Pixel Data (spatial join)
|   |-- 02_pool_pairs.R    Pixel Data x12    -> the analysis dataset
|   `-- 03_greg_crosswalk.R GREG shapefile  -> data/greg_crosswalk.dta
|-- code/
|   |-- 00_setup.do        shared preamble, called by every Stata script
|   |-- 01-17              Stata analysis, in master.do order
|   `-- 18-20              R analysis, in master.R order
|-- data/
|   |-- .gitignore
|   `-- README.md
`-- results/               created when the scripts run
```

`code/00_setup.do` is the shared preamble: it loads the data, drops exact
duplicate rows, builds `cluster_id`, applies the nearest-border rule (a DHS
cluster within range of several borders is kept only against its closest
one), applies sample restrictions and constructs the outcome variables.
Every Stata script begins with `do "code/00_setup.do"` and then applies its
own filters or global overrides.

---

## Script-to-output map

### Stata (`master.do`)

| Script | Output | Description |
|--------|--------|-------------|
| `01_main_reducedform.do` | `results_withedumargins_reducedform.tex` | Main reduced-form estimates, 13 outcomes x 3 bandwidths |
| `02_first_stage.do` | `results_firststage.tex` | Effect of treatment on child marriage, matched to each outcome sample |
| `03_heterogeneity_urban_rural.do` | `results_urban_reducedform.tex`, `results_rural_reducedform.tex` | Heterogeneity by urban vs rural residence |
| `04_robustness_migration.do` | `results_migration_reducedform.tex` | Exclude women who migrated in the same period as their marriage |
| `05_robustness_no_namibia.do` | `results_nonamibia_reducedform.tex` | Exclude the 4 Namibia border pairs |
| `06_robustness_no_niger.do` | `results_noniger_reducedform.tex` | Exclude the Benin-Niger pair |
| `07_robustness_no_pixel_controls.do` | `results_nopixel_reducedform.tex` | Drop pixel controls and ethnicity FE |
| `08_robustness_rd_2ndorder.do` | `results_rd2ndorder_reducedform.tex` | 2nd-order RD polynomial |
| `09_robustness_rd_3rdorder.do` | `results_rd3rdorder_reducedform.tex` | 3rd-order RD polynomial |
| `10_robustness_leave_one_pair_out.do` | `leave_one_pair_out_reducedform.{tex,csv}` | Drop each border pair in turn (100 km) |
| `11_lee_bounds.do` | `lee_bounds.{tex,csv}` | Lee (2009) bounds for conditional outcomes |
| `12_greg_ethnicity.do` | `results_greg_reducedform.{tex,csv}` | GREG ethnic map instead of Murdock |
| `13_twoway_cluster.do` | `twoway_cluster_comparison.{tex,csv}` | Two-way ethnicity x country clustering |
| `14_romano_wolf.do` | `romano_wolf.{tex,csv}` | Multiple-hypothesis correction |
| `15_summary_stats.do` | `summary_stats_{mean,median}.{tex,csv}` | Descriptives on the full country DHS samples |
| `16_balance_test.do` | `balance_table.tex` | RD covariate balance |
| `17_graphs.do` | 14 `*_rd_plot.png` | RD plots for child marriage and all 13 outcomes |

### R (`master.R`)

| Script | Output | Description |
|--------|--------|-------------|
| `18_conley_se_robustness.R` | `conley_reduced_form_results.{tex,csv}` | Conley spatial HAC SEs at 50, 100, 200 km |
| `19_lehner_conley.R` | `conley_50_table.tex`, `conley_lehner_table.tex`, `conley_50_and_lehner_results.csv`, `conley_lehner_covariogram_diagnostics.csv` | Conley 50 km and Lehner data-driven cutoff |
| `20_sample_map.R` | `sample_map.png`, `sample_map_split_homelands.png` | Sample map, Anderson (2018) Figure 1 style |

The GREG crosswalk that `12_greg_ethnicity.do` consumes is a data-construction
step, not an analysis one: it lives in `build/03_greg_crosswalk.R`.

---

## Known issues

- **The data build is reconstructed, not original.** `build/00_select_dhs.R`
  (raw DHS to per-pair files) was reverse-engineered from its outputs after
  the original interactive Step 0 was lost. It reproduces all twelve existing
  pair files exactly, but that shows it matches what was done, not that what
  was done was right. The legal-gap years behind the cohort restriction come
  from a spreadsheet rather than from the WORLD legal-age panel directly.
  See `build/README.md`.

- **`15_summary_stats.do` and `20_sample_map.R` read outside the package**
  (`../Pre Spatial Merge/`, `../Robustness Data/`). They need the wider
  project tree, not just this folder.
