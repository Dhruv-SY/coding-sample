# Data Construction (Merge Code)

This folder holds the code that builds the analysis dataset from raw
inputs. It is **not** run by `master.do` or `master.R` — the analysis
ships with `data/20yr_threshold_12pairs.dta` already built. Run these
only to rebuild that file from scratch.

Both scripts expect the working directory to be the package root (the
folder containing `master.do`), not this folder. All paths are relative
to the package root, and the raw inputs live one level up in the project
tree (`../Women/`, `../Pre Spatial Merge/`, `../Pixel Data/`,
`../Pixel Controls/`, `../Robustness Data/`, `../Ethno Controls/`).

---

## The chain

```text
../Women/<country>/<round>/            raw DHS individual recode (IR) + GPS (GE)
        |
        |   00_select_dhs.R
        v
../Pre Spatial Merge/<pair>.dta        one file per border pair, DHS vars kept,
        |                              country/year tagged, cohorts restricted
        |   01_build_pairs.R
        v
../Pixel Data/<pair>_pixel.dta         + distance to border, pixel controls,
        |                              ethnic homeland, bride price; |rd| <= 200km
        |   02_pool_pairs.R
        v
data/20yr_threshold_12pairs.dta        the 12-pair analysis dataset
        |
        |   code/00_setup.do
        v
                                       estimation sample (dedup, nearest-border
                                       rule, outcome construction)
```

---

## 00_select_dhs.R  (Step 0)

Turns the raw DHS country folders in `../Women/` into one file per border
pair in `../Pre Spatial Merge/`. For each country in each pair it reads
every survey round's individual recode (IR `.DTA`), keeps the analysis
variables, tags `country_name` and `survey_year`, merges cluster GPS from
that round's GE shapefile, stacks the two countries and applies the cohort
restriction.

Requires: `haven`, `dplyr`, `sf`. Takes roughly half an hour — countries
appearing in several pairs are re-read for each.

### Reconstructed, not original

The original Step 0 was done interactively and no script survived. This
one was reverse-engineered from its outputs: every rule below was
recovered by comparing the existing `../Pre Spatial Merge/*.dta` files
against the raw DHS rounds, not taken from documentation. See
"Verification" below for how far that reconstruction was checked.

### Variables kept

The 81 columns in the existing pair files (`KEEP_VARS` in the script),
including the full domestic violence module `d101a`-`d130a`. Rounds that
do not carry a variable get it as missing, so the stacked file has the
same columns throughout. DHS value labels are stripped on read: rounds
label the same variable differently and stacking would otherwise fail or
silently keep one round's labels.

### Survey rounds

All rounds available under `../Women/<country>/` are used. `survey_year`
is the first year in the round's folder name, so `BJ_2011-12_DHS_...`
becomes 2011. No round is excluded by hand — rounds too early to contain
any eligible woman simply empty out under the cohort restriction. Burkina
Faso 1998 disappears from Benin-Burkina Faso this way, because a gap
starting in 2004 leaves nobody aged 15 or over in a 1998 survey.

### Cohort restriction

A woman is kept only if her marriage-risk window — from the control
country's minimum legal age (`ctrl_age`) up to 18 — overlaps the years
the legal gap was actually in force:

| Bound | Condition | Meaning |
|-------|-----------|---------|
| Lower | `v011 > cmc(gap_start, 1) - 19*12` | had not yet turned 19 when the gap opened, so still at or below the top of the risk window |
| Upper | `v011 <= cmc(gap_end, 12) - ctrl_age*12` | reached `ctrl_age` on or before the gap closed, so entered the risk window while the laws still differed |

`gap_start`, `gap_end` and `ctrl_age` are in the `PAIRS` table at the top
of the script, taken from `../border_pairs_15yr_threshold.xlsx`.

The upper bound binds only where a survey runs past the end of the gap:
Mali 2023, Zambia 2024 and Burkina Faso 2021. In the other nine pairs the
newest survey is the binding constraint instead.

Note this is an **overlap** condition, not containment: a woman's window need
only touch the gap. It therefore admits partially-exposed cohorts -- someone
who was 17 when the gap opened passed ages 15 and 16 under identical laws on
both sides. Of the pooled file, 85.7% of women are fully exposed, 10.6%
partially and 3.7% not at all. The estimation sample restricts to the fully
exposed via the `fully_exposed` flag; see `02_pool_pairs.R`.

Note the lower bound uses 19, not 18. That is what the data shows, exactly
and in all twelve pairs, so it is reproduced here — but it means a woman
who was 18 when the gap opened is kept, which is one year more generous
than a literal reading of a risk window that tops out at 18.

### Under-18 never-married women are dropped

`v012 < 18 & v511 == .` — a woman under 18 who has never married may still
marry before 18, so her outcome is censored. `code/00_setup.do` applies the
identical restriction downstream, making this belt-and-braces, but it is
what the original did and it keeps the pair files internally consistent.

### Clusters with no GPS fix are kept

DHS writes `0, 0` for clusters it could not locate. The original build kept
those rows, coordinates and all, and so does this script. They are removed
one step later in `01_build_pairs.R`: a point at 0,0 sits in the Gulf of
Guinea, falls inside no country polygon, and drops out there. Reproducing
this matters — Benin 2017 alone has 15 such clusters carrying 240 women.

### Safety guard

`OUT_DIR` is the folder the rest of the build reads from, so the script
refuses to run when pair files are already there. Set `OVERWRITE <- TRUE`
to replace them, or point `OUT_DIR` somewhere else and diff the result
against the originals.

### Verification

Re-running this script reproduces the twelve existing
`../Pre Spatial Merge/*.dta` files. Each rebuilt pair was compared against
the original row for row and column for column after sorting on
`country_name, survey_year, v001, v011, v012, v005, v021`.

One immaterial difference: for the three pairs whose DHS rounds carry
neither `v542` nor `v861` (Benin-Niger, Namibia-South Africa,
Namibia-Zimbabwe), this script emits them as all-missing columns, giving
81 columns where the original file has 79. Every shared column is
identical. It changes nothing downstream — `01_build_pairs.R` creates both
as missing when absent, and `02_pool_pairs.R` fills them with missing when
stacking, so the pooled analysis dataset is the same either way.

What this does *not* establish: that the original Step 0 was correct, only
that this script matches it. The `gap_start`, `gap_end` and `ctrl_age`
values come from a spreadsheet, not from the WORLD legal-age panel
directly, and have not been re-derived from source. The Rwanda minimum-age
question flagged earlier in the project is unresolved and would change
`gap_start` for Rwanda-Tanzania if it turns out the law did not change
when the thesis assumes.

---

## 01_build_pairs.R

Takes each pair file in `../Pre Spatial Merge/` and adds the spatial layer:

- Country assignment and signed distance to border (`rd`), computed at the
  **individual GPS point**, not the pixel centroid
- PRIO-GRID pixel controls (nightlights, population density, agricultural
  suitability, malaria ecology, capital distance, wetland share, elevation,
  travel time, rainfall, temperature, distance to coast)
- Murdock (1959) ethnic homeland (`eth_fe`), via
  `../Robustness Data/Murdock_shapefile_v2_2024/`
- D-PLACE bride price (EA006) joined onto the homeland, plus a
  `treated x bride_wealth` interaction
- Distances to colonial railways, explorer routes and 1924 missions;
  ruggedness
- Exports at a 200 km bandwidth (`bw_export`)

`ONLY_FILES` at the top restricts the loop to the 12 estimation pairs.
Clear it to rebuild all 29 pairs present in `../Pre Spatial Merge/`.

Downloads GADM country boundaries on first run, cached in `../.gadm_cache`.

Requires: `sf`, `dplyr`, `haven`, `geodata`, `readr`.

## 02_pool_pairs.R

Appends the 12 per-pair files into `data/20yr_threshold_12pairs.dta`.
A plain row append — no deduplication, no bandwidth filter, no variable
selection; all of that is downstream in `code/00_setup.do`.

The 12 pairs are those whose female minimum-legal-marriage-age gap held
for at least 20 years. Panels at 15-year (22 pairs) and 10-year (29 pairs)
thresholds were built and explored but are not the estimation sample;
they are in `../Other/`.

Verified: running this script reproduces the shipped
`data/20yr_threshold_12pairs.dta` exactly — 139,686 rows, all 139
columns identical. The pair order in `PAIRS` is load-bearing for that
byte-level match only; it has no effect on any estimate.

Three pairs (Benin-Niger, Namibia-South Africa, Namibia-Zimbabwe) come
from DHS rounds lacking `v861` and `v542`, so those two columns append as
missing for them. Neither is used in any regression in `code/`. Check
coverage before adding either to a specification.

### The fully_exposed flag

`02_pool_pairs.R` also computes `fully_exposed`, which marks women whose
entire at-risk window -- from the control country's minimum legal age to the
treated country's (15-18, 16-18, 17-18, and 15-21 for Rwanda-Tanzania) -- fell
inside the years the gap was in force. These are the cohorts that were subject
to the differential regime for the whole of their exposure.

It is a flag only: no rows are dropped in the build. `code/00_setup.do` and
all four R scripts subset on it, so the unrestricted sample remains available
for comparison. Computing it once here, rather than in each analysis script,
keeps the Stata and R samples identical by construction.

Requires: `haven`, `dplyr`.

---

## 03_greg_crosswalk.R

Maps each DHS cluster coordinate onto a GREG (Atlas Narodov Mira) ethnic group,
writing `data/greg_crosswalk.dta`. Consumed by `code/12_greg_ethnicity.do` for
the alternative-ethnic-map robustness check, which uses GREG in place of the
Murdock homelands.

It lives here rather than in `code/` because it produces a dataset, not a
result. That also removes the only ordering dependency between `master.do` and
`master.R`, which can now be run in either order.

Reads `../Robustness Data/GREG/GREG.shp`. Requires: `haven`, `dplyr`, `sf`.
