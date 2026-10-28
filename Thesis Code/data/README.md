# Data

## Files

| File | Description |
|------|-------------|
| `20yr_threshold_12pairs.dta` | The analysis dataset. 139,686 women x 139 variables, pooled across the 12 border pairs, exported at a 200 km bandwidth. Every script in `code/` reads this. |
| `greg_crosswalk.dta` | GREG (Atlas Narodov Mira) ethnic-group crosswalk for the alternative-ethnic-map robustness check. Built by `build/03_greg_crosswalk.R`; consumed by `code/12_greg_ethnicity.do`. |

## Provenance

`20yr_threshold_12pairs.dta` is derived from DHS microdata with GPS-located
clusters. It is rebuilt by `build/00_select_dhs.R`, `build/01_build_pairs.R` and
`build/02_pool_pairs.R` — see `build/README.md` for the full chain and for
the one upstream step that is not scripted.

## Redistribution

DHS microdata is licensed and may not be redistributed. If this package is
shared outside the project, **remove `20yr_threshold_12pairs.dta`** and
leave this README in its place: a replicator with their own DHS license
can rebuild it from the raw country files.
