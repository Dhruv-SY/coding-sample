# ==============================================================================
# 03_greg_crosswalk.R
# Assign each DHS cluster coordinate its GREG ethnic group (Weidmann, Rod &
# Cederman 2010; digitized Soviet Atlas Narodov Mira), as an ALTERNATIVE
# ethnic-homeland map to Murdock (1959), following Anderson (2018)'s GREG
# robustness check. Output: data/greg_crosswalk.dta with one row per unique
# cluster coordinate and the primary GREG group name (nearest polygon used
# for coordinates falling outside all GREG polygons, e.g. coastal slivers).
# ==============================================================================

suppressMessages({
  library(sf)
  library(dplyr)
  library(haven)
})

sf_use_s2(FALSE)

GREG_SHP <- "../Robustness Data/GREG/GREG.shp"
DATA     <- "data/20yr_threshold_12pairs.dta"

d <- read_dta(DATA)
coords <- d %>% distinct(gps_lat, gps_lon) %>% filter(!is.na(gps_lat), !is.na(gps_lon))
message("Unique coordinates: ", nrow(coords))

pts <- st_as_sf(coords, coords = c("gps_lon", "gps_lat"), crs = 4326, remove = FALSE)

greg <- st_read(GREG_SHP, quiet = TRUE) %>% st_make_valid()

joined <- st_join(pts, greg %>% select(G1SHORTNAM), join = st_within, left = TRUE)
# a point on a polygon boundary can match twice; keep the first match
joined <- joined %>% group_by(gps_lat, gps_lon) %>% slice(1) %>% ungroup()

n_miss <- sum(is.na(joined$G1SHORTNAM))
message("Coordinates outside all GREG polygons: ", n_miss, " (nearest-polygon fill)")
if (n_miss > 0) {
  miss <- joined %>% filter(is.na(G1SHORTNAM))
  nearest <- st_nearest_feature(miss, greg)
  joined$G1SHORTNAM[is.na(joined$G1SHORTNAM)] <- greg$G1SHORTNAM[nearest]
}

out <- joined %>%
  st_drop_geometry() %>%
  transmute(gps_lat, gps_lon, greg_fe = as.character(G1SHORTNAM))

message("Distinct GREG groups in sample: ", n_distinct(out$greg_fe))
write_dta(out, "data/greg_crosswalk.dta")
message("Saved data/greg_crosswalk.dta (", nrow(out), " rows)")
