# =====================================================
# GENERALIZED PIPELINE — all border pairs
#  Based on "Angola and Namibia.R" but:
#   - generalized to loop over every Pre Spatial Merge/*.dta file
#   - distance-to-border and country assignment computed at the
#     INDIVIDUAL point level (not pixel centroid)
#   - adds bride-price (D-PLACE EA006) joined onto Murdock ethnic
#     homelands, plus a treated x bride-wealth interaction
#   - keeps every DHS variable already in the Pre Spatial Merge
#     files (including the full DV module), rather than a fixed
#     hardcoded keep-list
# =====================================================
library(sf)
library(dplyr)
library(haven)
library(geodata)
library(readr)

# -----------------------------------------------------
# 0) PATHS
# -----------------------------------------------------
DHS_DIR     <- "../Pre Spatial Merge"
OUT_DIR     <- "../Pixel Data"
PRIO_GPKG   <- "../Pixel Controls/prio_grid_polygons.gpkg"
PRIO_CSV    <- "../Pixel Controls/prio_pixel_controls.csv"
MURDOCK_SHP <- "../Robustness Data/Murdock_shapefile_v2_2024/Murdock_Map_2020.shp"
DPLACE_DIR  <- "../Ethno Controls/cldf"
GADM_CACHE  <- "../.gadm_cache"

ROBUST_DIR    <- "../Robustness Data"
MISSIONS_SHP  <- file.path(ROBUST_DIR, "DATA-22Religious-Conversion-roome_map_1924/Missions.shp")
RAILWAYS_SHP  <- file.path(ROBUST_DIR, "Colonial_railways/colonial_railways.shp")
ROUTES_SHP    <- file.path(ROBUST_DIR, "Colonial_railways/_extracted/Explorer_Routes_Final.shp")
RUGGED_DTA    <- file.path(ROBUST_DIR, "rugged_data/rugged_data.dta")

bw_export <- 200

dir.create(OUT_DIR, showWarnings = FALSE)
dir.create(GADM_CACHE, showWarnings = FALSE)

# only process these pair files if set (for quick testing); NULL = all
ONLY_FILES <- c(
  "Benin-Burkina Faso.dta", "Benin-Niger.dta", "Burundi-Tanzania.dta",
  "Cote d'Ivoire-Burkina Faso.dta", "Cote d'Ivoire-Mali.dta", "Ghana-Burkina Faso.dta",
  "Namibia-Angola.dta", "Namibia-South Africa.dta", "Namibia-Zambia.dta",
  "Namibia-Zimbabwe.dta", "Rwanda-Tanzania.dta", "Uganda-Tanzania.dta"
)

# -----------------------------------------------------
# 1) country name (as used in Pre Spatial Merge filenames) -> ISO3
# -----------------------------------------------------
iso3_map <- c(
  "Cote d'Ivoire" = "CIV",
  "Mali" = "MLI",
  "Rwanda" = "RWA",
  "Tanzania" = "TZA",
  "Burundi" = "BDI",
  "Burkina Faso" = "BFA",
  "Uganda" = "UGA",
  "Namibia" = "NAM",
  "South Africa" = "ZAF",
  "Angola" = "AGO",
  "Zambia" = "ZMB",
  "Zimbabwe" = "ZWE",
  "Ghana" = "GHA",
  "Benin" = "BEN",
  "Niger" = "NER",
  "Nigeria" = "NGA",
  "The Gambia" = "GMB",
  "Senegal" = "SEN",
  "Cameroon" = "CMR",
  "Kenya" = "KEN",
  "Togo" = "TGO",
  "Liberia" = "LBR",
  "Guinea" = "GIN",
  "Democratic Republic of the Congo" = "COD",
  "Ethiopia" = "ETH"
)

parse_pair <- function(fname) {
  base <- tools::file_path_sans_ext(basename(fname))
  matched <- names(iso3_map)[sapply(names(iso3_map), function(cn) startsWith(base, cn))]
  treat <- matched[which.max(nchar(matched))]
  ctrl <- substring(base, nchar(treat) + 2)
  list(treat = treat, ctrl = ctrl)
}

# -----------------------------------------------------
# 2) GADM cache (avoid re-downloading a country used in multiple pairs)
# -----------------------------------------------------
gadm_cache <- new.env()
get_country_sf <- function(iso3) {
  if (!exists(iso3, envir = gadm_cache)) {
    g <- geodata::gadm(country = iso3, level = 0, path = GADM_CACHE)
    s <- st_as_sf(g) %>% st_make_valid()
    assign(iso3, s, envir = gadm_cache)
  }
  get(iso3, envir = gadm_cache)
}

# -----------------------------------------------------
# 3) PRIO pixel controls (load once)
# -----------------------------------------------------
cat("Loading PRIO grid + controls...\n")
prio_grid     <- st_read(PRIO_GPKG, quiet = TRUE)
prio_controls <- read_csv(PRIO_CSV, show_col_types = FALSE)
prio_full     <- prio_grid %>% dplyr::select(gid) %>% left_join(prio_controls, by = "gid")
# population_density_mean and nightlights_mean from this CSV are replaced below
# with values built from the raw 05deg_yearly rasters (the CSV's population
# variable turned out to be a units bug: raw population COUNT per cell used
# directly as if it were density, giving impossible values up to ~22M/km^2).
prio_cols     <- setdiff(names(prio_controls), c("gid", "population_density_mean", "nightlights_mean"))

# -----------------------------------------------------
# 3b) Corrected population density, nightlights, elevation, travel time
#     (raw global 0.5-degree rasters, extracted directly at each woman's
#     own coordinate rather than through the 12.5km grid)
# -----------------------------------------------------
cat("Loading raw population/nightlights/elevation/travel-time rasters...\n")
RASTER_DIR <- "../05deg_yearly"

pop_raster <- terra::unwrap(readRDS(file.path(RASTER_DIR, "ghsl_population_grid.rds")))
pop_2020   <- pop_raster[["2020-12-31"]]
cell_area_km2 <- terra::cellSize(pop_2020, unit = "km")
popdens_raster <- pop_2020 / cell_area_km2
names(popdens_raster) <- "population_density_mean"

# nightlights: harmonized DMSP-OLS/simVIIRS 2022 layer at its native
# ~1km (30 arc-second) resolution; different grid from the 0.5-degree
# rasters, so it is extracted separately below rather than stacked
nl_2022 <- terra::rast("../Pixel Controls/Harmonized_DN_NTL_2022_simVIIRS.tif")
names(nl_2022) <- "nightlights_mean"

elev_raster  <- terra::unwrap(readRDS(file.path(RASTER_DIR, "ruggedterrain_elevation_mean.rds")))
names(elev_raster) <- "elevation_mean"

travel_raster <- terra::unwrap(readRDS(file.path(RASTER_DIR, "traveltime_mean.rds")))
names(travel_raster) <- "traveltime_mean"

# rainfall/temperature: long-run climatological normal (mean across all
# available years), matching the "long-run average" character of the other
# geographic controls (malaria ecology, agri suitability) rather than a
# survey-year-specific weather snapshot
cru_pre_raster <- terra::unwrap(readRDS(file.path(RASTER_DIR, "cru_pre.rds")))
rainfall_raster <- terra::mean(cru_pre_raster, na.rm = TRUE)
names(rainfall_raster) <- "rainfall_mean"

cru_tmp_raster <- terra::unwrap(readRDS(file.path(RASTER_DIR, "cru_tmp.rds")))
temperature_raster <- terra::mean(cru_tmp_raster, na.rm = TRUE)
names(temperature_raster) <- "temperature_mean"

# distance to coast: from the Natural Earth land/sea mask (TRUE = land,
# FALSE = ocean), compute each land cell's distance to the nearest ocean cell
landsea_raster <- terra::unwrap(readRDS(file.path(RASTER_DIR, "naturalearth_cover.rds")))
ocean_mask <- terra::ifel(landsea_raster == 0, 1, NA)
dist_coast_raster <- terra::distance(ocean_mask) / 1000
names(dist_coast_raster) <- "dist_coast_km"

raw_control_stack <- c(popdens_raster, elev_raster, travel_raster,
                        rainfall_raster, temperature_raster, dist_coast_raster)

# -----------------------------------------------------
# 4) Murdock ethnic homelands + bride-price crosswalk (build once)
# -----------------------------------------------------
cat("Loading Murdock shapefile + building bride-price crosswalk...\n")
murdock <- st_read(MURDOCK_SHP, quiet = TRUE)

ea_societies <- read_csv(file.path(DPLACE_DIR, "societies.csv"), show_col_types = FALSE)
ea_data      <- read_csv(file.path(DPLACE_DIR, "data.csv"), show_col_types = FALSE)

ea006 <- ea_data %>%
  filter(Var_ID == "EA006") %>%
  transmute(Soc_ID, bride_price_code = Value)

ea_points_df <- ea006 %>%
  inner_join(ea_societies %>% dplyr::select(ID, Latitude, Longitude), by = c("Soc_ID" = "ID")) %>%
  filter(!is.na(Latitude), !is.na(Longitude))

ea_points <- st_as_sf(ea_points_df, coords = c("Longitude", "Latitude"), crs = 4326)
ea_points <- st_transform(ea_points, st_crs(murdock))

murdock_joined <- st_join(murdock, ea_points, join = st_contains, left = TRUE)

bride_xwalk <- murdock_joined %>%
  st_drop_geometry() %>%
  group_by(OBJECTID) %>%
  summarise(bride_price_code = dplyr::first(na.omit(bride_price_code)), .groups = "drop")

murdock <- murdock %>% left_join(bride_xwalk, by = "OBJECTID")
murdock$bride_price_source <- ifelse(is.na(murdock$bride_price_code), NA_character_, "direct")

n_direct <- sum(!is.na(murdock$bride_price_code))
cat(sprintf("Bride-price code matched directly for %d / %d Murdock ethnic homelands\n", n_direct, nrow(murdock)))

# ---- nearest-neighbor fill for polygons with no EA point inside them ----
# assigns each unmatched polygon its nearest matched neighbor's bride-price code
# (neighboring ethnic groups often share similar marriage-payment institutions;
#  imputed vs directly-observed is kept distinguishable via bride_price_source)
has_code    <- !is.na(murdock$bride_price_code)
donors      <- murdock[has_code, ]
recipients  <- murdock[!has_code, ]

if (nrow(recipients) > 0) {
  nearest_idx <- st_nearest_feature(recipients, donors)
  recipients$bride_price_code   <- donors$bride_price_code[nearest_idx]
  recipients$bride_price_source <- "imputed_nearest_neighbor"
  murdock <- rbind(donors, recipients)
  murdock <- murdock[order(murdock$OBJECTID), ]
}

n_matched <- sum(!is.na(murdock$bride_price_code))
cat(sprintf("Bride-price code available for %d / %d Murdock ethnic homelands after nearest-neighbor fill (%d direct, %d imputed)\n",
            n_matched, nrow(murdock), n_direct, n_matched - n_direct))

# -----------------------------------------------------
# 4b) Missions, colonial railways, explorer trade routes,
#     and country-level controls (Anderson (2018) "additional
#     controls" robustness set) -- load once
# -----------------------------------------------------
cat("Loading missions, railways, trade routes, and country-level controls...\n")
missions_all <- st_read(MISSIONS_SHP, quiet = TRUE)
protestant_missions <- missions_all %>% filter(LAYER %in% c("Mission - Protestant", "Mission - BFBS"))
catholic_missions   <- missions_all %>% filter(LAYER == "Mission - Catholic")
railways     <- st_read(RAILWAYS_SHP, quiet = TRUE)
trade_routes <- st_read(ROUTES_SHP, quiet = TRUE)

rugged_country <- read_dta(RUGGED_DTA) %>%
  transmute(
    iso3 = isocode,
    gdppc_2000 = rgdppc_2000,
    country_lat = lat,
    country_lon = lon,
    ruggedness = rugged,
    ruggedness_popw = rugged_popw,
    dist_coast_ctry = dist_coast,
    region_n = africa_region_n,
    region_s = africa_region_s,
    region_w = africa_region_w,
    region_e = africa_region_e,
    region_c = africa_region_c
  )

# -----------------------------------------------------
# 5) Per-pair processing function
# -----------------------------------------------------
process_pair <- function(dhs_path, treat_name, ctrl_name, out_path) {

  treat_iso <- iso3_map[[treat_name]]
  ctrl_iso  <- iso3_map[[ctrl_name]]

  dhs <- read_dta(dhs_path)

  dhs_pts_ll <- dhs %>%
    filter(!is.na(gps_lon), !is.na(gps_lat)) %>%
    st_as_sf(coords = c("gps_lon", "gps_lat"), crs = 4326, remove = FALSE)

  # corrected population density (GHSL count / actual cell area, 2020),
  # elevation, and travel-time-to-market, extracted directly at each
  # woman's own coordinate (not routed through the grid); nightlights
  # (harmonized 2022, ~1km) extracted separately from its finer grid
  raw_vals <- terra::extract(raw_control_stack, terra::vect(dhs_pts_ll))
  nl_vals  <- terra::extract(nl_2022, terra::vect(dhs_pts_ll))
  dhs_pts_ll$population_density_mean <- raw_vals$population_density_mean
  dhs_pts_ll$nightlights_mean        <- nl_vals$nightlights_mean
  dhs_pts_ll$elevation_mean          <- raw_vals$elevation_mean
  dhs_pts_ll$traveltime_mean         <- raw_vals$traveltime_mean
  dhs_pts_ll$rainfall_mean           <- raw_vals$rainfall_mean
  dhs_pts_ll$temperature_mean        <- raw_vals$temperature_mean
  dhs_pts_ll$dist_coast_km           <- raw_vals$dist_coast_km

  # child marriage outcome (v511 = age at first union).
  # Missing v511 means never-married. The under-18-and-never-married right-
  # censoring case was already dropped upstream (Pre Spatial Merge cohort
  # filter), so everyone remaining with v511 missing is an adult who never
  # married -- a settled "did not marry as a child" outcome, i.e. 0, not NA.
  dhs_pts_ll$child_marriage <- ifelse(is.na(dhs_pts_ll$v511), 0L,
                                       as.integer(dhs_pts_ll$v511 < 18))

  if ("v861" %in% names(dhs_pts_ll)) {
    dhs_pts_ll$hiv_pos <- case_when(
      dhs_pts_ll$v861 == 1 ~ 1,
      dhs_pts_ll$v861 == 2 ~ 0,
      TRUE ~ NA_real_
    )
  } else {
    dhs_pts_ll$hiv_pos <- NA_real_
  }

  to_num <- function(x) suppressWarnings(as.numeric(x))
  dhs_pts_ll <- dhs_pts_ll %>%
    mutate(across(any_of(c("v012","v190","v136","v133","v025","v001")), to_num))

  # ---- country polygons ----
  t_sf <- get_country_sf(treat_iso)
  c_sf <- get_country_sf(ctrl_iso)
  both_ll <- rbind(t_sf, c_sf)

  # local equal-area projection centered on this pair's extent
  bb <- st_bbox(both_ll)
  lon0 <- unname((bb["xmin"] + bb["xmax"]) / 2)
  lat0 <- unname((bb["ymin"] + bb["ymax"]) / 2)
  spread <- unname((bb["ymax"] - bb["ymin"]) / 4)
  crs_aea <- sprintf(
    "+proj=aea +lat_1=%.3f +lat_2=%.3f +lat_0=%.3f +lon_0=%.3f +datum=WGS84 +units=m +no_defs",
    lat0 + spread, lat0 - spread, lat0, lon0
  )

  both_m <- st_transform(both_ll, crs_aea)
  dhs_m  <- st_transform(dhs_pts_ll, crs_aea)

  treat_m <- st_union(both_m[1, ])
  ctrl_m  <- st_union(both_m[2, ])

  border <- st_intersection(st_boundary(treat_m), st_boundary(ctrl_m))
  border <- st_collection_extract(border, "LINESTRING")
  border <- st_union(border)

  pair_poly <- st_union(treat_m, ctrl_m)

  # ---- grid (pixel_id only, for PRIO join) ----
  grid <- st_make_grid(pair_poly, cellsize = 12500, square = TRUE)
  grid_sf <- st_sf(pixel_id = seq_along(grid), geometry = grid)
  grid_pair <- grid_sf[st_intersects(grid_sf, pair_poly, sparse = FALSE), ]

  # ---- individual-level distance & country side ----
  # country_side is assigned from the respondent's KNOWN country_name (which
  # DHS survey round she was actually interviewed in), not from geometric
  # point-in-polygon containment. DHS displaces GPS coordinates by 0-10km
  # (further in rural clusters) for privacy, which can push a respondent's
  # displaced point across the border into the neighboring country's polygon
  # even though her true country of residence is unambiguous from the survey
  # itself -- using st_within() here would silently misclassify exactly the
  # points closest to the border, which is where the whole RD design's
  # identifying variation comes from.
  dhs_m$dist_km <- as.numeric(st_distance(dhs_m, border)) / 1000
  dhs_m$country_side <- dplyr::case_when(
    dhs_m$country_name == treat_name ~ "Treat",
    dhs_m$country_name == ctrl_name  ~ "Control",
    TRUE ~ NA_character_
  )
  dhs_m$rd <- ifelse(dhs_m$country_side == "Treat", dhs_m$dist_km, -dhs_m$dist_km)

  # ---- distance to nearest mission / colonial railway / explorer trade route ----
  # (Anderson (2018) "additional controls" robustness set, sourced from
  #  Nunn & Wantchekon (2011) and Roome (1924))
  protestant_m <- st_transform(protestant_missions, crs_aea)
  catholic_m   <- st_transform(catholic_missions, crs_aea)
  railways_m   <- st_transform(railways, crs_aea)
  routes_m     <- st_transform(trade_routes, crs_aea)

  nearest_prot <- st_nearest_feature(dhs_m, protestant_m)
  dhs_m$dist_protestant_mission_km <- as.numeric(st_distance(dhs_m, protestant_m[nearest_prot, ], by_element = TRUE)) / 1000

  nearest_cath <- st_nearest_feature(dhs_m, catholic_m)
  dhs_m$dist_catholic_mission_km <- as.numeric(st_distance(dhs_m, catholic_m[nearest_cath, ], by_element = TRUE)) / 1000

  nearest_rail <- st_nearest_feature(dhs_m, railways_m)
  dhs_m$dist_colonial_railway_km <- as.numeric(st_distance(dhs_m, railways_m[nearest_rail, ], by_element = TRUE)) / 1000

  nearest_route <- st_nearest_feature(dhs_m, routes_m)
  dhs_m$dist_trade_route_km <- as.numeric(st_distance(dhs_m, routes_m[nearest_route, ], by_element = TRUE)) / 1000

  # ---- country-level controls (Nunn & Puga 2012 ruggedness replication data) ----
  dhs_m$iso3 <- unname(iso3_map[dhs_m$country_name])
  dhs_m <- dhs_m %>% left_join(rugged_country, by = "iso3")

  dhs_pix <- dhs_m %>%
    mutate(
      treated = as.integer(rd > 0),
      rd1 = rd,
      rd2 = rd^2,
      trd1 = treated * rd,
      trd2 = treated * (rd^2)
    )

  dhs_pix <- st_join(dhs_pix, grid_pair, join = st_intersects, left = TRUE) %>%
    filter(!is.na(pixel_id))

  # ---- PRIO controls ----
  prio_t <- st_transform(prio_full, st_crs(grid_pair))
  grid_pts <- st_point_on_surface(grid_pair) %>% dplyr::select(pixel_id)
  grid_pts_prio <- st_join(grid_pts, prio_t, join = st_within, left = TRUE)
  grid_prio_df <- grid_pts_prio %>%
    st_drop_geometry() %>%
    dplyr::select(pixel_id, any_of(prio_cols))
  dhs_pix <- dhs_pix %>% left_join(grid_prio_df, by = "pixel_id")

  # ---- Murdock ethnicity + bride price ----
  murdock_t <- st_transform(murdock, st_crs(dhs_pix))
  dhs_pix <- st_join(dhs_pix, murdock_t, join = st_within, left = TRUE)

  dhs_pix$eth_fe <- as.character(dhs_pix$NAME)
  dhs_pix$bride_wealth <- case_when(
    dhs_pix$bride_price_code %in% c("Bride-wealth", "Token bride-wealth") ~ 1,
    dhs_pix$bride_price_code %in% c("Bride-service", "Gift exchange", "Woman exchange",
                                     "Insignificant", "Dowry") ~ 0,
    TRUE ~ NA_real_
  )
  dhs_pix$treated_x_bridewealth <- dhs_pix$treated * dhs_pix$bride_wealth

  # ---- bandwidth filter ----
  dhs_out <- dhs_pix %>% filter(abs(rd) <= bw_export)

  d <- st_drop_geometry(dhs_out)
  d <- d[, !sapply(d, is.list), drop = FALSE]
  d$country_pair <- paste(treat_name, ctrl_name, sep = "-")
  d$country_side <- as.character(d$country_side)
  d$bride_price_code <- as.character(d$bride_price_code)
  d$bride_price_source <- as.character(d$bride_price_source)

  write_dta(d, out_path)

  cat(sprintf(
    "  -> %s : %d rows (treat=%s n=%d, ctrl=%s n=%d), bride-price matched=%.0f%% (direct=%.0f%%, imputed=%.0f%%)\n",
    basename(out_path), nrow(d), treat_name, sum(d$country_side == "Treat"),
    ctrl_name, sum(d$country_side == "Control"),
    100 * mean(!is.na(d$bride_price_code)),
    100 * mean(d$bride_price_source == "direct", na.rm = TRUE),
    100 * mean(d$bride_price_source == "imputed_nearest_neighbor", na.rm = TRUE)
  ))
}

# -----------------------------------------------------
# 6) Run over every pair
# -----------------------------------------------------
pair_files <- list.files(DHS_DIR, pattern = "\\.dta$", full.names = TRUE)
if (!is.null(ONLY_FILES)) pair_files <- pair_files[basename(pair_files) %in% ONLY_FILES]

for (f in pair_files) {
  p <- parse_pair(f)
  out_name <- paste0(tools::file_path_sans_ext(basename(f)), "_pixel.dta")
  out_path <- file.path(OUT_DIR, out_name)
  cat(sprintf("=== %s - %s ===\n", p$treat, p$ctrl))
  tryCatch({
    process_pair(f, p$treat, p$ctrl, out_path)
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
  })
}

cat("\nDone.\n")
