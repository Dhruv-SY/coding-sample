# ==============================================================================
# 20_sample_map.R
# Sample map in the style of Anderson (2018) Figure 1: Africa with
# Murdock (1959) ethnic homeland boundaries overlaid, the 15 sample
# countries shaded by treatment status (stricter vs laxer female
# minimum-marriage-age law), and the 12 study borders highlighted.
# Output: results/sample_map.png
# ==============================================================================

suppressMessages({
  library(sf)
  library(dplyr)
  library(geodata)
  library(ggplot2)
})

GADM_CACHE  <- "../.gadm_cache"
MURDOCK_SHP <- "../Robustness Data/Murdock_shapefile_v2_2024/Murdock_Map_2020.shp"

sf_use_s2(FALSE)

# ------------------------------------------------------------------------------
# 1. Sample countries and roles (female-law treatment definition)
# ------------------------------------------------------------------------------

treat_countries <- c("CIV", "RWA", "BDI", "UGA", "NAM", "GHA", "BEN")
ctrl_countries  <- c("MLI", "TZA", "BFA", "ZAF", "AGO", "ZMB", "ZWE", "NER")

pairs <- list(
  c("CIV", "MLI"), c("RWA", "TZA"), c("BDI", "TZA"), c("CIV", "BFA"),
  c("UGA", "TZA"), c("NAM", "ZAF"), c("NAM", "AGO"), c("NAM", "ZMB"),
  c("NAM", "ZWE"), c("GHA", "BFA"), c("BEN", "NER"), c("BEN", "BFA")
)

# ------------------------------------------------------------------------------
# 2. Africa backdrop (coarse world boundaries cropped to the continent)
# ------------------------------------------------------------------------------

message("Loading world boundaries...")
w <- geodata::world(resolution = 3, path = GADM_CACHE)
world_sf <- st_as_sf(w) %>% st_make_valid()

african_iso <- c(
  "DZA","AGO","BEN","BWA","BFA","BDI","CMR","CPV","CAF","TCD","COM","COG",
  "COD","CIV","DJI","EGY","GNQ","ERI","SWZ","ETH","GAB","GMB","GHA","GIN",
  "GNB","KEN","LSO","LBR","LBY","MDG","MWI","MLI","MRT","MUS","MAR","MOZ",
  "NAM","NER","NGA","RWA","STP","SEN","SYC","SLE","SOM","ZAF","SSD","SDN",
  "TZA","TGO","TUN","UGA","ZMB","ZWE","ESH"
)
africa <- world_sf %>% filter(GID_0 %in% african_iso)

# ------------------------------------------------------------------------------
# 3. Sample country polygons (full-resolution GADM, already cached)
# ------------------------------------------------------------------------------

message("Loading sample-country polygons...")
get_country <- function(iso3) {
  g <- geodata::gadm(country = iso3, level = 0, path = GADM_CACHE)
  st_as_sf(g) %>% st_make_valid() %>% mutate(iso3 = iso3)
}

sample_iso <- c(treat_countries, ctrl_countries)
sample_sf <- do.call(rbind, lapply(sample_iso, get_country)) %>%
  mutate(role = ifelse(iso3 %in% treat_countries,
                       "Treated (stricter female marriage law)",
                       "Control (laxer female marriage law)"))

# ------------------------------------------------------------------------------
# 4. Study borders (shared boundary of each pair)
# ------------------------------------------------------------------------------

message("Computing study borders...")
border_list <- lapply(pairs, function(p) {
  a <- sample_sf %>% filter(iso3 == p[1])
  b <- sample_sf %>% filter(iso3 == p[2])
  br <- suppressWarnings(st_intersection(st_boundary(a), st_boundary(b)))
  br <- suppressWarnings(st_collection_extract(br, "LINESTRING"))
  br <- st_union(br)
  st_sf(pair = paste(p, collapse = "-"), geometry = br)
})
borders <- do.call(rbind, border_list)

# ------------------------------------------------------------------------------
# 5. Murdock ethnic homelands (boundary lines only)
# ------------------------------------------------------------------------------

message("Loading Murdock homelands...")
murdock <- st_read(MURDOCK_SHP, quiet = TRUE) %>% st_make_valid()

# ------------------------------------------------------------------------------
# 6. Plot
# ------------------------------------------------------------------------------

message("Drawing map...")
p <- ggplot() +
  geom_sf(data = africa, fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_sf(data = sample_sf, aes(fill = role), color = "grey40", linewidth = 0.25, alpha = 0.85) +
  geom_sf(data = murdock, fill = NA, color = "grey55", linewidth = 0.08) +
  geom_sf(data = borders, color = "black", linewidth = 0.9) +
  scale_fill_manual(values = c(
    "Treated (stricter female marriage law)" = "#4477AA",
    "Control (laxer female marriage law)"    = "#B2555B"
  )) +
  coord_sf(xlim = c(-18, 51.5), ylim = c(-35.5, 37.5), expand = FALSE) +
  labs(fill = NULL,
       caption = paste0("Thin grey lines: Murdock (1959) ethnic homeland boundaries. ",
                        "Thick black lines: the 12 study borders.")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.caption = element_text(size = 8, hjust = 0.5)
  )

ggsave("results/sample_map.png", p, width = 8, height = 8.5, dpi = 300)
message("Saved results/sample_map.png")

# ------------------------------------------------------------------------------
# 7. Panel B: split ethnic homelands (Murdock polygons straddling a study
#    border) highlighted -- the groups driving the ethnicity-FE identification
# ------------------------------------------------------------------------------

message("Identifying split homelands...")
border_union <- st_union(borders)
split_idx <- suppressWarnings(st_intersects(murdock, border_union, sparse = FALSE))[, 1]
split_homelands <- murdock[split_idx, ]
message("Split homelands: ", nrow(split_homelands), " of ", nrow(murdock))

p2 <- ggplot() +
  geom_sf(data = africa, fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_sf(data = sample_sf, aes(fill = role), color = "grey40", linewidth = 0.25, alpha = 0.45) +
  geom_sf(data = murdock, fill = NA, color = "grey60", linewidth = 0.08) +
  geom_sf(data = split_homelands, fill = "#F0C808", color = "#8C6D00",
          linewidth = 0.2, alpha = 0.75) +
  geom_sf(data = borders, color = "black", linewidth = 0.9) +
  scale_fill_manual(values = c(
    "Treated (stricter female marriage law)" = "#4477AA",
    "Control (laxer female marriage law)"    = "#B2555B"
  )) +
  coord_sf(xlim = c(-18, 51.5), ylim = c(-35.5, 37.5), expand = FALSE) +
  labs(fill = NULL,
       caption = paste0("Gold: Murdock (1959) ethnic homelands split by one of the 12 ",
                        "study borders (thick black lines). Thin grey lines: all other ",
                        "homeland boundaries.")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.caption = element_text(size = 8, hjust = 0.5)
  )

ggsave("results/sample_map_split_homelands.png", p2, width = 8, height = 8.5, dpi = 300)
message("Saved results/sample_map_split_homelands.png")
