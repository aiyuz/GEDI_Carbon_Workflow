## ------------------------------------------------------------
## 03_gedi_overlay.R
## Purpose:
##   (1) Search & download GEDI L4A / L2A / L2B data for XSBN
##   (2) Read GEDI footprints and apply quality filters
##   (3) Overlay GEDI with RF classification (rf_prediction.tif)
##   (4) Summarize & plot GEDI metrics by vegetation class
##   (5) Produce comparison plots (bamboo vs. trees):
##       - canopy height (RH98, RH profile)
##       - AGBD
##       - PAI
##       - vegetation volume (PAVD)
##       - footprint overlays
## ------------------------------------------------------------

## 0. Libraries & basic setup --------------------------------------------

# install.packages(c(
#   "GEDI4R", "rGEDI", "data.table", "sf", "terra",
#   "ggplot2", "dplyr", "stringr", "cowplot",
#   "tidyr", "ggridges", "ggspatial", "tidyterra", "rhdf5"
# ))

library(GEDI4R)     # devtools::install_github("VangiElia/GEDI4R")
library(rGEDI)      # for gedifinder, gediDownload
library(data.table)
library(sf)
library(stringr)
library(dplyr)
library(tidyr)
library(rhdf5)
library(terra)
library(ggplot2)
library(ggridges)
library(purrr)
library(cowplot)
library(ggspatial)  # annotation_scale, annotation_north_arrow
library(tidyterra)  # geom_spatraster

set.seed(42)

## Set your project root directory if needed
## setwd("/path/to/your/project")

# GEDI data directories (relative to project root)
data_dir_L4  <- "GEDI/GEDI_footprint/GEDI_L4A"
data_dir_2A  <- "GEDI/GEDI_footprint/GEDI02_A"
data_dir_2B  <- "GEDI/GEDI_footprint/GEDI02_B"

dir.create(data_dir_L4, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir_2A, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir_2B, recursive = TRUE, showWarnings = FALSE)

# RF classification result from 02_classifier
predicted_raster <- rast("rf_prediction.tif")

# Class names and colors (consistent across scripts)
classes <- c("bamboo","construction","farms","plantation","rubber","trees","water")

colors <- c(
  bamboo       = "#00CD10",
  trees        = "#ac03c9",
  farm         = "#E69F00",
  plantation   = "#F0E442",
  rubber       = "#CC79A7",
  construction = "#666666",
  water        = "#1E90FF"
)

# Read RAT once and standardize column names
rat <- terra::cats(predicted_raster)[[1]]
names(rat) <- c("ID", "category")  # ID = 1..7, category = class name

## ------------------------------------------------------------
## 1. Search & download GEDI data (L4A / L2A / L2B) for XSBN
## ------------------------------------------------------------

# ---- Inputs -------------------------------------------------

# Bounding box for GEDI search: c(west, south, east, north)
xsbn_bbox  <- c(100.4, 21.3, 101.0, 22.1)  # adjust if needed
gedi_year  <- 2021
date_range <- c(
  paste0(gedi_year, "-01-01"),
  paste0(gedi_year, "-12-31")
)

# Make sure your directories exist (assumes you defined these earlier)
# data_dir_L4, data_dir_2A, data_dir_2B
dir.create(data_dir_L4, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir_2A, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir_2B, recursive = TRUE, showWarnings = FALSE)

# ---- Parse bbox for packages --------------------------------
# xsbn_bbox = c(west, south, east, north)
west  <- xsbn_bbox[1]
south <- xsbn_bbox[2]
east  <- xsbn_bbox[3]
north <- xsbn_bbox[4]

# rGEDI uses upper-left / lower-right
ul_lat <- north
ul_lon <- west
lr_lat <- south
lr_lon <- east

# =============================================================
# 1.1 L4A (GEDI04_A.002)  -- use GEDI4R
# =============================================================

message("Searching & downloading GEDI04_A.002 (L4A) for XSBN...")

gedi_files_4A <- character(0)

if (requireNamespace("GEDI4R", quietly = TRUE)) {
  
  # GEDI4R L4A downloader expects coords + from/to
  # It handles search + download in one step.
  try({
    gedi_files_4A <- GEDI4R::l4_download(
      ul_lat = ul_lat,
      lr_lat = lr_lat,
      ul_lon = ul_lon,
      lr_lon = lr_lon,
      from   = date_range[1],
      to     = date_range[2],
      outdir = data_dir_L4,
      just_path = FALSE
    )
  }, silent = TRUE)
  
} else {
  message("Package 'GEDI4R' not installed.")
  message("Install with: remotes::install_github('VangiElia/GEDI4R')")
}

# Safety fallback: list whatever exists even if GEDI4R returns NULL
if (length(gedi_files_4A) == 0 || all(is.na(gedi_files_4A))) {
  gedi_files_4A <- list.files(
    data_dir_L4,
    pattern    = "GEDI04_A_.*\\.h5$",
    full.names = TRUE
  )
}

if (length(gedi_files_4A) == 0) {
  message("No L4A files found for this bbox/year (or download failed).")
}


# =============================================================
# 1.2 L2A (GEDI02_A.002)  -- use rGEDI
# =============================================================

message("Searching & downloading GEDI02_A.002 (L2A) for XSBN...")

urls_2A <- character(0)
gedi_files_2A <- character(0)

if (requireNamespace("rGEDI", quietly = TRUE)) {
  
  # rGEDI expects base product name + version argument
  urls_2A <- rGEDI::gedifinder(
    product   = "GEDI02_A",
    ul_lat    = ul_lat,
    ul_lon    = ul_lon,
    lr_lat    = lr_lat,
    lr_lon    = lr_lon,
    version   = "002",
    daterange = date_range
  )
  
  if (length(urls_2A) > 0) {
    rGEDI::gediDownload(urls_2A, outdir = data_dir_2A)
  } else {
    message("No L2A URLs returned for this bbox/year.")
  }
  
} else {
  message("Package 'rGEDI' not installed.")
}

gedi_files_2A <- list.files(
  data_dir_2A,
  pattern    = "GEDI02_A_.*\\.h5$",
  full.names = TRUE
)

# =============================================================
# 1.3 L2B (GEDI02_B.002)  -- use rGEDI
# =============================================================

message("Searching & downloading GEDI02_B.002 (L2B) for XSBN...")

urls_2B <- character(0)
gedi_files_2B <- character(0)

if (requireNamespace("rGEDI", quietly = TRUE)) {
  
  urls_2B <- rGEDI::gedifinder(
    product   = "GEDI02_B",
    ul_lat    = ul_lat,
    ul_lon    = ul_lon,
    lr_lat    = lr_lat,
    lr_lon    = lr_lon,
    version   = "002",
    daterange = date_range
  )
  
  if (length(urls_2B) > 0) {
    rGEDI::gediDownload(urls_2B, outdir = data_dir_2B)
  } else {
    message("No L2B URLs returned for this bbox/year.")
  }
  
} else {
  message("Package 'rGEDI' not installed.")
}

gedi_files_2B <- list.files(
  data_dir_2B,
  pattern    = "GEDI02_B_.*\\.h5$",
  full.names = TRUE
)

# ---- Summary messages ---------------------------------------

message("L4A files: ", length(gedi_files_4A))
message("L2A files: ", length(gedi_files_2A))
message("L2B files: ", length(gedi_files_2B))

# =============================================================
# 1.4 Alternatively, if files are downloaded locally
# skip this if you are downloading from 1.1-1.3
# =============================================================
# ------------------------------------------------------------
# Example project structure for readers:
# ./GEDI/GEDI_footprint/GEDI_L4A
# ./GEDI/GEDI_footprint/GEDI02_A
# ./GEDI/GEDI_footprint/GEDI02_B
# ------------------------------------------------------------

data_dir_L4 <- "./GEDI/GEDI_footprint/GEDI_L4A"
data_dir_2A <- "./GEDI/GEDI_footprint/GEDI02_A"
data_dir_2B <- "./GEDI/GEDI_footprint/GEDI02_B"
gedi_year<-2021

# List all files
gedi_files_4A <- list_gedi_files(data_dir_L4, "GEDI04_A")
gedi_files_2A <- list_gedi_files(data_dir_2A, "GEDI02_A")
gedi_files_2B <- list_gedi_files(data_dir_2B, "GEDI02_B")

# Subset by year
L4A_2021 <- list_gedi_files(data_dir_L4, "GEDI04_A", 2021)
L2A_2021 <- list_gedi_files(data_dir_2A, "GEDI02_A", 2021)
L2B_2021 <- list_gedi_files(data_dir_2B, "GEDI02_B", 2021)

# Quick sanity messages (optional)
message("L4A 2021 files: ", length(L4A_2021))
message("L2A 2021 files: ", length(L2A_2021))
message("L2B 2021 files: ", length(L2B_2021))

## ------------------------------------------------------------
## 2. Level 4A: quality filtering + overlay + AGBD summary
## ------------------------------------------------------------

L4A_2021 <- gedi_files_4A[str_detect(gedi_files_4A, paste0("GEDI04_A_", gedi_year))]

# Quick merge to check files (path not used later but kept for sanity check)
gediL4_path <- l4_getmulti(L4A_2021, merge = TRUE)

# Extra variables to add
l4_extra_cols <- c(
  "land_cover_data/leaf_off_flag",
  "agbd_pi_lower",
  "agbd_pi_upper",
  "agbd"
)

gediL4 <- l4_getmulti(
  L4A_2021,
  add_col = l4_extra_cols,
  source  = TRUE
)

# Quality filters following good-practice principles
gediL4 <- gediL4 %>%
  filter(
    l4_quality_flag == 1,
    l2_quality_flag == 1,
    degrade_flag == 0,
    sensitivity   >= 0.9,
    agbd > 0,
    agbd != -9999
  )

# Spatial clip to Xishuangbanna core area (adjust as needed)
b_box <- c(100.4, 100.5, 21.53, 21.55) # c(lon_min, lon_max, lat_min, lat_max)
gediL4 <- l4_clip(gediL4, clip = b_box)

# Convert to sf point object
gedi_sf <- st_as_sf(
  gediL4,
  coords = c("lon_lowestmode", "lat_lowestmode"),
  crs    = 4326
)

# Extract raster values (class IDs) at footprint locations
vals <- terra::extract(predicted_raster, gedi_sf)
gedi_sf$land_class_id <- as.numeric(vals[, 2])

# Drop footprints outside the classified raster
gedi_sf <- gedi_sf[!is.na(gedi_sf$land_class_id), ]

# Map ID → class name using RAT
idx <- match(gedi_sf$land_class_id, rat$ID)
gedi_sf$land_class <- factor(
  rat$category[idx],
  levels = rat$category
)

# Keep vegetation classes only; drop construction & water
gedi_sf <- gedi_sf |>
  dplyr::filter(
    !is.na(land_class),
    !land_class %in% c("construction", "water"),
    !is.na(agbd)
  ) |>
  droplevels()

# Summary AGBD by class
agb_summary <- gedi_sf |>
  dplyr::group_by(land_class) |>
  dplyr::summarise(
    mean_agbd   = mean(agbd,   na.rm = TRUE),
    median_agbd = median(agbd, na.rm = TRUE),
    min_agbd    = min(agbd,    na.rm = TRUE),
    max_agbd    = max(agbd,    na.rm = TRUE),
    count       = dplyr::n(),
    .groups     = "drop"
  )
print(agb_summary)

## ------------------------------------------------------------
## 3. Level 2A & 2B: quality filtering + overlay + metrics
## ------------------------------------------------------------

### 3.1 Level 2A (height) ------------------------------------

L2A_2021 <- gedi_files_2A[str_detect(gedi_files_2A, paste0("GEDI02_A_", gedi_year))]

gedi_l2a_list <- lapply(L2A_2021, readLevel2A)
level2AM      <- lapply(gedi_l2a_list, getLevel2AM)
level2AM_df   <- rbindlist(level2AM, use.names = TRUE, fill = TRUE)

# L2A quality filtering
level2AM_df <- level2AM_df %>%
  filter(
    quality_flag == 1,
    degrade_flag == 0,
    sensitivity  >= 0.9
  )

# Convert to sf and overlay land cover
level2AM_sf <- st_as_sf(
  level2AM_df,
  coords = c("lon_lowestmode", "lat_lowestmode"),
  crs    = 4326
)

vals_L2A <- terra::extract(predicted_raster, level2AM_sf)
level2AM_df$land_class_id <- as.numeric(vals_L2A[, 2])
level2AM_df <- level2AM_df[!is.na(level2AM_df$land_class_id), ]

# Map ID → class names using RAT
idx_2a <- match(level2AM_df$land_class_id, rat$ID)
level2AM_df$land_class <- factor(
  rat$category[idx_2a],
  levels = rat$category
)

# Remove non-vegetation & missing heights
level2AM_df <- level2AM_df %>%
  filter(
    !is.na(land_class),
    !land_class %in% c("construction", "water"),
    !is.na(rh100)
  ) %>%
  droplevels()

# IQR-based filtering of extreme rh100 outliers (optional)
Q1        <- quantile(level2AM_df$rh100, 0.25, na.rm = TRUE)
Q3        <- quantile(level2AM_df$rh100, 0.75, na.rm = TRUE)
IQR_value <- Q3 - Q1
lower_b   <- Q1 - 1.5 * IQR_value
upper_b   <- Q3 + 1.5 * IQR_value

level2AM_df <- level2AM_df %>%
  filter(rh100 >= lower_b & rh100 <= upper_b) %>%
  filter(!is.na(rh98), rh98 >= 0)

# Summary RH98 by class
maxH_summary <- level2AM_df %>%
  dplyr::group_by(land_class) %>%
  dplyr::summarise(
    mean_rh98   = mean(rh98),
    median_rh98 = median(rh98),
    min_rh98    = min(rh98),
    max_rh98    = max(rh98),
    count       = dplyr::n(),
    .groups     = "drop"
  )
print(maxH_summary)
#write.csv(level2AM_df,"level2AM_df.csv")
### 3.2 Level 2B (BVPM: PAI) ---------------------------------

L2B_2021 <- gedi_files_2B[str_detect(gedi_files_2B, paste0("GEDI02_B_", gedi_year))]

# Read only files that have the 'beam' group
gedi_l2b_list <- lapply(L2B_2021, function(file) {
  file_structure <- h5ls(file)
  if ("beam" %in% file_structure$name) {
    readLevel2B(file)
  } else {
    message("Skipping ", file, " - no 'beam' group found.")
    return(NULL)
  }
})

gedi_l2b_list <- Filter(Negate(is.null), gedi_l2b_list)

# Extract BVPM (PAI profiles)
level2BVPM <- lapply(gedi_l2b_list, function(data) {
  tryCatch({
    getLevel2BVPM(data)
  }, error = function(e) {
    message("Skipping file due to error: ", e$message)
    return(NULL)
  })
})

level2BVPM <- Filter(Negate(is.null), level2BVPM)
level2BVPM_df <- data.table::rbindlist(level2BVPM, use.names = TRUE, fill = TRUE)

# Quality filters for L2B BVPM
if ("quality_flag" %in% names(level2BVPM_df)) {
  level2BVPM_df <- level2BVPM_df |> dplyr::filter(quality_flag == 1)
}
if ("degrade_flag" %in% names(level2BVPM_df)) {
  level2BVPM_df <- level2BVPM_df |> dplyr::filter(degrade_flag == 0)
}

# Replace missing sentinel and drop negative PAI
level2BVPM_df$pai[level2BVPM_df$pai == -9999] <- NA
level2BVPM_df <- level2BVPM_df |>
  dplyr::filter(!is.na(pai), pai >= 0)

#write.csv(level2BVPM_df, "level2BVPM_df.csv")

# Overlay land cover (last bin coordinates)
level2BVPM_sf <- sf::st_as_sf(
  level2BVPM_df,
  coords = c("longitude_lastbin", "latitude_lastbin"),
  crs    = "EPSG:4326"
)

vals_L2B <- terra::extract(predicted_raster, level2BVPM_sf)
level2BVPM_sf$land_class_id <- as.numeric(vals_L2B[, 2])
level2BVPM_sf <- level2BVPM_sf[!is.na(level2BVPM_sf$land_class_id), ]

idx_2b <- match(level2BVPM_sf$land_class_id, rat$ID)
level2BVPM_sf$land_class <- factor(
  rat$category[idx_2b],
  levels = rat$category
)

# keep only vegetation classes
level2BVPM_sf <- level2BVPM_sf |>
  dplyr::filter(
    !is.na(land_class),
    !land_class %in% c("construction", "water")
  ) |>
  droplevels()

# Summary PAI by class
PAI_summary <- level2BVPM_sf |>
  dplyr::group_by(land_class) |>
  dplyr::summarise(
    mean_pai   = mean(pai,   na.rm = TRUE),
    median_pai = median(pai, na.rm = TRUE),
    min_pai    = min(pai,    na.rm = TRUE),
    max_pai    = max(pai,    na.rm = TRUE),
    count      = dplyr::n(),
    .groups    = "drop"
  )
print(PAI_summary)

### 3.3 Level 2B (PAVD profiles → vegetation volume) ----------

level2BVPAD <- lapply(gedi_l2b_list, function(data) {
  tryCatch({
    getLevel2BPAVDProfile(data)
  }, error = function(e) {
    message("Skipping file in PAVD due to error: ", e$message)
    return(NULL)
  })
})

level2BVPAD <- Filter(Negate(is.null), level2BVPAD)
level2BVPAD_df <- data.table::rbindlist(level2BVPAD, use.names = TRUE, fill = TRUE)

# Replace sentinel missing values
level2BVPAD_df[level2BVPAD_df == -9999] <- NA

# Optional quality filters
if ("quality_flag" %in% names(level2BVPAD_df)) {
  level2BVPAD_df <- level2BVPAD_df |>
    dplyr::filter(quality_flag == 1)
}
if ("degrade_flag" %in% names(level2BVPAD_df)) {
  level2BVPAD_df <- level2BVPAD_df |>
    dplyr::filter(degrade_flag == 0)
}

write.csv(level2BVPAD_df,"level2BVPAD_df")

# Convert to sf for overlay with the RF classification
level2BVPAD_sf <- sf::st_as_sf(
  level2BVPAD_df,
  coords = c("lon_lowestmode", "lat_lowestmode"),
  crs    = "EPSG:4326"
)

vals_BVPAD <- terra::extract(predicted_raster, level2BVPAD_sf)
level2BVPAD_sf$land_class_id <- as.numeric(vals_BVPAD[, 2])
level2BVPAD_sf <- level2BVPAD_sf[!is.na(level2BVPAD_sf$land_class_id), ]

idx_pavd <- match(level2BVPAD_sf$land_class_id, rat$ID)
level2BVPAD_sf$land_class <- factor(
  rat$category[idx_pavd],
  levels = rat$category
)

# Keep vegetation classes only
level2BVPAD_sf <- level2BVPAD_sf |>
  dplyr::filter(
    !is.na(land_class),
    !land_class %in% c("construction", "water")
  ) |>
  droplevels()

# Identify all PAVD columns (height bins)
pavd_columns <- grep("^pavd_z", names(level2BVPAD_sf), value = TRUE)

# Remove negative PAVD values
level2BVPAD_sf <- level2BVPAD_sf |>
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(pavd_columns),
      ~ ifelse(. < 0, NA_real_, .)
    )
  )

# Total vegetation volume per footprint = sum of PAVD across height bins
level2BVPAD_sf <- level2BVPAD_sf |>
  dplyr::rowwise() |>
  dplyr::mutate(
    total_veg_volume = sum(dplyr::c_across(dplyr::all_of(pavd_columns)), na.rm = TRUE)
  ) |>
  dplyr::ungroup() |>
  dplyr::filter(!is.na(total_veg_volume), total_veg_volume > 0)

# Summary by vegetation class
veg_volume_summary <- level2BVPAD_sf |>
  dplyr::group_by(land_class) |>
  dplyr::summarise(
    mean_veg_volume   = mean(total_veg_volume,   na.rm = TRUE),
    median_veg_volume = median(total_veg_volume, na.rm = TRUE),
    min_veg_volume    = min(total_veg_volume,    na.rm = TRUE),
    max_veg_volume    = max(total_veg_volume,    na.rm = TRUE),
    count             = dplyr::n(),
    .groups           = "drop"
  )
print(veg_volume_summary)

## ------------------------------------------------------------
## 4. Combine metrics into a long table (for later use)
## ------------------------------------------------------------

selected_classes <- c("farms","plantation","bamboo","rubber","trees")

# Canopy height (L2A)
df_height <- level2AM_df %>%
  filter(
    land_class %in% selected_classes,
    !is.na(rh100)
  ) %>%
  transmute(
    land_class,
    value  = rh100,
    metric = "Canopy Height (m)"
  )

# AGBD (L4A)
df_agbd <- gedi_sf %>%
  filter(
    land_class %in% selected_classes,
    !is.na(agbd)
  ) %>%
  transmute(
    land_class,
    value  = agbd,
    metric = "Aboveground Biomass Density (Mg/ha)"
  )

# PAI (L2B BVPM)
df_pai <- level2BVPM_sf %>%
  filter(
    land_class %in% selected_classes,
    !is.na(pai)
  ) %>%
  transmute(
    land_class,
    value  = pai,
    metric = "Plant Area Index (m²/m²)"
  )

# Total vegetation volume (L2B PAVD)
df_volume <- level2BVPAD_sf %>%
  filter(
    land_class %in% selected_classes,
    !is.na(total_veg_volume)
  ) %>%
  transmute(
    land_class,
    value  = total_veg_volume,
    metric = "Total Vegetation Volume (m³/m²)"
  )

all_metrics <- bind_rows(df_height, df_agbd, df_pai, df_volume)

all_metrics$land_class <- factor(
  all_metrics$land_class,
  levels = selected_classes
)

all_metrics$metric <- factor(
  all_metrics$metric,
  levels = c(
    "Aboveground Biomass Density (Mg/ha)",
    "Canopy Height (m)",
    "Plant Area Index (m²/m²)",
    "Total Vegetation Volume (m³/m²)"
  )
)

# Save for later analyses / plotting
write.csv(all_metrics, "GEDI_all_metrics_by_class.csv", row.names = FALSE)
## ------------------------------------------------------------
## 5. One-row boxplot: bamboo vs trees across four metrics
## ------------------------------------------------------------

bt_metrics <- all_metrics %>%
  dplyr::filter(land_class %in% c("bamboo", "trees"))

cols2 <- c(
  bamboo = "#00CD10",  # green
  trees  = "#ac03c9"   # purple
)

p_bt_box <- ggplot(
  bt_metrics,
  aes(x = value,
      y = land_class,
      fill  = land_class,
      color = land_class)
) +
  geom_boxplot(
    outlier.shape = NA,
    width         = 0.5,
    alpha         = 0.7,
    linewidth     = 0.6
  ) +
  scale_fill_manual(values = cols2, name = "Class") +
  scale_color_manual(values = cols2, guide = "none") +
  facet_wrap(
    ~ metric,
    nrow   = 4,
    scales = "free_x"
  ) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "GEDI-derived canopy structure metrics"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 20, hjust = 1),
    legend.position  = "bottom",
    strip.text       = element_text(face = "bold")
  )

p_bt_box
# ggsave("Fig3_box_bamboo_vs_trees.png", p_bt_box,
#        width = 10, height = 3.5, dpi = 600)

## ------------------------------------------------------------
## 6. Summary statistics: bamboo vs trees for four metrics
## ------------------------------------------------------------

stats_summary <- all_metrics %>%
  filter(land_class %in% c("bamboo", "trees")) %>%
  group_by(metric, land_class) %>%
  summarise(
    n          = sum(!is.na(value)),
    mean_value = mean(value, na.rm = TRUE),
    median     = median(value, na.rm = TRUE),
    sd         = sd(value, na.rm = TRUE),
    min_value  = min(value,na.rm = TRUE),
    max_value  = max(value, na.rm = TRUE),
    .groups    = "drop"
  )

print(stats_summary)


## ------------------------------------------------------------
## 7. RF × GEDI overlay maps (full + zoom)
## ------------------------------------------------------------

# predicted_raster: 1–7 (1=bamboo, 6=trees)
predicted_fac <- as.factor(predicted_raster)

# Bamboo = 1, trees = 6, others = white
rf_pal <- c(
  "1" = "darkred",  # bamboo
  "2" = "white",    # construction
  "3" = "white",    # farms
  "4" = "white",    # plantation
  "5" = "white",    # rubber
  "6" = "black",    # trees
  "7" = "white"     # water
)

# GEDI bamboo / trees with valid AGBD
gedi_bt <- gedi_sf %>%
  dplyr::filter(
    land_class %in% c("bamboo", "trees"),
    !is.na(agbd)
  ) %>%
  droplevels()

# Helper to draw overlay
make_gedi_overlay <- function(pred_rast, gedi_bt, rf_pal, title, bbox_sf = NULL) {
  p <- ggplot() +
    geom_spatraster(data = pred_rast) +
    scale_fill_manual(
      name     = "Class (RF)",
      values   = rf_pal,
      breaks   = c("1", "6"),
      labels   = c("bamboo", "trees"),
      na.value = "white"
    )
  
  if (!is.null(bbox_sf)) {
    p <- p +
      geom_sf(
        data      = bbox_sf,
        fill      = NA,
        color     = "black",
        linewidth = 0.6
      )
  }
  
  p +
    geom_sf(
      data  = gedi_bt,
      aes(color = agbd),
      size  = 1.2,
      alpha = 0.8
    ) +
    scale_color_viridis_c(
      name      = "AGBD (Mg/ha)",
      option    = "plasma",
      direction = -1
    ) +
    annotation_scale(
      location    = "bl",
      style       = "bar",
      width_hint  = 0.2
    ) +
    annotation_north_arrow(
      location    = "tl",
      which_north = "true",
      style       = north_arrow_fancy_orienteering
    ) +
    coord_sf(crs = terra::crs(pred_rast)) +
    labs(
      title = title,
      x     = NULL,
      y     = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position  = "right"
    )
}

# 7.1 Full map
p_gedi_full <- make_gedi_overlay(
  pred_rast = predicted_fac,
  gedi_bt   = gedi_bt,
  rf_pal    = rf_pal,
  title     = "GEDI AGBD footprints over RF bamboo / trees (Xishuangbanna)"
)

# 7.2 Zoom around GEDI beams
gedi_bt_proj <- sf::st_transform(gedi_bt, crs = terra::crs(predicted_fac))
bbox_raw     <- sf::st_bbox(gedi_bt_proj)
bbox_sfc     <- sf::st_as_sfc(bbox_raw)
bbox_buf     <- sf::st_buffer(bbox_sfc, dist = 500)   # 500 m buffer
bbox_vect    <- terra::vect(bbox_buf)

predicted_masked <- terra::mask(predicted_fac, bbox_vect)
predicted_carved <- terra::crop(predicted_masked, bbox_vect)

p_gedi_zoom <- make_gedi_overlay(
  pred_rast = predicted_carved,
  gedi_bt   = gedi_bt,
  rf_pal    = rf_pal,
  title     = "GEDI AGBD footprints over RF bamboo / trees (zoomed)",
  bbox_sf   = bbox_buf
)

p_gedi_full
p_gedi_zoom

## ------------------------------------------------------------
## 8. Additional structural difference plots
## ------------------------------------------------------------

### 8.1 AGBD vs local tree cover (bamboo + trees) ------------

df_tc <- gedi_sf %>%
  dplyr::filter(
    land_class %in% c("bamboo", "trees"),
    !is.na(tree_cover),
    !is.na(agbd)
  ) %>%
  dplyr::filter(
    tree_cover >= 0, tree_cover <= 100,
    agbd > 0,
    agbd < quantile(agbd, 0.99, na.rm = TRUE)
  )

p_tc_agbd <- ggplot(df_tc,
                    aes(x = tree_cover, y = agbd, color = land_class)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = cols2, name = "Class") +
  labs(
    x     = "Local tree cover (%) defined by GEDI",
    y     = "Aboveground biomass density (Mg/ha)",
    title = "AGBD vs. local tree cover (GEDI L4A footprints)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black")
  )

p_tc_agbd

### 8.1.1 Show that bamboo has systematically lower tree cover + stats test

# df_tc should include:
#   - land_class (factor with levels including "bamboo", "trees")
#   - tree_cover (0-100)
#   - agbd (>0)
# Optional: ensure order
df_tc <- df_tc %>%
  mutate(
    land_class = factor(land_class, levels = c("bamboo", "trees"))
  )

# ------------------------------------------------------------
# (A) ANCOVA: tree_cover ~ land_class * agbd
#     Tests whether bamboo has lower tree cover after
#     accounting for AGBD
# ------------------------------------------------------------

m_tc <- lm(tree_cover ~ land_class * agbd, data = df_tc)
a_tc <- anova(m_tc)

p_land <- a_tc["land_class", "Pr(>F)"]
p_agbd <- a_tc["agbd", "Pr(>F)"]
p_int  <- a_tc["land_class:agbd", "Pr(>F)"]

lab_ancova <- paste0(
  "ANCOVA: tree_cover ~ class * AGBD",
  "\nclass p = ", signif(p_land, 3),
  "\nAGBD p = ", signif(p_agbd, 3),
  "\ninteraction p = ", signif(p_int, 3)
)

line_table <- df_tc %>%
  group_by(land_class) %>%
  summarise(
    intercept = coef(lm(agbd ~ tree_cover))[1],
    slope     = coef(lm(agbd ~ tree_cover))[2],
    .groups   = "drop"
  )

line_table # slopes and intercepts for lm fits

# ------------------------------------------------------------
# Figure 8a: Scatter + class-specific linear fits + ANCOVA text
# ------------------------------------------------------------

p_tc_agbd <- ggplot(
  df_tc,
  aes(x = tree_cover, y = agbd, color = land_class)
) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  scale_color_manual(values = cols2, name = "Class") +
  labs(
    x     = "Local tree cover (%)",
    y     = "Aboveground biomass density (Mg/ha)",
    title = "AGBD vs. local tree cover (GEDI L4A footprints)"
  ) +
  #  annotate(
  #    "text",
  #    x = Inf, y = Inf,
  #    label = lab_ancova,
  #    hjust = 1.05, vjust = 1.15,
  #    size = 4, color = "black"
  #  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black")
  )

p_tc_agbd

# ------------------------------------------------------------
# (B) Residual-based robustness test
#     Fit pooled model tree_cover ~ agbd,
#     then test whether bamboo residuals are more negative
# ------------------------------------------------------------

m0 <- lm(tree_cover ~ agbd, data = df_tc)

df_resid <- df_tc %>%
  mutate(resid_tc = resid(m0))

p_t <- t.test(resid_tc ~ land_class, data = df_resid)$p.value
p_w <- wilcox.test(resid_tc ~ land_class, data = df_resid)$p.value

lab_resid <- paste0(
  "Residuals from pooled model: tree_cover ~ AGBD",
  "\nWelch t p = ", signif(p_t, 3),
  "\nWilcoxon p = ", signif(p_w, 3)
)

# ------------------------------------------------------------
# Figure 8b: Residual violin + box + jitter + p-values
# ------------------------------------------------------------

p_resid <- ggplot(
  df_resid,
  aes(y = land_class, x = resid_tc, fill = land_class, color = land_class)
) +
  geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.5, linewidth = 1) +
  geom_jitter(width = 0.1, alpha = 0.7, size = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7) +
  scale_fill_manual(values = cols2, name = "Class") +
  scale_color_manual(values = cols2, guide = "none") +
  labs(
    x     = NULL,
    y     = "Residual tree cover (%)",
    title = "Bamboo shows systematically lower\nGEDI tree cover at comparable AGBD"
  ) +
  #  annotate(
  #    "text",
  #    x = 1.5,
  #    y = max(df_resid$resid_tc, na.rm = TRUE),
  #    label = lab_resid,
  #    vjust = 1.2,
  #    size  = 4,
  #    color = "black"
  #  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black"),
    legend.position  = "none"
  )

p_resid

### 8.2 RH percentile profiles (bamboo vs trees) ------------- 

# Define RH columns: 0–100
rh_cols_all <- grep("^rh[0-9]+$", names(level2AM_df), value = TRUE)
rh_perc_all <- as.numeric(sub("rh", "", rh_cols_all))
keep_idx    <- rh_perc_all >= 0 & rh_perc_all <= 100
rh_cols     <- rh_cols_all[keep_idx]

max_n <- 5000  # max footprints per class

# Remove negative RH values
level2AM_df <- level2AM_df %>%
  mutate(across(all_of(rh_cols), ~ ifelse(. < 0, NA_real_, .)))

# Wide format & subsample
profiles_wide <- level2AM_df %>%
  filter(
    land_class %in% c("bamboo", "trees"),
    !is.na(rh100),
    rh100 >= 0
  ) %>%
  group_by(land_class) %>%
  group_modify(~ dplyr::slice_sample(.x, n = min(nrow(.x), max_n))) %>%
  ungroup() %>%
  mutate(foot_id = row_number()) %>%
  dplyr::select(foot_id, land_class, dplyr::all_of(rh_cols))

# Long format, drop degenerate profiles
profiles_long <- profiles_wide %>%
  pivot_longer(
    cols      = all_of(rh_cols),
    names_to  = "rh_bin",
    values_to = "height"
  ) %>%
  mutate(
    rh = as.numeric(sub("rh", "", rh_bin))
  ) %>%
  filter(!is.na(height), height >= 0) %>%
  group_by(land_class, foot_id) %>%
  filter(max(height, na.rm = TRUE) > 0.5) %>%  # drop tiny/no-canopy footprints
  ungroup()

p_rh_profile <- ggplot() +
  geom_line(
    data  = profiles_long,
    aes(x = rh, y = height,
        group = interaction(land_class, foot_id),
        color = land_class),
    alpha     = 0.05,
    linewidth = 0.3
  ) +
  scale_color_manual(values = cols2, guide = "none") +
  labs(
    x     = "Relative Height Percentile (RH)",
    y     = "Height (m)",
    title = "GEDI RH profiles (RH0–RH100)\nAll bamboo & tree footprints"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black")
  )

p_rh_profile

# calculate each land_class's mean + IQR（25–75%）with each RH percentile
profile_summ <- profiles_long %>%
  group_by(land_class, rh) %>%
  summarise(
    med_h = median(height, na.rm = TRUE),
    q25   = quantile(height, 0.25, na.rm = TRUE),
    q75   = quantile(height, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

p_rh_profile_median <- ggplot(profile_summ, aes(x = rh)) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = land_class),
              alpha = 0.5, color = NA) +
  geom_line(aes(y = med_h, color = land_class),
            linewidth = 1.2) +
  scale_color_manual(values = cols2, guide = "none") +
  scale_fill_manual(values = cols2, guide = "none") +
  labs(
    x = "Relative Height Percentile (RH)",
    y = "Height (m)",
    title = "GEDI RH profiles (RH0–RH100)",
    subtitle = "Solid line: median profile; shaded band: interquartile range (25–75%)"
  ) +
  coord_cartesian(ylim = c(0, 50)) +
  scale_y_continuous(breaks = seq(0, 50, by = 10))+
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black")
  )

p_rh_profile_median

### 8.3 PAVD vertical distribution (pyramid: bamboo vs trees) -----------

# Use previously defined pavd_columns, parse height bins
pavd_bounds       <- gsub("pavd_z|m", "", pavd_columns)     # "0_5", "5_10", ...
pavd_bounds_split <- strsplit(pavd_bounds, "_", fixed = TRUE)

lower_z  <- as.numeric(vapply(pavd_bounds_split, `[`, 1, FUN.VALUE = "0"))
upper_z  <- as.numeric(vapply(pavd_bounds_split, `[`, 2, FUN.VALUE = "0"))
bin_depth <- upper_z - lower_z
z_mid    <- (lower_z + upper_z) / 2                        # bin mid-height (m)

# Long format for bamboo & trees
pavd_long <- level2BVPAD_sf |>
  sf::st_drop_geometry() |>
  dplyr::filter(land_class %in% c("bamboo", "trees")) |>
  dplyr::select(land_class, dplyr::all_of(pavd_columns)) |>
  tidyr::pivot_longer(
    cols      = dplyr::all_of(pavd_columns),
    names_to  = "z_bin",
    values_to = "pavd"
  ) |>
  dplyr::mutate(
    height_mid = z_mid[match(z_bin, pavd_columns)],
    depth_m    = bin_depth[match(z_bin, pavd_columns)],
    pavd       = ifelse(pavd < 0, NA_real_, pavd)
  ) |>
  dplyr::filter(!is.na(pavd))

pavd_hist <- pavd_long |>
  dplyr::group_by(land_class, height_mid) |>
  dplyr::summarise(
    mean_pavd = mean(pavd, na.rm = TRUE),
    .groups   = "drop"
  ) |>
  dplyr::filter(height_mid <= 77.5) |>
  dplyr::mutate(
    height_bin = factor(
      height_mid,
      levels = sort(unique(height_mid))
    ),
    mean_pavd_signed = dplyr::case_when(
      land_class == "bamboo" ~ -mean_pavd,
      land_class == "trees"  ~  mean_pavd,
      TRUE                   ~  NA_real_
    )
  )

p_pavd_pyramid <- ggplot(
  pavd_hist,
  aes(x = mean_pavd_signed, y = height_bin, fill = land_class)
) +
  geom_col(alpha = 0.8) +
  # geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = cols2, guide = "none") +
  scale_x_continuous(
    name   = "Mean PAVD per height bin (m²/m³)",
    labels = function(x) abs(x)
  ) +
  labs(
    y     = "Height bin (m)",
    title = "Vertical distribution of PAVD\nBamboo (left) vs trees (right)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid   = element_blank(),
    axis.title.y = element_text(),
    axis.title.x = element_text(),
    axis.text.y  = element_text(),
    axis.text.x  = element_text(),
    plot.title   = element_text(face = "bold")
  )

p_pavd_pyramid

### 8.4 Vertical distribution of PAI contribution (pyramid) -----------

pai_long <- level2BVPAD_sf |>
  sf::st_drop_geometry() |>
  filter(land_class %in% c("bamboo", "trees")) |>
  select(land_class, all_of(pavd_columns)) |>
  pivot_longer(
    cols      = all_of(pavd_columns),
    names_to  = "z_bin",
    values_to = "pavd"
  ) |>
  mutate(
    height_mid = z_mid[match(z_bin, pavd_columns)],
    depth_m    = bin_depth[match(z_bin, pavd_columns)],
    pavd       = ifelse(pavd < 0, NA_real_, pavd),
    pai_contrib = pavd * depth_m
  ) |>
  filter(!is.na(pai_contrib))

pai_hist <- pai_long |>
  group_by(land_class, height_mid) |>
  summarise(
    mean_pai = mean(pai_contrib, na.rm = TRUE),
    .groups  = "drop"
  ) |>
  filter(height_mid <= 77.5) |>
  mutate(
    height_bin = factor(
      height_mid,
      levels = sort(unique(height_mid))
    ),
    mean_pai_signed = dplyr::case_when(
      land_class == "bamboo" ~ -mean_pai,
      land_class == "trees"  ~  mean_pai,
      TRUE                   ~  NA_real_
    )
  )

p_pai_pyramid <- ggplot(
  pai_hist,
  aes(x = mean_pai_signed, y = height_bin, fill = land_class)
) +
  geom_col(alpha = 0.8) +
  scale_fill_manual(values = cols2, guide = "none") +
  scale_x_continuous(
    name   = "Mean PAI contribution per height bin (m²/m²)",
    labels = function(x) abs(x)
  ) +
  labs(
    y     = "Height bin (m)",
    title = "Vertical distribution of PAI contributions\nBamboo (left) vs trees (right)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid   = element_blank(),
    axis.title.y = element_text(),
    axis.title.x = element_text(),
    axis.text.y  = element_text(),
    axis.text.x  = element_text(),
    plot.title   = element_text(face = "bold")
  )

p_pai_pyramid
