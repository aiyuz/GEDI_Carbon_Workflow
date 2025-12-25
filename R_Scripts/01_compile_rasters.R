## 01_compile_rasters_xsbn.R
## Purpose:
##   - Read seasonal Sentinel-2 composites for Xishuangbanna (from GEE) plus
##     canopy height, assemble them into a single feature stack (xsbn_final)
##     for one season (currently Oct–Dec).
##   - Load training polygons and inspect their spectral/seasonal signatures.
##   - Produce the following Supplementary Figures:
##       * Figure S1: DTM + training polygons (label distribution in terrain)
##       * Figure S2: Seasonal mean reflectance by category × band (heatmap)
##       * Figure S3: Spectral band distributions by class (ridge densities)
##   - Save xsbn_final.tif for downstream classification.
##
## Inputs (local folders expected):
##   - MergedSeasons/*.tif : 4 seasonal Sentinel-2 stacks (JanMar, AprJun, JulSep, OctDec)
##   - Terrain/DTM_Xishuangbanna.tif
##   - Terrain/ETH_GlobalCanopyHeight_2020_10m_Xishuangbanna.tif
##   - Terrain/ETH_GlobalCanopyHeightSD_2020_10m_Xishuangbanna.tif
##   - training_updated/*.kml : KML polygons for each land-cover class
##
## Outputs:
##   - xsbn_final.tif
##   - ref_curve_data.csv
##   - training_data.rds
##   - Figure_S1_DTM_training_labels.png
##   - Figure_S2_seasonal_reflectance_heatmap.png
##   - Figure_S3_spectral_band_distributions.png

# setwd("/yourworkdirectory")
source("00_functions.R")   # vegetation indices, glcm wrapper, etc.

library(sp)
library(sf)
library(rgee)   
library(glcm)
library(terra)
library(tidyr)
library(dplyr)
library(purrr)
library(raster)
library(elevatr)
library(ggplot2)
library(viridis)
library(rasterVis)
library(ggridges)
library(FactoMineR)   # For PCA
library(factoextra)   # For Visualization)
library(ggspatial)

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Step 1: gather raw data together into tiffs
#   - Read 4 seasonal Sentinel-2 composites (MergedSeasons)
#   - Read global canopy height + SD, reproject, crop, and append as layers
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

tif_dir   <- "MergedSeasons"
tif_files <- list.files(tif_dir, pattern = "\\.tif$", full.names = TRUE)
seasons   <- c("JanMar","AprJun","JulSep","OctDec")

# Read rasters into a named list: merged_rasters$JanMar, $AprJun, ...
merged_rasters <- setNames(lapply(tif_files, rast), seasons)

# Add the canopy height model (exported from GEE to Terrain/)
canopy    <- rast("Terrain/ETH_GlobalCanopyHeight_2020_10m_Xishuangbanna.tif")
canopy_sd <- rast("Terrain/ETH_GlobalCanopyHeightSD_2020_10m_Xishuangbanna.tif")

# Reproject if necessary to match the first seasonal raster
if (!compareCRS(canopy, merged_rasters[[1]])) {
  canopy    <- project(canopy,    merged_rasters[[1]])
  canopy_sd <- project(canopy_sd, merged_rasters[[1]])
}

# Crop the CHM to the extent of AOI (XSBN) and align resolution
canopy_cropped    <- crop(canopy,    merged_rasters[[1]]$B2)
canopy_sd_cropped <- crop(canopy_sd, merged_rasters[[1]]$B2)

canopy_aligned <- resample(
  canopy_cropped,
  merged_rasters[[1]]$B2,
  method = "bilinear"
)

canopy_sd_aligned <- resample(
  canopy_sd_cropped,
  merged_rasters[[1]]$B2,
  method = "bilinear"
)

xsbn_canopy    <- mask(canopy_aligned,    merged_rasters[[1]]$B2)
xsbn_canopy_sd <- mask(canopy_sd_aligned, merged_rasters[[1]]$B2)

# Append canopy and canopy_sd to each seasonal stack
for (i in 1:4) {
  merged_rasters[[i]]$canopy    <- xsbn_canopy
  merged_rasters[[i]]$canopy_sd <- xsbn_canopy_sd
}

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Step 2: load training labels (polygons) and terrain
#   - Used for Figure S1 (DTM + labels) and later extraction
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# Read in label polygons
kml_files <- list.files("training_updated", 
                        pattern = "\\.kml$", full.names = TRUE)
feature_list <- lapply(kml_files, st_read)
category_names <- gsub(".*\\/|\\.kml$", "", kml_files)  # Extract label names
names(feature_list) <- category_names

# Color scheme for labels (must match later classification)
colors <- c(
  bamboo       = "#00CD10",  # neon green (fixed)
  trees        = "#ac03c9",  # purple-ish
  farm         = "#E69F00",  # orange
  plantation   = "#F0E442",  # yellow
  rubber       = "#CC79A7",  # magenta / pink-purple
  construction = "#666666",  # dark grey
  water        = "#1E90FF"   # blue
)

# DTM for Figure S1 background
dtm    <- rast("Terrain/DTM_Xishuangbanna.tif")
dtm_df <- as.data.frame(dtm, xy = TRUE)
colnames(dtm_df) <- c("x", "y", "elevation")

# Stack feature_list into one sf, adding a 'category' column
features_sf <- do.call(rbind, lapply(names(feature_list), function(cat) {
  sf_obj <- feature_list[[cat]]
  sf_obj$category <- cat
  sf_obj
}))

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Figure S1: DEM + training labels (label distribution across terrain)
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

fig_S1_dem_labels <- ggplot() +
  geom_raster(
    data = dtm_df,
    aes(x = x, y = y, fill = elevation)
  ) + 
  scale_fill_gradient(
    name = "Elevation (m)", 
    low  = "black",
    high = "white"
  ) + 
  geom_sf(
    data  = features_sf, 
    aes(color = category),
    fill  = NA,
    size  = 0.8
  ) +
  scale_color_manual(
    name   = "Category",
    values = colors
  ) +
  ggspatial::annotation_scale(
    location  = "bl",
    style     = "bar",
    width_hint = 0.2,
    bar_cols  = c("grey60","white"),
    text_cex  = 0.7
  ) +
  coord_sf() +
  theme_minimal() +
  theme(
    legend.position  = "right",
    panel.grid.major = element_line(color = "transparent")
  )

ggsave(
  filename = "Figure_S1_DTM_training_labels.png",
  plot     = fig_S1_dem_labels,
  width    = 7,
  height   = 7,
  dpi      = 300
)

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Step 3: Seasonal reflectance profiles by class (Figure S2)
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

band_mapping <- c(
  B2  = "Blue",
  B3  = "Green",
  B4  = "Red",
  B5  = "Red Edge 1",
  B6  = "Red Edge 2",
  B7  = "Red Edge 3",
  B8  = "NIR",
  B11 = "SWIR 1",
  B12 = "SWIR 2"
)

spectral_data_list <- vector("list", 4)
names(spectral_data_list) <- paste0("spectral_data_", seasons)
ref_curve_data    <- data.frame()
season_names      <- c("January-March","April-June","July-September","October-December")

for (i in 1:4) {
  spectral_features <- lapply(names(feature_list), function(name) {
    polygons  <- feature_list[[name]]
    extracted <- terra::extract(merged_rasters[[seasons[i]]], 
                                polygons, df = TRUE)
    extracted$category <- name
    extracted
  })
  spectral_data_list[[i]] <- do.call(rbind, spectral_features)
  
  season_long <- spectral_data_list[[i]] %>%
    pivot_longer(
      cols      = starts_with("B"), 
      names_to  = "band", 
      values_to = "reflectance"
    ) %>%
    filter(!category %in% c("water", "construction")) %>%
    mutate(band = recode(band, !!!band_mapping)) %>%
    filter(!band %in% c("NDVI","MBI","SI","BI","BPCI")) %>%
    filter(!is.na(reflectance))
  
  summary <- season_long %>%
    group_by(category, band) %>%
    summarise(
      mean_reflectance = mean(reflectance, na.rm = TRUE),
      .groups          = "drop"
    )
  summary$season <- season_names[i]
  ref_curve_data <- rbind(summary, ref_curve_data)
}

season_order <- c("January-March", "April-June", 
                  "July-September", "October-December")

ref_curve_data <- ref_curve_data %>%
  mutate(season = factor(season, levels = season_order)) %>%
  filter(!category %in% c("water", "construction"))

write.csv(ref_curve_data, "ref_curve_data.csv", row.names = FALSE)

fig_S2_reflectance <- ggplot(
  ref_curve_data, 
  aes(x = band, y = category, fill = mean_reflectance)
) +
  geom_tile(color = "white") +
  scale_x_discrete(limits = band_mapping) +
  scale_fill_viridis_c(
    name   = "Reflectance",
    option = "inferno"
  ) +
  facet_wrap(~season, ncol = 1) +
  labs(
    title = "Seasonal Mean Reflectance by Category × Band",
    x     = "Spectral Band",
    y     = "Category"
  ) +
  theme_minimal() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "Figure_S2_seasonal_reflectance_heatmap.png",
  plot     = fig_S2_reflectance,
  width    = 7,
  height   = 9,
  dpi      = 300
)

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Step 4: Choose Oct–Dec stack, compute indices + texture, save xsbn_final
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

xsbn_rast <- merged_rasters$OctDec

# Vegetation indices
xsbn <- getBPCI(xsbn_rast)
xsbn <- getNDVI(xsbn)

# GLCM texture on B8
selected_bands   <- "B8"
texture_features <- c("homogeneity", "dissimilarity")
window_size      <- c(5, 5)  # window size used in the paper

glcm_results <- lapply(selected_bands, function(band) {
  message("Computing GLCM for ", band, " with 5×5 window")
  raster_band <- subset(xsbn, band)
  compute_glcm(raster_band)
})

glcm_stack <- rast(stack(glcm_results))
glcm_layer_names <- paste0(
  rep(selected_bands, each = length(texture_features)), "_",
  rep(texture_features, times = length(selected_bands)), "_5x5"
)
names(glcm_stack) <- glcm_layer_names

# Final feature stack for Oct–Dec
xsbn_final <- c(xsbn, glcm_stack)
writeRaster(xsbn_final, "xsbn_final.tif", overwrite = TRUE)

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Step 5: Extract training_data and build Figure S3 (spectral band distributions)
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# Convert polygons to SpatVector
features_vect <- terra::vect(features_sf)

# Extract raster values for each polygon (from xsbn_final)
training_data <- terra::extract(xsbn_final, features_vect, df = TRUE, mask = TRUE)

# Use 'category' from features_sf as the class label
training_data$label <- features_sf$category[training_data$ID]

training_data <- training_data %>%
  dplyr::select(-ID) %>%
  na.omit()

training_data$label <- factor(
  training_data$label,
  levels = c("bamboo", "trees", "plantation", "rubber",
             "farm", "construction", "water")
)

saveRDS(training_data, "training_data.rds")

# Figure S3: spectral band distributions by class (density ridges)
band_cols <- intersect(
  names(training_data),
  c("B2","B3","B4","B5","B6","B7","B8","B11","B12")
)

bands_long <- training_data %>%
  dplyr::select(all_of(c("label", band_cols))) %>%
  pivot_longer(
    cols      = all_of(band_cols),
    names_to  = "band",
    values_to = "value"
  ) %>%
  mutate(
    band = recode(band, !!!band_mapping)
  ) %>%
  filter(!is.na(value))

fig_S3_bands <- ggplot(
  bands_long,
  aes(x = value, y = label, fill = label)
) +
  geom_density_ridges(
    alpha          = 0.8,
    color          = "black",
    scale          = 1.1,
    size           = 0.3,
    rel_min_height = 0.01
  ) +
  scale_fill_viridis_d(option = "inferno", direction = -1) +
  facet_wrap(
    ~ band,
    scales = "free_x",
    ncol   = 3
  ) +
  labs(
    title = "Spectral band distributions by class",
    x     = "Reflectance",
    y     = NULL,
    fill  = "Class"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text      = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  filename = "Figure_S3_spectral_band_distributions.png",
  plot     = fig_S3_bands,
  width    = 8,
  height   = 8,
  dpi      = 300
)
