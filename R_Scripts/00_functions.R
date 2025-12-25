## 0. Utility functions and pre-installation of packages
##    Purpose: Load packages and define helper functions used throughout the
##    GEDI–Sentinel bamboo mapping and carbon-bias analysis pipeline.
##
##    - Computes vegetation indices (NDVI, SI, BI, BPCI, MBI) from Sentinel-2 bands.
##    - Provides a convenience wrapper (process_raster) to add all indices.
##    - Splits rasters into tiles for memory–efficient texture computation.
##    - Samples classified pixels for visual validation (exports KML).
##    - Wraps GLCM computation for a given band.
##
##    Inputs:  None directly (functions operate on rasters / classification outputs).
##    Outputs: In-memory functions (no files written except KML from sample_class()).
##
##    Notes:
##      * Source this script before running classification, validation, or
##        carbon-simulation scripts.
##      * compute_glcm() assumes that `window_size` and `texture_features`
##        are defined in the global environment.

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Packages
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

library(sp)
library(sf)
library(rgee)       
library(glcm)
library(caret)
library(terra)
library(torch)
library(tidyr)
library(dplyr)
library(purrr)
library(raster)
library(elevatr)
library(stringr)
library(ggplot2)
library(viridis)
library(ggridges)
library(bigmemory)
library(rasterVis)
library(neuralnet)
library(FactoMineR)   # PCA
library(factoextra)   # PCA visualisation
library(torchvision)
library(randomForest)

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 1. Vegetation indices (from Huang et al. / Goswami et al.)
#   Assumes Sentinel-2 style bands: B4 (red), B5, B6, B8 (NIR), B11 (SWIR).
#   Each function takes a SpatRaster / Raster* stack and appends a new layer.
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

getNDVI <- function(raster) {
  raster$NDVI <- (raster$B8 - raster$B4) / (raster$B8 + raster$B4)
  return(raster)
}

getSI <- function(raster) {
  raster$SI <- (raster$B8 - raster$B11) / (raster$B8 + raster$B11)
  return(raster)
}

getBI <- function(raster) {
  raster$BI <- (raster$NDVI - raster$SI) / (raster$NDVI + raster$SI)
  return(raster)
}

getBPCI <- function(raster) {
  raster$BPCI <- (raster$B5 / raster$B8 - raster$B8 / raster$B11)
  return(raster)
}

getMBI <- function(raster) {
  raster$MBI <- (raster$B11 - raster$B6) / (raster$B11 + raster$B6)
  return(raster)
}

# Convenience wrapper: add all indices in one call
process_raster <- function(raster) {
  raster <- getNDVI(raster)
  raster <- getSI(raster)
  raster <- getBI(raster)
  raster <- getBPCI(raster)
  raster <- getMBI(raster)
  return(raster)
}

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 2. Tile a raster into a regular grid (for texture / batch processing)
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

splitRasterToTiles <- function(r, ncol_tiles, nrow_tiles) {
  # Full extent of the raster
  r_ext <- ext(r)
  
  # Width and height of each tile
  tile_width  <- (r_ext[2] - r_ext[1]) / ncol_tiles
  tile_height <- (r_ext[4] - r_ext[3]) / nrow_tiles
  
  tiles   <- list()
  counter <- 1
  
  for (i in 0:(ncol_tiles - 1)) {
    for (j in 0:(nrow_tiles - 1)) {
      tile_ext <- ext(
        r_ext[1] + i * tile_width,
        r_ext[1] + (i + 1) * tile_width,
        r_ext[3] + j * tile_height,
        r_ext[3] + (j + 1) * tile_height
      )
      tiles[[counter]] <- crop(r, tile_ext)
      counter <- counter + 1
    }
  }
  return(tiles)
}

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 3. Sample classified pixels for visual validation (exports KML)
#    - pred_int_raster must exist in the global environment
#    - class_value: numeric ID in raster
#    - class_name: character label (used in output filename)
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

sample_class <- function(class_value, class_name, n_samples = 100) {
  all_cells   <- which(values(pred_int_raster) == class_value)
  n_available <- length(all_cells)
  
  if (n_available == 0L) {
    message("No pixels for '", class_name, "'. Skipping.")
    return(NULL)
  }
  
  sample_size <- if (n_available < n_samples) {
    message("Only ", n_available, " pixels for '", class_name, 
            "'. Sampling all instead of ", n_samples, ".")
    n_available
  } else {
    n_samples
  }
  
  set.seed(42)
  chosen_cells <- sample(all_cells, sample_size)
  xy_mat       <- xyFromCell(pred_int_raster, chosen_cells)
  
  pts_df  <- data.frame(category = class_name,
                        x = xy_mat[, 1],
                        y = xy_mat[, 2])
  pts_sf  <- st_as_sf(pts_df, coords = c("x", "y"),
                      crs = crs(pred_int_raster))
  pts_wgs84 <- st_transform(pts_sf, "EPSG:4326")
  
  fname <- paste0(class_name, "_sample.kml")
  st_write(pts_wgs84, fname,
           driver    = "KML",
           delete_dsn = TRUE,
           quiet      = FALSE)
  
  message("✓ Saved '", fname, "' with ", sample_size, " points.")
  return(pts_wgs84)
}

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 4. Wrapper for GLCM texture computation
#    - Assumes `window_size` and `texture_features` are defined globally
#      (e.g., window_size <- c(5,5); texture_features <- c("homogeneity","dissimilarity",...))
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

compute_glcm <- function(raster_band) {
  # Convert SpatRaster → RasterLayer for glcm()
  raster_band_raster <- raster(raster_band)
  
  glcm_textures <- glcm(
    raster_band_raster,
    window     = window_size,
    statistics = texture_features
  )
  return(glcm_textures)
}

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 5. Helper for listing local GEDI files by product + year
#    - Assumes users have downloaded GEDI .h5 files locally
#    - `dir` points to a folder containing one GEDI product
#      (e.g., GEDI04_A, GEDI02_A, GEDI02_B files)
#    - `product` should be a base shortname string:
#      "GEDI04_A", "GEDI02_A", "GEDI02_B"
#    - If `year` is provided, filters files that start with:
#      <product>_<year>
#    - Returns a character vector of full file paths
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

list_gedi_files <- function(dir, product, year = NULL) {
  stopifnot(is.character(dir), length(dir) == 1)
  stopifnot(is.character(product), length(product) == 1)
  
  if (!dir.exists(dir)) {
    message("Directory not found: ", dir)
    return(character(0))
  }
  
  pattern <- paste0("^", product, "_.*\\.h5$")
  
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  
  if (!is.null(year)) {
    year  <- as.character(year)
    files <- files[stringr::str_detect(
      basename(files),
      paste0("^", product, "_", year)
    )]
  }
  
  return(files)
}


