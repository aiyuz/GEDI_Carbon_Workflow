## ============================================================
## 06_bamboo_carbon_pipeline_revised.R
## Purpose:
##   (A1) Fit RH50 ~ RH98 for trees and bamboo (from GEDI L2A metrics)
##   (A2) Attach canopy-derived RH98 (= canopy height) to GEDI footprints
##   (A3) Predict RH50 using the fitted models; apply Kellner et al. EBT allometry
##   (A4) Compute footprint-level carbon_est and residuals vs GEDI L4A (carbon_obs)
##   (A5) Classify bamboo pixels into forest vs fallow (3x3 neighborhood rule)
##   (A6) Build pixel-level carbon_est_rast (EBT-based)
##   (A7) Save outputs to disk for later scripts (avoid terra pointer issues)
##
## Inputs (must exist before running):
##   - level2AM_df   : from 03_gedi_overlay.R (contains rh50, rh98, land_class, shot_number)
##   - gedi_sf       : GEDI footprint sf with 'land_class' and L4A 'agbd'
##   - predicted_raster : RF classification raster (1=bamboo, 6=trees)
##   - xsbn_final.tif : contains a canopy layer named 'canopy' (10 m height)
##
## Outputs (written to disk):
##   - data/gedi_sf_proj.rds
##   - data/bamboo_type_raster.tif
##   - data/carbon_est_rast.tif
##   - data/rh_models.rds (coefficients + helper)
## ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
})

## ---------------------------
## CONFIG
## ---------------------------
dir.create("data", showWarnings = FALSE)

PATH_XSBN_FINAL   <- "xsbn_final.tif"
OUT_GEDI_RDS      <- file.path("data", "gedi_sf_proj.rds")
OUT_BAMBOO_TYPE   <- file.path("data", "bamboo_type_raster.tif")
OUT_CARBON_RAST   <- file.path("data", "carbon_est_rast.tif")
OUT_RH_MODELS_RDS <- file.path("data", "rh_models.rds")

CLASS_BAMBOO <- 1
CLASS_TREES  <- 6
NB_THRESHOLD <- 5  # 3x3 window; >=5 neighbors (incl. center) => "forest"

## ---------------------------
## Global functions
## ---------------------------

# Kellner et al. (2023) EBT global allometry: AGBD (Mg/ha)
agbd_ebt_global <- function(rh50, rh98) {
  rh50 <- as.numeric(rh50)
  rh98 <- as.numeric(rh98)
  term <- -104.965 + 6.802 * sqrt(rh50 + 100) + 3.955 * sqrt(rh98 + 100)
  1.113 * (term^2)
}

fit_rh_model <- function(df, class_name) {
  stopifnot(all(c("rh50", "rh98") %in% names(df)))
  df2 <- df %>% filter(!is.na(rh50), !is.na(rh98))
  if (nrow(df2) < 30) {
    warning("RH model for ", class_name, " has low sample size: n = ", nrow(df2))
  }
  m <- lm(rh50 ~ rh98, data = df2)
  coef(m)
}

predict_rh50 <- function(rh98, coef_vec) {
  as.numeric(coef_vec[1] + coef_vec[2] * rh98)
}

## ============================================================
## A1) Fit RH50~RH98 models for trees and bamboo
## ============================================================
stopifnot(exists("level2AM_df"))

bt_rh <- level2AM_df %>%
  filter(land_class %in% c("bamboo", "trees")) %>%
  dplyr::select(shot_number, land_class, rh50, rh98)

coef_tree <- fit_rh_model(bt_rh %>% filter(land_class == "trees"),  "trees")
coef_bam  <- fit_rh_model(bt_rh %>% filter(land_class == "bamboo"), "bamboo")

saveRDS(
  list(coef_tree = coef_tree, coef_bamboo = coef_bam),
  OUT_RH_MODELS_RDS
)

### report stats:
library(dplyr)
library(broom)

# Fit model and extract stats
rh_model_stats <- function(df, class_label){
  m <- lm(rh50 ~ rh98, data = df)
  
  s <- summary(m)
  rmse <- sqrt(mean(residuals(m)^2, na.rm = TRUE))
  
  tibble(
    class = class_label,
    model = "RH50 = β0 + β1·RH98",
    beta0 = coef(m)[1],
    beta1 = coef(m)[2],
    n     = nobs(m),
    r2    = s$r.squared,
    rmse  = rmse,
    se_beta0 = s$coefficients[1,2],
    se_beta1 = s$coefficients[2,2]
  )
}

tab2_panelA <- bind_rows(
  rh_model_stats(filter(bt_rh, land_class=="trees"),  "Trees"),
  rh_model_stats(filter(bt_rh, land_class=="bamboo"), "Bamboo")
) %>%
  mutate(
    beta0 = round(beta0, 3),
    beta1 = round(beta1, 3),
    r2    = round(r2, 3),
    rmse  = round(rmse, 2),
    se_beta0 = round(se_beta0, 3),
    se_beta1 = round(se_beta1, 3)
  )

tab2_panelA

## ============================================================
## A2) Attach canopy-derived RH98 to GEDI footprints
## ============================================================
stopifnot(file.exists(PATH_XSBN_FINAL))
stopifnot(exists("gedi_sf"))

xsbn_final  <- rast(PATH_XSBN_FINAL)
stopifnot("canopy" %in% names(xsbn_final))
canopy_rast <- xsbn_final[["canopy"]]
crs_canopy  <- crs(canopy_rast)

# project footprints to canopy CRS
gedi_sf_proj <- st_transform(gedi_sf, crs = crs_canopy)

# extract canopy height at footprints
canopy_vals <- terra::extract(canopy_rast, vect(gedi_sf_proj))
gedi_sf_proj$canopy_10m <- canopy_vals[, 2]
gedi_sf_proj <- gedi_sf_proj %>% filter(!is.na(canopy_10m))

# approximate RH98 by canopy height
gedi_sf_proj$rh98 <- gedi_sf_proj$canopy_10m

## ============================================================
## A3) Predict RH50 by land_class and compute EBT AGBD + carbon
## ============================================================
stopifnot("land_class" %in% names(gedi_sf_proj))
stopifnot("agbd" %in% names(gedi_sf_proj))  # GEDI L4A AGBD

tree_idx   <- gedi_sf_proj$land_class == "trees"
bamboo_idx <- gedi_sf_proj$land_class == "bamboo"

gedi_sf_proj$rh50 <- NA_real_
gedi_sf_proj$rh50[tree_idx]   <- predict_rh50(gedi_sf_proj$rh98[tree_idx],   coef_tree)
gedi_sf_proj$rh50[bamboo_idx] <- predict_rh50(gedi_sf_proj$rh98[bamboo_idx], coef_bam)

# EBT allometry
gedi_sf_proj$est_agbd <- agbd_ebt_global(gedi_sf_proj$rh50, gedi_sf_proj$rh98)

# carbon (Mg C / ha)
gedi_sf_proj$carbon_obs <- gedi_sf_proj$agbd     * 0.5   # GEDI L4A-based carbon
gedi_sf_proj$carbon_est <- gedi_sf_proj$est_agbd * 0.5   # canopy->EBT estimate

# residuals (estimate - GEDI L4A)
gedi_sf_proj$residual_agbd   <- gedi_sf_proj$est_agbd   - gedi_sf_proj$agbd
gedi_sf_proj$residual_carbon <- gedi_sf_proj$carbon_est - gedi_sf_proj$carbon_obs

## ============================================================
## A4) Optional: quick bias summary (kept lightweight)
##      (If you want figures/diagnostics, do them in Script B)
## ============================================================
bias_summary <- gedi_sf_proj %>%
  st_drop_geometry() %>%
  filter(land_class %in% c("trees", "bamboo")) %>%
  group_by(land_class) %>%
  summarise(
    n = sum(is.finite(residual_carbon)),
    mean_resid_C = mean(residual_carbon, na.rm = TRUE),
    sd_resid_C   = sd(residual_carbon, na.rm = TRUE),
    rmse_resid_C = sqrt(mean(residual_carbon^2, na.rm = TRUE)),
    .groups = "drop"
  )
bias_summary

# Optional figure: residuals by RH98 bin for trees
gedi_sf_proj$height_bin <- cut(gedi_sf_proj$rh98,
                               breaks = seq(0, 60, by = 10),
                               include.lowest = TRUE)

ggplot(gedi_sf_proj[tree_idx, ], aes(x = height_bin, y = residual_agbd)) +
  geom_boxplot(width = 0.5,
               fill  = scales::alpha("#ac03c9", 0.4),
               color = "#ac03c9",
               outlier.shape = NA) +
  geom_hline(yintercept = 0, color = "black",
             linetype = "dashed", linewidth = 0.8) +
  labs(
    title = "Residual Bias of GEDI EBT Model by RH98 Bin (Trees)",
    x     = "RH98 Bin (m)",
    y     = "Residual AGBD (Mg/ha; est - GEDI L4A)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )

ggplot(gedi_sf_proj[bamboo_idx, ], aes(x = height_bin, y = residual_agbd)) +
  geom_boxplot(width = 0.5,
               fill  = scales::alpha("#00CD10", 0.4),
               color = "#00CD10",
               outlier.shape = NA) +
  geom_hline(yintercept = 0, color = "black",
             linetype = "dashed", linewidth = 0.8) +
  labs(
    title = "Residual Bias of GEDI EBT Model by RH98 Bin (Bamboo)",
    x     = "RH98 Bin (m)",
    y     = "Residual AGBD (Mg/ha; est - GEDI L4A)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )


## ============================================================
## A5) Classify bamboo forest vs fallow on raster grid
## ============================================================
stopifnot(exists("predicted_raster"))
predicted_raster <- rast(predicted_raster)  # ensure SpatRaster

# mask: vegetation = bamboo or trees => 1; else NA
bt_mask <- classify(
  predicted_raster,
  rbind(c(CLASS_BAMBOO, 1),
        c(CLASS_TREES,  1)),
  others = NA
)

w3 <- matrix(1, 3, 3)
neighbor_sum <- focal(
  bt_mask,
  w = w3,
  fun = sum,
  na.policy = "omit",
  na.rm = TRUE
)

bamboo_mask        <- predicted_raster == CLASS_BAMBOO
bamboo_forest_mask <- bamboo_mask & (neighbor_sum >= NB_THRESHOLD)
bamboo_fallow_mask <- bamboo_mask & (neighbor_sum <  NB_THRESHOLD)

bamboo_type_raster <- rast(predicted_raster)
values(bamboo_type_raster) <- NA
bamboo_type_raster[bamboo_forest_mask] <- 1  # forest
bamboo_type_raster[bamboo_fallow_mask] <- 2  # fallow

# attach bamboo_type to bamboo footprints only
btype_vals <- terra::extract(bamboo_type_raster, vect(gedi_sf_proj))[, 2]
gedi_sf_proj$bamboo_type <- NA_character_
gedi_sf_proj$bamboo_type[bamboo_idx & btype_vals == 1] <- "forest"
gedi_sf_proj$bamboo_type[bamboo_idx & btype_vals == 2] <- "fallow"

## ============================================================
## A6) Build pixel-level carbon_est_rast (EBT-based)
## ============================================================
rh98_rast <- canopy_rast

rh50_rast <- ifel(
  predicted_raster == CLASS_TREES,
  predict_rh50(rh98_rast, coef_tree),
  ifel(
    predicted_raster == CLASS_BAMBOO,
    predict_rh50(rh98_rast, coef_bam),
    NA
  )
)

agbd_est_rast <- lapp(
  c(rh50_rast, rh98_rast),
  fun = function(x1, x2) agbd_ebt_global(x1, x2)
)
carbon_est_rast <- agbd_est_rast * 0.5

## ============================================================
## A7) Save outputs (critical for stability across sessions)
## ============================================================
saveRDS(gedi_sf_proj, OUT_GEDI_RDS)
writeRaster(bamboo_type_raster, OUT_BAMBOO_TYPE, overwrite = TRUE)
writeRaster(carbon_est_rast,    OUT_CARBON_RAST, overwrite = TRUE)

message("Saved: ", OUT_GEDI_RDS)
message("Saved: ", OUT_BAMBOO_TYPE)
message("Saved: ", OUT_CARBON_RAST)
message("Saved: ", OUT_RH_MODELS_RDS)