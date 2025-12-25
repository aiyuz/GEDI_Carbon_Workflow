## ============================================================
## 07_bamboo_bias_figures.R
## Purpose:
##   (B1) Compare bamboo carbon (GEDI L4A) vs Yuen (2017) with CI & Welch test
##   (B2) Quantify misclassification-aware bamboo bias via Monte Carlo resampling
##   (B3) Combine MC bias with EBT simplification bias -> mean bamboo bias raster
##   (B4) Map mean % bias relative to EBT-estimated carbon (carbon_est_rast)
##   (B5) Regional bamboo bias maps and zonal summaries
##   (B6) Summary plots of EBT-based carbon and areas (prior vs posterior)
##
## Inputs (read from disk):
##   - data/gedi_sf_proj.rds
##   - data/bamboo_type_raster.tif
##   - data/carbon_est_rast.tif
##   - rf_prediction.tif
##   - (from 05) post_bayes + sample_true_class() available in session
##     OR loaded from disk if saved in 05
##   - Xishuangbanna shapefile:
##       Xishuangbanna_Polygon/Xishuangbanna_Polygon.shp
##
## Outputs:
##   - Figures (e.g. for Figure 4 panels and related plots)
##   - Optional rasters (e.g. data/bias_mean_rast.tif)
##   - Summary tables for Table 2 (panel C) and B7
## ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(viridisLite)
})

dir.create("figs", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

## ---------------------------
## CONFIG
## ---------------------------
IN_GEDI_RDS    <- file.path("data", "gedi_sf_proj.rds")
IN_BAMBOO_TYPE <- file.path("data", "bamboo_type_raster.tif")
IN_CARBON_RAST <- file.path("data", "carbon_est_rast.tif")

# Colors for bamboo contexts
pal_bamboo <- c(
  "Forest bamboo" = "#00CD10",
  "Fallow bamboo" = "#795C34"
)

# Yuen (2017) empirical summary values (Mg C / ha)
bamboo_fallow_mean <- 14.7
bamboo_fallow_sd   <- 14.1
bamboo_fallow_max  <- 56.4
n_fallow_emp       <- 35

bamboo_forest_mean <- 27.5
bamboo_forest_sd   <- 43.1
bamboo_forest_max  <- 162
n_forest_emp       <- 24

# Monte Carlo settings for misclassification-aware bias
RUN_MC        <- TRUE   # Set FALSE if mc_bamboo.rds is already saved
N_ITER_BAMBOO <- 1000
OUT_MC_RDS    <- file.path("data", "mc_bamboo.rds")

# Map settings
EPS       <- 1e-6
PBIAS_CAP <- 100   # Clamp % bias for visualization

## ---------------------------
## Load inputs
## ---------------------------
stopifnot(file.exists(IN_GEDI_RDS))
stopifnot(file.exists(IN_BAMBOO_TYPE))
stopifnot(file.exists(IN_CARBON_RAST))

gedi_sf_proj       <- readRDS(IN_GEDI_RDS)
bamboo_type_raster <- rast(IN_BAMBOO_TYPE)
carbon_est_rast    <- rast(IN_CARBON_RAST)
prediction_raster  <- rast("rf_prediction.tif")

## ============================================================
## B1) Empirical vs GEDI L4A: CI plot (outlier-filtered GEDI)
##     - Compare GEDI L4A bamboo carbon to Yuen (2017)
##     - Compute t-based 95% CIs (GEDI & empirical)
##     - Welch’s t-tests (GEDI > empirical)
## ============================================================

compute_ci_summary <- function(x_L4, mean_emp, sd_emp, n_emp, group_name) {
  # Empirical CI
  emp_se <- sd_emp / sqrt(n_emp)
  emp_t  <- qt(0.975, df = n_emp - 1)
  emp_ci <- mean_emp + c(-1, 1) * emp_t * emp_se
  
  # GEDI CI (x_L4 already cleaned); robust n
  x <- x_L4[is.finite(x_L4)]
  n_L4 <- length(x)
  mean_L4 <- mean(x, na.rm = TRUE)
  sd_L4   <- sd(x,   na.rm = TRUE)
  
  if (n_L4 >= 2) {
    se_L4 <- sd_L4 / sqrt(n_L4)
    t_L4  <- qt(0.975, df = n_L4 - 1)
    ci_L4 <- mean_L4 + c(-1, 1) * t_L4 * se_L4
  } else {
    ci_L4 <- c(NA_real_, NA_real_)
  }
  
  tibble(
    bamboo_group = group_name,
    source       = c("Empirical", "GEDI L4A"),
    mean         = c(mean_emp,   mean_L4),
    ci_low       = c(emp_ci[1],  ci_L4[1]),
    ci_high      = c(emp_ci[2],  ci_L4[2]),
    n            = c(n_emp,      n_L4)
  )
}

# Build bamboo footprint data frame from GEDI L4A carbon (carbon_obs)
bamboo_df <- gedi_sf_proj %>%
  st_drop_geometry() %>%
  filter(land_class == "bamboo", !is.na(bamboo_type)) %>%
  select(shot_number, bamboo_type, carbon_obs)

bamboo_fallow_df <- bamboo_df %>% filter(bamboo_type == "fallow")
bamboo_forest_df <- bamboo_df %>% filter(bamboo_type == "forest")

# Tukey outlier rule per group (used to remove extreme GEDI values)
flag_outliers_tukey <- function(x) {
  q   <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  lower <- q[1] - 1.5 * iqr
  upper <- q[2] + 1.5 * iqr
  x < lower | x > upper
}

bamboo_fallow_df_plot <- bamboo_fallow_df %>%
  mutate(
    bamboo_group = "Fallow bamboo",
    is_outlier   = flag_outliers_tukey(carbon_obs)
  )

bamboo_forest_df_plot <- bamboo_forest_df %>%
  mutate(
    bamboo_group = "Forest bamboo",
    is_outlier   = flag_outliers_tukey(carbon_obs)
  )

# Cleaned GEDI L4A values (non-outliers only)
fallow_L4_clean <- bamboo_fallow_df_plot$carbon_obs[!bamboo_fallow_df_plot$is_outlier]
forest_L4_clean <- bamboo_forest_df_plot$carbon_obs[!bamboo_forest_df_plot$is_outlier]

# CI summaries for empirical vs GEDI
fallow_summary <- compute_ci_summary(
  fallow_L4_clean,
  bamboo_fallow_mean, bamboo_fallow_sd, n_fallow_emp,
  "Fallow bamboo"
)
forest_summary <- compute_ci_summary(
  forest_L4_clean,
  bamboo_forest_mean, bamboo_forest_sd, n_forest_emp,
  "Forest bamboo"
)

summary_all <- bind_rows(fallow_summary, forest_summary)
bamboo_all_df_plot <- bind_rows(bamboo_fallow_df_plot, bamboo_forest_df_plot)

y_max <- max(summary_all$ci_high, bamboo_fallow_max, bamboo_forest_max, na.rm = TRUE) * 1.10

# Welch’s t-test from summary statistics (GEDI vs empirical)
welch_from_summary <- function(mean_x, sd_x, n_x, mean_y, sd_y, n_y,
                               alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  
  diff_mean <- mean_x - mean_y
  se_diff   <- sqrt(sd_x^2 / n_x + sd_y^2 / n_y)
  
  df <- se_diff^4 / (
    (sd_x^2 / n_x)^2 / (n_x - 1) +
      (sd_y^2 / n_y)^2 / (n_y - 1)
  )
  
  t_stat <- diff_mean / se_diff
  
  p_val <- switch(
    alternative,
    "two.sided" = 2 * (1 - pt(abs(t_stat), df = df)),
    "greater"   = 1 - pt(t_stat, df = df),
    "less"      = pt(t_stat, df = df)
  )
  
  tibble(diff_mean = diff_mean, t = t_stat, df = df, p_value = p_val, alternative = alternative)
}

# Tests: GEDI > empirical
fallow_test <- welch_from_summary(
  mean_x = mean(fallow_L4_clean, na.rm = TRUE),
  sd_x   = sd(fallow_L4_clean,   na.rm = TRUE),
  n_x    = sum(is.finite(fallow_L4_clean)),
  mean_y = bamboo_fallow_mean,
  sd_y   = bamboo_fallow_sd,
  n_y    = n_fallow_emp,
  alternative = "greater"
)

forest_test <- welch_from_summary(
  mean_x = mean(forest_L4_clean, na.rm = TRUE),
  sd_x   = sd(forest_L4_clean,   na.rm = TRUE),
  n_x    = sum(is.finite(forest_L4_clean)),
  mean_y = bamboo_forest_mean,
  sd_y   = bamboo_forest_sd,
  n_y    = n_forest_emp,
  alternative = "greater"
)

# Labels for n per group/source
n_labels_all <- summary_all %>%
  transmute(bamboo_group, source, label = paste0("n = ", n))

# Labels for p-values per bamboo group
p_labels <- bind_rows(
  tibble(
    bamboo_group = "Fallow bamboo",
    source       = "GEDI L4A",
    p_value      = forest_test$p_value * 0 + fallow_test$p_value
  ),
  tibble(
    bamboo_group = "Forest bamboo",
    source       = "GEDI L4A",
    p_value      = forest_test$p_value
  )
) %>%
  mutate(label = paste0("p = ", signif(p_value, 2)))

# Reference maxima from Yuen (2017)
max_df <- tibble(
  bamboo_group = c("Fallow bamboo", "Forest bamboo"),
  max_c        = c(bamboo_fallow_max, bamboo_forest_max)
) %>%
  mutate(
    max_label = paste0("ref max = ", round(max_c, 1), " Mg C/ha"),
    label_y   = max_c + 0.05 * y_max
  )

# Panel B1: GEDI vs empirical CIs and tests
p_bamboo_ci <- ggplot(summary_all, aes(x = source, y = mean)) +
  # Footprints (non-outliers only)
  geom_jitter(
    data  = bamboo_all_df_plot %>% filter(!is_outlier),
    aes(x = "GEDI L4A", y = carbon_obs, color = bamboo_group),
    width = 0.25, alpha = 0.4, size = 2.5
  ) +
  # Empirical CI
  geom_errorbar(
    data = summary_all %>% filter(source == "Empirical"),
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.12, linewidth = 0.5, color = "black"
  ) +
  # GEDI CI
  geom_errorbar(
    data = summary_all %>% filter(source == "GEDI L4A"),
    aes(ymin = ci_low, ymax = ci_high, color = bamboo_group),
    width = 0.12, linewidth = 0.5
  ) +
  # Means
  geom_point(data = summary_all %>% filter(source == "Empirical"), size = 4, color = "black") +
  geom_point(data = summary_all %>% filter(source == "GEDI L4A"), aes(color = bamboo_group), size = 4) +
  # Empirical maxima
  geom_hline(data = max_df, aes(yintercept = max_c),
             linewidth = 0.7, linetype = "dashed", color = "black") +
  geom_text(data = max_df, aes(x = 1.2, y = label_y, label = max_label),
            inherit.aes = FALSE, size = 3.5, color = "black") +
  # n labels
  geom_text(data = n_labels_all, aes(x = source, y = 0.3, label = label),
            vjust = 1.8, size = 3.5, color = "black") +
  # p-value labels
  geom_text(data = p_labels, aes(x = source, y = y_max * 0.8, label = label),
            vjust = 1.0, size = 4, color = "darkred") +
  coord_cartesian(ylim = c(0, y_max), clip = "off") +
  facet_wrap(~ bamboo_group, nrow = 1) +
  scale_color_manual(values = pal_bamboo, guide = "none") +
  labs(
    x = NULL,
    y = "Carbon density (Mg C / ha)",
    title = "Bamboo Aboveground Carbon Density",
    subtitle = "GEDI L4A footprints (non-outliers) vs Yuen (2017) empirical means ± 95% CI"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin        = ggplot2::margin(t = 10, r = 10, b = 25, l = 10, unit = "pt")
  )

print(p_bamboo_ci)

## ============================================================
## B2) MC: Misclassification-aware bias (footprint level)
##     - Use posterior class probabilities (post_bayes) and
##       sample_true_class() to simulate “true” classes
##     - Compute bias per iteration (GEDI L4A – empirical mean)
##     - Summarize distribution and visualize
## ============================================================

# EBT simplification bias (mean residual from canopy->EBT vs GEDI L4A) for bamboo
bamboo_idx <- gedi_sf_proj$land_class == "bamboo"
ebt_bias_bamboo_mean <- mean(gedi_sf_proj$residual_carbon[bamboo_idx], na.rm = TRUE)
ebt_bias_bamboo_sd   <- sd(gedi_sf_proj$residual_carbon[bamboo_idx],   na.rm = TRUE)

message("EBT simplification bias (bamboo): mean=", round(ebt_bias_bamboo_mean, 3),
        " sd=", round(ebt_bias_bamboo_sd, 3))

# Build bamboo_L4 table for MC using non-outlier GEDI L4A values
bamboo_L4 <- bind_rows(
  bamboo_fallow_df_plot %>%
    filter(!is_outlier) %>%
    mutate(land_class = "bamboo", bamboo_type = "fallow"),
  bamboo_forest_df_plot %>%
    filter(!is_outlier) %>%
    mutate(land_class = "bamboo", bamboo_type = "forest")
) %>%
  select(shot_number, land_class, bamboo_type, carbon_obs)

# MC helper for a single iteration
run_one_mc_bamboo <- function(bamboo_L4, post_mat, mean_emp_fallow, mean_emp_forest) {
  
  # Requires: sample_true_class(class_label, post_mat)
  stopifnot(exists("sample_true_class"))
  stopifnot(is.matrix(post_mat) || is.data.frame(post_mat))
  
  df_draw <- bamboo_L4 %>%
    mutate(
      land_class = as.character(land_class),
      draw_class = vapply(
        land_class,
        sample_true_class,
        FUN.VALUE = character(1),
        post_mat  = post_mat
      )
    )
  
  df_true_bamboo <- df_draw %>% filter(draw_class == "bamboo")
  
  # If too few bamboo footprints are retained, skip this iteration
  if (nrow(df_true_bamboo) < 5) {
    return(tibble(
      fallow_bias = NA_real_,
      forest_bias = NA_real_,
      n_fallow    = NA_integer_,
      n_forest    = NA_integer_
    ))
  }
  
  fallow_vals <- df_true_bamboo$carbon_obs[df_true_bamboo$bamboo_type == "fallow"]
  forest_vals <- df_true_bamboo$carbon_obs[df_true_bamboo$bamboo_type == "forest"]
  
  fallow_bias <- if (sum(is.finite(fallow_vals)) >= 2) mean(fallow_vals, na.rm = TRUE) - mean_emp_fallow else NA_real_
  forest_bias <- if (sum(is.finite(forest_vals)) >= 2) mean(forest_vals, na.rm = TRUE) - mean_emp_forest else NA_real_
  
  tibble(
    fallow_bias = fallow_bias,
    forest_bias = forest_bias,
    n_fallow    = sum(is.finite(fallow_vals)),
    n_forest    = sum(is.finite(forest_vals))
  )
}

if (RUN_MC) {
  stopifnot(exists("post_bayes"))  # from 05_monte_carlo.R
  mc_bamboo <- bind_rows(lapply(seq_len(N_ITER_BAMBOO), function(i) {
    if (i %% 100 == 0 || i == 1) message("Running bamboo MC iteration: ", i, " / ", N_ITER_BAMBOO)
    run_one_mc_bamboo(
      bamboo_L4       = bamboo_L4,
      post_mat        = post_bayes,
      mean_emp_fallow = bamboo_fallow_mean,
      mean_emp_forest = bamboo_forest_mean
    ) %>% mutate(iter = i)
  }))
  saveRDS(mc_bamboo, OUT_MC_RDS)
} else {
  stopifnot(file.exists(OUT_MC_RDS))
  mc_bamboo <- readRDS(OUT_MC_RDS)
}

# Long format for plotting MC biases
mc_long_bias <- mc_bamboo %>%
  select(iter, fallow_bias, forest_bias) %>%
  pivot_longer(cols = c(fallow_bias, forest_bias),
               names_to = "group", values_to = "bias") %>%
  mutate(group = recode(group,
                        fallow_bias = "Fallow bamboo",
                        forest_bias = "Forest bamboo"))

# Naive (no misclassification) biases for reference
raw_bias_df <- tibble(
  group    = c("Fallow bamboo", "Forest bamboo"),
  raw_bias = c(
    mean(fallow_L4_clean, na.rm = TRUE) - bamboo_fallow_mean,
    mean(forest_L4_clean, na.rm = TRUE) - bamboo_forest_mean
  )
)

# Panel B2: MC distributions of bias vs naive bias
p_bamboo_mc_bias <- ggplot(
  mc_long_bias,
  aes(x = "GEDI L4A", y = bias, color = group, fill = group)
) +
  geom_jitter(
    width = 0.25,
    alpha = 0.4,
    size  = 2.5
  ) +
  geom_boxplot(
    width         = 0.25,
    outlier.shape = NA,
    alpha         = 1,
    linewidth     = 0.8,
    color         = "black"
  ) +
  geom_point(
    data        = raw_bias_df,
    aes(x = "GEDI L4A", y = raw_bias),
    inherit.aes = FALSE,
    shape       = 17,
    size        = 3,
    color       = "black"
  ) +
  facet_wrap(~ group, nrow = 1) +
  coord_cartesian(xlim = c(0.5, 1.5)) +
  scale_color_manual(values = pal_bamboo, guide = "none") +
  scale_fill_manual(values  = pal_bamboo, guide = "none") +
  labs(
    x = NULL,
    y = expression("Monte Carlo–simulated bias (GEDI L4A - empirical mean; Mg C ha"^{-1}*")"),
    title = "Misclassification-aware Carbon Bias",
    subtitle = paste0(
      "Colored boxes/jitter: posterior-resampled biases (n = ", N_ITER_BAMBOO,
      "); black triangles: naive biases (no misclassification correction)"
    )
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    strip.text         = element_text(face = "bold"),
    plot.margin        = ggplot2::margin(t = 10, r = 10, b = 25, l = 10, unit = "pt")
  )

print(p_bamboo_mc_bias)

## ============================================================
## B3) Summarize MC results for Table 2 (Panel C)
##     - Summaries of MC-only bias
##     - Total bias (MC + EBT discrepancy)
##     - Add naive (no correction) for comparison
## ============================================================

summ_mc <- function(x) {
  x <- x[is.finite(x)]
  tibble(
    mean  = mean(x, na.rm = TRUE),
    med   = quantile(x, 0.50,  na.rm = TRUE),
    q2.5  = quantile(x, 0.025, na.rm = TRUE),
    q97.5 = quantile(x, 0.975, na.rm = TRUE)
  )
}

# (C1) Misclassification-only bias: (GEDI L4A - empirical mean)
mc_panelC_misclass <- bind_rows(
  summ_mc(mc_bamboo$fallow_bias) %>%
    mutate(bamboo_context = "Fallow bamboo", component = "Misclassification-only"),
  summ_mc(mc_bamboo$forest_bias) %>%
    mutate(bamboo_context = "Forest bamboo", component = "Misclassification-only")
) %>%
  select(component, bamboo_context, mean, med, q2.5, q97.5)

# Naive (no misclassification correction) biases
raw_panelC <- raw_bias_df %>%
  rename(mean = raw_bias) %>%
  mutate(
    component = "Naive (no correction)",
    med = NA_real_, q2.5 = NA_real_, q97.5 = NA_real_
  ) %>%
  select(component, bamboo_context = group, mean, med, q2.5, q97.5)

# (C2) Total bias used for mapping: add EBT re-estimation discrepancy (mean residual)
mc_bamboo_total <- mc_bamboo %>%
  mutate(
    fallow_bias_total = fallow_bias + ebt_bias_bamboo_mean,
    forest_bias_total = forest_bias + ebt_bias_bamboo_mean
  )

mc_panelC_total <- bind_rows(
  summ_mc(mc_bamboo_total$fallow_bias_total) %>%
    mutate(bamboo_context = "Fallow bamboo",
           component      = "Total (MC + EBT discrepancy)"),
  summ_mc(mc_bamboo_total$forest_bias_total) %>%
    mutate(bamboo_context = "Forest bamboo",
           component      = "Total (MC + EBT discrepancy)")
) %>%
  select(component, bamboo_context, mean, med, q2.5, q97.5)

# Combine and format for reporting (Table 2, panel C)
table2_panelC <- bind_rows(
  raw_panelC,
  mc_panelC_misclass,
  mc_panelC_total
) %>%
  mutate(across(c(mean, med, q2.5, q97.5), ~ round(.x, 2)))

print(table2_panelC)

# Optional: summarize effective sample sizes per iteration (diagnostic)
n_diag <- mc_bamboo %>%
  summarise(
    n_fallow_med = median(n_fallow, na.rm = TRUE),
    n_fallow_q25 = quantile(n_fallow, 0.25, na.rm = TRUE),
    n_fallow_q75 = quantile(n_fallow, 0.75, na.rm = TRUE),
    n_forest_med = median(n_forest, na.rm = TRUE),
    n_forest_q25 = quantile(n_forest, 0.25, na.rm = TRUE),
    n_forest_q75 = quantile(n_forest, 0.75, na.rm = TRUE)
  )
print(n_diag)

## ============================================================
## B4) Combine MC bias + EBT simplification bias => mean bias raster
##     - Compute mean total bias for fallow/forest bamboo
##     - Map those biases back to bamboo_type_raster
## ============================================================

mc_bamboo_total <- mc_bamboo %>%
  mutate(
    fallow_bias_total = fallow_bias + ebt_bias_bamboo_mean,
    forest_bias_total = forest_bias + ebt_bias_bamboo_mean
  )

bias_mean_tbl <- mc_bamboo_total %>%
  summarise(
    fallow_bias_mean = mean(fallow_bias_total, na.rm = TRUE),
    forest_bias_mean = mean(forest_bias_total, na.rm = TRUE)
  )

print(bias_mean_tbl)

b_fallow_mean_tot <- bias_mean_tbl$fallow_bias_mean
b_forest_mean_tot <- bias_mean_tbl$forest_bias_mean

# Mean total bias raster for bamboo pixels (Mg C ha^-1)
bias_mean_rast <- ifel(
  bamboo_type_raster == 1, b_forest_mean_tot,
  ifel(bamboo_type_raster == 2, b_fallow_mean_tot, NA)
)

# Optional save for reuse
# writeRaster(bias_mean_rast, file.path("data", "bias_mean_rast.tif"), overwrite = TRUE)

## ============================================================
## B5) Map percent bias relative to EBT-estimated carbon
##     - percent_bias = (mean bias / carbon_est_rast) * 100
## ============================================================

percent_bias <- ifel(
  carbon_est_rast <= EPS,
  NA,
  (bias_mean_rast / carbon_est_rast) * 100
)

percent_bias_plot <- clamp(percent_bias, lower = 0, upper = PBIAS_CAP)

## Load the Xishuangbanna polygon shapefile
xbn <- vect("Xishuangbanna_Polygon/Xishuangbanna_Polygon.shp")
xbn <- project(xbn, crs(percent_bias_plot))

## ============================================================
## B6) Region-wide bamboo bias maps and zonal summaries
##     Uses:
##       - bias_mean_rast        (from B4)
##       - percent_bias_plot     (from B5)
##       - bamboo_type_raster    (1 = forest, 2 = fallow)
##       - carbon_est_rast       (EBT carbon, Mg C ha^-1)
##       - bamboo_forest_mean / bamboo_fallow_mean (Yuen means)
## ============================================================

# Quick map of EBT-based carbon
plot(carbon_est_rast, col = viridis(5))
plot(xbn, add = TRUE, border = "black", lwd = 1)

# Set plotting parameters (and keep old par)
par_old <- par(
  bg       = "transparent",
  fg       = "black",
  col.axis = "black",
  col.lab  = "black",
  col.main = "black"
)

## ------------------------------------------------------------
## B6.1  Map: mean carbon bias as % of EBT-estimated value
## ------------------------------------------------------------
plot(
  percent_bias_plot,
  col   = turbo(5),
  colNA = "transparent",
  axes  = FALSE,
  box   = FALSE,
  main  = "Mean Carbon Bias as % GEDI-Estimated Value"
)
plot(xbn, add = TRUE, border = "black", lwd = 1)

## ------------------------------------------------------------
## B6.2  Map: bamboo in forests vs fallow (empirical mean carbon)
## ------------------------------------------------------------

empirical_rast <- rast(bamboo_type_raster)
values(empirical_rast) <- NA_real_

# 1 = forest, 2 = fallow
empirical_rast[bamboo_type_raster == 1] <- bamboo_forest_mean
empirical_rast[bamboo_type_raster == 2] <- bamboo_fallow_mean

plot(
  empirical_rast,
  col   = c("#795C34", "#00CD10"),  # fallow / forest
  colNA = "transparent",
  axes  = FALSE,
  box   = FALSE,
  main  = "Bamboo in Forests vs Fallow (Empirical mean C)"
)
plot(xbn, add = TRUE, border = "black", lwd = 1)

## ------------------------------------------------------------
## B6.3  Map: pixel-wise bias (EBT carbon - empirical mean)
##       and zonal mean bias by bamboo type
## ------------------------------------------------------------

bias_rast <- carbon_est_rast - empirical_rast  # Mg C ha^-1

plot(
  bias_rast,
  col   = viridis(5),
  colNA = "transparent",
  axes  = FALSE,
  box   = FALSE,
  main  = "Carbon Bias (EBT - empirical mean, Mg C ha^-1)"
)
plot(xbn, add = TRUE, border = "black", lwd = 1)

# Zonal mean bias for forest vs fallow bamboo
z_bias <- zonal(
  bias_rast,
  bamboo_type_raster,
  fun   = "mean",
  na.rm = TRUE
)
z_bias$type <- ifelse(z_bias$class == 1, "forest", "fallow")
z_bias <- z_bias[order(z_bias$class), ]

print(z_bias)

# Restore plotting parameters
par(par_old)

## ============================================================
## B7) Summary plots: EBT carbon estimates & areas (prior vs posterior)
##     - B7a: EBT-based carbon distributions for trees / forest / fallow bamboo
##     - B7b: Areas of trees / forest / fallow bamboo (ha), prior vs posterior
## ============================================================

## ---------------------------
## B7a. Carbon estimates (EBT) for trees / forest / fallow bamboo
## ---------------------------

# 3-class raster: 1 = Trees, 2 = Forest bamboo, 3 = Fallow bamboo
class3 <- rast(prediction_raster)
values(class3) <- NA_integer_

# Trees: RF class = 6
class3[prediction_raster == 6] <- 1
# Forest bamboo: bamboo_type_raster == 1
class3[bamboo_type_raster == 1] <- 2
# Fallow bamboo: bamboo_type_raster == 2
class3[bamboo_type_raster == 2] <- 3

class_labels <- c("Trees", "Forest bamboo", "Fallow bamboo")

# (1) Sample size (n) per class3 class
freq_class3 <- freq(class3)
# freq_class3: columns "value" (class id) and "count"

n_df <- as.data.frame(freq_class3) |>
  transmute(
    class_id = value,
    n        = count,
    class    = factor(class_labels[class_id], levels = class_labels)
  )

# (2) Mean and SD calculated using zonal statistics on the raster
z_mean <- zonal(carbon_est_rast, class3, fun = "mean", na.rm = TRUE)
z_sum  <- zonal(carbon_est_rast,        class3, fun = "sum", na.rm = TRUE)
z_sum2 <- zonal(carbon_est_rast^2,      class3, fun = "sum", na.rm = TRUE)

# (3) Merge mean, sum, sum2, n -> compute SD (avoid full in-memory raster -> df)
stats_df <- data.frame(
  zone = z_mean[, 1],
  mean = z_mean[, 2],
  sum  = z_sum[,  2],
  sum2 = z_sum2[, 2]
)

stats_df$n <- freq_class3[match(stats_df$zone, freq_class3[, "value"]), "count"]
stats_df <- stats_df |>
  mutate(
    var   = pmax((sum2 / n) - (mean^2), 0),
    sd    = sqrt(var),
    class = factor(class_labels[zone], levels = class_labels)
  ) |>
  select(class, mean, sd, n)

stats_df

# Sample a manageable subset of pixels for plotting / quantiles
set.seed(123)
stack_3 <- c(carbon_est_rast, class3)
names(stack_3) <- c("carbon", "class_id")
n_target <- 5000
n_target <- min(n_target, ncell(stack_3))

samp <- spatSample(
  stack_3,
  size         = n_target,
  method       = "random",
  na.rm        = TRUE
)

samp <- samp %>%
  filter(!is.na(class_id)) %>%
  mutate(class = factor(class_labels[class_id],
                        levels = class_labels))

# Quantiles of EBT-based carbon per class (from sample)
q_stats <- samp %>%
  group_by(class) %>%
  summarise(
    q2.5   = quantile(carbon, 0.025, na.rm = TRUE),
    q25    = quantile(carbon, 0.25,  na.rm = TRUE),
    median = quantile(carbon, 0.50,  na.rm = TRUE),
    q75    = quantile(carbon, 0.75,  na.rm = TRUE),
    q97.5  = quantile(carbon, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

carbon_stats_B7 <- stats_df %>%
  left_join(q_stats, by = "class")

print(carbon_stats_B7)

# Panel B7a: EBT-based carbon distributions (boxplot + jitter)
p_B7_carbon <- ggplot(samp, aes(x = class, y = carbon)) +
  geom_jitter(
    aes(color = class),
    width = 0.4,
    size  = 2.0,
    alpha = 0.2
  ) +
  geom_boxplot(
    aes(fill = class),
    outlier.shape = NA,
    alpha = 0.8,
    linewidth = 0.6
  ) +
  labs(
    x = NULL,
    y = expression("Carbon density (EBT model, Mg C ha"^{-1}*")"),
    title = "EBT-based carbon estimates"
  ) +
  scale_fill_manual(
    values = c(
      "Trees"          = "#ac03c9",
      "Forest bamboo"  = "#00CD10",
      "Fallow bamboo"  = "#795C34"
    ),
    guide = "none"
  ) +
  scale_color_manual(
    values = c(
      "Trees"          = "#ac03c9",
      "Forest bamboo"  = "#00CD10",
      "Fallow bamboo"  = "#795C34"
    ),
    guide = "none"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 15, hjust = 1)
  )

print(p_B7_carbon)
