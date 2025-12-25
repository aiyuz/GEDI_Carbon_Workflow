## ------------------------------------------------------------
## 05_monte_carlo.R
## Misclassification-aware robustness check (ID-free)
##   - Uses confusion matrix + raster-derived priors
##   - Builds posterior P(actual | predicted)
##   - Re-samples "true" class for each ROW once per iteration
##   - Repeats global bamboo vs trees tests for 4 structure metrics
##   - Compares raw vs corrected effects
##   - Visualizes MC effect distributions
## ------------------------------------------------------------

library(terra)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)

set.seed(42)

## -------------------------
## 0) Class order (must match CM + RF map)
## -------------------------
classes <- c("bamboo", "construction", "farm",
             "plantation", "rubber", "trees", "water")

## -------------------------
## 1) Confusion matrix
##    rows = actual (true), cols = predicted (RF)
## -------------------------
cm <- matrix(
  c(
    92,0,0,3,5,0,0,
    1,76,14,7,2,0,0,
    1,3,88,8,0,0,0,
    2,0,7,78,9,4,0,
    4,0,0,5,76,15,0,
    6,0,0,0,2,92,0,
    0,5,2,0,0,1,92
  ),
  nrow = 7, byrow = FALSE,
  dimnames = list(actual = classes, predicted = classes)
)

## Likelihood: P(predicted | actual)
lik <- prop.table(cm, margin = 1)
lik <- lik[classes, classes]


## -------------------------
## 2) Prior from RF raster
##    Computes P(actual) using predicted class pixel frequencies
## -------------------------
get_prior_from_raster <- function(predicted_raster, classes) {
  
  freq_tbl <- terra::freq(predicted_raster)
  freq_tbl <- as.data.frame(freq_tbl)
  
  # terra::freq() may return:
  # (layer, value, count) OR (value, count)
  if (all(c("layer","value","count") %in% names(freq_tbl))) {
    freq_tbl <- freq_tbl %>% dplyr::select(value, count)
  } else if (all(c("value","count") %in% names(freq_tbl))) {
    freq_tbl <- freq_tbl %>% dplyr::select(value, count)
  } else {
    stop("Unexpected terra::freq() output structure.")
  }
  
  freq_tbl <- freq_tbl %>% dplyr::filter(!is.na(value))
  
  # Case A: value is class name (character/factor)
  if (is.character(freq_tbl$value) || is.factor(freq_tbl$value)) {
    
    prior_df <- freq_tbl %>%
      mutate(
        category = as.character(value),
        n_pix    = as.numeric(count)
      ) %>%
      group_by(category) %>%
      summarise(n_pix = sum(n_pix), .groups = "drop")
    
  } else {
    # Case B: numeric class IDs
    prior_df <- freq_tbl %>%
      mutate(
        class_id = as.integer(value),
        n_pix    = as.numeric(count),
        category = classes[class_id]
      ) %>%
      filter(!is.na(category)) %>%
      group_by(category) %>%
      summarise(n_pix = sum(n_pix), .groups = "drop")
  }
  
  prior_df <- prior_df %>%
    filter(category %in% classes) %>%
    mutate(category = factor(category, levels = classes)) %>%
    arrange(category) %>%
    mutate(prior = n_pix / sum(n_pix))
  
  prior_vec <- setNames(prior_df$prior, as.character(prior_df$category))
  prior_vec <- prior_vec[classes]
  
  list(prior_df = prior_df, prior_vec = prior_vec)
}

# Requires your RF map object:
# predicted_raster <- rast("rf_prediction.tif")
pri_obj   <- get_prior_from_raster(predicted_raster, classes)
prior_df  <- pri_obj$prior_df
prior_vec <- pri_obj$prior_vec

print(prior_df)


## -------------------------
## 3) Posterior matrix
##    P(actual | predicted) ∝ P(pred | actual) * P(actual)
## -------------------------
post_bayes <- sapply(classes, function(pred) {
  unnorm <- lik[, pred] * prior_vec
  unnorm / sum(unnorm)
})

post_bayes <- as.matrix(post_bayes)
rownames(post_bayes) <- classes
colnames(post_bayes) <- classes


## -------------------------
## 4) Row-level sampler
## -------------------------
sample_true_class <- function(pred_class, post_mat) {
  if (is.na(pred_class) || !(pred_class %in% colnames(post_mat))) {
    return(NA_character_)
  }
  probs <- post_mat[, pred_class]
  sample(rownames(post_mat), size = 1, prob = probs)
}


## -------------------------
## 5) One MC iteration (ID-free)
##    Each ROW draws ONCE based on its predicted class
## -------------------------
run_one_mc_rowwise <- function(df_long, post_mat) {
  
  df_draw <- df_long %>%
    mutate(
      land_class = as.character(land_class),
      draw_class = vapply(
        land_class,
        sample_true_class,
        FUN.VALUE = character(1),
        post_mat = post_mat
      )
    )
  
  # Only keep posterior bamboo vs trees for comparison
  df_bt <- df_draw %>%
    filter(draw_class %in% c("bamboo","trees")) %>%
    mutate(draw_class = factor(draw_class, levels = c("bamboo","trees")))
  
  if (nrow(df_bt) < 10) {
    return(tibble(
      metric = unique(df_long$metric),
      n_bamboo = NA_integer_,
      n_trees  = NA_integer_,
      median_bamboo = NA_real_,
      median_trees  = NA_real_,
      delta_median  = NA_real_,
      p_wilcox = NA_real_
    ))
  }
  
  df_bt %>%
    group_by(metric) %>%
    summarise(
      n_bamboo = sum(draw_class == "bamboo"),
      n_trees  = sum(draw_class == "trees"),
      median_bamboo = median(value[draw_class == "bamboo"], na.rm = TRUE),
      median_trees  = median(value[draw_class == "trees"],  na.rm = TRUE),
      delta_median  = median_trees - median_bamboo,
      p_wilcox = ifelse(
        n_bamboo > 5 & n_trees > 5,
        wilcox.test(value ~ draw_class)$p.value,
        NA_real_
      ),
      .groups = "drop"
    )
}


## -------------------------
## 6) Prepare global metrics table
##    all_metrics should include:
##      land_class, value, metric
## -------------------------
stopifnot(exists("all_metrics"))

target_metrics <- c(
  "Canopy Height (m)",
  "Aboveground Biomass Density (Mg/ha)",
  "Plant Area Index (m²/m²)",
  "Total Vegetation Volume (m³/m²)"
)

bt_global <- all_metrics %>%
  filter(metric %in% target_metrics) %>%
  filter(!is.na(value)) %>%
  filter(land_class %in% classes)


## -------------------------
## 7) Monte Carlo runs
## -------------------------
n_iter <- 1000

mc_res <- bind_rows(lapply(seq_len(n_iter), function(i) {
  if (i %% 10 == 0 || i == 1) {
    message("Running MC iteration: ", i, " / ", n_iter)
  }
  
  run_one_mc_rowwise(bt_global, post_bayes) %>%
    mutate(iter = i)
}))


## -------------------------
## 8) MC summary (effect + robustness)
## -------------------------
mc_summary_full <- mc_res %>%
  group_by(metric) %>%
  summarise(
    delta_med_2.5  = quantile(delta_median, 0.025, na.rm = TRUE),
    delta_med_50   = quantile(delta_median, 0.50,  na.rm = TRUE),
    delta_med_97.5 = quantile(delta_median, 0.975, na.rm = TRUE),
    prop_p_lt_0.05 = mean(p_wilcox < 0.05, na.rm = TRUE),
    n_bamboo_med   = round(median(n_bamboo, na.rm = TRUE)),
    n_trees_med    = round(median(n_trees,  na.rm = TRUE)),
    .groups = "drop"
  )

print(mc_summary_full)


## ------------------------------------------------------------
## 9) RAW (no correction) bamboo vs trees
## ------------------------------------------------------------
raw_bt <- bt_global %>%
  filter(land_class %in% c("bamboo","trees")) %>%
  mutate(land_class = factor(land_class, levels = c("bamboo","trees")))

raw_summary <- raw_bt %>%
  group_by(metric) %>%
  summarise(
    n_bamboo = sum(land_class == "bamboo"),
    n_trees  = sum(land_class == "trees"),
    median_bamboo = median(value[land_class == "bamboo"], na.rm = TRUE),
    median_trees  = median(value[land_class == "trees"],  na.rm = TRUE),
    delta_median_raw = median_trees - median_bamboo,
    p_wilcox_raw = wilcox.test(value ~ land_class)$p.value,
    .groups = "drop"
  )

raw_summary


## ------------------------------------------------------------
## 10) RAW vs MC comparison table
## ------------------------------------------------------------
compare_table <- raw_summary %>%
  left_join(mc_summary_full, by = "metric") %>%
  mutate(
    delta_shift = delta_med_50 - delta_median_raw
  )

compare_table


## ------------------------------------------------------------
## 11) Visualization: MC delta distributions
##     + overlay raw delta
## ------------------------------------------------------------

plot_mc_box <- mc_res %>%
  filter(!is.na(delta_median)) %>%
  mutate(metric = factor(metric, levels = target_metrics))

raw_pts <- raw_summary %>%
  transmute(metric, delta_median_raw) %>%
  mutate(metric = factor(metric, levels = target_metrics))

p_mc_box <- ggplot(plot_mc_box, aes(x = metric, y = delta_median)) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    alpha = 0.6,
    linewidth = 0.6
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.15,
    size = 0.9
  ) +
  geom_point(
    data = raw_pts,
    aes(x = metric, y = delta_median_raw),
    inherit.aes = FALSE,
    shape = 17,   # triangle
    size = 3
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 18)) +
  labs(
    x = NULL,
    y = "Median difference (trees − bamboo)",
    title = "Misclassification-corrected effects (Monte Carlo)",
    subtitle = paste0(
      "Box/jitter = posterior resampling (n=", n_iter,
      "); triangle = raw estimate"
    )
  ) +
  facet_wrap(~metric,scales="free_y")+
  theme_bw (base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 11)
  )

p_mc_box

## optional save
# write.csv(as.data.frame(post_bayes), "MC_posterior_P_actual_given_pred.csv")
# write.csv(compare_table, "MC_raw_vs_corrected_compare.csv", row.names = FALSE)
# ggsave("Fig_MC_box_raw_overlay.png", p_mc_box, width = 8, height = 4.5, dpi = 600)
