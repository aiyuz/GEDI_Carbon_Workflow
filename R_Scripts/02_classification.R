## 02_classifier.R
##
## Purpose:
##   - Use the compiled feature stack (xsbn_final.tif) and polygon labels
##     (training_data.rds) to:
##       1. Explore class separability in feature space (Figure 2a).
##       2. Train and compare multiple classifiers (RF, 1D-CNN, ANN, SVM).
##       3. Apply the best-performing classifier (RF) to map XSBN at 10 m.
##       4. Derive class areas and correct them via Bayesian updating
##          using an external confusion matrix.
##       5. Visualize final classification and performance summaries.
##
## Inputs (produced by 00_functions.R + 01_compiling_rasters.R):
##   - xsbn_final.tif        : feature stack for Oct–Dec (bands + indices + texture)
##   - training_data.rds     : per-pixel training data with 'label' column
##   - (optional) rf_model.csv: RF per-class performance table (for Figure on performance)
##
## Outputs:
##   - Figure_2a_feature_space_ridges.png : Ridge plot of key predictors by class
##   - Figure 2b: rf_prediction.tif              : 10-m RF classification map
##   - Figure 2c: performance chart of sensitivity and specificity
##   - KML samples (class-wise) for visual validation
##   - Area tables (observed vs corrected) per class
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 0. Libraries & Setup
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# Helper functions (NDVI/SI/BI/BPCI/MBI, GLCM helper, sample_class, etc.)
source("00_functions.R")

library(terra)          # rasters
library(sf)             # vector labels if needed
library(dplyr)          # data wrangling
library(tidyr)          # pivoting data
library(ggplot2)        # plotting
library(caret)          # train/test split + SVM & ANN via caret
library(randomForest)   # RF
library(e1071)          # SVM backend
library(keras)          # CNN + ANN
library(tensorflow)
library(patchwork)
library(scales)
library(ggspatial)
library(tidyterra)
library(ggridges)
library(tibble)

# If you haven't configured TensorFlow/Keras, run once in R:
#   install_keras()

set.seed(42)  # Reproducibility

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 1. Load data: feature stack & training data
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# Feature stack (Oct–Dec Sentinel-2 + indices + GLCM + canopy height)
xsbn_final    <- rast("xsbn_final.tif")

# Training data extracted under polygons in 01_compiling_rasters.R
training_data <- readRDS("training_data.rds")

# Consistent color scheme for classes (must match class names & order)
colors <- c(
  bamboo       = "#00CD10",  # neon green
  trees        = "#ac03c9",  # purple-ish for forest
  farm         = "#E69F00",  # warm orange
  plantation   = "#F0E442",  # yellow
  rubber       = "#CC79A7",  # magenta/pink
  construction = "#666666",  # dark grey
  water        = "#1E90FF"   # bright blue
)

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 2. Feature-space ridge plot (Figure 2a)
#    Goal: show that training classes occupy distinct regions of
#          feature space for four key predictors:
#          - B8 homogeneity (5×5)
#          - B8 dissimilarity (5×5)
#          - BPCI
#          - Canopy height
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# 2.1 Features to include ----------------------------------------------
feat_keep <- c(
  "B8_homogeneity_5x5",
  "B8_dissimilarity_5x5",
  "BPCI",
  "canopy"
)

# 2.2 Order labels for plotting ----------------------------------------
training_data <- training_data |>
  dplyr::mutate(
    label = factor(
      label,
      levels = c("bamboo", "trees", "plantation", "rubber",
                 "farm", "construction", "water")
    )
  )

# 2.3 Long format + nice feature labels -------------------------------
ridge_df <- training_data |>
  tidyr::pivot_longer(
    cols      = all_of(feat_keep),
    names_to  = "feature",
    values_to = "value"
  ) |>
  dplyr::mutate(
    feature = factor(
      feature,
      levels = feat_keep,
      labels = c(
        "B8 homogeneity (5×5)",
        "B8 dissimilarity (5×5)",
        "BPCI",
        "Canopy height"
      )
    )
  )

# 2.4 Standardize (z-score) within each feature ------------------------
ridge_df_std <- ridge_df |>
  dplyr::group_by(feature) |>
  dplyr::mutate(
    value_std = as.numeric(scale(value))
  ) |>
  dplyr::ungroup()

# 2.5 Color scheme for ridges (must match label levels) ----------------
cols_use <- c(
  "#00CD10",  # bamboo
  "purple",   # trees
  "#FFE5CC",  # plantation
  "#FFF68F",  # rubber
  "orange",   # farm
  "darkred",  # construction
  "#1E90FF"   # water
)

# 2.6 Ridge plot (Figure 2a) -------------------------------------------
fig_2a_ridges <- ggplot(
  ridge_df_std,
  aes(x = value_std, y = label, fill = label)
) +
  ggridges::geom_density_ridges(
    alpha          = 0.8,
    scale          = 1.1,
    color          = "black",
    rel_min_height = 0.01
  ) +
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    color      = "grey30",
    linewidth  = 0.3
  ) +
  scale_fill_manual(
    values = cols_use,
    breaks = levels(ridge_df_std$label)
  ) +
  facet_wrap(
    ~ feature,
    nrow   = 1,
    scales = "fixed"
  ) +
  labs(
    title = "Standardized distributions of key predictors by class (Figure 2a)",
    x     = "Standardized value (z-score)",
    y     = NULL,
    fill  = "Class"
  ) +
  scale_x_continuous(
    breaks = seq(-10, 10, by = 2)
  ) +
  coord_flip() +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid      = element_blank(),
    strip.text      = element_text(face = "bold"),
    legend.position = "bottom"
  )

# print to screen (optional)
fig_2a_ridges

# save to file
ggsave(
  filename = "Figure_2a_feature_space_ridges.png",
  plot     = fig_2a_ridges,
  width    = 8,
  height   = 4,
  dpi      = 300
)

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 3. Train / validation split
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# Ensure label is factor
training_data$label <- as.factor(training_data$label)

set.seed(42)
train_index <- createDataPartition(training_data$label, p = 0.7, list = FALSE)
train_set   <- training_data[train_index, ]
valid_set   <- training_data[-train_index, ]

# For convenience: feature matrices & label vectors
x_train <- train_set %>% dplyr::select(-label)
y_train <- train_set$label

x_valid <- valid_set %>% dplyr::select(-label)
y_valid <- valid_set$label

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 4. RANDOM FOREST (RF)
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

rf_model <- randomForest(
  label ~ .,
  data       = train_set,
  ntree      = 100,
  mtry       = 4,
  importance = TRUE
)

# Predict on validation set
valid_pred_rf <- predict(rf_model, valid_set)

# Confusion matrix + accuracy
rf_conf_mat  <- confusionMatrix(valid_pred_rf, y_valid)
rf_accuracy  <- rf_conf_mat$overall["Accuracy"] |> as.numeric()

cat("RF accuracy:", round(rf_accuracy, 4), "\n")

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 5. CONVOLUTIONAL NEURAL NETWORK (1D-CNN)
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# 5.1 One-hot labels ----------------------------------------------------
train_labels <- as.numeric(as.factor(train_set$label)) - 1
valid_labels <- as.numeric(as.factor(valid_set$label)) - 1

y_train_cnn <- tf$keras$utils$to_categorical(train_labels)
y_valid_cnn <- tf$keras$utils$to_categorical(valid_labels)

num_classes <- dim(y_train_cnn)[2]

# 5.2 Feature matrices (scaled) ----------------------------------------
x_train_mat <- as.matrix(x_train)
x_valid_mat <- as.matrix(x_valid)

x_train_mat <- scale(x_train_mat)
x_valid_mat <- scale(
  x_valid_mat,
  center = attr(x_train_mat, "scaled:center"),
  scale  = attr(x_train_mat, "scaled:scale")
)

# 5.3 Reshape for 1D CNN ------------------------------------------------
num_features <- as.integer(ncol(x_train_mat))

x_train_reshaped <- array_reshape(
  x_train_mat,
  c(as.integer(nrow(x_train_mat)), num_features, 1L)
)
x_valid_reshaped <- array_reshape(
  x_valid_mat,
  c(as.integer(nrow(x_valid_mat)), num_features, 1L)
)
input_shape <- c(num_features, 1L)

# 5.4 Build CNN ---------------------------------------------------------
cnn_model <- keras_model_sequential()
cnn_model$add(
  layer_conv_1d(
    filters     = 32L,
    kernel_size = 3L,
    activation  = "relu",
    input_shape = input_shape
  )
)
cnn_model$add(layer_max_pooling_1d(pool_size = 2L))
cnn_model$add(layer_dropout(rate = 0.3))
cnn_model$add(
  layer_conv_1d(
    filters     = 64L,
    kernel_size = 3L,
    activation  = "relu"
  )
)
cnn_model$add(layer_max_pooling_1d(pool_size = 2L))
cnn_model$add(layer_dropout(rate = 0.3))
cnn_model$add(layer_flatten())
cnn_model$add(layer_dense(units = 256L, activation = "relu"))
cnn_model$add(layer_dropout(rate = 0.5))
cnn_model$add(layer_dense(units = num_classes, activation = "softmax"))

cnn_model$compile(
  loss      = "categorical_crossentropy",
  optimizer = optimizer_adam(learning_rate = 0.0005),
  metrics   = list("accuracy")
)

early_stop <- callback_early_stopping(monitor = "val_loss", patience = 5)
lr_schedule <- function(epoch, lr) {
  if (epoch > 10) return(lr * 0.5)
  lr
}
lr_callback <- callback_learning_rate_scheduler(schedule = lr_schedule)

cnn_history <- cnn_model$fit(
  x_train_reshaped, y_train_cnn,
  epochs          = as.integer(50),
  batch_size      = as.integer(32),
  validation_data = list(x_valid_reshaped, y_valid_cnn),
  callbacks       = list(early_stop, lr_callback)
)

cnn_eval     <- cnn_model$evaluate(x_valid_reshaped, y_valid_cnn)
cnn_loss     <- cnn_eval[1]
cnn_accuracy <- cnn_eval[2]
cat("CNN accuracy:", round(cnn_accuracy, 4), "\n")

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 6. Simple feed-forward ANN (Dense)
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

ann_model <- keras_model_sequential()
ann_model$add(
  layer_dense(
    units       = 64,
    activation  = "relu",
    input_shape = c(num_features)
  )
)
ann_model$add(layer_dropout(rate = 0.3))
ann_model$add(layer_dense(units = 32, activation = "relu"))
ann_model$add(layer_dropout(rate = 0.3))
ann_model$add(layer_dense(units = num_classes, activation = "softmax"))

ann_model$compile(
  optimizer = optimizer_adam(learning_rate = 0.001),
  loss      = "categorical_crossentropy",
  metrics   = list("accuracy")
)

ann_history <- ann_model$fit(
  x               = x_train_mat,
  y               = y_train_cnn,
  epochs          = as.integer(50),
  batch_size      = as.integer(32),
  validation_data = list(x_valid_mat, y_valid_cnn),
  callbacks       = list(early_stop)
)

ann_eval     <- ann_model$evaluate(x_valid_mat, y_valid_cnn)
ann_loss     <- ann_eval[1]
ann_accuracy <- ann_eval[2]
cat("ANN accuracy:", round(ann_accuracy, 4), "\n")

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 7. SVM (linear) via caret
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

ctrl <- trainControl(method = "none")

svm_model <- caret::train(
  label ~ .,
  data       = train_set,
  method     = "svmLinear",
  preProcess = c("center", "scale"),
  trControl  = ctrl
)

valid_pred_svm    <- predict(svm_model, valid_set)
all_levels        <- levels(train_set$label)
valid_pred_factor <- factor(valid_pred_svm, levels = all_levels)
y_ref_factor      <- factor(valid_set$label, levels = all_levels)
svm_conf_mat      <- confusionMatrix(
  data      = valid_pred_factor,
  reference = y_ref_factor
)
svm_accuracy <- svm_conf_mat$overall["Accuracy"]
cat("SVM accuracy:", round(svm_accuracy, 4), "\n")

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 8. Collect & plot model performance (optional figure)
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

results_df <- tibble(
  model    = c("Random Forest", "CNN (1D)", "ANN (Dense)", "SVM (Linear)"),
  accuracy = c(rf_accuracy, cnn_accuracy, ann_accuracy, svm_accuracy)
)

p_models <- ggplot(results_df, aes(x = model, y = accuracy, fill = model)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(
    aes(label = scales::percent(accuracy, accuracy = 1)),
    vjust = -0.5,
    size  = 3
  ) +
  labs(
    title = "Comparison of classifier accuracy on validation set",
    x     = "Model",
    y     = "Accuracy"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  theme_minimal() +
  theme(
    axis.text.x      = element_text(angle = 20, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# p_models  # print if needed

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 9. Final RF prediction over XSBN
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

predicted_raster <- predict(xsbn_final, rf_model, progress = 1)

# save classification raster
writeRaster(predicted_raster, "rf_prediction.tif", overwrite = TRUE)

# If needed to reload:
# predicted_raster <- rast("rf_prediction.tif")

# Quick base plot (optional)
# plot(predicted_raster)

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 10. Sample prediction results for visual validation (export KMLs)
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

rat <- data.frame(
  ID       = 1:7,
  category = c("bamboo",
               "construction",
               "farms",
               "plantation",
               "rubber",
               "trees",
               "water"),
  stringsAsFactors = FALSE
)

class_raster_fac <- as.factor(predicted_raster)
levels(class_raster_fac) <- rat
pred_int_raster <- class_raster_fac  # used by sample_class in 00_functions.R

# Sample up to 100 pixels per class for KML export
for (i in 1:7) {
  sample_class(rat$ID[i], rat$category[i], n_samples = 100)
}

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 11. Confusion matrix & Bayesian area correction
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Data below were collected by visually inspecting Maxar maps in 2020-2022
# Confusion matrix C: columns = predicted, rows = “true”
C <- matrix(
  c(
    # predicted: bamboo, construction, farms, plantation, rubber, trees, water
    # true = bamboo
    92, 1, 1, 2, 4, 6, 0,
    # true = construction
    0, 76, 3, 0, 0, 0, 5,
    # true = farms
    0, 14, 88, 7, 0, 0, 2,
    # true = plantation
    3, 7, 8, 78, 5, 0, 0,
    # true = rubber
    5, 2, 0, 9, 76, 2, 0,
    # true = trees
    0, 0, 0, 4, 15, 92, 1,
    # true = water
    0, 0, 0, 0, 0, 0, 92
  ),
  nrow = 7,
  byrow = TRUE
)

# 11.1 Compute per-cell area in m²
area_rast <- cellSize(predicted_raster)  # area at each pixel

# 11.2 Sum areas by integer class code
area_by_code <- terra::zonal(
  x   = area_rast,
  z   = predicted_raster,
  fun = "sum",
  na.rm = TRUE
)

# Make sure columns are named consistently
# typically: first column = zone, second = sum
colnames(area_by_code) <- c("class_id", "area_m2")

area_by_code <- as_tibble(area_by_code) |>
  dplyr::mutate(class_id = as.integer(class_id))

# 11.3 Map integer IDs to class names via raster attribute table (RAT)
cat_tbl <- cats(predicted_raster)[[1]]
colnames(cat_tbl) <- c("class_id", "class")

area_df <- dplyr::left_join(area_by_code, cat_tbl, by = "class_id") |>
  dplyr::mutate(
    area_ha  = area_m2 / 1e4,
    area_km2 = area_m2 / 1e6
  ) |>
  dplyr::select(class, area_m2, area_ha, area_km2)

# In case class names are mis-ordered, enforce them:
area_df$class <- factor(
  area_df$class,
  levels = c("bamboo", "construction", "farms",
             "plantation", "rubber", "trees", "water")
)

# 11.4 Bayesian updating: redistribute predicted areas into “true” areas
P_true_given_pred <- C / colSums(C)  # columns sum to 1

rownames(P_true_given_pred) <- levels(area_df$class)
colnames(P_true_given_pred) <- levels(area_df$class)

area_pred <- setNames(area_df$area_km2, area_df$class)
area_true_est <- as.numeric(P_true_given_pred %*% area_pred)

area_true_df <- data.frame(
  class        = levels(area_df$class),
  est_area_km2 = area_true_est,
  row.names    = NULL
)

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 12. Visualization of classification map & area corrections
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# 12.1 Classification map (prediction_map) ------------------------------
fact_rast <- as.factor(predicted_raster)
levels(fact_rast) <- rat  # ensure ID→category mapping

prediction_map <- ggplot() +
  geom_spatraster(data = fact_rast) +
  scale_fill_manual(
    name    = "Class",
    values  = colors,
    na.value = "transparent",
    drop    = FALSE
  ) +
  annotation_scale(
    location    = "bl",
    style       = "bar",
    width_hint  = 0.2,
    text_cex    = 1,
    pad_x       = unit(0.5, "cm"),
    pad_y       = unit(0.5, "cm")
  ) +
  annotation_north_arrow(
    location      = "tl",
    which_north   = "true",
    pad_x         = unit(0.5, "cm"),
    pad_y         = unit(0.5, "cm"),
    style         = north_arrow_fancy_orienteering
  ) +
  coord_sf(crs = crs(predicted_raster)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text        = element_text(size = 16),
    axis.title       = element_text(size = 18),
    legend.title     = element_text(size = 16),
    legend.text      = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# 12.2 Observed vs corrected area by class (p_area) --------------------
area_comp_df <- cbind(area_true_df, obs_area_km2 = area_df$area_km2)
colnames(area_comp_df) <- c("class", "corrected_area_km2", "obs_area_km2")

area_comp_df <- area_comp_df |>
  arrange(corrected_area_km2) |>
  mutate(class = factor(class, levels = class))

plot_df <- area_comp_df |>
  tidyr::pivot_longer(
    cols      = c(corrected_area_km2, obs_area_km2),
    names_to  = "type",
    values_to = "area"
  ) |>
  dplyr::mutate(
    type = recode(
      type,
      corrected_area_km2 = "True",
      obs_area_km2       = "Observed"
    ),
    fill_group = ifelse(type == "Observed", "Observed", as.character(class)),
    type       = factor(type, levels = c("True", "Observed"))
  )

fill_vals2 <- c(colors, Observed = alpha("black", 0.2))

p_area <- ggplot(
  plot_df,
  aes(
    y     = class,
    x     = area,
    fill  = fill_group,
    group = type
  )
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(
    name   = "",
    values = fill_vals2,
    drop   = FALSE
  ) +
  labs(
    title = "Observed vs. corrected area by class",
    x     = "Area (km²)",
    y     = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    axis.text.y      = element_text(size = 14),
    axis.text.x      = element_text(size = 12),
    axis.title.x     = element_text(size = 14),
    legend.position  = "right"
  )

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 13. Model vs visual validation: Sensitivity & Specificity (p_perf)
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# Requires rf_model.csv exported separately (per-class metrics from caret or custom)
rf_model_perf <- read.csv("rf_model.csv")
# overwrite class names with area_true_df$class to ensure consistent order
rf_model_perf$X <- area_true_df$class

perf_df <- rf_model_perf[, c(1:3, 6, 12)]
colnames(perf_df) <- c("class", "Sensitivity", "Specificity",
                       "Precision", "Balanced.Accuracy")

df_model <- perf_df |>
  tibble::as_tibble() |>
  dplyr::select(class, Sensitivity, Specificity) |>
  tidyr::pivot_longer(
    cols      = c(Sensitivity, Specificity),
    names_to  = "metric",
    values_to = "value"
  ) |>
  dplyr::mutate(
    source = "Model",
    class  = factor(class, levels = perf_df$class)
  )

# Manually defined confusion matrix for visual validation
cm <- matrix(
  c(92,0,0,3,5,0,0,
    1,76,14,7,2,0,0,
    1,3,88,8,0,0,0,
    2,0,7,78,9,4,0,
    4,0,0,5,76,15,0,
    6,0,0,0,2,92,0,
    0,5,2,0,0,1,92),
  nrow = 7, byrow = FALSE,
  dimnames = list(
    actual    = c("bamboo","construction","farms","plantation","rubber","trees","water"),
    predicted = c("bamboo","construction","farms","plantation","rubber","trees","water")
  )
)

total <- sum(cm)

metrics <- lapply(rownames(cm), function(cls) {
  TP <- cm[cls, cls]
  FN <- sum(cm[cls, ]) - TP
  FP <- sum(cm[, cls]) - TP
  TN <- total - TP - FN - FP
  data.frame(
    class             = cls,
    TP                = TP,
    FN                = FN,
    FP                = FP,
    TN                = TN,
    Sensitivity       = TP / (TP + FN),
    Specificity       = TN / (TN + FP),
    FPR               = FP / (FP + TN),
    FNR               = FN / (FN + TP),
    UsersAccuracy     = TP / (TP + FP),
    ProducersAccuracy = TP / (TP + FN),
    row.names         = NULL
  )
}) %>% bind_rows()

df_valid <- metrics |>
  tibble::as_tibble() |>
  dplyr::select(class, Sensitivity, Specificity) |>
  tidyr::pivot_longer(
    cols      = c(Sensitivity, Specificity),
    names_to  = "metric",
    values_to = "value"
  ) |>
  dplyr::mutate(
    source = "Validation",
    class  = factor(class, levels = perf_df$class)
  )

comp_perf_df <- bind_rows(df_model, df_valid) |>
  dplyr::mutate(
    class  = factor(class, levels = unique(class)),
    metric = factor(metric, levels = c("Sensitivity","Specificity")),
    source = factor(source, levels = c("Model","Validation"))
  )

p_perf <- ggplot() +
  geom_point(
    data  = dplyr::filter(comp_perf_df, source == "Model"),
    aes(x = value, y = class, shape = metric),
    color  = "grey",
    size   = 4.5,
    stroke = 1
  ) +
  geom_point(
    data = dplyr::filter(comp_perf_df, source == "Validation"),
    aes(x = value, y = class, shape = metric, color = class),
    size = 5
  ) +
  scale_shape_manual(
    name   = "Metric",
    values = c(Sensitivity = 16, Specificity = 17)
  ) +
  scale_color_manual(
    name   = "Class",
    values = colors
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0.7, 1.0),
    name   = "Rate"
  ) +
  facet_wrap(~ metric, nrow = 2) +
  labs(
    title = "Sensitivity & Specificity: Model vs. visual validation",
    y     = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_line(color = "grey80"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 12),
    axis.text.x        = element_text(size = 12),
    axis.title.x       = element_text(size = 14),
    legend.position    = "right"
  )

# End of 02_classifier.R
