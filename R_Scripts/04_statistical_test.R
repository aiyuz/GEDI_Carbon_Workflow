## ------------------------------------------------------------
## 04_statistical_test.R
## Purpose:
##   Summarize and test bamboo vs trees structural/biomass differences
##   based on outputs produced in 03_gedi_overlay.R.
##
## Main outputs:
##   (1) Wilcoxon tests for four GEDI-derived metrics (bamboo vs trees)
##   (2) Per-height-bin Wilcoxon tests for PAVD and PAI contribution
##   (3) Tree cover discrepancy tests:
##       - ANCOVA: tree_cover ~ land_class * agbd
##       - Residual test from pooled model
##   (4) Optional: slope/intercept table for Fig 8A per-class fits
##
## Required objects from 03_gedi_overlay:
##   - bt_metrics
##   - df_tc
##   - pavd_long
##   - pai_long
##   - pavd_columns, z_mid, bin_depth (if you rebuild long tables here)
##
## Packages
## ------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)

## ------------------------------------------------------------
## 1) Helper: prepare bamboo vs trees metric dataframe
## ------------------------------------------------------------

prep_bt_metric <- function(bt_metrics, metric_name) {
  bt_metrics %>%
    dplyr::filter(metric == metric_name) %>%
    dplyr::mutate(
      land_class = factor(land_class, levels = c("bamboo", "trees"))
    )
}

## ------------------------------------------------------------
## 2) Helper: run Wilcoxon test + short summary
## ------------------------------------------------------------

run_wilcox_bt <- function(df, label) {
  w <- wilcox.test(value ~ land_class, data = df)
  
  tibble::tibble(
    metric     = label,
    n_bamboo   = sum(df$land_class == "bamboo", na.rm = TRUE),
    n_trees    = sum(df$land_class == "trees",  na.rm = TRUE),
    median_bamboo = median(df$value[df$land_class == "bamboo"], na.rm = TRUE),
    median_trees  = median(df$value[df$land_class == "trees"],  na.rm = TRUE),
    p_wilcox   = w$p.value
  )
}

## ------------------------------------------------------------
## 3) (Global) Bamboo vs trees: 4 GEDI metrics
## ------------------------------------------------------------

bt_agbd <- prep_bt_metric(bt_metrics, "Aboveground Biomass Density (Mg/ha)")
bt_height <- prep_bt_metric(bt_metrics, "Canopy Height (m)")
bt_pai <- prep_bt_metric(bt_metrics, "Plant Area Index (m²/m²)")
bt_volume <- prep_bt_metric(bt_metrics, "Total Vegetation Volume (m³/m²)") 

global_tests <- dplyr::bind_rows(
  run_wilcox_bt(bt_agbd,   "AGBD"),
  run_wilcox_bt(bt_height, "Height"),
  run_wilcox_bt(bt_pai,    "PAI"),
  run_wilcox_bt(bt_volume, "Total vegetation volume")
) %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_wilcox, method = "BH"))

global_tests

## ------------------------------------------------------------
## 4) (Vertical) Per-bin tests for PAVD and PAI contribution
## ------------------------------------------------------------

# ---- PAVD per-bin tests ----
pavd_bin_tests <- pavd_long %>%
  dplyr::filter(
    land_class %in% c("bamboo","trees"),
    !is.na(pavd), pavd >= 0
  ) %>%
  dplyr::group_by(height_mid) %>%
  dplyr::summarise(
    p_wilcox = wilcox.test(pavd ~ land_class)$p.value,
    .groups  = "drop"
  ) %>%
  dplyr::mutate(p_adj = p.adjust(p_wilcox, method = "BH"))

pavd_bin_tests


# ---- PAI contribution per-bin tests ----
pai_bin_tests <- pai_long %>%
  dplyr::filter(
    land_class %in% c("bamboo","trees"),
    !is.na(pai_contrib), pai_contrib >= 0
  ) %>%
  dplyr::group_by(height_mid) %>%
  dplyr::summarise(
    p_wilcox = wilcox.test(pai_contrib ~ land_class)$p.value,
    .groups  = "drop"
  ) %>%
  dplyr::mutate(p_adj = p.adjust(p_wilcox, method = "BH"))

pai_bin_tests
## ------------------------------------------------------------
## 5) Tree cover discrepancy tests (Fig 8A/8B logic)
## ------------------------------------------------------------
## df_tc should include:
##   - land_class (bamboo / trees)
##   - tree_cover (0-100)
##   - agbd (>0)
# ensure factor order
df_tc <- df_tc %>%
  mutate(land_class = factor(land_class, levels = c("bamboo", "trees")))

# ------------------------------------------------------------
# (A) ANCOVA aligned with Fig 8A
#     Response = AGBD, predictor = tree_cover, class interaction
# ------------------------------------------------------------

m_tc_fig8a <- lm(agbd ~ land_class * tree_cover, data = df_tc)
a_tc_fig8a <- anova(m_tc_fig8a)

p_land <- a_tc_fig8a["land_class", "Pr(>F)"]
p_tc   <- a_tc_fig8a["tree_cover", "Pr(>F)"]
p_int  <- a_tc_fig8a["land_class:tree_cover", "Pr(>F)"]

lab_ancova <- paste0(
  "ANCOVA: AGBD ~ class * tree_cover",
  "\nclass p = ", signif(p_land, 3),
  "\ntree cover p = ", signif(p_tc, 3),
  "\ninteraction p = ", signif(p_int, 3)
)

a_tc_fig8a
summary(m_tc_fig8a)

# ------------------------------------------------------------
# (B) Extract equations for the two class-specific lines
#     (consistent with Fig 8A direction)
# ------------------------------------------------------------

b <- coef(m_tc_fig8a)

eq_bamboo <- paste0(
  "Bamboo: AGBD = ",
  round(b["(Intercept)"], 3), " + ",
  round(b["tree_cover"], 3), " × tree_cover"
)

eq_trees <- paste0(
  "Trees: AGBD = ",
  round(b["(Intercept)"] + b["land_classtrees"], 3), " + ",
  round(b["tree_cover"] + b["land_classtrees:tree_cover"], 3),
  " × tree_cover"
)

eq_bamboo
eq_trees


