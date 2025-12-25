# GEDI Carbon Workflow (GEDI + Sentinel-2) — Xishuangbanna Bamboo Carbon Bias

This repository contains scripts supporting a workflow that integrates GEDI Level 2A, 2B, and 4A products and Sentinel-2 (spectral and texture features) to classify bamboo-dominated canopies and quantify carbon bias associated with grassy-tree misclassification in Xishuangbanna.

## Repository structure

- `GEE/`  
  Google Earth Engine scripts for Sentinel-2 preprocessing and raster compositing/feature engineering.

- `R_Scripts/`  
  R scripts for compiling rasters, classification post-processing, GEDI overlay, statistical tests, Monte Carlo uncertainty propagation, and bias mapping.

## Data policy

This repository distributes **code and documentation only**. Large raster exports and intermediate outputs are not tracked.

Recommended local folders:
```bash
mkdir -p data/raw data/intermediate data/derived outputs figures 


