/************************************************************
 * Terrain rasters for Xishuangbanna
 *
 * Purpose:
 *   - Export four terrain-related rasters for Xishuangbanna:
 *       (1) DTM (COPERNICUS GLO30 DEM)
 *       (2) DSM (JAXA ALOS AW3D30)
 *       (3) 10 m canopy height (ETH global canopy height)
 *       (4) 10 m canopy height standard deviation
 *   - All exports are clipped to the study area and written
 *     to Google Drive folder "Terrain" as GeoTIFFs.
 *
 * Requirements:
 *   1. Upload your study-area shapefile as an Earth Engine
 *      asset (ee.FeatureCollection).
 *   2. Replace the asset ID below (STUDY_AREA_ASSET_ID) with
 *      your own asset path.
 *   3. ETH global canopy height must be accessible via the
 *      public GEE assets used below (users/nlang/...).
 *
 * Outputs (in Google Drive → "Terrain"):
 *   - DTM_Xishuangbanna.tif
 *   - DSM_Xishuangbanna.tif
 *   - ETH_GlobalCanopyHeight_2020_10m_Xishuangbanna.tif
 *   - ETH_GlobalCanopyHeightSD_2020_10m_Xishuangbanna.tif
 *
 * After export:
 *   - Download these four GeoTIFFs and place them into a
 *     local folder called "Terrain" next to your R scripts,
 *     so that 01_compile_rasters_xsbn.R can find:
 *       Terrain/DTM_Xishuangbanna.tif
 *       Terrain/DSM_Xishuangbanna.tif
 *       Terrain/ETH_GlobalCanopyHeight_2020_10m_Xishuangbanna.tif
 *       Terrain/ETH_GlobalCanopyHeightSD_2020_10m_Xishuangbanna.tif
 ************************************************************/

// ------------------------
// 1. Import Study Area
// ------------------------

// TODO: replace with your own asset if different
// e.g. 'projects/ee-aiyuzheng/assets/Xishuangbanna_Polygon'
var studyArea = ee.FeatureCollection('STUDY_AREA_ASSET_ID');

Map.centerObject(studyArea, 8);
Map.addLayer(studyArea, {color: 'blue'}, 'Study Area');

// ------------------------
// 2. Load DTM (GLO30) and DSM (AW3D30)
// ------------------------

// Digital Terrain Model (ground elevation)
var dtm = ee.ImageCollection('COPERNICUS/DEM/GLO30')
  .filterBounds(studyArea)
  .select('DEM')
  .mosaic()
  .clip(studyArea)
  .rename('DTM');

// Digital Surface Model (top of canopy/buildings)
var dsm = ee.ImageCollection('JAXA/ALOS/AW3D30/V4_1')
  .filterBounds(studyArea)
  .select('DSM')
  .mosaic()
  .clip(studyArea)
  .rename('DSM');

// Optional quick-look visualization
var dtmVis = {min: 0, max: 2500, palette: ['blue', 'green', 'brown']};
var dsmVis = {min: 0, max: 2500, palette: ['white', 'gray', 'black']};

Map.addLayer(dtm, dtmVis, 'DTM (GLO30)');
Map.addLayer(dsm, dsmVis, 'DSM (AW3D30)');

// ------------------------
// 3. Load ETH global canopy height (mean + SD)
//    See: https://langnico.github.io/globalcanopyheight/
// ------------------------

// These are public GEE image assets provided by the dataset authors.
// If the IDs ever change, update them here according to the dataset docs.
var canopyMean = ee.Image('users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1')
  .clip(studyArea);

var canopySD = ee.Image('users/nlang/ETH_GlobalCanopyHeightSD_2020_10m_v1')
  .clip(studyArea);

// Simple visualization for sanity-check
var chVis = {min: 0, max: 40, palette: ['#f7fbff', '#6baed6', '#08306b']};
Map.addLayer(canopyMean, chVis, 'Canopy height (10 m)');

// ------------------------
// 4. Export to Google Drive ("Terrain" folder)
// ------------------------

// Note:
//  - CRS choices here are pragmatic; in your R workflow you
//    reproject all terrain rasters to match the Sentinel-2
//    composites (EPSG:4326).
//  - You can change CRS if you prefer, as long as you handle
//    it consistently in R.

// 4.1 DTM (GLO30, ~30 m)
Export.image.toDrive({
  image: dtm,
  description: 'DTM_Xishuangbanna',
  folder: 'Terrain',
  fileNamePrefix: 'DTM_Xishuangbanna',
  region: studyArea.geometry(),
  scale: 30,
  crs: 'EPSG:32647',   // UTM 47N (same as your original script)
  maxPixels: 1e13
});

// 4.2 DSM (AW3D30, ~30 m)
Export.image.toDrive({
  image: dsm,
  description: 'DSM_Xishuangbanna',
  folder: 'Terrain',
  fileNamePrefix: 'DSM_Xishuangbanna',
  region: studyArea.geometry(),
  scale: 30,
  crs: 'EPSG:32647',
  maxPixels: 1e13
});

// 4.3 Canopy height mean (10 m)
Export.image.toDrive({
  image: canopyMean,
  description: 'ETH_GlobalCanopyHeight_2020_10m_Xishuangbanna',
  folder: 'Terrain',
  fileNamePrefix: 'ETH_GlobalCanopyHeight_2020_10m_Xishuangbanna',
  region: studyArea.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// 4.4 Canopy height standard deviation (10 m)
Export.image.toDrive({
  image: canopySD,
  description: 'ETH_GlobalCanopyHeightSD_2020_10m_Xishuangbanna',
  folder: 'Terrain',
  fileNamePrefix: 'ETH_GlobalCanopyHeightSD_2020_10m_Xishuangbanna',
  region: studyArea.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});
