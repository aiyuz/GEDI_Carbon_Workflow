/************************************************************
 * Seasonal Sentinel-2 composites for Xishuangbanna
 *
 * Purpose:
 *   - Generate four seasonal Sentinel-2 surface reflectance
 *     composites (Jan–Mar, Apr–Jun, Jul–Sep, Oct–Dec) for
 *     the Xishuangbanna study area.
 *   - Apply QA60 and s2cloudless-based cloud masking.
 *   - Export each seasonal composite as a 10 m GeoTIFF to
 *     Google Drive in the folder "MergedSeasons".
 *
 * Requirements:
 *   1. Upload your Xishuangbanna (or other) study-area shapefile
 *      from your local working directory as an Earth Engine asset
 *      (ee.FeatureCollection).
 *   2. Replace the asset ID below (STUDY_AREA_ASSET_ID) with
 *      your own asset path.
 *
 * Outputs (in Google Drive → "MergedSeasons"):
 *   - Sentinel2_Composite_2021_JanMar.tif
 *   - Sentinel2_Composite_2021_AprJun.tif
 *   - Sentinel2_Composite_2021_JulSep.tif
 *   - Sentinel2_Composite_2021_OctDec.tif
 ************************************************************/

// ------------------------
// 1. Import Study Area
// ------------------------

// TODO: Replace this asset ID with your own study-area asset.
// Example after upload: 'users/your_username/Xishuangbanna_Polygon'
var studyArea = ee.FeatureCollection('STUDY_AREA_ASSET_ID');

// Center the map and display the polygon
Map.centerObject(studyArea, 8);
Map.addLayer(studyArea, {color: 'blue'}, 'Study Area');

// ------------------------
// 2. Define Seasons and Bands
// ------------------------

// Year can be changed if desired
var year = 2021;

// Four seasonal windows for this year
var seasons = [
  {name: 'JanMar', start: year + '-01-01', end: year + '-03-31'},
  {name: 'AprJun', start: year + '-04-01', end: year + '-06-30'},
  {name: 'JulSep', start: year + '-07-01', end: year + '-09-30'},
  {name: 'OctDec', start: year + '-10-01', end: year + '-12-31'}
];

// Sentinel-2 bands used in the R workflow (B2–B12 subset)
var bands = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B11', 'B12'];

// ------------------------
// 3. Cloud-masking function
// ------------------------

function cloudMaskS2(image) {
  // 3.1 QA60 bitmask (clouds and cirrus)
  var qa = image.select('QA60');
  var cloudBitMask  = 1 << 10;  // clouds
  var cirrusBitMask = 1 << 11;  // cirrus

  var qaMask = qa.bitwiseAnd(cloudBitMask).eq(0)
                 .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  // 3.2 s2cloudless cloud probability (0–100)
  var cloudProb = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
      .filter(ee.Filter.eq('system:index', image.get('system:index')))
      .first()
      .select('probability');

  // Cloud probability threshold (can be tuned)
  var cloudThreshold = 20;
  var s2cloudlessMask = cloudProb.lt(cloudThreshold);

  // 3.3 Combine QA60 and s2cloudless masks
  var finalMask = qaMask.and(s2cloudlessMask);

  // 3.4 Apply mask, scale reflectance (DN / 10000 → 0–1), select bands
  return image
    .updateMask(finalMask)
    .divide(10000)
    .select(bands)
    .copyProperties(image, ['system:time_start']);
}

// ------------------------
// 4. Seasonal compositing loop
// ------------------------

seasons.forEach(function(season) {
  print('Processing season: ' + season.name);

  // 4.1 Filter Sentinel-2 TOA (harmonized) to find low-cloud scenes
  var s2_toa = ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
      .filterDate(season.start, season.end)
      .filterBounds(studyArea)
      .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 5));

  print('TOA images for ' + season.name + ' (cloudy_pixel_percentage < 5):');
  print(s2_toa.size());

  // 4.2 Match TOA indices to SR_HARMONIZED for consistent imagery
  var indices = s2_toa.aggregate_array('system:index');

  var s2_sr = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
      .filter(ee.Filter.inList('system:index', indices));

  print('SR images for ' + season.name + ' after matching TOA indices:');
  print(s2_sr.size());

  // 4.3 Apply cloud masking
  var s2_masked = s2_sr.map(cloudMaskS2);

  print('Masked images for ' + season.name + ':');
  print(s2_masked.size());

  // 4.4 Median composite over the season, clipped to the study area
  var composite = s2_masked.median().clip(studyArea);

  print('Composite for ' + season.name + ':');
  print(composite);

  // 4.5 Add composite to the map (true color)
  var visParams = {
    bands: ['B4', 'B3', 'B2'],
    min: 0,
    max: 0.3,
    gamma: 1.4
  };
  Map.addLayer(
    composite,
    visParams,
    season.name + ' Composite (' + year + ')'
  );

  // 4.6 Export composite to Google Drive
  // Files will appear under Drive → "MergedSeasons"
  Export.image.toDrive({
    image: composite,
    description: 'Sentinel2_Composite_' + year + '_' + season.name,
    folder: 'MergedSeasons',  // Drive folder name
    fileNamePrefix: 'Sentinel2_Composite_' + year + '_' + season.name,
    region: studyArea.geometry(),
    scale: 10,                // 10 m resolution
    crs: 'EPSG:4326',         // can be changed to a projected CRS if desired
    maxPixels: 1e13
  });

});
