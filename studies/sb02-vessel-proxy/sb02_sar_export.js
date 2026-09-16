/*******************************************************************************
 * sb02_sar_export.js
 * Sentinel-1 per-scene vessel-detection table for SanctSound SB02
 * (Stellwagen Bank NMS, Massachusetts Bay; 42.470793 N, -70.24294 W)
 *
 * -----------------------------------------------------------------------------
 * WHAT THIS SCRIPT DOES
 * -----------------------------------------------------------------------------
 * It builds ONE ROW PER SENTINEL-1 ACQUISITION over the SB02 hydrophone and
 * exports it to Drive as a CSV. Each row carries:
 *   - acquisition metadata (id, UTC timestamp, orbit pass, relative orbit,
 *     platform),
 *   - vessel-detection statistics (target count, total bright-target area,
 *     largest single target) computed independently inside each of five
 *     concentric search discs: 2, 5, 10, 20 and 40 km,
 *   - a per-radius footprint-coverage fraction and a full-coverage flag, so
 *     scenes that clip the edge of a disc (and therefore UNDERCOUNT) can be
 *     identified and excluded offline,
 *   - ERA5 10 m wind (u, v, speed, direction) at the acquisition hour,
 *   - the detection parameter values used, echoed into every row, so a
 *     parameter sweep produces self-describing CSVs.
 *
 * -----------------------------------------------------------------------------
 * WHAT THIS SCRIPT DOES *NOT* DO
 * -----------------------------------------------------------------------------
 * It computes NO statistics, NO correlations, NO regressions, NO skill scores,
 * and makes NO claim whatsoever about whether SAR vessel counts are a good
 * proxy for acoustic ship detections. It does not touch, read, or join the
 * 15,255 acoustic ship-detection intervals. This script is a DATA EXTRACTION
 * step and nothing else. Every inferential step happens OFFLINE, downstream,
 * against this CSV. Nothing printed by this script should be read as a result.
 *
 * -----------------------------------------------------------------------------
 * ES5 ONLY - THIS IS A HARD CONSTRAINT, NOT A STYLE PREFERENCE
 * -----------------------------------------------------------------------------
 * The Earth Engine Code Editor sandbox parses ES5 only. It does NOT accept
 * const, let, arrow functions, template literals, Object.assign, spread,
 * for...of, class, Promise, Map/Set, Array.prototype.includes/find, or
 * String.prototype.repeat/startsWith/endsWith. Code using any of those fails at
 * PARSE time with an opaque message and nothing runs. This has broken this
 * project twice in production. Everything below is `var` and `function()`.
 * If you edit this file, keep it that way.
 *
 * -----------------------------------------------------------------------------
 * DETECTION PARAMETERS - SUPERSEDED NOTICE (v2)
 * -----------------------------------------------------------------------------
 * This section originally said every value in PARAMS was a placeholder pending
 * a parameter specification. THAT IS NO LONGER TRUE and the sentence is
 * corrected rather than deleted, because a file carrying two opposite claims
 * about its own parameters is worse than either claim alone.
 * The specification landed (sar_params.md) and its values were applied in v2:
 * VH-only, no speckle pre-filter, annulus-median background de-biased by m(L),
 * ENL-derived threshold at k=8 (+6.82 dB), 150/400 m CFAR geometry, 3 px
 * minimum target, 10 m analysis scale, WorldCover land mask at 1000 m buffer.
 * The function SIGNATURE and RETURN CONTRACT remain fixed; the body is still
 * swappable. See the RECONCILED v2 annotations inline for each prior value.
 *
 * -----------------------------------------------------------------------------
 * COMPUTE COST - READ BEFORE PRESSING RUN
 * -----------------------------------------------------------------------------
 * The 40 km disc is ~5.03e9 m2. At the v2 analysis scale of 10 m that is
 * ~50.3 million pixels PER SCENE, times 278 scenes - 4x the original 20 m
 * draft. (The v1 text here said 20 m and ~12.6 M px; that was superseded when
 * the spec moved the scale to native 10 m so the 3 px / 300 m2 minimum target
 * means what it says.) This is a genuinely large job. Interactive prints are therefore deliberately cheap:
 * the previewed feature is built over CONFIG.PREVIEW_RADII_M only (2 km by
 * default). The Export task is NOT subject to the interactive timeout, but it
 * may still take hours. If the export fails with "computation timed out",
 * split it by year via CONFIG.DATE_START / CONFIG.DATE_END and concatenate the
 * CSVs offline, or coarsen PARAMS.analysisScaleM.
 ******************************************************************************/


/* =============================================================================
 * RECONCILED v2 - the detection parameters in this file are NO LONGER the
 * first draft's placeholders. They were replaced with the values from the SAR
 * parameter specification (sar_params.md) and each swap is annotated inline
 * with its previous value. Changed: polarisation (VV+VH -> VH only), speckle
 * pre-filter (on -> off), background estimator (sample mean+k*sd -> annulus
 * median de-biased by m(L)), threshold (k=5 on sample sd -> k=8 on ENL-derived
 * sigma, +6.82 dB), CFAR geometry (90/300 m -> 150/400 m), minimum target
 * (2 px -> 3 px), analysis scale (20 m -> 10 m native), land mask (Hansen 30 m
 * -> ESA WorldCover 10 m, buffer 300 m -> 1000 m).
 * VERIFIED OUTSIDE EARTH ENGINE: node --check, an ES5-subset scan with strings
 * and comments stripped, the annulus kernel weight count (4,316), the threshold
 * multiplier (4.8139 = +6.82 dB) and the 3 px = 300 m2 minimum target.
 * NOT VERIFIED, and unverifiable in this environment: anything requiring Earth
 * Engine to execute. A clean syntax and ES5 pass means the script is
 * well-formed, NOT that the detector works.
 *
 * 1. CONFIG - every tunable lives here, nothing is buried in the body
 * ========================================================================== */

var CONFIG = {

  // ---- site ---------------------------------------------------------------
  SITE_NAME: 'SanctSound_SB02',
  SITE_LAT:  42.470793,
  SITE_LON:  -70.24294,

  // ---- temporal window ----------------------------------------------------
  // Verified: exactly 278 IW scenes, all dual-pol [VV, VH], relative orbits
  // 142 / 40 / 62, ASCENDING 93 + DESCENDING 185.
  DATE_START: '2018-11-12',
  DATE_END:   '2021-12-16',
  // ee.Filter.date() treats the end as EXCLUSIVE. The verified count of 278 was
  // produced with the half-open interval, so this stays false. If your sanity
  // print comes back 279 (or whatever count), flip this to true and re-check -
  // do not "fix" it by editing the expected number.
  DATE_END_INCLUSIVE: false,

  // ---- search radii, METRES ----------------------------------------------
  RADII_M: [2000, 5000, 10000, 20000, 40000],

  // ---- collection ---------------------------------------------------------
  S1_COLLECTION: 'COPERNICUS/S1_GRD',   // NOTE: this collection is in dB (log
                                        // scale). COPERNICUS/S1_GRD_FLOAT is
                                        // the linear-power variant.
  INSTRUMENT_MODE: 'IW',
  REQUIRE_DUAL_POL: true,               // adds listContains VV *and* VH filters.
                                        // Verified not to change the count of
                                        // 278 (all scenes are dual-pol), but it
                                        // makes the assumption explicit and
                                        // fail-loud if the archive changes.

  // Scene selection geometry. 'point' reproduces the verified 278. Using the
  // 40 km disc instead would ADD scenes that touch the disc edge but miss the
  // hydrophone, raising the count above 278 - correct for some purposes, but
  // it is not the verified set, so it is not the default.
  FILTER_BOUNDS: 'point',               // 'point' | 'maxdisc'

  // ---- land / harbour exclusion ------------------------------------------
  // RECONCILED v2: ESA WorldCover v200, 10 m, band 'Map'; class 80 = permanent
  // water. Replaces Hansen datamask (30 m) - the parameter spec calls for a
  // 10 m marine-capable source, and Hansen's grid is coarser than the SAR grid
  // it is masking. See the LAND MASKING section for the no-data handling, which
  // keeps the original script's defensive "water = NOT land" property.
  // ---- v3: persistent-target mask -----------------------------------------
  // Modes:
  //   'off'     - no persistence masking (v2 behaviour)
  //   'compute' - build the masks inline. CORRECT but evaluates the detector
  //               over every scene of every relative orbit IN ADDITION to the
  //               per-scene work. Expect this to be too slow for one export.
  //   'asset'   - load masks previously written by the STAGE A block below.
  //               This is the intended production path.
  PERSISTENCE_MODE: 'compute',
  EXPORT_PERSISTENCE_ASSETS: true,   // STAGE A; see section 13
  PERSISTENCE_ASSET_PREFIX: 'users/CHANGE_ME/sb02_persist_orbit_',
  RELATIVE_ORBITS: [142, 40, 62],       // verified: 93 / 92 / 93 scenes
  // ---- v3: per-detection long table ---------------------------------------
  // The spec's primary output. One row per DETECTION, not per scene, carrying
  // its own range from the hydrophone - so all five radii are derived offline
  // by filtering one table, and radius differences are pure geometry rather
  // than five independent detector runs. Also the only form that supports
  // azimuth-ambiguity flagging, which needs pairwise geometry between
  // detections and cannot be done from per-scene counts.
  EXPORT_DETECTIONS: true,
  DETECTION_MAX_PER_SCENE: 5000,        // guard against a pathological scene
  LAND_SOURCE: 'ESA/WorldCover/v200',
  LAND_BAND:   'Map',                   // class 80 = permanent water bodies
  LAND_WATER_CLASS: 80,
  LAND_BUFFER_M: 1000,                  // spec: 1000 m, was 300 m
  MASK_SCALE_M: 30,                     // scale at which water area is tallied

  // Extra hard-exclusion polygons for harbour structures, piers, breakwaters
  // and moored/rafted vessels, which an adaptive threshold will happily call
  // "targets". Each entry is [west, south, east, north] in degrees.
  //
  // GEOMETRY CHECK, because the brief's framing is slightly off and it matters:
  // from (42.470793, -70.24294), 40 km reaches Gloucester Harbor / Annisquam
  // (~37 km) and the whole Cape Ann shore, which IS inside the disc. Boston
  // Inner Harbor is ~67 km away, Nahant/Lynn ~56 km, Marblehead ~51 km,
  // Scituate ~49 km - ALL OUTSIDE the 40 km disc. So Cape Ann is the harbour
  // problem here; Boston is not, at these radii. The Boston boxes are listed
  // and left enabled anyway: they cost nothing (they simply never intersect)
  // and they keep the list correct if anyone widens RADII_M later.
  EXCLUDE_BOXES: [
    [-70.700, 42.575, -70.615, 42.680],   // Gloucester Hbr + Annisquam River
    [-70.640, 42.645, -70.590, 42.685],   // Rockport / Pigeon Cove
    [-70.920, 42.495, -70.820, 42.570],   // Salem / Beverly / Marblehead
    [-71.070, 42.310, -70.940, 42.420]    // Boston Inner Harbor + Chelsea
  ],
  APPLY_EXCLUDE_BOXES: true,

  // ---- ERA5 wind ----------------------------------------------------------
  // MUST be ECMWF/ERA5/HOURLY, NOT ECMWF/ERA5_LAND/HOURLY. ERA5-LAND is masked
  // over ocean, so at an offshore hydrophone every single wind value would come
  // back null. This is the exact trap the brief hints at.
  ERA5_COLLECTION: 'ECMWF/ERA5/HOURLY',
  ERA5_U_BAND: 'u_component_of_wind_10m',
  ERA5_V_BAND: 'v_component_of_wind_10m',
  ERA5_MAX_DIFF_MILLIS: 3600000,  // +/- 1 h; saveBest keeps the NEAREST hour
  ERA5_SCALE_M: 27830,            // ERA5 native ~0.25 deg
  // Joining an hourly asset per scene is not awkward here: ee.Join.saveBest
  // with an ee.Filter.maxDifference on system:time_start does it in one
  // server-side pass, no getInfo, no client loop. ERA5 is hourly-instantaneous
  // at the top of the hour while S1 acquires mid-hour, so what is emitted IS
  // the NEAREST-HOUR value, never an interpolation. era5_dt_min records the
  // signed-magnitude offset in minutes so you can audit that offline.
  // If ECMWF/ERA5/HOURLY's ingested range ever stops short of 2021-12, rows
  // will carry era5_matched = 0; swap ERA5_COLLECTION rather than patching
  // around it.

  // ---- output -------------------------------------------------------------
  PARAM_SET_ID: 'p001',                 // bump this for every parameter sweep
  DRIVE_FOLDER: 'GEE_SB02_SAR',
  EXPORT_DESCRIPTION: 'SB02_S1_vessel_detections',
  FILE_PREFIX: 'sb02_s1_vessel_detections',

  // ---- run control --------------------------------------------------------
  INCLUDE_WATER_AREA: true,   // per-radius unmasked water area (for density
                              // normalisation offline). Computed ONCE, not per
                              // scene, then attached as a constant column.
  INCLUDE_PARAM_COLUMNS: true,// echo PARAMS into every row (sweep provenance)
  PRINT_FIRST_FEATURE: true,  // cheap: uses PREVIEW_RADII_M only
  PREVIEW_RADII_M: [2000],    // radii used for the interactive preview row
  SHOW_MAP: true,
  TEST_LIMIT: 0               // >0 exports only the first N scenes (smoke test).
                              // Sanity prints always use the FULL collection.
};


/* -----------------------------------------------------------------------------
 * PARAMS - detection tunables. v2: these are the SPECIFICATION values, not
 * placeholders (see the superseded notice in the header).
 * Sweep by editing this object and bumping CONFIG.PARAM_SET_ID; the values are
 * echoed into the CSV so a row always knows how it was made.
 * -------------------------------------------------------------------------- */

var PARAMS = {
  // ==========================================================================
  // RECONCILED v2 - these values now come from the SAR parameter specification
  // (sar_params.md), not from the first draft's placeholders. Every change is
  // annotated with the previous value so the swap is auditable.
  // ==========================================================================

  // VH ONLY. Spec: VH has lower Bragg sea clutter than VV, so target-to-clutter
  // is higher. WAS ['VV','VH'] combined with OR, which roughly doubles the
  // false-alarm rate for little sensitivity gain and has no principled
  // threshold. VV is still READ, but only as a wind/diagnostic covariate.
  detectionPols: ['VH'],
  polCombine: 'first',
  useLinearPower: true,          // CFAR maths is only valid in linear power

  // NO SPECKLE PRE-FILTER. WAS true/30 m. Spec's argument, and it is correct:
  // the CFAR annulus IS the speckle averaging, and pre-smoothing (a) cuts
  // small-target peak amplitude by 5-10 dB, which removes exactly the vessels
  // at the edge of detectability, and (b) breaks the sd = mean/sqrt(L) relation
  // that the ENL-based threshold below depends on. Kept switchable for a
  // sensitivity run only.
  speckleApply: false,
  speckleWindowM: 30,

  // CFAR geometry, spec values. WAS 90 m / 300 m.
  // Guard 150 m: must exclude the target's own energy plus azimuth smearing;
  // a ~300 m LOA ship centred in the cell reaches +/-150 m.
  // Outer 400 m: enough independent looks that the background mean is precise
  // relative to the threshold offset.
  useGuardRing: true,
  bgGuardRadiusM: 150,
  bgOuterRadiusM: 400,

  // THRESHOLD. Spec form: T = mu_bg * (1 + k/sqrt(L)), NOT mean + k*sd.
  // WAS thresholdK 5.0 against a sample sd. The spec derives sigma from the
  // known ENL instead of estimating it per-pixel, because a sample sd computed
  // over an annulus that contains a second vessel is inflated, which raises the
  // threshold and hides both - the convoy self-masking case.
  // k = 8.0, L = 4.4 -> T = 4.81 * mu_bg = +6.8 dB above local background.
  thresholdK: 8.0,
  enl: 4.4,                      // ESA nominal for IW GRDH - MEASURE IT (spec 2.4)
  // Median-to-mean de-bias factor for a unit-mean Gamma(L). m(4.4) = 0.9253,
  // so mu_bg = median / 0.9253. Recompute if enl is re-measured.
  medianToMean: 0.9253,
  // Robust background: annulus MEDIAN (spec) vs annulus MEAN (cheaper).
  // Median needs an explicit annulus kernel and is materially more expensive;
  // the guard ring already excludes the target's own energy, so the mean path
  // is a defensible fallback if the export will not complete. If you switch
  // this to false, SAY SO when reporting results - it changes the estimator.
  useMedianBackground: true,
  minBackgroundPixels: 30,

  minBackscatterDb: { VV: -15.0, VH: -22.0 },   // VH floor at the IW NESZ

  // TARGET SIZE. Spec: 3 connected pixels at 10 m = 300 m2, ~20-30 m LOA in
  // calm seas. WAS 2 px at a 20 m analysis scale = 800 m2, which is a ~28 m
  // cell and would have quietly raised the minimum detectable vessel.
  // These two MUST move together - see analysisScaleM below.
  minTargetPixels: 3,
  eightConnected: true,
  maxComponentSizePx: 1024,
  edgeErodeM: 250,

  // ANALYSIS SCALE. WAS 20 m for tractability. Moved to the native 10 m so the
  // 3 px / 300 m2 minimum means what the spec says it means. THIS IS THE MAIN
  // COMPUTE RISK: it is 4x the pixels of the 20 m draft, on an export already
  // flagged as possibly needing to be split by year.
  // DOCUMENTED FALLBACK, if the export will not complete: set analysisScaleM
  // to 20 AND minTargetPixels to 2, and report the minimum detectable vessel
  // as ~800 m2 rather than ~300 m2. Do not change one without the other.
  // ---- v3: persistence ----------------------------------------------------
  // A pixel detected in >= this fraction of a relative orbit's scenes is a
  // FIXED object, not a vessel: navigation buoys, wrecks, the SanctSound
  // mooring itself, unmasked rock. Data-driven, so no prior list is needed.
  //
  // PER RELATIVE ORBIT, NOT POOLED. Viewing geometry differs between orbits,
  // so a fixed target's apparent position and detectability differ too;
  // pooling would smear each object across three positions.
  //
  // THE RISK THIS CARRIES, stated plainly: SB02 sits IN a shipping lane. If
  // traffic is channelled tightly enough that one 10 m pixel is occupied in
  // >=20% of scenes, this mask would delete real signal at exactly the site
  // the study is about. That is why CHECK 10 prints the persistence histogram:
  // adopt 0.20 only if the histogram is BIMODAL (a fixed-object mode near 1.0,
  // clearly separated from a traffic mode near 0). If it is unimodal, raise
  // the threshold or disable persistence and mask known objects by hand.
  persistenceThreshold: 0.20,
  persistenceDilateM: 30,
  // ---- v3: azimuth ambiguity (parameters are ECHOED, not applied here) -----
  // Flagging happens OFFLINE on the detection table - it needs pairwise range
  // and bearing between detections, which is trivial in a table and painful in
  // image space. These live here so the CSV is self-describing.
  ambigMinKm: 3.5,
  ambigMaxKm: 7.5,
  ambigParentDb: 15.0,
  analysisScaleM: 10,
  tileScale: 4,
  maxPixels: 1e10,
  geomMaxErrorM: 10,
  fullCoverageMinFrac: 0.999
};


/* =============================================================================
 * 2. DERIVED CONSTANTS (client-side; no server calls, no getInfo)
 * ========================================================================== */

var SITE_POINT = ee.Geometry.Point([CONFIG.SITE_LON, CONFIG.SITE_LAT]);

// Radius tags: 2000 -> '2km'. Used to build column names.
var RADIUS_TAGS = [];
var i;
for (i = 0; i < CONFIG.RADII_M.length; i++) {
  RADIUS_TAGS.push((CONFIG.RADII_M[i] / 1000) + 'km');
}
var MAX_RADIUS_M = CONFIG.RADII_M[CONFIG.RADII_M.length - 1];

// Harbour exclusion geometry (client-side union of the CONFIG boxes).
var EXCLUDE_GEOM = null;
if (CONFIG.APPLY_EXCLUDE_BOXES && CONFIG.EXCLUDE_BOXES.length > 0) {
  var excludeFeats = [];
  for (i = 0; i < CONFIG.EXCLUDE_BOXES.length; i++) {
    // Default (geodesic) rectangles, to match the geodesic buffered discs they
    // are subtracted from - mixing geodesic and planar geometries in
    // Geometry.difference() is a reliable way to get an obscure error.
    excludeFeats.push(ee.Feature(ee.Geometry.Rectangle(CONFIG.EXCLUDE_BOXES[i])));
  }
  EXCLUDE_GEOM = ee.FeatureCollection(excludeFeats).geometry(PARAMS.geomMaxErrorM);
}

// One disc per radius, with harbour boxes cut out. These are the AOIs handed to
// detectVessels(). The raster water mask is applied separately, inside the
// detection image - geometry alone cannot follow a coastline.
var DISCS = [];      // pure discs, used for the coverage fraction denominator
var AOIS = [];       // discs minus harbour boxes, used for the reductions
for (i = 0; i < CONFIG.RADII_M.length; i++) {
  var disc = SITE_POINT.buffer(CONFIG.RADII_M[i]);
  DISCS.push(disc);
  if (EXCLUDE_GEOM === null) {
    AOIS.push(disc);
  } else {
    AOIS.push(disc.difference(EXCLUDE_GEOM, PARAMS.geomMaxErrorM));
  }
}


/* =============================================================================
 * 3. LAND MASKING
 * ========================================================================== */
/*
 * RECONCILED v2. Source: ESA WorldCover v200, band 'Map', 10 m; class 80 is
 * permanent water. The parameter spec selected this over Hansen datamask (30 m)
 * because the mask should not be coarser than the SAR grid it masks.
 *
 * THE NO-DATA TRAP, AND WHY THE NEGATION IS KEPT. The obvious encoding is
 * land = Map.neq(80). That is WRONG here in a way that fails silently: ESA
 * WorldCover is a LAND cover product, and its offshore footprint is not
 * guaranteed to extend across open ocean. Wherever WorldCover has no data,
 * neq(80) evaluates TRUE - so every unmapped ocean pixel becomes "land", the
 * whole search area is masked away, and the export returns zero vessels at
 * every scene while looking perfectly healthy.
 *
 * So land is asserted only where WorldCover ACTUALLY HAS DATA and that data is
 * not water; unmapped pixels fall through to water. Water is then still defined
 * as NOT(dilated land), preserving the original script's defensive property:
 * the worst case is that an unmapped patch gets searched, never that the ocean
 * silently disappears.
 *
 * Agent-flagged as unverified: WorldCover's exact ocean-footprint edge
 * behaviour. This construction makes that uncertainty harmless rather than
 * fatal, but CHECK 7 (water area per radius) must still be read before the
 * export - a water area far below the geometric disc area means the mask is
 * eating ocean.
 *
 * The land class is then dilated seaward by LAND_BUFFER_M, which removes the
 * coastal pixels where a 30 m coastline and a 10-20 m SAR grid disagree, plus
 * shore-attached piers and breakwaters. Marinas and moored fleets that sit
 * clear of the dilation are handled by CONFIG.EXCLUDE_BOXES above.
 *
 * Alternative sources if Hansen is ever unavailable or too coarse:
 *   ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017')  - vector, simplified
 *   ee.Image('MODIS/006/MOD44W/2015_01_01')        - 250 m, too coarse here
 *   ee.Image('NOAA/NGDC/ETOPO1').select('bedrock').lt(0) - bathymetric, crude
 */

var landRaw = ee.Image(CONFIG.LAND_SOURCE).select(CONFIG.LAND_BAND);
// hasData=1 where WorldCover mapped the pixel at all; 0 where it did not.
var landHasData = landRaw.mask().gt(0);
// land = mapped AND not the permanent-water class. Unmapped -> 0 -> water.
var landBinary = landHasData.and(landRaw.neq(CONFIG.LAND_WATER_CLASS)).unmask(0);
var landDilated = landBinary.focal_max({
  radius: CONFIG.LAND_BUFFER_M,
  kernelType: 'circle',
  units: 'meters'
});
var WATER_MASK = landDilated.not().rename('water');   // 1 = searchable water

// Per-radius searchable water area. Computed ONCE here, not inside the per-scene
// map, so it is a single node in the graph rather than 278 repeats. Used offline
// to convert counts into densities (vessels per km2 of searchable water).
var WATER_AREA_M2 = [];
if (CONFIG.INCLUDE_WATER_AREA) {
  for (i = 0; i < CONFIG.RADII_M.length; i++) {
    var wa = ee.Image.pixelArea().updateMask(WATER_MASK).reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: AOIS[i],
      scale: CONFIG.MASK_SCALE_M,
      maxPixels: PARAMS.maxPixels,
      tileScale: PARAMS.tileScale
    }).get('area');
    WATER_AREA_M2.push(ee.Number(ee.Algorithms.If(wa, wa, 0)));
  }
}


/* =============================================================================
 * 4. PLUGGABLE DETECTION FUNCTION
 * =============================================================================
 *
 *  +-----------------------------------------------------------------------+
 *  |  *** v2: SPECIFICATION VALUES APPLIED - CONTRACT STILL FIXED ***      |
 *  |                                                                       |
 *  |  Everything between this banner and the "END PLUGGABLE BLOCK" banner  |
 *  |  is a reasonable default standing in for the detection specification  |
 *  |  being written in a parallel workstream. Every constant it uses lives |
 *  |  in PARAMS above. The SIGNATURE and the RETURN CONTRACT below are     |
 *  |  fixed and are relied on by the per-scene mapper; the body is not.    |
 *  |  Swap the body, keep the contract, and nothing downstream changes.    |
 *  +-----------------------------------------------------------------------+
 *
 * detectVessels(img, aoi, params)
 *   img    - one ee.Image from COPERNICUS/S1_GRD (dual-pol, dB)
 *   aoi    - ee.Geometry to search within (already land-masked upstream)
 *   params - the PARAMS object from CONFIG
 * returns an ee.Dictionary with AT MINIMUM:
 *   {count: <number of detected vessel targets>,
 *    area_m2: <total bright-target area>,
 *    max_target_m2: <largest single target>}
 *
 * NOTE ON `aoi` AND LAND: the AOIs handed in are already cut against the
 * harbour exclusion polygons, but a polygon cannot follow a coastline. The
 * raster WATER_MASK is therefore ALSO applied inside the detection image below.
 * Both are applied; neither alone is sufficient.
 *
 * Default method:
 *   1. select polarisation, optionally dB -> linear power
 *   2. boxcar speckle filter
 *   3. CFAR-style adaptive threshold: pixel > localMean + k * localStdDev,
 *      where localMean/localStdDev are taken over an ANNULUS (guard ring), so a
 *      bright target cannot raise its own detection threshold
 *   4. absolute dB floor per polarisation
 *   5. OR/AND across polarisations, then AND with the water mask and the
 *      eroded valid-data mask
 *   6. connected-component labelling; components smaller than minTargetPixels
 *      are discarded (this is what rejects single-pixel speckle)
 *   7. count = countDistinctNonNull(labels); area = sum(pixelArea);
 *      max_target_m2 = max(componentPixelCount * pixelArea)
 */

// --- internal helper: annulus (guard-ring) background statistics -------------
// Returns a 2-band image ['mean','sd'] of the background statistics of a
// single-band image, estimated over the annulus between guardR and outerR.
// Implemented from ring sums/counts because ee.Kernel has no arithmetic:
//   annulus = outer circle MINUS inner circle, for sum, count and sum-of-squares
function ringStats(band1, params) {
  var kOut = ee.Kernel.circle({
    radius: params.bgOuterRadiusM, units: 'meters', normalize: false
  });
  var kIn = ee.Kernel.circle({
    radius: params.bgGuardRadiusM, units: 'meters', normalize: false
  });

  var x = band1.rename('b');
  var x2 = x.multiply(x).rename('b');

  var sumOut = x.reduceNeighborhood({
    reducer: ee.Reducer.sum(), kernel: kOut, skipMasked: true}).rename('b');
  var cntOut = x.reduceNeighborhood({
    reducer: ee.Reducer.count(), kernel: kOut, skipMasked: true}).rename('b');
  var sqOut = x2.reduceNeighborhood({
    reducer: ee.Reducer.sum(), kernel: kOut, skipMasked: true}).rename('b');

  var sumAnn, cntAnn, sqAnn;

  if (params.useGuardRing) {
    var sumIn = x.reduceNeighborhood({
      reducer: ee.Reducer.sum(), kernel: kIn, skipMasked: true}).rename('b');
    var cntIn = x.reduceNeighborhood({
      reducer: ee.Reducer.count(), kernel: kIn, skipMasked: true}).rename('b');
    var sqIn = x2.reduceNeighborhood({
      reducer: ee.Reducer.sum(), kernel: kIn, skipMasked: true}).rename('b');
    sumAnn = sumOut.subtract(sumIn);
    cntAnn = cntOut.subtract(cntIn);
    sqAnn = sqOut.subtract(sqIn);
  } else {
    // Cheaper path: plain disc neighbourhood, no guard cells. Bright targets
    // then contribute to their own background estimate and the effective
    // threshold rises - this is the classic CFAR self-masking problem.
    sumAnn = sumOut;
    cntAnn = cntOut;
    sqAnn = sqOut;
  }

  // Degenerate neighbourhoods (scene edge, a narrow lead between masked land)
  // are MASKED OUT, not clamped. Clamping them would leave a mean of ~0 and an
  // sd of ~0, i.e. a threshold of ~0, and in linear power every pixel clears
  // that - a ring of phantom "vessels" around every coastline and scene edge.
  var enoughBg = cntAnn.gte(params.minBackgroundPixels);
  cntAnn = cntAnn.max(1);   // keeps the division finite; enoughBg does the work

  // RECONCILED v2: this path now returns only the background MEAN estimate.
  // The sample sd is no longer used - see annulusMu() and buildDetectionImage()
  // for why sigma is derived from the ENL instead.
  var mean = sumAnn.divide(cntAnn).updateMask(enoughBg).rename('mu');
  return mean;
}

/* ---------------------------------------------------------------------------
 * RECONCILED v2: annulus kernel + robust background mean.
 *
 * WHY NOT JUST REUSE THE SUMS ABOVE. The guard-ring trick (outer sums minus
 * inner sums) works for sums because sums are additive. Medians are not: you
 * cannot subtract the median of the inner disc from the median of the outer
 * disc. A median over an annulus therefore needs the annulus to exist as an
 * explicit kernel, which is what annulusKernel() builds.
 *
 * WHY A MEDIAN AT ALL. The spec's argument, and it is the right one for this
 * site: SB02 sits in a shipping lane, so a second vessel inside the background
 * annulus is common, not exotic. A sample mean (or worse, a sample sd) computed
 * over an annulus containing another ship is inflated, the threshold rises, and
 * BOTH ships are hidden. That is classic CFAR self-masking, and a guard ring
 * does not fix it - the guard only excludes the target's OWN energy. A median
 * is insensitive to a minority of bright pixels in the annulus.
 *
 * THE DE-BIAS. The median of a unit-mean Gamma(L) is not 1, so a raw median
 * underestimates the mean. mu = median / m(L), with m(4.4) = 0.9253 from the
 * spec. If enl is re-measured, medianToMean must be recomputed with it.
 *
 * COST. This is the expensive part of the whole script: an 81x81 kernel at
 * 10 m evaluated per pixel. MEASURED, by replaying annulusKernel() outside
 * Earth Engine: 4,316 non-zero weights (the spec text said ~5,600, which is
 * 1.30x too high). The conclusion is unchanged - the fractional standard error
 * of the background mean is 1/sqrt(4316 * 4.4) = 0.73%, i.e. 0.031 dB, still
 * negligible against a 6.82 dB threshold offset. PARAMS.useMedianBackground
 * = false falls back to the cheap additive mean path above.
 * ------------------------------------------------------------------------- */
function annulusKernel(outerM, innerM, scaleM) {
  var outerPx = Math.round(outerM / scaleM);
  var innerPx = Math.round(innerM / scaleM);
  var size = 2 * outerPx + 1;
  var w = [];
  var r, c, dy, dx, d;
  for (r = 0; r < size; r++) {
    var row = [];
    for (c = 0; c < size; c++) {
      dy = r - outerPx;
      dx = c - outerPx;
      d = Math.sqrt(dx * dx + dy * dy);
      row.push((d > innerPx && d <= outerPx) ? 1 : 0);
    }
    w.push(row);
  }
  return ee.Kernel.fixed(size, size, w, outerPx, outerPx, false);
}

function annulusMu(work, params) {
  if (!params.useMedianBackground) {
    return ringStats(work, params).select('mu');
  }
  var kAnn = annulusKernel(params.bgOuterRadiusM,
                           params.bgGuardRadiusM,
                           params.analysisScaleM);
  var med = work.reduceNeighborhood({
    reducer: ee.Reducer.median(), kernel: kAnn, skipMasked: true
  }).rename('b');
  // Independently require enough background samples, exactly as the mean path
  // does, so a pixel in a narrow lead between masked land is dropped rather
  // than thresholded against three pixels of sea.
  var cnt = work.reduceNeighborhood({
    reducer: ee.Reducer.count(), kernel: kAnn, skipMasked: true
  }).rename('b');
  var enoughBg = cnt.gte(params.minBackgroundPixels);
  return med.divide(params.medianToMean).updateMask(enoughBg).rename('mu');
}

// --- internal helper: the detection image ------------------------------------
// Split out from detectVessels() ON PURPOSE: it depends only on (img, params),
// so calling detectVessels() five times (once per radius) produces five
// IDENTICAL sub-graphs here, which Earth Engine deduplicates. Only the five
// reduceRegion calls differ. Do not inline this back into detectVessels().
// Returns a 2-band image: 'target' (1 on accepted targets, masked elsewhere)
// and 'npix' (pixel count of the component each target pixel belongs to).
function buildDetectionImage(img, params) {
  img = ee.Image(img);

  var pols = params.detectionPols;
  var combined = null;
  var p, polName, dbBand, work, mu, thr, hit;

  for (p = 0; p < pols.length; p++) {
    polName = pols[p];
    dbBand = img.select(polName).rename('b');

    // dB -> linear power
    work = params.useLinearPower
      ? ee.Image(10).pow(dbBand.divide(10)).rename('b')
      : dbBand;

    // speckle pre-filter
    if (params.speckleApply) {
      work = work.focal_mean({
        radius: params.speckleWindowM, kernelType: 'circle', units: 'meters'
      }).rename('b');
    }

    // RECONCILED v2 THRESHOLD.
    //   was:  T = sample_mean + k * sample_sd        (k = 5, sd estimated)
    //   now:  T = mu_bg * (1 + k / sqrt(L))          (k = 8, sigma from ENL)
    // For speckle in linear power, sd = mean / sqrt(L) exactly, so the ENL form
    // is the SAME statistic with a variance that is known rather than estimated
    // from the same annulus that may contain another vessel. At k=8, L=4.4 the
    // multiplier is 4.81, i.e. +6.8 dB above local background.
    mu = annulusMu(work, params);
    thr = mu.multiply(1 + params.thresholdK / Math.sqrt(params.enl)).rename('b');

    hit = work.gt(thr);
    // absolute floor, always applied on the ORIGINAL dB band
    hit = hit.and(dbBand.gt(params.minBackscatterDb[polName]));
    hit = hit.rename('b');

    if (combined === null) {
      combined = hit;
    } else if (params.polCombine === 'and') {
      combined = combined.and(hit).rename('b');
    } else if (params.polCombine === 'first') {
      // keep the first polarisation only; later ones are deliberately ignored
    } else {
      combined = combined.or(hit).rename('b');   // default 'or'
    }
  }

  // valid-data mask of the scene, optionally eroded to kill GRD border noise
  var validMask = img.select(pols[0]).mask().rename('b');
  if (params.edgeErodeM > 0) {
    validMask = validMask.focal_min({
      radius: params.edgeErodeM, kernelType: 'circle', units: 'meters'
    }).rename('b');
  }

  // water + valid data
  var raw = combined.and(WATER_MASK.rename('b')).and(validMask).rename('b');

  // connected-component size filter: this is what rejects single-pixel speckle
  var npix = raw.selfMask().connectedPixelCount({
    maxSize: params.maxComponentSizePx,
    eightConnected: params.eightConnected
  }).rename('npix');

  var target = raw.selfMask()
                  .updateMask(npix.gte(params.minTargetPixels))
                  .rename('target');

  return target.addBands(npix.updateMask(target));
}

// --- the contracted entry point ---------------------------------------------
function detectVessels(img, aoi, params) {
  var det = buildDetectionImage(img, params);
  // v3 PERSISTENCE IS APPLIED HERE, NOT INSIDE buildDetectionImage().
  // That placement is load-bearing: buildPersistence() CALLS
  // buildDetectionImage() over every scene of an orbit, so applying the mask
  // inside it would make the detector depend on itself - infinite recursion in
  // the graph. Applying it one level up breaks the cycle and keeps
  // buildDetectionImage() a pure per-scene function.
  var persist = persistenceFor(img);
  var target = det.select('target').updateMask(persist.not());
  var npix = det.select('npix').updateMask(persist.not());

  // Unique component labels. connectedComponents() is the documented GEE idiom
  // for object counting; countDistinctNonNull over the label band gives the
  // number of distinct targets in the AOI.
  // (An algebraically equivalent alternative is sum(1/npix) over target pixels,
  // since each component of n pixels contributes n * 1/n = 1. It is cheaper but
  // returns a float that can drift at tile seams, so the label count is used.)
  var labels = target.connectedComponents({
    connectedness: params.eightConnected ? ee.Kernel.square(1) : ee.Kernel.plus(1),
    maxSize: params.maxComponentSizePx
  }).select('labels');

  var areaImg = ee.Image.pixelArea().updateMask(target).rename('area');
  var tAreaImg = npix.multiply(ee.Image.pixelArea()).updateMask(target).rename('tarea');

  var common = {
    geometry: aoi,
    scale: params.analysisScaleM,
    maxPixels: params.maxPixels,
    tileScale: params.tileScale,
    bestEffort: false
  };

  // Three reductions rather than one combined reducer: combining reducers with
  // sharedInputs:false silently depends on band ORDER for its output naming,
  // which is exactly the kind of fragility that should not sit in an export
  // that takes hours. Clarity wins here.
  var cnt = labels.reduceRegion({
    reducer: ee.Reducer.countDistinctNonNull(),
    geometry: common.geometry, scale: common.scale,
    maxPixels: common.maxPixels, tileScale: common.tileScale,
    bestEffort: common.bestEffort
  }).get('labels');

  var area = areaImg.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: common.geometry, scale: common.scale,
    maxPixels: common.maxPixels, tileScale: common.tileScale,
    bestEffort: common.bestEffort
  }).get('area');

  var maxT = tAreaImg.reduceRegion({
    reducer: ee.Reducer.max(),
    geometry: common.geometry, scale: common.scale,
    maxPixels: common.maxPixels, tileScale: common.tileScale,
    bestEffort: common.bestEffort
  }).get('tarea');

  // A scene with no targets reduces to null, not 0. Normalise to 0 so the CSV
  // never carries an empty cell that offline code could read as "missing".
  return ee.Dictionary({
    count:         ee.Number(ee.Algorithms.If(cnt, cnt, 0)),
    area_m2:       ee.Number(ee.Algorithms.If(area, area, 0)),
    max_target_m2: ee.Number(ee.Algorithms.If(maxT, maxT, 0))
  });
}

/* ---------------------------------------------------------------------------
 * v3: PERSISTENT-TARGET MASK
 *
 * buildPersistence(relOrbit) returns a 1/0 image: 1 where a pixel was detected
 * in at least PARAMS.persistenceThreshold of that relative orbit's scenes.
 *
 * CIRCULARITY WARNING. The mask is built BY the detector, using the current
 * PARAMS. Change the threshold, the polarisation or the analysis scale and the
 * mask changes too. A persistence asset is therefore only valid for the
 * parameter set that produced it - the asset name carries param_set_id for
 * exactly this reason, and reusing a mask across parameter sets is a silent
 * error, not a caught one.
 * ------------------------------------------------------------------------- */
function buildPersistence(relOrbit, params) {
  var col = SCENES.filter(ee.Filter.eq('relativeOrbitNumber_start', relOrbit));
  var n = ee.Number(col.size()).max(1);
  var hits = ee.ImageCollection(col.map(function (img) {
    return buildDetectionImage(ee.Image(img), params)
             .select('target').unmask(0).rename('h');
  })).sum().rename('h');
  var frac = hits.divide(n).rename('frac');
  var persist = frac.gte(params.persistenceThreshold).rename('persist');
  if (params.persistenceDilateM > 0) {
    persist = persist.focal_max({
      radius: params.persistenceDilateM, kernelType: 'circle', units: 'meters'
    }).rename('persist');
  }
  return persist.addBands(frac);
}

// One mask per relative orbit, built ONCE here rather than inside the per-scene
// map. In 'asset' mode these are loaded instead of recomputed, which is the
// whole point of STAGE A.
var PERSIST = {};
var PERSIST_FRAC = {};
(function () {
  var i, ro;
  for (i = 0; i < CONFIG.RELATIVE_ORBITS.length; i++) {
    ro = CONFIG.RELATIVE_ORBITS[i];
    if (CONFIG.PERSISTENCE_MODE === 'asset') {
      PERSIST[ro] = ee.Image(CONFIG.PERSISTENCE_ASSET_PREFIX + ro).select('persist');
      PERSIST_FRAC[ro] = ee.Image(CONFIG.PERSISTENCE_ASSET_PREFIX + ro).select('frac');
    } else if (CONFIG.PERSISTENCE_MODE === 'compute') {
      var pb = buildPersistence(ro, PARAMS);
      PERSIST[ro] = pb.select('persist');
      PERSIST_FRAC[ro] = pb.select('frac');
    } else {
      PERSIST[ro] = ee.Image.constant(0).rename('persist');
      PERSIST_FRAC[ro] = ee.Image.constant(0).rename('frac');
    }
  }
})();

// Look the mask up from an image's own orbit. ee.Dictionary keys must be
// strings; the relative orbit arrives as a number.
function persistenceFor(img) {
  if (CONFIG.PERSISTENCE_MODE === 'off') return ee.Image.constant(0).rename('persist');
  var d = ee.Dictionary({});
  var i, ro;
  for (i = 0; i < CONFIG.RELATIVE_ORBITS.length; i++) {
    ro = CONFIG.RELATIVE_ORBITS[i];
    d = d.set(ee.String(ee.Number(ro).format('%d')), PERSIST[ro]);
  }
  var key = ee.String(ee.Number(ee.Image(img).get('relativeOrbitNumber_start')).format('%d'));
  // A scene from an orbit not in RELATIVE_ORBITS falls back to "not persistent"
  // rather than throwing, so an unexpected orbit degrades to v2 behaviour.
  return ee.Image(ee.Algorithms.If(d.contains(key), d.get(key),
                                   ee.Image.constant(0))).rename('persist');
}

/* ---------------------------------------------------------------------------
 * v3: PER-DETECTION VECTORS
 *
 * One feature per detected target: centroid, size, peak backscatter, and its
 * RANGE from the hydrophone. All five radii are then a filter on range_m in
 * the offline join, so the radius sweep measures geometry rather than five
 * separate detector runs - and pairwise ambiguity flagging becomes possible.
 *
 * reduceToVectors applies its reducer to every band AFTER the label band.
 * npix is constant within a component, so max(npix) IS npix and the area
 * follows as npix * scale^2 - which avoids needing a sum reducer here.
 * ------------------------------------------------------------------------- */
function detectionVectors(img, aoi, params) {
  var det = buildDetectionImage(ee.Image(img), params);
  // Same recursion-avoidance as detectVessels(): mask applied here, not inside.
  var persist = persistenceFor(img);
  var target = det.select('target').updateMask(persist.not());
  var npix = det.select('npix').updateMask(persist.not());
  var vh = ee.Image(img).select(params.detectionPols[0]).rename('peak_db');

  var labels = target.connectedComponents({
    connectedness: params.eightConnected ? ee.Kernel.square(1) : ee.Kernel.plus(1),
    maxSize: params.maxComponentSizePx
  }).select('labels');

  var stack = labels.rename('label')
                    .addBands(vh.updateMask(target))
                    .addBands(npix.rename('npix').updateMask(target));

  var fc = stack.reduceToVectors({
    geometry: aoi,
    scale: params.analysisScaleM,
    geometryType: 'centroid',
    labelProperty: 'label',
    reducer: ee.Reducer.max(),
    maxPixels: params.maxPixels,
    tileScale: params.tileScale,
    bestEffort: false
  });

  return fc.limit(CONFIG.DETECTION_MAX_PER_SCENE);
}

/*  +----------------------- END PLUGGABLE BLOCK -------------------------+  */


/* =============================================================================
 * 5. SCENE COLLECTION
 * ========================================================================== */

var startDate = ee.Date(CONFIG.DATE_START);
var endDate = CONFIG.DATE_END_INCLUSIVE
  ? ee.Date(CONFIG.DATE_END).advance(1, 'day')
  : ee.Date(CONFIG.DATE_END);

var boundsGeom = (CONFIG.FILTER_BOUNDS === 'maxdisc') ? DISCS[DISCS.length - 1]
                                                      : SITE_POINT;

// s1Base is the verified set: bounds + date + IW, nothing else. It is kept so
// the dual-pol assumption can be CHECKED rather than asserted (CHECK 5 below).
var s1Base = ee.ImageCollection(CONFIG.S1_COLLECTION)
  .filterBounds(boundsGeom)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.eq('instrumentMode', CONFIG.INSTRUMENT_MODE));

var s1 = s1Base;
if (CONFIG.REQUIRE_DUAL_POL) {
  s1 = s1
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'));
}

s1 = s1.sort('system:time_start');


/* =============================================================================
 * 6. ERA5 WIND JOIN (nearest hour, one server-side pass)
 * ========================================================================== */

var era5 = ee.ImageCollection(CONFIG.ERA5_COLLECTION)
  .filterDate(startDate.advance(-1, 'day'), endDate.advance(1, 'day'))
  .filterBounds(SITE_POINT)
  .select([CONFIG.ERA5_U_BAND, CONFIG.ERA5_V_BAND], ['u10', 'v10']);

var windFilter = ee.Filter.maxDifference({
  difference: CONFIG.ERA5_MAX_DIFF_MILLIS,
  leftField: 'system:time_start',
  rightField: 'system:time_start'
});

var windJoin = ee.Join.saveBest({
  matchKey: 'era5_img',
  measureKey: 'era5_dt_millis'
});

var s1Joined = ee.ImageCollection(windJoin.apply(s1, era5, windFilter));

// Fallback image used when a scene finds no ERA5 match. It is a VALID image, so
// the arithmetic below never sees a null; the era5_matched flag (0) is what
// tells you to drop the row's wind offline. Do not read wind_speed_ms = 0 as
// "calm" without checking era5_matched.
var ERA5_DUMMY = ee.Image.constant([0, 0]).rename(['u10', 'v10']).toFloat();

function windDict(img) {
  // Presence is tested on the property NAME LIST, not by leaning on
  // ee.Algorithms.If's coercion of an object to "true". That coercion does work,
  // but a list membership test is unambiguous and costs nothing.
  var hasEra5 = ee.List(img.propertyNames()).contains('era5_img');
  var m = img.get('era5_img');

  var matched = ee.Number(ee.Algorithms.If(hasEra5, 1, 0));
  // Note the arithmetic is OUTSIDE the If. ee.Algorithms.If does not reliably
  // short-circuit, so a branch containing ee.Number(null).divide(...) can fail
  // the whole scene. Select a raw value first, then do the maths.
  var dtMillis = ee.Number(
    ee.Algorithms.If(hasEra5, img.get('era5_dt_millis'), -60000));
  var dtMin = dtMillis.divide(60000);

  var w = ee.Image(ee.Algorithms.If(hasEra5, m, ERA5_DUMMY));
  var u = w.select('u10');
  var v = w.select('v10');
  var spd = u.hypot(v).rename('wind_speed_ms');
  // Meteorological convention: direction the wind is coming FROM, degrees
  // clockwise from true north. from-north = 0, from-east = 90.
  var dir = u.atan2(v).multiply(180.0 / Math.PI).add(180).mod(360)
             .rename('wind_dir_deg');

  var vals = u.rename('wind_u10_ms')
              .addBands(v.rename('wind_v10_ms'))
              .addBands(spd)
              .addBands(dir)
              .reduceRegion({
                reducer: ee.Reducer.mean(),
                geometry: SITE_POINT,
                scale: CONFIG.ERA5_SCALE_M,
                maxPixels: 1e6,
                bestEffort: true
              });

  return ee.Dictionary({
    era5_matched:  matched,
    era5_dt_min:   dtMin,
    wind_u10_ms:   ee.Number(ee.Algorithms.If(vals.get('wind_u10_ms'),
                                              vals.get('wind_u10_ms'), 0)),
    wind_v10_ms:   ee.Number(ee.Algorithms.If(vals.get('wind_v10_ms'),
                                              vals.get('wind_v10_ms'), 0)),
    wind_speed_ms: ee.Number(ee.Algorithms.If(vals.get('wind_speed_ms'),
                                              vals.get('wind_speed_ms'), 0)),
    wind_dir_deg:  ee.Number(ee.Algorithms.If(vals.get('wind_dir_deg'),
                                              vals.get('wind_dir_deg'), 0))
  });
}


/* =============================================================================
 * 7. PER-SCENE FEATURE
 * ========================================================================== */
/*
 * radii/tags/aois are passed in so the interactive preview can use a cheap
 * subset (CONFIG.PREVIEW_RADII_M) while the export uses all five. The loop over
 * radii is a CLIENT-SIDE loop over five CONSTANT numbers, which is fine - there
 * is no client-side iteration over collection elements anywhere in this script,
 * and no getInfo() inside any loop.
 */
function makeSceneFeature(img, radiiM, tags, aois, discs) {
  img = ee.Image(img);

  var t = ee.Date(img.get('system:time_start'));
  var footprint = img.geometry();

  var props = {};
  props.site = CONFIG.SITE_NAME;
  props.param_set_id = CONFIG.PARAM_SET_ID;
  props.scene_id = img.get('system:index');
  props.timestamp_iso = t.format("YYYY-MM-dd'T'HH:mm:ss'Z'", 'UTC');
  props.timestamp_millis = img.get('system:time_start');
  props.date_utc = t.format('YYYY-MM-dd', 'UTC');
  props.hour_utc = t.get('hour');
  props.orbit_pass = img.get('orbitProperties_pass');
  props.relative_orbit = img.get('relativeOrbitNumber_start');
  props.orbit_number_start = img.get('orbitNumber_start');
  // Platform from the scene id ('S1A_IW_...' / 'S1B_IW_...'). Taken from the id
  // rather than the 'platform_number' property because the id prefix is present
  // on every archived scene and needs no assumption about metadata naming.
  props.platform = ee.String(img.get('system:index')).slice(0, 3);
  props.instrument_mode = img.get('instrumentMode');

  // --- ERA5 wind ---
  var wd = windDict(img);
  props.era5_matched = wd.get('era5_matched');
  props.era5_dt_min = wd.get('era5_dt_min');
  props.wind_u10_ms = wd.get('wind_u10_ms');
  props.wind_v10_ms = wd.get('wind_v10_ms');
  props.wind_speed_ms = wd.get('wind_speed_ms');
  props.wind_dir_deg = wd.get('wind_dir_deg');

  // --- detection + coverage, per radius ---
  var k, tag, d, covFrac;
  for (k = 0; k < radiiM.length; k++) {
    tag = tags[k];

    d = detectVessels(img, aois[k], PARAMS);
    props['vessels_' + tag] = d.get('count');
    props['area_m2_' + tag] = d.get('area_m2');
    props['max_target_m2_' + tag] = d.get('max_target_m2');

    // Partial-coverage detection. A scene that clips the edge of a disc sees
    // less water and therefore UNDERCOUNTS; such rows must be excludable
    // offline. cov_frac is the fraction of the disc inside the scene footprint;
    // full_cov is 1 when that clears PARAMS.fullCoverageMinFrac.
    // This is a FOOTPRINT-GEOMETRY test. It does not see interior no-data
    // (border noise, dropped bursts). A stricter alternative is to count
    // unmasked VV pixels inside the disc against the water-pixel total, at the
    // cost of one more reduceRegion per radius per scene.
    covFrac = discs[k].intersection(footprint, PARAMS.geomMaxErrorM)
                      .area(PARAMS.geomMaxErrorM)
                      .divide(discs[k].area(PARAMS.geomMaxErrorM));
    props['cov_frac_' + tag] = covFrac;
    props['full_cov_' + tag] = covFrac.gte(PARAMS.fullCoverageMinFrac);

    if (CONFIG.INCLUDE_WATER_AREA) {
      // Constant across scenes; attached per row for convenient offline
      // density normalisation. Indexed against the FULL radius list.
      props['water_area_m2_' + tag] = WATER_AREA_M2[indexOfRadius(radiiM[k])];
    }
  }

  // --- parameter provenance, so a sweep's CSVs are self-describing ---
  if (CONFIG.INCLUDE_PARAM_COLUMNS) {
    props.p_threshold_k = PARAMS.thresholdK;
    props.p_min_target_px = PARAMS.minTargetPixels;
    props.p_analysis_scale_m = PARAMS.analysisScaleM;
    props.p_speckle_window_m = PARAMS.speckleApply ? PARAMS.speckleWindowM : 0;
    props.p_bg_outer_m = PARAMS.bgOuterRadiusM;
    props.p_bg_guard_m = PARAMS.useGuardRing ? PARAMS.bgGuardRadiusM : 0;
    props.p_min_bg_px = PARAMS.minBackgroundPixels;
    props.p_edge_erode_m = PARAMS.edgeErodeM;
    props.p_pol_combine = PARAMS.polCombine;
    props.p_land_buffer_m = CONFIG.LAND_BUFFER_M;
    // RECONCILED v2 additions - the estimator is now part of the parameter set,
    // so a CSV cannot be mistaken for one produced by the sample-sd threshold.
    props.p_enl = PARAMS.enl;
    props.p_median_to_mean = PARAMS.medianToMean;
    props.p_use_median_bg = PARAMS.useMedianBackground ? 1 : 0;
    props.p_speckle_applied = PARAMS.speckleApply ? 1 : 0;
    props.p_detection_pols = PARAMS.detectionPols.join('+');
    props.p_threshold_db = 10 * Math.log10(1 + PARAMS.thresholdK / Math.sqrt(PARAMS.enl));
    props.p_land_source = CONFIG.LAND_SOURCE;
  }

  // Null geometry: the rows are joined offline by time, not by location, and a
  // geometry column would bloat the CSV for no benefit.
  return ee.Feature(null, props);
}

// Client-side lookup of a radius's position in the FULL radius list.
// Written as an explicit loop because Array.prototype.indexOf on numbers is
// fine in ES5 but indexOf/includes habits are what break this codebase; the
// loop is unambiguous.
function indexOfRadius(r) {
  var j;
  for (j = 0; j < CONFIG.RADII_M.length; j++) {
    if (CONFIG.RADII_M[j] === r) { return j; }
  }
  return 0;
}


/* =============================================================================
 * 8. BUILD THE OUTPUT TABLE
 * ========================================================================== */

var exportSource = s1Joined;
if (CONFIG.TEST_LIMIT > 0) {
  exportSource = s1Joined.limit(CONFIG.TEST_LIMIT);
  print('*** TEST_LIMIT ACTIVE: exporting only the first ' + CONFIG.TEST_LIMIT +
        ' scenes. This is a smoke test, NOT the study table. ***');
}

var features = ee.FeatureCollection(exportSource.map(function (img) {
  return makeSceneFeature(img, CONFIG.RADII_M, RADIUS_TAGS, AOIS, DISCS);
}));


/* =============================================================================
 * 9. DETERMINISTIC COLUMN ORDER
 * ========================================================================== */

var SELECTORS = [
  'site',
  'param_set_id',
  'scene_id',
  'timestamp_iso',
  'timestamp_millis',
  'date_utc',
  'hour_utc',
  'orbit_pass',
  'relative_orbit',
  'orbit_number_start',
  'platform',
  'instrument_mode',
  'era5_matched',
  'era5_dt_min',
  'wind_u10_ms',
  'wind_v10_ms',
  'wind_speed_ms',
  'wind_dir_deg'
];

for (i = 0; i < RADIUS_TAGS.length; i++) {
  SELECTORS.push('vessels_' + RADIUS_TAGS[i]);
  SELECTORS.push('area_m2_' + RADIUS_TAGS[i]);
  SELECTORS.push('max_target_m2_' + RADIUS_TAGS[i]);
  SELECTORS.push('cov_frac_' + RADIUS_TAGS[i]);
  SELECTORS.push('full_cov_' + RADIUS_TAGS[i]);
  if (CONFIG.INCLUDE_WATER_AREA) {
    SELECTORS.push('water_area_m2_' + RADIUS_TAGS[i]);
  }
}

if (CONFIG.INCLUDE_PARAM_COLUMNS) {
  SELECTORS.push('p_threshold_k');
  SELECTORS.push('p_min_target_px');
  SELECTORS.push('p_analysis_scale_m');
  SELECTORS.push('p_speckle_window_m');
  SELECTORS.push('p_bg_outer_m');
  SELECTORS.push('p_bg_guard_m');
  SELECTORS.push('p_min_bg_px');
  SELECTORS.push('p_edge_erode_m');
  SELECTORS.push('p_pol_combine');
  SELECTORS.push('p_land_buffer_m', 'p_enl', 'p_median_to_mean', 'p_use_median_bg', 'p_speckle_applied', 'p_detection_pols', 'p_threshold_db', 'p_land_source');
}


/* =============================================================================
 * 10. SANITY PRINTS - check these BEFORE starting the export
 * ========================================================================== */

print('=== SB02 Sentinel-1 vessel-detection export: pre-flight checks ===');
print('Site: ' + CONFIG.SITE_NAME + '  (' + CONFIG.SITE_LAT + ', ' +
      CONFIG.SITE_LON + ')');
print('Window: ' + CONFIG.DATE_START + ' .. ' + CONFIG.DATE_END +
      (CONFIG.DATE_END_INCLUSIVE ? ' (end inclusive)' : ' (end EXCLUSIVE)'));
print('Radii (m): ' + CONFIG.RADII_M.join(', ') +
      '   (largest = ' + MAX_RADIUS_M + ' m)');
print('Parameter set: ' + CONFIG.PARAM_SET_ID);

print('CHECK 1 - scene count (EXPECT 278):', s1.size());
print('CHECK 2 - by orbit pass (EXPECT ASCENDING 93, DESCENDING 185):',
      s1.aggregate_histogram('orbitProperties_pass'));
print('CHECK 3 - by relative orbit (EXPECT keys 142, 40, 62):',
      s1.aggregate_histogram('relativeOrbitNumber_start'));
// v3 FIX: was aggregate_histogram('platform_number'), which prints an empty
// dictionary if that property is absent or named differently - a silent blank
// rather than a caught error. Derived from the scene id prefix instead, which
// is the same source both output tables use.
print('CHECK 4 - by platform, from the scene id prefix (S1A / S1B split):',
      s1.map(function (im) {
        return im.set('platform_id',
                      ee.String(im.get('system:index')).slice(0, 3));
      }).aggregate_histogram('platform_id'));
// Dual-pol is CHECKED, not assumed: if these two numbers differ, some scene in
// the window is single-pol and the verified 278 no longer describes this set.
print('CHECK 5a - scenes before the dual-pol filter (EXPECT 278):', s1Base.size());
print('CHECK 5b - scenes after the dual-pol filter (EXPECT the same 278):',
      s1.size());
print('CHECK 6 - scenes surviving the ERA5 join (should equal CHECK 1):',
      s1Joined.size());

if (CONFIG.INCLUDE_WATER_AREA) {
  print('CHECK 7 - searchable water area per radius, m2 ' +
        '(land + harbour boxes removed):', ee.List(WATER_AREA_M2));
}

// The preview row. Built over CONFIG.PREVIEW_RADII_M only, so this print stays
// inside the interactive compute budget; the exported rows use all five radii.
// If this still times out, set CONFIG.PRINT_FIRST_FEATURE = false - the Export
// task does not share the interactive timeout.
/* v3 FIX - CHECK 8 IS NOT CHEAP ANY MORE IN 'compute' MODE, AND THAT WAS A
 * REGRESSION I INTRODUCED. The preview was designed to stay inside the
 * interactive budget by using PREVIEW_RADII_M (2 km) only. Persistence broke
 * that: makeSceneFeature -> detectVessels -> persistenceFor, and in 'compute'
 * mode that evaluates the DETECTOR OVER ALL ~93 SCENES of the orbit before it
 * can print one row. The radius no longer bounds the cost.
 * Worse, it would have failed FIRST, before CHECK 10 - the histogram that
 * decides whether the persistence threshold is defensible at all - ever
 * printed. So the preview is skipped in 'compute' mode and says why.
 * CHECK 1-7 and 9 do not touch persistence and still print normally. */
if (CONFIG.PRINT_FIRST_FEATURE && CONFIG.PERSISTENCE_MODE === 'compute') {
  print('CHECK 8 - SKIPPED. PERSISTENCE_MODE is "compute", so building one ' +
        'preview row would first evaluate the detector over every scene in ' +
        'the orbit and would very likely time out before CHECK 10 printed. ' +
        'Run STAGE A, switch PERSISTENCE_MODE to "asset", and CHECK 8 becomes ' +
        'cheap again. To see it now without persistence, set ' +
        'PERSISTENCE_MODE to "off".');
}
if (CONFIG.PRINT_FIRST_FEATURE && CONFIG.PERSISTENCE_MODE !== 'compute') {
  var previewTags = [];
  var previewAois = [];
  var previewDiscs = [];
  var idx;
  for (i = 0; i < CONFIG.PREVIEW_RADII_M.length; i++) {
    idx = indexOfRadius(CONFIG.PREVIEW_RADII_M[i]);
    previewTags.push(RADIUS_TAGS[idx]);
    previewAois.push(AOIS[idx]);
    previewDiscs.push(DISCS[idx]);
  }
  print('CHECK 8 - first scene, preview row at radii ' +
        CONFIG.PREVIEW_RADII_M.join(', ') + ' m ' +
        '(the exported rows carry all ' + CONFIG.RADII_M.length + ' radii):',
        makeSceneFeature(ee.Image(s1Joined.first()),
                         CONFIG.PREVIEW_RADII_M, previewTags,
                         previewAois, previewDiscs));
}

/* CHECK 10 - v3 PERSISTENCE DIAGNOSTIC.
 * READ THIS BEFORE TRUSTING persistenceThreshold. The 0.20 default is only
 * defensible if the per-pixel detection-frequency histogram is BIMODAL: a
 * fixed-object mode pushed up near 1.0, clearly separated from a traffic mode
 * decaying toward 0. If it is unimodal, there is no clean cut and 0.20 is
 * arbitrary - raise it, or set PERSISTENCE_MODE to 'off' and mask known fixed
 * objects by hand. SB02 sits IN a shipping lane, so the failure mode that
 * matters is masking a tightly-channelled traffic pixel as if it were a buoy. */
if (CONFIG.PERSISTENCE_MODE !== 'off') {
  for (i = 0; i < CONFIG.RELATIVE_ORBITS.length; i++) {
    var roD = CONFIG.RELATIVE_ORBITS[i];
    print('CHECK 10 - persistence fraction histogram, orbit ' + roD +
          ' (want BIMODAL; a spike near 1.0 is fixed objects)',
      PERSIST_FRAC[roD].updateMask(WATER_MASK).reduceRegion({
        reducer: ee.Reducer.histogram(20, 0.05),
        geometry: AOIS[AOIS.length - 1],
        scale: CONFIG.MASK_SCALE_M,
        maxPixels: PARAMS.maxPixels,
        tileScale: PARAMS.tileScale
      }));
    print('CHECK 10b - water area masked as persistent, orbit ' + roD + ' (m2)',
      ee.Image.pixelArea().updateMask(PERSIST[roD]).updateMask(WATER_MASK)
        .reduceRegion({
          reducer: ee.Reducer.sum(),
          geometry: AOIS[AOIS.length - 1],
          scale: CONFIG.MASK_SCALE_M,
          maxPixels: PARAMS.maxPixels,
          tileScale: PARAMS.tileScale
        }));
  }
}

print('CHECK 9 - exported column order (' + SELECTORS.length + ' columns):',
      SELECTORS);

print('REMINDER: this script produces a TABLE ONLY. It computes no statistics ' +
      'and makes no claim about proxy skill. All inference happens offline.');
print('NOTE: PARAMS carries the v2 specification values (VH-only, annulus-median ' +
      'background, k=8 on ENL-derived sigma = +6.82 dB, 3 px at 10 m). They are ' +
      'no longer placeholders. PARAMS.enl = 4.4 is the ESA nominal and is the ' +
      'one value that still wants measuring on calm scenes.');


/* =============================================================================
 * 11. MAP LAYERS (visual QA; tiles compute only for the visible viewport)
 * ========================================================================== */

if (CONFIG.SHOW_MAP) {
  Map.centerObject(SITE_POINT, 9);
  Map.addLayer(WATER_MASK.selfMask(), {palette: ['0b3d91']}, 'searchable water', false);
  for (i = 0; i < DISCS.length; i++) {
    Map.addLayer(ee.Image().byte().paint(ee.FeatureCollection([
        ee.Feature(DISCS[i])]), 1, 2),
      {palette: ['ffff00']}, 'disc ' + RADIUS_TAGS[i], i === DISCS.length - 1);
  }
  if (EXCLUDE_GEOM !== null) {
    Map.addLayer(ee.Image().byte().paint(
        ee.FeatureCollection([ee.Feature(EXCLUDE_GEOM)]), 1, 2),
      {palette: ['ff0000']}, 'harbour exclusions', true);
  }
  Map.addLayer(SITE_POINT, {color: 'ff00ff'}, 'SB02 hydrophone', true);

  var previewImg = ee.Image(s1Joined.first());
  Map.addLayer(previewImg.select('VV'), {min: -25, max: 0}, 'first scene VV', false);
  Map.addLayer(buildDetectionImage(previewImg, PARAMS).select('target'),
    {palette: ['ff8800']}, 'first scene detections', false);
}


/* =============================================================================
 * 12. EXPORT
 * ========================================================================== */

Export.table.toDrive({
  collection: features,
  description: CONFIG.EXPORT_DESCRIPTION + '_' + CONFIG.PARAM_SET_ID,
  folder: CONFIG.DRIVE_FOLDER,
  fileNamePrefix: CONFIG.FILE_PREFIX + '_' +
                  CONFIG.DATE_START + '_' + CONFIG.DATE_END + '_' +
                  CONFIG.PARAM_SET_ID,
  fileFormat: 'CSV',
  selectors: SELECTORS
});


/* =============================================================================
 * 13. v3 EXPORTS - STAGE A (persistence assets) AND THE DETECTION TABLE
 * ==========================================================================
 *
 * RUN ORDER, and it is not optional:
 *
 *   1. Set PERSISTENCE_MODE = 'compute' and run ONLY the STAGE A exports.
 *      Read CHECK 10 first - if the histogram is not bimodal, stop and
 *      reconsider the threshold rather than exporting a mask you cannot defend.
 *   2. Point PERSISTENCE_ASSET_PREFIX at the written assets and set
 *      PERSISTENCE_MODE = 'asset'.
 *   3. Run the per-scene and detection exports.
 *
 * Doing step 3 with PERSISTENCE_MODE = 'compute' is legal and produces the
 * same numbers, but it re-evaluates the detector over every scene of every
 * orbit inside every scene's own computation. On 278 scenes at 10 m that is
 * very unlikely to finish.
 *
 * The asset name carries PARAM_SET_ID because a persistence mask is only valid
 * for the parameters that built it - see the circularity warning above.
 */
if (CONFIG.PERSISTENCE_MODE === 'compute' && CONFIG.EXPORT_PERSISTENCE_ASSETS) {
  for (i = 0; i < CONFIG.RELATIVE_ORBITS.length; i++) {
    var roE = CONFIG.RELATIVE_ORBITS[i];
    Export.image.toAsset({
      image: buildPersistence(roE, PARAMS).clip(AOIS[AOIS.length - 1]),
      description: 'sb02_persist_orbit_' + roE + '_' + CONFIG.PARAM_SET_ID,
      assetId: CONFIG.PERSISTENCE_ASSET_PREFIX + roE + '_' + CONFIG.PARAM_SET_ID,
      region: AOIS[AOIS.length - 1],
      scale: PARAMS.analysisScaleM,
      maxPixels: PARAMS.maxPixels
    });
  }
}

/* The per-detection long table. One row per detected target, carrying its own
 * range from the hydrophone, so offline code derives every radius by filtering
 * range_m - and azimuth-ambiguity flagging (which needs pairwise geometry
 * between detections) becomes possible at all. */
if (CONFIG.EXPORT_DETECTIONS) {
  var detFeatures = ee.FeatureCollection(s1Joined.map(function (img) {
    img = ee.Image(img);
    var sid = ee.String(img.get('system:index'));
    var tms = ee.Number(img.get('system:time_start'));
    var ro  = img.get('relativeOrbitNumber_start');
    var pas = img.get('orbitProperties_pass');
    // v3 FIX: was img.get('platform_number'). The per-scene table derives
    // platform from the scene id prefix precisely because that property's
    // presence is an assumption, and the two tables disagreeing on the same
    // field is worse than either choice. Both now use the id.
    var plat = ee.String(img.get('system:index')).slice(0, 3);
    return detectionVectors(img, AOIS[AOIS.length - 1], PARAMS)
      .map(function (f) {
        f = ee.Feature(f);
        var pt = f.geometry();
        var npx = ee.Number(f.get('npix'));
        return f.set({
          site: CONFIG.SITE_NAME,
          param_set_id: CONFIG.PARAM_SET_ID,
          scene_id: sid,
          timestamp_millis: tms,
          timestamp_iso: ee.Date(tms).format('YYYY-MM-dd\'T\'HH:mm:ss\'Z\''),
          relative_orbit: ro,
          orbit_pass: pas,
          platform: plat,
          lon: pt.coordinates().get(0),
          lat: pt.coordinates().get(1),
          // Range from the hydrophone. THIS is what makes the radius sweep a
          // filter rather than five detector runs.
          range_m: pt.distance(SITE_POINT, PARAMS.geomMaxErrorM),
          npix: npx,
          area_m2: npx.multiply(PARAMS.analysisScaleM * PARAMS.analysisScaleM),
          peak_db: f.get('peak_db')
        });
      });
  })).flatten();

  Export.table.toDrive({
    collection: detFeatures,
    description: CONFIG.EXPORT_DESCRIPTION + '_detections_' + CONFIG.PARAM_SET_ID,
    folder: CONFIG.DRIVE_FOLDER,
    fileNamePrefix: CONFIG.FILE_PREFIX + '_detections_' +
                    CONFIG.DATE_START + '_' + CONFIG.DATE_END + '_' +
                    CONFIG.PARAM_SET_ID,
    fileFormat: 'CSV',
    selectors: ['site', 'param_set_id', 'scene_id', 'timestamp_iso',
                'timestamp_millis', 'relative_orbit', 'orbit_pass', 'platform',
                'lon', 'lat', 'range_m', 'npix', 'area_m2', 'peak_db']
  });
}
