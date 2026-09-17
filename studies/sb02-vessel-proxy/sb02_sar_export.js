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
  // PROTOTYPE PROFILE (p002), and what to restore for the full run.
  // STAGE A did not deliver: two 8-10 hour orbit tasks left the assets empty,
  // so persistence masking is unavailable and the study was blocked on it.
  // It should not have been. A persistent fixed object - a buoy, a wreck, the
  // SanctSound mooring - contributes a roughly CONSTANT count to every scene of
  // an orbit. Constant offsets do not change a rank correlation or an AUC at
  // all, and move a logistic fit's intercept rather than its discrimination.
  // For the question this study asks - does SAR vessel presence predict
  // ACOUSTIC ship presence - persistence masking is close to irrelevant. It
  // matters for absolute counts and per-km2 densities, which the prototype does
  // not claim.
  // NOT irrelevant, and stated rather than buried: a fixed object's
  // DETECTABILITY varies with sea state, and wind drives acoustic noise too, so
  // there is a confound path from wind to both sides. That is what the ERA5
  // columns are for, and proxy_skill.py conditions on them rather than trusting
  // a raw correlation.
  // Restore for the full run: RADII_M to the five radii, PERSISTENCE_MODE to
  // 'asset' once STAGE A completes, and bump PARAM_SET_ID.
  //   full:  RADII_M: [2000, 5000, 10000, 20000, 40000]
  RADII_M: [2000, 5000, 10000],

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
  //   full:  PERSISTENCE_MODE: 'asset'  (after STAGE A)
  PERSISTENCE_MODE: 'off',
  EXPORT_PERSISTENCE_ASSETS: true,   // STAGE A; see section 13. Only has an
                                     // effect in 'compute' mode, so the
                                     // prototype queues no asset tasks.
  // VERIFY THIS ROOT BEFORE RUNNING STAGE A - it is the one pre-flight check
  // no CHECK in this file can make. The export dialog validates the asset root
  // only when you press RUN on the task, so a root that does not exist fails
  // there and nowhere earlier: 'projects/ee-lathavis' was a guess and returned
  // "Resource projects/ee-lathavis could not be found." Read the real one off
  // the Code Editor's Assets tab (left panel) - it lists every root the account
  // can write to, as CLOUD ASSETS (projects/<id>/assets/) and LEGACY ASSETS
  // (users/<name>/). This account has stem-marine-engine and stem-earth-engine
  // as cloud projects plus the legacy users/lathavis2000; the marine project is
  // the one that matches this study.
  // AND CHANGE IT HERE, NOT IN THE EXPORT DIALOG. The dialog sets where this
  // one task writes. The READ path at PERSISTENCE_MODE 'asset' builds its id
  // from this same string, so fixing only the dialog leaves STAGE A writing to
  // one place and step 2 reading from another - the exact split PR #19 fixed,
  // and it costs a whole STAGE A run to rediscover.
  PERSISTENCE_ASSET_PREFIX: 'projects/stem-marine-engine/assets/sb02_persist_orbit_',
  RELATIVE_ORBITS: [142, 40, 62],
  // PROVISIONAL, from an exploratory query whose window and filters were not
  // recorded alongside it. A live run of the configuration in this file gave
  // 184 scenes over 2 relative orbits. Treat these as the number to re-measure,
  // not the number to trust; CHECK 1b reports the comparison either way.
  EXPECTED_SCENES: 278,
  EXPECTED_ASC: 93,
  EXPECTED_DESC: 185,
  // Radius of the bounds filter CHECK 3c probes with. NOT tied to RADII_M: the
  // 278 hypothesis is about a 40 km bounds filter, and RADII_M is 10 km under
  // the p002 prototype, so deriving this from MAX_RADIUS_M would quietly test
  // the wrong thing whenever the radius list changes.
  CHECK3C_BOUNDS_RADIUS_M: 40000,
  // CHECK 10 reduces at analysisScaleM (10 m). Over the full 40 km disc in
  // 'compute' mode that exceeded the Code Editor's memory limit on a live run.
  // The documented remedy is to shrink the REGION, never to coarsen the scale -
  // so here is the lever, defaulted to a radius that fits. Raise it once
  // PERSISTENCE_MODE is 'asset' and the detector is no longer being re-evaluated
  // per scene. The radius used is printed in the CHECK 10 label.
  CHECK10_RADIUS_M: 10000,
  // 'compute' mode samples rather than censuses - see the CHECK 10 block. 10k
  // pixels is far more than bimodality needs, and it bounds the OUTPUT.
  CHECK10_SAMPLE_N: 10000,
  // WHAT SAMPLING DOES NOT DO, corrected. An earlier note here claimed the
  // sampled path was "~300x cheaper than the disc". It is not: ee.Image.sample
  // bounds how many pixels come BACK, not how many are computed - the image is
  // still evaluated tile by tile across the whole region, so every pixel of it
  // still pays for its 81x81 annulus on every scene of the orbit. That is why
  // a 10 km sampled spot check still returned "User memory limit exceeded".
  // The input has to be bounded on BOTH axes, so it is, explicitly:
  //   region  - a small disc, independent of CHECK10_RADIUS_M (which is the
  //             asset-mode region, where the detector is not re-evaluated);
  //   scenes  - the first N acquisitions of the orbit, not all ~92.
  // THIRD SETTING, AND THE LAST ONE TO GUESS AT. 2 km x 20 scenes was an
  // arithmetic estimate (~1.1e10 kernel operations against the ~1.2e12 of the
  // exhaustive 10 km version) and a live run rejected it - not for memory this
  // time, but "Computation timed out", on both non-empty orbits. The memory
  // fixes held; what is left is the ~5 minute INTERACTIVE wall clock, and an
  // annulus MEDIAN over 4,316 weights is a selection per pixel per scene, which
  // is slower per operation than the kernel count suggests.
  // 1 km x 12 scenes is ~6.6x cheaper again. 1 km is chosen over a wider disc
  // on purpose: the SanctSound mooring sits AT the site, so the pixels most
  // likely to carry a genuine near-1.0 fixed-object mode are the ones closest
  // to it. 12 scenes is chosen over 8 because persistence is a question about
  // repetition in TIME, and 8 repeats quantise frac to steps of 0.125, which is
  // too coarse to see a shape in.
  // IF IT TIMES OUT AGAIN, DO NOT TUNE A FOURTH TIME. A spot check that will
  // not run is a missing gate, not a blocker: go to STAGE A. Section 13 sets
  // out why exporting before the threshold is settled is recoverable - frac is
  // exported raw, so a different threshold is a re-read, not a re-export.
  // Neither lever touches analysisScaleM: coarsening the scale would measure a
  // different detector (see the CHECK 10 block).
  CHECK10_SPOT_RADIUS_M: 1000,
  // The first N scenes by acquisition time, NOT a random draw - EE has no
  // cheap deterministic subsample of an ImageCollection. For one relative
  // orbit that is N consecutive 12-day repeats (12 -> ~5 months), enough to
  // show whether a near-1.0 mode exists at all. It is not the orbit's
  // persistence and the label says so. Note that it also quantises frac to
  // steps of 1/N, so the 0.20 threshold cannot be READ off this check - 12
  // scenes put it between 2/12 and 3/12. Bimodality is the question here;
  // the threshold is settled in 'asset' mode.
  CHECK10_MAX_SCENES: 12,
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

  // HOW THE 1000 m LAND BUFFER IS GROWN. 'distance' (default) uses
  // fastDistanceTransform; 'focal' uses the original focal_max. They select
  // the SAME set of pixels - focal_max(r) is exactly "within r of land" - but
  // focal_max sizes its kernel in PIXELS at the REQUEST scale, so at the 10 m
  // analysis scale it is a radius-100 px disc, ~31,400 weights per output
  // pixel, evaluated wherever WATER_MASK is read: the detector (every scene),
  // CHECK 7, CHECK 10 and CHECK 10b. That made the land mask the single
  // largest cost multiplier in the graph, and it is not where the science is.
  // Keep 'focal' only to reproduce a pre-v3.1 run.
  LAND_DILATE_METHOD: 'distance',       // 'distance' | 'focal'
  // Search radius for the distance transform, in PIXELS of the request
  // projection. It must exceed LAND_BUFFER_M at every scale WATER_MASK is read
  // at, or an ocean pixel beyond the search radius reports a clamped distance:
  // 256 px is 2560 m at 10 m and 7680 m at MASK_SCALE_M, both well past 1000 m.
  LAND_DILATE_NEIGHBOURHOOD_PX: 256,

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
  ERA5_MAX_DIFF_MILLIS: 3600000,  // +/- 1 h; the NEAREST hour inside it wins
  ERA5_SCALE_M: 27830,            // ERA5 native ~0.25 deg
  // WHY THIS IS NO LONGER AN ee.Join.saveBest, and it is a fix, not a taste.
  // The join was applied to the WHOLE filtered ERA5 collection: ~37 months of
  // hourly imagery is ~27,000 right-hand elements, matched against 184 scenes,
  // and the result had to be materialised before anything downstream could be
  // counted. On a live run CHECK 6 - which was only s1Joined.size() - returned
  // "User memory limit exceeded". That collection feeds the per-scene table,
  // the detection table AND buildPersistence, so the same cost sat under every
  // export, not just the print that exposed it.
  // The replacement is era5At(): per scene, filterDate a +/-1 h window (at most
  // three hourly images), sort by |dt| and take the first. It is still one
  // server-side pass with no getInfo and no client loop, it reads the same
  // NEAREST value the join did, and it touches ~3 images per scene instead of
  // 27,000 for the collection. It also honours ERA5_MAX_DIFF_MILLIS for any
  // collection, not just an hourly one.
  // era5_dt_min is now SIGNED (ERA5 hour minus acquisition time, minutes).
  // saveBest's measureKey was an ABSOLUTE difference, so the old comment here
  // promising a signed value was wrong about its own output.
  // If ECMWF/ERA5/HOURLY's ingested range ever stops short of 2021-12, rows
  // will carry era5_matched = 0; swap ERA5_COLLECTION rather than patching
  // around it.

  // ---- output -------------------------------------------------------------
  // p002 = the prototype profile above. The DETECTION parameters are identical
  // to p001; what differs is PERSISTENCE_MODE and the radius list, and both
  // change the numbers in the table, so they get their own id. The mode is also
  // echoed into every row now - see p_persistence_mode - because a table built
  // with persistence off and one built with it on are not comparable and the
  // CSV could not previously tell you which it was.
  PARAM_SET_ID: 'p002',                 // bump this for every parameter sweep
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

// ESA/WorldCover/v200 is an IMAGE COLLECTION, not an Image - the STAC record
// says gee:type: image_collection. ee.Image(<collection id>) builds a valid
// lazy node, so nothing failed at parse or graph build; it failed on
// EVALUATION with "Image.load: ... is not an image". Because WATER_MASK feeds
// CHECK 7, CHECK 10/10b, the detector and all three exports, the symptom was
// CHECK 1-6 printing fine, CHECK 7 throwing, and every export task failing
// later with an opaque asset-load error.
//   .mosaic() rather than .first(): mosaic is correct whether the collection
//   holds one global image or many tiles, while .first() is only correct in
//   the first case and would silently cover part of the AOI in the second.
//   The STAC record does not state the tiling, so take the option that does
//   not depend on knowing it.
var landRaw = ee.ImageCollection(CONFIG.LAND_SOURCE).select(CONFIG.LAND_BAND).mosaic();
// hasData=1 where WorldCover mapped the pixel at all; 0 where it did not.
var landHasData = landRaw.mask().gt(0);
// land = mapped AND not the permanent-water class. Unmapped -> 0 -> water.
var landBinary = landHasData.and(landRaw.neq(CONFIG.LAND_WATER_CLASS)).unmask(0);

// GROWING THE 1000 m BUFFER - THE COST FIX, same set of pixels either way.
// focal_max's kernel is sized in PIXELS at the REQUEST scale, so a 1000 m
// circle at the 10 m analysis scale is radius 100 px, ~31,400 weights per
// output pixel. WATER_MASK is read by the detector on every scene and by
// CHECK 7, CHECK 10 and CHECK 10b, so that factor multiplied the whole graph:
// the CHECK 10 histogram for the EMPTY orbit 142, whose frac image is a
// constant, still had to evaluate it over 3.1e6 pixels of a 10 km disc.
// fastDistanceTransform is a linear-time transform returning the SQUARED
// distance in pixels to the nearest non-zero (land) pixel. "Within 1000 m of
// land" is then one comparison, and it is EXACTLY what focal_max(1000 m)
// selects - a dilation by a disc of radius r is the set at distance <= r.
var landDilated;
if (CONFIG.LAND_DILATE_METHOD === 'focal') {
  landDilated = landBinary.focal_max({
    radius: CONFIG.LAND_BUFFER_M,
    kernelType: 'circle',
    units: 'meters'
  });
} else {
  // Squared distance in pixels -> distance in pixels -> metres. pixelArea() is
  // the area of the pixel in the REQUEST projection, so its square root is the
  // pixel side in metres at whatever scale the consumer asked for, which is
  // the same projection the transform counted pixels in.
  var landDistM = landBinary
    .fastDistanceTransform({
      neighborhood: CONFIG.LAND_DILATE_NEIGHBOURHOOD_PX,
      units: 'pixels',
      metric: 'squared_euclidean'
    })
    .sqrt()
    .multiply(ee.Image.pixelArea().sqrt());
  // Beyond the search radius the transform returns a CLAMPED distance. That is
  // safe here only because the clamp (2560 m at 10 m, 7680 m at MASK_SCALE_M)
  // is larger than LAND_BUFFER_M - see LAND_DILATE_NEIGHBOURHOOD_PX. If either
  // the buffer grows or the neighbourhood shrinks, open ocean starts reporting
  // a distance below the buffer and the mask eats the search area.
  landDilated = landDistM.lte(CONFIG.LAND_BUFFER_M);
}
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
function buildPersistence(relOrbit, params, maxScenes) {
  // s1Scenes, not s1: every consumer downstream (the per-scene table, the
  // detection table, the CHECK 8 preview) maps over s1Scenes, so the mask is
  // built from exactly the scene set the detections come from. Under the old
  // ee.Join.saveBest those could differ - a scene that found no ERA5 match was
  // dropped from the joined collection - and this line is what kept the mask
  // and the detections in step. The join is gone (section 6), so s1Scenes IS
  // s1 and they cannot diverge at all; the reference is kept rather than
  // rewritten to s1 so that reinstating a filtering join cannot silently
  // reintroduce the mismatch.
  var col = s1Scenes.filter(ee.Filter.eq('relativeOrbitNumber_start', relOrbit));
  // The CHECK 10 spot check bounds the scene count as well as the region -
  // see CONFIG.CHECK10_MAX_SCENES. maxScenes is the FIRST n acquisitions, and
  // n below is then that subset's size, so frac stays a true fraction of what
  // was actually looked at. Omitted (STAGE A, the exports) means every scene.
  if (maxScenes && maxScenes > 0) {
    col = ee.ImageCollection(col.limit(maxScenes, 'system:time_start'));
  }
  var n = ee.Number(col.size()).max(1);

  // An orbit with NO scenes in the configured window used to crash here:
  // ee.ImageCollection([]).sum() is an image with ZERO bands, and dividing it
  // by a number raised
  //   "Image.divide: If one image has no bands, the other must also have no
  //    bands. Got 0 and 1."
  // which is exactly what a live run produced for orbit 142 while the other two
  // orbits ran. RELATIVE_ORBITS is a fixed client list, so it can name an orbit
  // the date window does not contain, and the whole diagnostic then died on the
  // first such orbit rather than reporting it.
  //
  // Merging one all-zero image guarantees .sum() yields a banded image. It adds
  // nothing to the hit count, and n is the TRUE scene count, so a non-empty
  // orbit is arithmetically unchanged. An empty orbit now yields frac = 0
  // everywhere, i.e. "nothing is persistent here", which is the correct reading
  // of no evidence - and CHECK 3b below says plainly that the orbit was empty,
  // so a zero mask is never mistaken for a measured one.
  var zeroHit = ee.Image.constant(0).rename('h').toFloat();
  var hits = ee.ImageCollection(col.map(function (img) {
    return buildDetectionImage(ee.Image(img), params)
             .select('target').unmask(0).rename('h').toFloat();
  })).merge(ee.ImageCollection([zeroHit])).sum().rename('h');
  var frac = hits.divide(n).rename('frac');
  var persist = frac.gte(params.persistenceThreshold).rename('persist');
  if (params.persistenceDilateM > 0) {
    persist = persist.focal_max({
      radius: params.persistenceDilateM, kernelType: 'circle', units: 'meters'
    }).rename('persist');
  }
  return persist.addBands(frac);
}

// The PERSIST tables are built further down, immediately before section 8,
// because buildPersistence() reads s1Scenes and that collection does not
// exist yet at this point in the file. Do not move them back up here: this
// block runs at module level, so it would read an undeclared variable and
// the script would die before CHECK 1 prints.

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
 * 6. ERA5 WIND LOOKUP (nearest hour, per scene, no collection-wide join)
 * ========================================================================== */
/*
 * WAS ee.Join.saveBest OVER THE WHOLE ERA5 COLLECTION, AND THAT IS WHAT BROKE
 * CHECK 6. The join matched 184 scenes against ~27,000 hourly ERA5 images and
 * had to materialise the result before the joined collection could even be
 * counted; a live run returned "User memory limit exceeded" on
 * `print(..., s1Joined.size())`, which is nothing but a count. The joined
 * collection is also what the per-scene table, the detection table and
 * buildPersistence() all map over, so the cost was under every export too.
 *
 * era5At() asks the same question per scene: filterDate a +/-ERA5_MAX_DIFF
 * window, which for an hourly collection is at most three images, sort by
 * |dt| and take the nearest. Same value as saveBest, same single server-side
 * pass, no getInfo and no client-side loop over collection elements - just
 * ~3 candidates per scene instead of 27,000 for the whole collection.
 */

var era5 = ee.ImageCollection(CONFIG.ERA5_COLLECTION)
  .filterDate(startDate.advance(-1, 'day'), endDate.advance(1, 'day'))
  .select([CONFIG.ERA5_U_BAND, CONFIG.ERA5_V_BAND], ['u10', 'v10']);

// The scene collection every table and the persistence masks are built from.
// It is now just s1: the wind is attached per scene by windDict() rather than
// by a join, so there is no second collection that could silently hold a
// different scene set. CHECK 6 therefore checks what actually matters - how
// many scenes FOUND an ERA5 hour - instead of re-counting the same images.
var s1Scenes = s1;

// Fallback image used when a scene finds no ERA5 match. It is a VALID image, so
// the arithmetic below never sees a null; the era5_matched flag (0) is what
// tells you to drop the row's wind offline. Do not read wind_speed_ms = 0 as
// "calm" without checking era5_matched.
var ERA5_DUMMY = ee.Image.constant([0, 0]).rename(['u10', 'v10']).toFloat();

// Nearest ERA5 image to one scene, within +/-CONFIG.ERA5_MAX_DIFF_MILLIS.
// Returns a CLIENT-side object holding server-side values, so the caller can
// pull them out without an ee.Dictionary round trip.
function era5At(img) {
  var t = ee.Number(ee.Image(img).get('system:time_start'));
  var tol = CONFIG.ERA5_MAX_DIFF_MILLIS;
  // filterDate's end is exclusive, so +tol+1 keeps a match exactly at +tol.
  var win = era5.filterDate(ee.Date(t.subtract(tol)), ee.Date(t.add(tol + 1)))
    .map(function (e) {
      e = ee.Image(e);
      var dt = ee.Number(e.get('system:time_start')).subtract(t);
      // dt is SIGNED (ERA5 hour minus acquisition). abs() is only the sort key
      // - it never reaches the CSV, so the sign survives into era5_dt_min.
      return e.set({ era5_dt_millis: dt, era5_abs_dt: dt.abs() });
    })
    .sort('era5_abs_dt');
  var has = win.size().gt(0);
  return {
    matched: ee.Number(has),
    // Arithmetic OUTSIDE the If, as below: pull a raw value first.
    dtMillis: ee.Number(ee.Algorithms.If(has, win.first().get('era5_dt_millis'),
                                         -60000)),
    image: ee.Image(ee.Algorithms.If(has, win.first(), ERA5_DUMMY))
  };
}

function windDict(img) {
  var e = era5At(img);
  var matched = e.matched;
  // Note the arithmetic is OUTSIDE the If. ee.Algorithms.If does not reliably
  // short-circuit, so a branch containing ee.Number(null).divide(...) can fail
  // the whole scene. Select a raw value first, then do the maths.
  var dtMin = e.dtMillis.divide(60000);

  var w = e.image;
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
    // Math.log10 is an ES2015 LIBRARY addition and the Code Editor sandbox does
    // not implement it - the same class as String.repeat (v10.149) and
    // Object.assign (v10.158), both of which broke this codebase in the browser
    // while passing node --check. This line is on the default path
    // (INCLUDE_PARAM_COLUMNS is true) and map() invokes the callback once at
    // graph-build time, so it threw at MODULE LEVEL before CHECK 1 printed.
    props.p_threshold_db =
      10 * Math.log(1 + PARAMS.thresholdK / Math.sqrt(PARAMS.enl)) / Math.LN10;
    props.p_land_source = CONFIG.LAND_SOURCE;
    // PERSISTENCE PROVENANCE. Without these two a persistence-off table and a
    // persistence-on one are byte-indistinguishable apart from their numbers,
    // which is exactly the kind of unfalsifiable provenance the rest of this
    // file exists to avoid. The threshold is echoed even when the mode is
    // 'off', where it is inert, so the column is never missing.
    props.p_persistence_mode = CONFIG.PERSISTENCE_MODE;
    props.p_persistence_threshold = PARAMS.persistenceThreshold;
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

// Built HERE, not beside buildPersistence(): this runs at module level and
// buildPersistence() reads s1Scenes, which is not assigned until the ERA5
// section above. Sitting next to the function it calls, it ran ~115 lines too
// early and the script threw before any CHECK could print. Everything that
// consumes PERSIST (CHECK 9, the preview, STAGE A) comes after this point.
// One mask per relative orbit, built ONCE here rather than inside the per-scene
// map. In 'asset' mode these are loaded instead of recomputed, which is the
// whole point of STAGE A.
var PERSIST = {};
var PERSIST_FRAC = {};
(function () {
  var i, ro, aid;
  for (i = 0; i < CONFIG.RELATIVE_ORBITS.length; i++) {
    ro = CONFIG.RELATIVE_ORBITS[i];
    if (CONFIG.PERSISTENCE_MODE === 'asset') {
      // The id MUST match what STAGE A writes, PARAM_SET_ID included. A mask
      // built under one parameter set is not valid for another - that is the
      // whole reason the asset name carries the id. Reading without the
      // suffix asked for an asset STAGE A never wrote.
      aid = CONFIG.PERSISTENCE_ASSET_PREFIX + ro + '_' + CONFIG.PARAM_SET_ID;
      PERSIST[ro] = ee.Image(aid).select('persist');
      PERSIST_FRAC[ro] = ee.Image(aid).select('frac');
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

var exportSource = s1Scenes;
if (CONFIG.TEST_LIMIT > 0) {
  exportSource = s1Scenes.limit(CONFIG.TEST_LIMIT);
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
  SELECTORS.push('p_persistence_mode', 'p_persistence_threshold');
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

/* THE EXPECTED COUNTS ARE PROVISIONAL. They came from an exploratory query
 * whose date window, bounds filter and polarisation filter were never written
 * down beside them, and a live run of THIS configuration returned 184 scenes
 * across 2 relative orbits, not 278 across 3. One of the two is wrong and the
 * file cannot tell you which, because the conditions behind 278 were not
 * recorded. So they are a stated expectation to be compared against, not an
 * assertion of fact, and CHECK 1b says MATCH or MISMATCH rather than leaving a
 * reader to notice. Confirm the window you want, re-measure, and update
 * EXPECTED_* together with the conditions - a count without its query is the
 * same kind of unverifiable claim this file exists to avoid. */
print('CHECK 1 - scene count (expected ' + CONFIG.EXPECTED_SCENES +
      ', PROVISIONAL - see note in source):', s1.size());
print('CHECK 1b - does the scene count match the stated expectation?',
      ee.Algorithms.If(s1.size().eq(CONFIG.EXPECTED_SCENES),
        'MATCH',
        ee.String('MISMATCH: this window yields ').cat(s1.size().format('%d'))
          .cat(', the recorded expectation is ')
          .cat(ee.Number(CONFIG.EXPECTED_SCENES).format('%d'))
          .cat('. Neither is authoritative until the query behind the ')
          .cat('expectation is restated.')));
print('CHECK 2 - by orbit pass (expected ASCENDING ' + CONFIG.EXPECTED_ASC +
      ', DESCENDING ' + CONFIG.EXPECTED_DESC + ', PROVISIONAL):',
      s1.aggregate_histogram('orbitProperties_pass'));
print('CHECK 3 - by relative orbit (configured: ' +
      CONFIG.RELATIVE_ORBITS.join(', ') + '):',
      s1.aggregate_histogram('relativeOrbitNumber_start'));
/* CHECK 3b. RELATIVE_ORBITS is a fixed client-side list, so it can name an
 * orbit the window does not contain. That used to surface as a crash inside
 * buildPersistence (a 0-band sum divided by a number) rather than as a fact
 * about the data. It is a fact about the data, so report it as one: an empty
 * orbit now yields an all-zero mask, and this line is what stops that zero
 * being read as "nothing persistent here" when it means "nothing looked at". */
print('CHECK 3b - configured orbits with NO scenes in this window ' +
      '(each yields an ALL-ZERO persistence mask, not a measured one):',
      ee.List(CONFIG.RELATIVE_ORBITS).removeAll(
        s1.aggregate_array('relativeOrbitNumber_start').distinct()
          .map(function (o) { return ee.Number(o).toInt(); })));
// v3 FIX: was aggregate_histogram('platform_number'), which prints an empty
// dictionary if that property is absent or named differently - a silent blank
// rather than a caught error. Derived from the scene id prefix instead, which
// is the same source both output tables use.
/* CHECK 3c-3f - THE 278 HYPOTHESIS, TESTED RATHER THAN ARGUED.
 *
 * THE ARITHMETIC THAT MOTIVATES IT. A live run of the shipped configuration
 * gives 184 scenes: 92 ASCENDING + 92 DESCENDING, orbits 40 and 62 at 92 each,
 * and orbit 142 absent. The recorded expectation is 278 = 93 ASC + 185 DESC.
 * A third, DESCENDING orbit of ~93 scenes closes both gaps at once:
 *     DESC  92 + 93 = 185   ASC  92 + 1 = 93   total 278
 * Orbit 142 is the obvious candidate - it can cover part of a 40 km disc around
 * the site without covering the hydrophone itself, which is exactly the set a
 * 'maxdisc' bounds filter takes in and a 'point' filter does not.
 *
 * WHY THIS IS A CHECK AND NOT AN INSTRUCTION TO FLIP FILTER_BOUNDS. Flipping
 * it and flipping it back is a footgun: left flipped, the EXPORT silently
 * carries a different scene set than the one CHECK 1 reported, and nothing
 * downstream would say so. This probe builds its own collection, inside an
 * IIFE so the identifier cannot be reached from the export path at all, and
 * CONFIG.FILTER_BOUNDS is untouched. It costs nothing to leave on: metadata
 * filters only, no pixels.
 *
 * READ 3f FIRST. If it says CONFIRMED, the recorded 278 was a maxdisc query and
 * EXPECTED_* can finally be restated WITH its conditions. If NOT CONFIRMED, the
 * hypothesis is dead and 278 came from somewhere else - a different window, a
 * different site, or a different product - and it should be retired rather than
 * carried around as an expectation nothing can reproduce.
 */
(function () {
  var probeBounds = SITE_POINT.buffer(CONFIG.CHECK3C_BOUNDS_RADIUS_M);
  var probe = ee.ImageCollection(CONFIG.S1_COLLECTION)
    .filterBounds(probeBounds)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.eq('instrumentMode', CONFIG.INSTRUMENT_MODE));
  if (CONFIG.REQUIRE_DUAL_POL) {
    probe = probe
      .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
      .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'));
  }
  var km = CONFIG.CHECK3C_BOUNDS_RADIUS_M / 1000;
  print('CHECK 3c - scene count under a ' + km + ' km BOUNDS filter instead of ' +
        'the hydrophone point (diagnostic only - FILTER_BOUNDS is unchanged and ' +
        'the export still uses "' + CONFIG.FILTER_BOUNDS + '"):', probe.size());
  print('CHECK 3d - that set by orbit pass (the hypothesis predicts ASCENDING ' +
        CONFIG.EXPECTED_ASC + ', DESCENDING ' + CONFIG.EXPECTED_DESC + '):',
        probe.aggregate_histogram('orbitProperties_pass'));
  print('CHECK 3e - that set by relative orbit (the hypothesis predicts orbit ' +
        '142 present, DESCENDING, with ~93 scenes):',
        probe.aggregate_histogram('relativeOrbitNumber_start'));
  print('CHECK 3f - VERDICT on the recorded ' + CONFIG.EXPECTED_SCENES + ':',
    ee.Algorithms.If(probe.size().eq(CONFIG.EXPECTED_SCENES),
      ee.String('CONFIRMED: a ').cat(ee.Number(km).format('%d'))
        .cat(' km bounds filter reproduces the recorded ')
        .cat(ee.Number(CONFIG.EXPECTED_SCENES).format('%d'))
        .cat(' exactly. The recorded figures were a maxdisc query, not a point ')
        .cat('query. Restate EXPECTED_* WITH that condition, or set ')
        .cat('FILTER_BOUNDS to "maxdisc" if that wider set is the one you want ')
        .cat('- they are different studies, not different spellings.'),
      ee.String('NOT CONFIRMED: this bounds filter yields ')
        .cat(probe.size().format('%d')).cat(', not ')
        .cat(ee.Number(CONFIG.EXPECTED_SCENES).format('%d'))
        .cat('. The orbit-142 hypothesis does not explain the recorded figures, ')
        .cat('so they came from some other query - a different window, site or ')
        .cat('product. Retire them rather than carry an expectation nothing ')
        .cat('reproduces. Compare CHECK 3d and 3e to see how it differs.')));
})();

print('CHECK 4 - by platform, from the scene id prefix (S1A / S1B split):',
      s1.map(function (im) {
        return im.set('platform_id',
                      ee.String(im.get('system:index')).slice(0, 3));
      }).aggregate_histogram('platform_id'));
// Dual-pol is CHECKED, not assumed: if these two numbers differ, some scene in
// the window is single-pol and the verified 278 no longer describes this set.
print('CHECK 5a - scenes before the dual-pol filter:', s1Base.size());
print('CHECK 5b - scenes after the dual-pol filter (5a and 5b must agree; ' +
      'if they differ, some scene in the window is single-pol):', s1.size());
// CHECK 6 USED TO BE s1Joined.size() - the size of the saveBest join - and on
// a live run that print itself returned "User memory limit exceeded" (see
// section 6). There is no join now, so re-counting the collection would be a
// tautology. The question the check was FOR still stands, so it is asked
// directly: how many scenes found an ERA5 hour inside ERA5_MAX_DIFF_MILLIS.
// Anything below CHECK 1 means those rows carry era5_matched = 0 and their
// wind columns are zeros, NOT calm.
print('CHECK 6 - scenes with an ERA5 hour within +/-' +
      (CONFIG.ERA5_MAX_DIFF_MILLIS / 60000) + ' min (should equal CHECK 1):',
      s1Scenes.map(function (im) {
        return ee.Image(im).set('era5_ok', era5At(im).matched);
      }).aggregate_sum('era5_ok'));

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
        makeSceneFeature(ee.Image(s1Scenes.first()),
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
/* CHECK 10 - THE PERSISTENCE DIAGNOSTIC, AND WHEN IT IS AFFORDABLE.
 *
 * Both reductions run at PARAMS.analysisScaleM, the scale the mask is applied
 * at. Reducing at MASK_SCALE_M (30 m) measured a different detector entirely -
 * annulusKernel() is a FIXED kernel sized in PIXELS with no .reproject(), so a
 * 30 m request inflates the 400/150 m annulus to ~1200/450 m and the 300 m2
 * minimum target to ~2,700 m2, suppressing the near-1.0 mode this check exists
 * to find. In 'asset' mode a coarse read is worse still: it serves a pyramid
 * level, and block-averaging frac is precisely what destroys bimodality.
 *
 * COST, MEASURED ON THREE LIVE RUNS. In 'compute' mode an exhaustive reduction
 * is not affordable at any useful region size. The detector runs per scene, and
 * every pixel needs an 81x81 annulus (4,316 weights):
 *
 *     10 km disc : 3.1e6 px x 92 scenes x 4,316 = 1.2e12 kernel ops
 *      2 km disc : 1.3e5 px x 92 scenes x 4,316 = 5.0e10
 *
 * The 40 km disc failed first, then the 10 km disc, both with "User memory
 * limit exceeded".
 *
 * THE THIRD RUN IS THE ONE THAT MATTERS, because it falsified the fix. Sampling
 * 10,000 px over the 10 km disc failed the SAME way, on all three orbits -
 * including orbit 142, which CHECK 3b reports as EMPTY and whose frac image is
 * therefore a constant. A constant image cannot cost 1.2e12 operations, so the
 * cost was never only the detector. Two things were wrong:
 *
 *   1. sample() bounds the OUTPUT, not the input. The claim that it was "~300x
 *      cheaper" was wrong: numPixels caps how many pixels come back, while the
 *      image is still evaluated tile by tile across the whole region. Every
 *      pixel still paid for its annulus on every scene.
 *   2. WATER_MASK, which both checks apply, was a focal_max with a 1000 m
 *      circle kernel - radius 100 px at the 10 m analysis scale, ~31,400
 *      weights per pixel - on top of everything else. That is what made the
 *      EMPTY orbit expensive, and it is fixed in the LAND MASKING section:
 *      fastDistanceTransform selects the identical set of pixels.
 *
 * So the two modes ask the question differently, and say which they used:
 *
 *   'compute' - a RANDOM SAMPLE of CHECK10_SAMPLE_N pixels, binned to 0.05,
 *      over CHECK10_SPOT_RADIUS_M and the FIRST CHECK10_MAX_SCENES scenes of
 *      the orbit. Bimodality is a question about a distribution's SHAPE, which
 *      a sample answers honestly, but the sample only helps if the INPUT is
 *      bounded too - hence the second cap. It is a spot check before committing
 *      to STAGE A, not the study-area histogram, and the label says so. Neither
 *      lever is analysisScaleM; coarsening that measures a different detector.
 *      CHECK 10b has no honest cheap form here and is skipped - see its note.
 *   'asset'   - exhaustive over the FULL disc, with CHECK 10b alongside it.
 *      Once STAGE A has written the asset the detector is no longer being
 *      re-evaluated, so this is an ordinary raster read: cheap, exact, and the
 *      one to base a decision on.
 *
 * WHY THAT ORDERING IS SOUND, and why the old "read CHECK 10 before exporting"
 * instruction was circular: STAGE A exports BOTH bands, persist AND frac.
 * persistenceThreshold is applied only inside buildPersistence(), so frac is
 * the raw per-pixel detection frequency, untouched by the threshold. Choosing a
 * different threshold is therefore a re-read of frac, NOT a re-export. Only a
 * change to the DETECTOR - k, the scale, the polarisation - invalidates the
 * asset. The threshold decision can safely be made after STAGE A, which is the
 * only point at which the evidence for it is affordable.
 */
if (CONFIG.PERSISTENCE_MODE !== 'off') {
  var CHECK10_AOI = SITE_POINT.buffer(CONFIG.CHECK10_RADIUS_M);
  var CHECK10_SPOT_AOI = SITE_POINT.buffer(CONFIG.CHECK10_SPOT_RADIUS_M);
  var c10km = (CONFIG.CHECK10_RADIUS_M / 1000).toFixed(1);
  var cSpotKm = (CONFIG.CHECK10_SPOT_RADIUS_M / 1000).toFixed(1);
  var fromAsset = (CONFIG.PERSISTENCE_MODE === 'asset');
  var emptyNote = ' An orbit listed in CHECK 3b is EMPTY - its histogram is ' +
                  'all zeros and says nothing about persistence.';
  for (i = 0; i < CONFIG.RELATIVE_ORBITS.length; i++) {
    var roD = CONFIG.RELATIVE_ORBITS[i];

    if (fromAsset) {
      var fracD = PERSIST_FRAC[roD].updateMask(WATER_MASK);
      print('CHECK 10 - persistence fraction histogram, orbit ' + roD + ' at ' +
            PARAMS.analysisScaleM + ' m, EXHAUSTIVE over the full ' +
            (MAX_RADIUS_M / 1000).toFixed(1) + ' km disc (want BIMODAL; a spike ' +
            'near 1.0 is fixed objects).' + emptyNote,
        fracD.reduceRegion({
          reducer: ee.Reducer.histogram(20, 0.05),
          geometry: AOIS[AOIS.length - 1],
          scale: PARAMS.analysisScaleM,
          maxPixels: PARAMS.maxPixels,
          tileScale: PARAMS.tileScale
        }));

      print('CHECK 10b - water area masked as persistent, orbit ' + roD + ' at ' +
            PARAMS.analysisScaleM + ' m over a ' + c10km + ' km disc (m2)',
        ee.Image.pixelArea().updateMask(PERSIST[roD]).updateMask(WATER_MASK)
          .reduceRegion({
            reducer: ee.Reducer.sum(),
            geometry: CHECK10_AOI,
            scale: PARAMS.analysisScaleM,
            maxPixels: PARAMS.maxPixels,
            tileScale: PARAMS.tileScale
          }));
    } else {
      // BOUNDED ON BOTH AXES - region AND scene count. The previous version
      // sampled 10k pixels over a 10 km disc and still failed with "User
      // memory limit exceeded", because sample() bounds the OUTPUT while the
      // image is still evaluated tile by tile across the whole region. So the
      // spot check gets its own small disc and its own scene cap, and a
      // persistence image built from that capped subset. Reaching for the
      // all-scenes frac image here would drag that whole graph back in and
      // silently undo the cap, so this branch does not name it at all - and
      // CI asserts the branch stays clear of it.
      var spotFrac = buildPersistence(roD, PARAMS, CONFIG.CHECK10_MAX_SCENES)
                       .select('frac').updateMask(WATER_MASK);
      // Bin on the IMAGE before sampling: aggregate_histogram over a continuous
      // float would return one bucket per distinct value, which is not a
      // histogram. 0.05 bins match the exhaustive reducer above, so the two
      // modes are directly comparable.
      print('CHECK 10 - persistence fraction histogram, orbit ' + roD + ' at ' +
            PARAMS.analysisScaleM + ' m, RANDOM SAMPLE of ' +
            CONFIG.CHECK10_SAMPLE_N + ' px over a ' + cSpotKm + ' km disc, from ' +
            'the FIRST ' + CONFIG.CHECK10_MAX_SCENES + ' scenes of the orbit - ' +
            'a SPOT CHECK, not the study-area histogram (want BIMODAL; a spike ' +
            'near 1.0 is fixed objects).' + emptyNote,
        spotFrac.multiply(20).floor().divide(20).rename('frac_bin').sample({
          region: CHECK10_SPOT_AOI,
          scale: PARAMS.analysisScaleM,
          numPixels: CONFIG.CHECK10_SAMPLE_N,
          seed: 42,
          dropNulls: true,
          tileScale: PARAMS.tileScale
        }).aggregate_histogram('frac_bin'));

      // CHECK 10b IS SKIPPED IN 'compute' MODE, and skipping it is the honest
      // option rather than the cheap one. The masked AREA is a property of the
      // mask the export will actually apply, which is built from ALL the
      // orbit's scenes; measuring it off the capped subset would put a number
      // on screen that no export ever uses, and measuring it off the full mask
      // is the reduction that returned "Earth Engine memory capacity exceeded"
      // on a live run. It costs nothing to wait: in 'asset' mode it is an
      // ordinary raster read, and STAGE A does not depend on it.
      print('CHECK 10b - SKIPPED for orbit ' + roD + '. PERSISTENCE_MODE is ' +
            '"compute", so the masked-area reduction would re-evaluate the ' +
            'detector over every scene in the orbit - that is the reduction ' +
            'that exceeded the memory limit. Run STAGE A, switch ' +
            'PERSISTENCE_MODE to "asset", and CHECK 10b becomes a raster read.');
    }
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

  var previewImg = ee.Image(s1Scenes.first());
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
 *   1. Set PERSISTENCE_MODE = 'compute'. Read CHECK 10's SAMPLED spot check. If
 *      it already looks clearly unimodal, stop and reconsider before spending a
 *      STAGE A run. Run ONLY the STAGE A exports - and run the task for EVERY
 *      orbit in RELATIVE_ORBITS, including one CHECK 3b reports as empty. See
 *      the note below; this is not the optional step it looks like.
 *   2. Point PERSISTENCE_ASSET_PREFIX at the written assets and set
 *      PERSISTENCE_MODE = 'asset'. Read CHECK 10 again - now EXHAUSTIVE over
 *      the full disc, and this is the histogram to base the threshold on.
 *   3. Run the per-scene and detection exports.
 *
 * WHY AN EMPTY ORBIT STILL NEEDS ITS STAGE A TASK. The loader above walks
 * CONFIG.RELATIVE_ORBITS, not the orbits that have scenes, and builds
 * ee.Image(PERSISTENCE_ASSET_PREFIX + orbit + '_' + PARAM_SET_ID) for each.
 * CHECK 10 then reduces every one. ee.Image() on an asset that was never
 * written is a LAZY node: it builds fine and throws "Image.load: asset not
 * found" on EVALUATION - the same failure class as the ee.Image-on-an-
 * ImageCollection bug in the LAND MASKING section, and it surfaces the same
 * way, with the early CHECKs printing clean first.
 * An empty orbit's task is cheap: no scenes means no detector runs, so it is
 * a 50 M-pixel write of zeros and nothing else. There is nothing to save by
 * skipping it, and a run of step 2 to lose by doing so.
 *
 * STEP 1 USED TO SAY "read CHECK 10 first" AS THOUGH THE FULL HISTOGRAM WERE
 * AVAILABLE THERE. It is not, and that instruction was circular. Live runs
 * returned "User memory limit exceeded" at 40 km, then at 10 km, and then
 * again at 10 km with sampling switched on - because an exhaustive reduction
 * in 'compute' mode is ~1.2e12 kernel operations and sampling bounds the
 * output rather than the input. The spot check in step 1 now caps the region
 * AND the scene count; see the CHECK 10 block.
 *
 * The circle breaks on a fact about the export rather than about the compute
 * budget: STAGE A writes BOTH bands, persist AND frac. persistenceThreshold is
 * applied only inside buildPersistence(), so frac is the raw per-pixel
 * detection frequency with the threshold nowhere in it. Choosing a different
 * threshold is a RE-READ of frac, not a re-export. Only a change to the
 * DETECTOR - k, analysisScaleM, detectionPols - invalidates the asset.
 *
 * So committing to STAGE A before the threshold is settled costs nothing that
 * cannot be recovered, and it is the only point at which the evidence for the
 * threshold becomes affordable. The sampled spot check in step 1 exists so the
 * step is not taken blind.
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
      // persist is BOOLEAN and frac is a fraction. Earth Engine's default
      // pyramiding policy is MEAN, which averages a 0/1 mask into a fraction at
      // every coarser overview: a 30 m cell holding one persistent 10 m pixel
      // reads 0.11, and any read above native scale silently returns something
      // that is neither the mask nor an honest summary of it. mode preserves a
      // boolean; frac is a mean by construction so mean is right for it.
      // Set BEFORE the assets are written - changing it afterwards means
      // re-running STAGE A.
      pyramidingPolicy: { persist: 'mode', frac: 'mean' },
      maxPixels: PARAMS.maxPixels
    });
  }
}

/* The per-detection long table. One row per detected target, carrying its own
 * range from the hydrophone, so offline code derives every radius by filtering
 * range_m - and azimuth-ambiguity flagging (which needs pairwise geometry
 * between detections) becomes possible at all. */
if (CONFIG.EXPORT_DETECTIONS) {
  var detFeatures = ee.FeatureCollection(s1Scenes.map(function (img) {
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
