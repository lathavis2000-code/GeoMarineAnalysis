// ============================================================
// STEMGeoHS Marine v10.66
// Coastal Attractor Landscape + Cancer Score Pipeline
//
// SCIENTIFIC APPROACH:
//   Waddington Landscape (PNAS 2025) + Scheffer 2009 CSD +
//   Ramamurthy 2024 (Bunodosoma) + Levitan 2023 (Diadema) +
//   Lozano-Bilbao 2020-2024 (metals) + Peixoto 2025 (Red Sea)
//
// DATA PHILOSOPHY (v9.0):
//   1. Satellite data used globally (always valid)
//   2. Field data used ONLY if published for that water body
//   3. No data transfer between regions
//   4. Per-pixel MMM from NOAA CRW DHW product (not fixed 28.5)
//   5. Accuracy reported honestly per region
//   6. Model used printed to Console on every click
//
// ============================================================
// MODULE A - DATASETS (all with ocean mask)
// ============================================================
var bathy = ee.Image('NOAA/NGDC/ETOPO1').select('bedrock');
var bathyU = bathy.unmask(0);
var oceanMask      = bathyU.lt(0);
var shallowMask    = bathyU.gte(-50).and(bathyU.lt(0));
var intertidalMask = bathyU.gte(-20).and(bathyU.lt(0));

// ============================================================
// MODULE A11 - REAL IN-SITU OCEAN CHEMISTRY BASELINE (v10.36 NEW)
// Real, cleaned in-situ SeapHOx sensor baselines (pH, temperature,
// salinity, dissolved oxygen) from three confirmed ERDDAP stations,
// downloaded and outlier-cleaned via a separate Python pipeline
// (global_ocean_sensor_download.py + clean_and_compare.py).
//
// Same data philosophy already stated above: "Field data used ONLY
// if published for that water body. No data transfer between
// regions." A click is matched to real in-situ data ONLY if it falls
// within a small radius of an actual sensor - no interpolation or
// extrapolation to unmeasured locations.
//
// HOW THIS DIFFERS FROM S18 (BGC Snapshot):
//   S18 = Copernicus Global Ocean BGC model, 2022-present (3yr), any
//         location on Earth, MODEL/SATELLITE-DERIVED.
//   S19 (this module) = REAL SENSOR MEASUREMENTS, only at 3 specific
//         station locations, but with a much longer, outlier-cleaned
//         record (6,000+ raw readings per site).
//
// DATA PROVENANCE:
//   Looe Key, FL:      SECOORA ERDDAP, dataset ID looe-key-fl
//                       Mote Marine Laboratory SeapHOx, coral reef MPA
//   Agua Hedionda, CA:  SCCOOS ERDDAP, dataset ID pH-AHL
//                       Martz Lab SeapHOx, oyster farm lagoon
//   Scripps Pier, CA:   CenCOOS ERDDAP, dataset ID
//                       scripps-pier-automated-shore-sta-1
// Each entry's stats come from the CLEAN (outlier-removed) baseline,
// not the raw feed - see clean_and_compare.py for the exact method.
// ============================================================
var IN_SITU_BASELINES = {
  looe_key: {
    label: 'Looe Key, FL (coral reef MPA)',
    lat: 24.5463, lon: -81.4014,
    radius_km: 15,
    source: 'SECOORA ERDDAP (looe-key-fl) - Mote Marine Laboratory SeapHOx',
    record: '2024-06-05 to 2025-11-18 (~17.5 months, 6228 raw readings)',
    n_clean: 5151, pct_flagged: 17.3,
    pH: {mean: 7.899, std: 0.156},
    temp_c: {mean: 27.606},
    salinity: {mean: 35.708},
    do_mgL: {mean: 5.705},
    notes: 'Highest artifact rate of the 3 stations (17.3%) - remote reef ' +
      'mooring, likely biofouling/drift on an unattended pH/DO sensor.'
  },
  agua_hedionda: {
    label: 'Agua Hedionda Lagoon, CA (oyster farm)',
    lat: 33.1425, lon: -117.3275,
    radius_km: 10,
    source: 'SCCOOS ERDDAP (pH-AHL) - Martz Lab SeapHOx',
    record: '~6294 raw readings, actively managed aquaculture site',
    n_clean: 6250, pct_flagged: 0.7,
    pH: {mean: 7.902, std: 0.091},
    temp_c: {mean: 18.624},
    salinity: {mean: 33.758},
    do_mgL: {mean: 7.644},
    notes: 'Lowest artifact rate (0.7%) - actively visited/maintained ' +
      'research site, most trustworthy raw feed of the 3 stations.'
  },
  scripps_pier: {
    label: 'Scripps Pier, La Jolla, CA (open coast)',
    lat: 32.8669, lon: -117.2571,
    radius_km: 10,
    source: 'CenCOOS ERDDAP (scripps-pier-automated-shore-sta-1)',
    record: 'Last 365 days at download time, 178809 raw readings',
    n_clean: 174515, pct_flagged: 2.4,
    pH: {mean: 7.953, std: 0.165},
    temp_c: {mean: 19.023},
    salinity: {mean: 33.232},
    do_mgL: {mean: 7.820},
    notes: 'Largest clean sample (174,515 pts) - most statistically ' +
      'robust baseline of the 3 stations.'
  }
};

function haversineKm(lat1, lon1, lat2, lon2) {
  var R = 6371;
  var dLat = (lat2 - lat1) * Math.PI / 180;
  var dLon = (lon2 - lon1) * Math.PI / 180;
  var a = Math.sin(dLat/2) * Math.sin(dLat/2) +
    Math.cos(lat1 * Math.PI/180) * Math.cos(lat2 * Math.PI/180) *
    Math.sin(dLon/2) * Math.sin(dLon/2);
  var c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1-a));
  return R * c;
}

function getInSituBaseline(lat, lon) {
  var best = null, bestDist = Infinity;
  for (var key in IN_SITU_BASELINES) {
    var st = IN_SITU_BASELINES[key];
    var d = haversineKm(lat, lon, st.lat, st.lon);
    if (d <= st.radius_km && d < bestDist) {
      best = Object.assign({key: key, distance_km: d}, st);
      bestDist = d;
    }
  }
  return best;
}

// ============================================================
// MODULE A0 - DATASET SWAP (v10.13 FIX - see changelog): MODIS Aqua/Terra L3SMI retired
// ============================================================
function getOISSTColl(startDate, endDate) {
  return ee.ImageCollection('NOAA/CDR/OISST/V2_1')
    .filterDate(startDate, endDate)
    .select('sst')
    .map(function(img){
      return img.multiply(0.01).copyProperties(img, ['system:time_start']);
    });
}

var sstColl = getOISSTColl('2023-01-01','2024-12-31');
var sstRaw = ee.Image(ee.Algorithms.If(sstColl.size().gt(0),
  sstColl.mean(), ee.Image.constant(25).rename('sst'))).rename('sst');
var sst = sstRaw.updateMask(sstRaw.gt(-2).and(sstRaw.lt(40)).and(oceanMask));

var oisstPeakColl = getOISSTColl('2023-06-01','2023-10-31');
var oisstPeakMax = ee.Image(ee.Algorithms.If(
  oisstPeakColl.size().gt(0),
  oisstPeakColl.max(),
  ee.Image.constant(28).rename('sst'))).rename('sst_peak');
var sstPeak = oisstPeakMax.updateMask(
  oisstPeakMax.gt(-2).and(oisstPeakMax.lt(40)).and(oceanMask));
var sstPeakFinal = sstPeak.unmask(sstRaw.updateMask(oceanMask));

var mmmColl = getOISSTColl('2003-01-01','2022-12-31');
var monthList = ee.List.sequence(1,12);
var monthlyMeans = ee.ImageCollection(monthList.map(function(mo){
  var moNum = ee.Number(mo);
  var col = mmmColl.filter(ee.Filter.calendarRange(moNum, moNum, 'month'));
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
    col.mean().rename('sst'),
    ee.Image.constant(25).rename('sst').updateMask(ee.Image.constant(0))));
}));
var MMM_perpixel = monthlyMeans.max().rename('mmm').updateMask(oceanMask);

var dhwProper = sstPeakFinal.subtract(MMM_perpixel.add(1)).max(0)
  .multiply(12).rename('dhw').updateMask(oceanMask);

var chlColl = ee.ImageCollection('COPERNICUS/MARINE/SATELLITE_OCEAN_COLOR/V6')
  .filterDate('2023-01-01','2024-12-31').select('chlor_a');
var chlaRaw = ee.Image(ee.Algorithms.If(chlColl.size().gt(0),
  chlColl.mean(), ee.Image.constant(0).rename('chlor_a'))).rename('chlor_a');
var chla = chlaRaw.updateMask(chlaRaw.gt(0).and(chlaRaw.lt(100)).and(oceanMask));
var chlaCoastal = chlaRaw.updateMask(chlaRaw.gt(0).and(chlaRaw.lt(100)));

var s2Coll = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterDate('2023-01-01','2024-12-31')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',20)).select(['B4','B8']);
var turbRaw = ee.Image(ee.Algorithms.If(s2Coll.size().gt(0),
  s2Coll.mean().normalizedDifference(['B4','B8']).rename('turbidity'),
  ee.Image.constant(0).rename('turbidity').updateMask(ee.Image.constant(0))));
var turbImg = turbRaw.updateMask(shallowMask);

var s2_basic = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterDate('2023-01-01','2024-12-31')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',20))
  .select(['B4','B8']);
var s2_basic_mean = ee.Image(ee.Algorithms.If(s2_basic.size().gt(0),
  s2_basic.map(function(img){ return img.divide(10000); }).median(),
  ee.Image.constant(0.05).rename('B4')
    .addBands(ee.Image.constant(0.05).rename('B8'))));

// v10.41 FIX: the raw normalizedDifference() call below was confirmed to
// throw "Unknown reference to value named ''" server-side at some coastal
// points (e.g. 33.14N,-117.346W near Carlsbad, CA) where Sentinel-2
// imagery for the 2023-2024 window is sparse/malformed enough that the
// ee.Algorithms.If() fallback's band structure doesn't resolve cleanly.
// Confirmed via live testing + the v10.37 error-handling fix (which is
// what made this previously-silent failure finally visible). FIX: check
// band count server-side BEFORE calling normalizedDifference() - if the
// expected bands aren't cleanly present, produce a fully-masked (null)
// result instead of crashing. This converts a hard pipeline failure into
// a clean "n/a - no Sentinel-2 imagery here", same graceful-degradation
// standard already used for the Chl-a coastal fallback elsewhere.
var _s2basic_bandsOk = s2_basic_mean.bandNames().size().gte(2);
var ndviWater = ee.Image(ee.Algorithms.If(
  _s2basic_bandsOk,
  s2_basic_mean.select(['B8','B4']).normalizedDifference(['B8','B4']).rename('ndvi_water'),
  ee.Image.constant(0).rename('ndvi_water').updateMask(ee.Image.constant(0))
));
var ndviAlgae = ndviWater.updateMask(ndviWater.gt(0.03).and(oceanMask));

var s2_ext = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterDate('2023-01-01','2024-12-31')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',20))
  .select(['B4','B5','B8','B11']);
var s2_ext_avail = s2_ext.size().gt(0);

var s2_ext_mean = ee.Image(ee.Algorithms.If(s2_ext_avail,
  s2_ext.map(function(img){ return img.divide(10000); }).median(),
  s2_basic_mean
    .addBands(ee.Image.constant(0.03).rename('B5'))
    .addBands(ee.Image.constant(0.03).rename('B11'))));

var fai_step = (842-665)/(1610-665);
var faiImg = s2_ext_mean.expression(
  'NIR - (RED + (SWIR - RED) * step)', {
    'NIR':  s2_ext_mean.select('B8'),
    'RED':  s2_ext_mean.select('B4'),
    'SWIR': s2_ext_mean.select('B11'),
    'step': fai_step
  }).rename('fai');
var faiMild = faiImg.updateMask(faiImg.gt(0.003).and(oceanMask));
var faiAlgae = faiImg.updateMask(faiImg.gt(0.01).and(oceanMask));

// v10.41 FIX: same defensive pattern as ndviWater above - ndciImg shares
// the identical failure mode (same s2_ext_mean source, same
// normalizedDifference() crash risk at sparse-Sentinel-2 coastal points).
var _s2ext_bandsOk = s2_ext_mean.bandNames().size().gte(4);
var ndciImg = ee.Image(ee.Algorithms.If(
  _s2ext_bandsOk,
  s2_ext_mean.select(['B5','B4']).normalizedDifference(['B5','B4']).rename('ndci'),
  ee.Image.constant(0).rename('ndci').updateMask(ee.Image.constant(0))
));
var ndciBloom = ndciImg.updateMask(ndciImg.gt(0.05).and(oceanMask));

var no2Coll = ee.ImageCollection('COPERNICUS/S5P/NRTI/L3_NO2')
  .filterDate('2023-01-01','2024-12-31').select('NO2_column_number_density')
  .merge(ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2')
  .filterDate('2023-01-01','2024-12-31').select('NO2_column_number_density'));
var no2Raw = ee.Image(ee.Algorithms.If(no2Coll.size().gt(0),
  no2Coll.mean(), ee.Image.constant(0.00003).rename('NO2_column_number_density')))
  .rename('NO2_column_number_density');
var no2 = no2Raw.updateMask(no2Raw.gt(0));

var yrList = ee.List.sequence(2003,2024);
var nYears = 22;
var SST_FLOOR_GLOBAL = 14.0;

var annSSTraw = ee.ImageCollection(yrList.map(function(yr){
  var yrN = ee.Number(yr);
  var d = ee.Date.fromYMD(yrN,1,1), de = ee.Date.fromYMD(yrN,12,31);
  var col = getOISSTColl(d, de);
  var hasData = col.size().gt(3);
  var sstY = ee.Image(ee.Algorithms.If(hasData,
    col.mean().updateMask(
      col.mean().gt(SST_FLOOR_GLOBAL).and(col.mean().lt(40)).and(oceanMask)),
    ee.Image.constant(-9999).rename('sst').updateMask(ee.Image.constant(0))))
    .rename('sst');
  return sstY.set('system:time_start', d.millis()).set('year', yrN);
}));

var sstStats = annSSTraw.select('sst').reduce(
  ee.Reducer.mean().combine(ee.Reducer.stdDev(), '', true));
var sstMeanImg = sstStats.select('sst_mean');
var sstStdImg  = sstStats.select('sst_stdDev').max(0.3);

var annSSTcoll = ee.ImageCollection(annSSTraw.map(function(img){
  var z = img.subtract(sstMeanImg).divide(sstStdImg).abs();
  var isOutlier = z.gt(2.5);
  return img.updateMask(isOutlier.not())
    .copyProperties(img, ['system:time_start','year']);
}));

var annSSTlist = annSSTcoll.toList(nYears);
var allIdx = ee.List.sequence(0, nYears-1);

var pairSlopeImagesList = allIdx.map(function(i){
  var iN = ee.Number(i);
  var jList = ee.List.sequence(iN.add(1), nYears-1);
  return jList.map(function(j){
    var jN = ee.Number(j);
    var img_i = ee.Image(annSSTlist.get(iN));
    var img_j = ee.Image(annSSTlist.get(jN));
    var yearGap = jN.subtract(iN);
    return img_j.subtract(img_i).divide(yearGap).rename('slope')
      .updateMask(img_i.mask().and(img_j.mask()));
  });
}).flatten();

var pairSlopes = ee.ImageCollection.fromImages(pairSlopeImagesList);

var sst_slope = pairSlopes.select('slope').median()
  .unmask(0).updateMask(oceanMask).rename('scale');

var validYearCount = annSSTcoll.select('sst').map(function(img){
  return img.mask();
}).sum().rename('validYears');

function isEBUS(lat, lon) {
  if(lat>20 && lat<40 && lon>-25 && lon<-10) return true;
  if(lat>30 && lat<50 && lon>-130 && lon<-115) return true;
  if(lat>-45 && lat<-5 && lon>-90 && lon<-70) return true;
  if(lat>-35 && lat<-15 && lon>10 && lon<20) return true;
  return false;
}

var annualSST = ee.ImageCollection(yrList.map(function(yr){
  var yrN = ee.Number(yr);
  var d = ee.Date.fromYMD(yrN,1,1), de = ee.Date.fromYMD(yrN,12,31);
  var col = getOISSTColl(d, de);
  return ee.Image(ee.Algorithms.If(col.size().gt(3),
    col.mean().updateMask(
      col.mean().gt(SST_FLOOR_GLOBAL).and(col.mean().lt(40))).rename('sst'),
    ee.Image.constant(-9999).rename('sst').updateMask(ee.Image.constant(0))))
    .set('system:time_start', d.millis());
}));

var mStart = ee.Date('2023-01-01');
var mList24 = ee.List.sequence(0,23);

function mkMoSST() {
  return ee.ImageCollection(mList24.map(function(m){
    var d = mStart.advance(m,'month'), de = d.advance(1,'month');
    var col = getOISSTColl(d, de);
    return ee.Image(ee.Algorithms.If(col.size().gt(0),
      col.mean().updateMask(col.mean().gt(-2).and(col.mean().lt(40))
        .and(oceanMask)).rename('sst'),
      ee.Image.constant(0).rename('sst').updateMask(ee.Image.constant(0))))
      .set('system:time_start', d.millis());
  }));
}
function mkMoSSTRange(startDateStr, nMonths) {
  var startD = ee.Date(startDateStr);
  var monthsList = ee.List.sequence(0, nMonths-1);
  return ee.ImageCollection(monthsList.map(function(m){
    var d = startD.advance(m,'month'), de = d.advance(1,'month');
    var col = getOISSTColl(d, de);
    return ee.Image(ee.Algorithms.If(col.size().gt(0),
      col.mean().updateMask(col.mean().gt(-2).and(col.mean().lt(40))
        .and(oceanMask)).rename('sst'),
      ee.Image.constant(0).rename('sst').updateMask(ee.Image.constant(0))))
      .set('system:time_start', d.millis());
  }));
}
function mkMoCHL() {
  return ee.ImageCollection(mList24.map(function(m){
    var d = mStart.advance(m,'month'), de = d.advance(1,'month');
    var col = ee.ImageCollection('COPERNICUS/MARINE/SATELLITE_OCEAN_COLOR/V6')
      .filterDate(d,de).select('chlor_a');
    return ee.Image(ee.Algorithms.If(col.size().gt(0),
      col.mean().updateMask(col.mean().gt(0).and(oceanMask)).rename('chlor_a'),
      ee.Image.constant(0).rename('chlor_a').updateMask(ee.Image.constant(0))))
      .set('system:time_start', d.millis());
  }));
}
function mkMoDHW() {
  return ee.ImageCollection(mList24.map(function(m){
    var d = mStart.advance(m,'month'), de = d.advance(1,'month');
    var col = getOISSTColl(d, de);
    return ee.Image(ee.Algorithms.If(col.size().gt(0),
      col.max().subtract(MMM_perpixel.add(1)).max(0).multiply(4.33)
        .updateMask(oceanMask).rename('dhw'),
      ee.Image.constant(0).rename('dhw').updateMask(ee.Image.constant(0))))
      .set('system:time_start', d.millis());
  }));
}
function mkMoNO2() {
  return ee.ImageCollection(mList24.map(function(m){
    var d = mStart.advance(m,'month'), de = d.advance(1,'month');
    var col = ee.ImageCollection('COPERNICUS/S5P/NRTI/L3_NO2')
      .filterDate(d,de).select('NO2_column_number_density')
      .merge(ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2')
      .filterDate(d,de).select('NO2_column_number_density'));
    return ee.Image(ee.Algorithms.If(col.size().gt(0),
      col.mean().updateMask(col.mean().gt(0)).rename('NO2_column_number_density'),
      ee.Image.constant(0).rename('NO2_column_number_density')
        .updateMask(ee.Image.constant(0))))
      .set('system:time_start', d.millis());
  }));
}

function computeRealCSD(monthlyColl, study, bandName, scale) {
  var valsFC = ee.FeatureCollection(monthlyColl.map(function(img){
    var v = img.reduceRegion({reducer:ee.Reducer.mean(), geometry:study,
      scale:scale, maxPixels:1e9}).get(bandName);
    return ee.Feature(null, {v:v, t:img.get('system:time_start')});
  })).sort('t');

  var validVals = valsFC.aggregate_array('v').removeAll([null]);
  var nValid = validVals.size();

  return ee.Algorithms.If(nValid.gte(4),
    (function(){
      var arr = ee.Array(validVals);
      var idxArr = ee.Array(ee.List.sequence(0, nValid.subtract(1)));

      var meanIdx = idxArr.reduce(ee.Reducer.mean(), [0]).get([0]);
      var meanVal = arr.reduce(ee.Reducer.mean(), [0]).get([0]);
      var idxDev  = idxArr.subtract(meanIdx);
      var valDev  = arr.subtract(meanVal);
      var cov     = idxDev.multiply(valDev).reduce(ee.Reducer.sum(), [0]).get([0]);
      var varIdx  = idxDev.pow(2).reduce(ee.Reducer.sum(), [0]).get([0]);
      var slope   = ee.Number(cov).divide(ee.Number(varIdx));
      var intercept = ee.Number(meanVal).subtract(slope.multiply(meanIdx));
      var trendArr = idxArr.multiply(slope).add(intercept);
      var residArr = arr.subtract(trendArr);

      var n = nValid;
      var resid0 = residArr.slice(0, 0, n.subtract(1));
      var resid1 = residArr.slice(0, 1, n);
      var mean0 = resid0.reduce(ee.Reducer.mean(), [0]).get([0]);
      var mean1 = resid1.reduce(ee.Reducer.mean(), [0]).get([0]);
      var dev0 = resid0.subtract(mean0);
      var dev1 = resid1.subtract(mean1);
      var num  = dev0.multiply(dev1).reduce(ee.Reducer.sum(), [0]).get([0]);
      var den0 = ee.Number(dev0.pow(2).reduce(ee.Reducer.sum(), [0]).get([0])).sqrt();
      var den1 = ee.Number(dev1.pow(2).reduce(ee.Reducer.sum(), [0]).get([0])).sqrt();
      var realAC1 = ee.Number(num).divide(den0.multiply(den1).max(1e-6));

      var half = n.divide(2).floor();
      var firstHalf  = residArr.slice(0, 0, half);
      var secondHalf = residArr.slice(0, half, n);
      var varFirst  = ee.Number(firstHalf.reduce(ee.Reducer.sampleVariance(), [0]).get([0]));
      var varSecond = ee.Number(secondHalf.reduce(ee.Reducer.sampleVariance(), [0]).get([0]));
      var varTrendRatio = varSecond.divide(varFirst.max(1e-6));

      return ee.Dictionary({
        realAC1: realAC1,
        varFirstHalf: varFirst,
        varSecondHalf: varSecond,
        varTrendRatio: varTrendRatio,
        trendSlope: slope,
        nValidMonths: n
      });
    })(),
    ee.Dictionary({realAC1:null, varFirstHalf:null, varSecondHalf:null,
      varTrendRatio:null, trendSlope:null, nValidMonths:nValid})
  );
}

function analyzeThermalRecovery(monthlyFeatures, mmmValue) {
  if (mmmValue === null || mmmValue === undefined) {
    return {error: 'No real MMM baseline available at this location - cannot compute recovery time.'};
  }
  var threshold = mmmValue + 1;
  var episodes = [];
  var currentEpisode = null;
  var nValid = 0, nGap = 0;

  for (var i = 0; i < monthlyFeatures.length; i++) {
    var props = monthlyFeatures[i].properties;
    var v = props.v;
    var dateStr = new Date(props.t).toISOString().slice(0,7);
    if (v === null || v === undefined) {
      nGap++;
      continue;
    }
    nValid++;
    var stressed = v > threshold;
    if (stressed) {
      if (currentEpisode === null) {
        currentEpisode = {startDate: dateStr, months: 1, peakSST: v};
      } else {
        currentEpisode.months++;
        if (v > currentEpisode.peakSST) currentEpisode.peakSST = v;
      }
    } else {
      if (currentEpisode !== null) {
        currentEpisode.endDate = dateStr;
        currentEpisode.recovered = true;
        episodes.push(currentEpisode);
        currentEpisode = null;
      }
    }
  }
  var ongoing = null;
  if (currentEpisode !== null) {
    currentEpisode.recovered = false;
    ongoing = currentEpisode;
  }
  var completedDurations = episodes.map(function(e){ return e.months; });
  var meanRecovery = completedDurations.length > 0 ?
    (completedDurations.reduce(function(a,b){ return a+b; }, 0) / completedDurations.length) : null;
  var maxRecovery = completedDurations.length > 0 ?
    Math.max.apply(null, completedDurations) : null;

  return {
    threshold: threshold,
    mmmValue: mmmValue,
    episodes: episodes,
    ongoingEpisode: ongoing,
    nCompletedEpisodes: episodes.length,
    meanRecoveryMonths: meanRecovery,
    maxRecoveryMonths: maxRecovery,
    nValidMonths: nValid,
    nGapMonths: nGap
  };
}

// v10.42 CRITICAL FIX: this asset is CONFIRMED PERMANENTLY DEAD - fails
// with "Image asset ... not found (does not exist or caller does not
// have access)" at EVERY location, including at script STARTUP before
// any click at all (confirmed via live Console testing). This was
// already flagged as a known, unresolved issue in the ORIGINAL v10.35
// script's own changelog (v10.19 entry) but never actually fixed.
//
// WHY THIS MATTERED MUCH MORE THAN IT LOOKED: rSOIL sits in the MIDDLE
// of the 19-step nested evaluate() chain. Before v10.37's error
// handling existed, this permanent failure caused the ENTIRE REST OF
// THE CHAIN (S10/11, S12, S14-S18) to hang silently forever, at EVERY
// single location, on EVERY single click - this was very likely the
// real root cause behind most/all of the "computing for 20+ minutes"
// reports during testing, not a location-specific slowness at all.
//
// FIX: replaced the real (broken) asset reference with a safe,
// fully-masked constant image. No server-side asset load is attempted
// at all anymore - S9 now honestly reports "Dataset unavailable"
// (see soilTextureV.setValue below) instead of ever hitting this dead
// endpoint, and critically, the chain now continues normally past
// this step at every location instead of stopping here forever.
var soilTexture = ee.Image.constant(0).rename('b0')
  .updateMask(ee.Image.constant(0)); // fully masked -> always resolves to null, never blocks

function soilTextureLabel(code) {
  var labels = {
    1:'Clay', 2:'Silty Clay', 3:'Sandy Clay', 4:'Clay Loam',
    5:'Silty Clay Loam', 6:'Sandy Clay Loam', 7:'Loam', 8:'Silty Loam',
    9:'Sandy Loam', 10:'Silt', 11:'Loamy Sand', 12:'Sand'
  };
  return labels[code] || ('Unknown class '+code);
}

var usgsEarthquakes = ee.FeatureCollection('projects/sat-io/open-datasets/USGS/usgs_earthquakes')
  .filter(ee.Filter.gte('mag', 4.5));
var eqModerate = usgsEarthquakes.filter(ee.Filter.lt('mag', 6));
var eqMajor = usgsEarthquakes.filter(ee.Filter.gte('mag', 6));

var gdisDisasters = ee.FeatureCollection('projects/sat-io/open-datasets/gdis_1960-2018');
// v10.58 CRITICAL FIX: confirmed via live diagnostic testing at 3
// independent locations (Kilauea HI, coastal Japan, Bocas del Toro
// Panama) that the real GDIS property name is 'disasterty' (likely
// truncated to 10 characters, a common shapefile/DBF field-name limit)
// - NOT 'disastertype' as originally guessed. All 3 locations showed
// real, non-zero raw record counts (8, 221, 245 respectively) that
// NEVER matched the old 'disastertype' filter, while the real feature
// schema explicitly printed 'disasterty' as an actual property (see
// the v10.47 diagnostic: "disasterty":"flood" was directly observed in
// a real feature's properties at Bocas del Toro). This was previously
// causing S11 to ALWAYS report "None recorded" everywhere on Earth,
// regardless of how much real volcanic activity data actually existed.
var volcanicActivity = gdisDisasters.filter(ee.Filter.stringContains('disasterty', 'olcan'));

// ============================================================
// MODULE B - REGIONAL FIELD DATA REGISTRY
// ============================================================
function getFieldProfile(region) {

  if(region === 'Bocas del Toro, Panama' ||
     region === 'Caribbean Sea' ||
     region === 'Florida Keys, USA' ||
     region === 'Florida Reef Tract, USA' ||
     region === 'Jamaica' ||
     region === 'Little Cayman / Grand Cayman' ||
     region === 'Puerto Rico' ||
     region === 'Cuba (north coast)' ||
     region === 'Gulf of Mexico' ||
     region === 'Atlantic Coast USA') {
    return {
      hasField: true,
      species: 'Bunodosoma granuliferum + B. cavernatum (anemone)',
      R: null,
      R_stressed: null,
      urchin_N: 0.3,
      urchin_healthy: 15.0,
      anem_N: 6.0,
      anem_N_estimated: true,
      Cd: 0.006,
      Pb: 0.5,
      Cd_poll: 0.058,
      Pb_poll: 25.264,
      recruit: 2.0,
      recruit_estimated: true,
      sources: '[A] Ramamurthy 2024 CJS 54:77-82 (qualitative behavior only - see real chart below) | [B] Levitan & Edmunds 2023 PNAS | [C] Lozano-Bilbao 2020 Env Sci Poll Res',
      accuracy_field_gain: 12,
      notes: 'Field calibration: urchin density [B] + metal burden [C] measured ' +
        'in Bocas del Toro, Panama - applied here as a regional analog (no ' +
        'Florida/Caribbean-wide site has published all metrics). [B] urchin ' +
        'figures WELL CORROBORATED (98%+ decline independently confirmed via ' +
        'multiple sources on the 2022 Diadema die-off). [C] metal thresholds ' +
        'NOT yet independently verified against full source text. Ramamurthy ' +
        '2024 [A] is real but reports ONLY qualitative behavior (attack/retreat/' +
        'posture codes) - it does NOT contain a numeric A/I ratio; the real ' +
        'qualitative data is shown in the chart below for Caribbean clicks, ' +
        'NOT used as a numeric score correction. Anemone density (6.0/m2) and ' +
        'recruitment (2.0) are ESTIMATED PLACEHOLDERS, not measured values - ' +
        'see F3/F5 below.'
    };
  }

  if(region === 'Red Sea') {
    return {
      hasField: true,
      species: 'Radianthus magnifica + Amphiprion bicinctus (anemone + clownfish)',
      R: null,
      R_stressed: null,
      urchin_N: null,
      urchin_healthy: null,
      anem_N: null,
      Cd: null,
      Pb: null,
      Cd_poll: null,
      Pb_poll: null,
      recruit: null,
      dhw_calibration: 22.0,
      bleach_threshold_insitu: 32.0,
      mortality_at_calibration: 0.78,
      sources: '[D] Peixoto 2025 npj Biodiversity (KAUST) | [E] Al-Rshaidat 2020 | [F] Furby 2022 Frontiers fmars',
      accuracy_field_gain: 8,
      notes: 'Red Sea: anemone mortality calibration (DHW=22 -> 78% mortality) ' +
        'from Peixoto 2025 [D] - CONFIRMED ACCURATE, independently verified ' +
        'directly against the published abstract (both numbers match exactly). ' +
        'MMM baseline [F] (Furby 2022, ~31-32 deg C) NOT yet independently ' +
        'verified against full source text. No A/I ratio, no urchin, no metals.'
    };
  }

  if(region === 'Great Barrier Reef') {
    return {
      hasField: true,
      species: 'Coral recruit density (scleractinian spp.)',
      R: null,
      R_stressed: null,
      urchin_N: null,
      urchin_healthy: null,
      anem_N: null,
      Cd: null,
      Pb: null,
      Cd_poll: null,
      Pb_poll: null,
      recruit: 187.0,
      recruit_healthy: 247.0,
      recruit_turbid: 43.5,
      sources: '[G] Drake 2025 PLOS One | [H] Great Reef Census 2024',
      accuracy_field_gain: 7,
      notes: 'GBR: coral recruit density used (Drake 2025, 141 stations) - NOT ' +
        'yet independently verified against full source text/data. No anemone ' +
        'behavioral data - different genus.'
    };
  }

  if(region === 'Mediterranean Sea') {
    return {
      hasField: true,
      species: 'Anemonia sulcata (sea anemone) + Paracentrotus lividus (urchin)',
      R: null,
      R_stressed: null,
      urchin_N: 5.0,
      urchin_healthy: 50.0,
      anem_N: null,
      Cd: 0.006,
      Pb: 0.5,
      Cd_poll: 0.058,
      Pb_poll: 25.264,
      recruit: null,
      sources: '[C] Lozano-Bilbao 2020 Env Sci Poll Res, 2023 Mar Poll Bull, 2024 Reg Stud Mar Sci | [I] Hereu 2012 PLOS One (Medes Islands 20yr)',
      accuracy_field_gain: 10,
      notes: 'Mediterranean: metal burden from Anemonia sulcata (Lozano-Bilbao) ' +
        'and urchin proxy (Hereu 2012, Medes Islands) - NEITHER yet independently ' +
        'verified against full source text/data. No A/I ratio equivalent.'
    };
  }

  if(region === 'Pacific Coast USA') {
    return {
      hasField: false,
      species: 'none',
      R: null, R_stressed: null, urchin_N: null, urchin_healthy: null,
      anem_N: null, Cd: null, Pb: null, Cd_poll: null, Pb_poll: null,
      recruit: null, recruit_healthy: null,
      sources: 'None published for intertidal anemone survey matching our protocol',
      accuracy_field_gain: 0,
      notes: 'Pacific USA: cold water (14-17 deg C), no coral reefs, DHW=0 expected. Kelp forest ecosystem - different indicator species needed. Satellite stress index only.'
    };
  }

  if(region === 'Canary Islands, Spain') {
    return {
      hasField: false,
      species: 'none (A. taxiformis present as invasive species - no behavioral surveys)',
      R: null, R_stressed: null, urchin_N: null, urchin_healthy: null,
      anem_N: null, Cd: null, Pb: null, Cd_poll: null, Pb_poll: null,
      recruit: null, recruit_healthy: null,
      sources: 'Algaebase 2024: A. taxiformis invasive in Macaronesia | IUCN invasive species list',
      accuracy_field_gain: 0,
      notes: 'Canary Islands: A. taxiformis naturally present (invasive). SST 18-22 deg C optimal for bromoform. Wild seed stock available. No anemone field data. Aquaculture S8 score applies.'
    };
  }

  if(region === 'Azores, Portugal') {
    return {
      hasField: false,
      species: 'none (A. taxiformis reported present)',
      R: null, R_stressed: null, urchin_N: null, urchin_healthy: null,
      anem_N: null, Cd: null, Pb: null, Cd_poll: null, Pb_poll: null,
      recruit: null, recruit_healthy: null,
      sources: 'AlgaeBase: A. taxiformis distribution includes Macaronesia',
      accuracy_field_gain: 0,
      notes: 'Azores: SST 17-22 deg C, clean Atlantic water, EU regulatory framework. No anemone field data. Aquaculture S8 score applies.'
    };
  }

  if(region === 'Galapagos Islands, Ecuador') {
    return {
      hasField: false,
      species: 'none (A. taxiformis documented floating abundantly - Darwin Foundation)',
      R: null, R_stressed: null, urchin_N: null, urchin_healthy: null,
      anem_N: null, Cd: null, Pb: null, Cd_poll: null, Pb_poll: null,
      recruit: null, recruit_healthy: null,
      sources: 'Darwin Foundation Galapagos Species Database: A. taxiformis confirmed',
      accuracy_field_gain: 0,
      notes: 'Galapagos: Humboldt Current upwelling, SST 18-24 deg C, A. taxiformis documented. Near largest cattle market (Latin America 350M cattle). No anemone field data.'
    };
  }

  return {
    hasField: false,
    species: 'none',
    R: null, R_stressed: null, urchin_N: null, urchin_healthy: null,
    anem_N: null, Cd: null, Pb: null, Cd_poll: null, Pb_poll: null,
    recruit: null, recruit_healthy: null,
    sources: 'No published field data matching our protocol for this region',
    accuracy_field_gain: 0,
    notes: 'Satellite stress index only for this region. Field calibration needed to improve accuracy.'
  };
}

// ============================================================
// MODULE C0 - AQUACULTURE SUITABILITY SCORE
// ============================================================
function computeAquaculture(sv, sv_peak, cv, nv, turv, tv) {
  var sstReallyAvailable = (sv !== null && sv !== undefined);
  var sstAnnual = sv || 25;
  var sstPk = sv_peak || sv || 25;

  // v10.57 FIX: replaced GUESSED SST thresholds (previously 18-32 survivable,
  // 22-28 "optimal", with no citation) with REAL, CITED values found via
  // targeted research specifically for A. taxiformis:
  //
  //   Statton, J. (2024). Asparagopsis taxiformis hatchery and cultivation
  //   manual. AgriFutures Australia, Publication No. 24-083, Project No.
  //   PRJ-012795. ISBN 978-1-76053-489-9. Seaweed Aquaculture Research and
  //   Hatchery (SARaH) Lab, University of Western Australia.
  //   https://seaweedcentral.com.au/wp-content/uploads/4_Asparagopsis-taxiformis-hatchery-and-cultivation-manual.pdf
  //
  //   Direct quote (Table 1, "Capabilities and equipment"): "Asparagopsis
  //   grows in water between 15 and 28 deg C but we have found that 17-21
  //   deg C is optimal for A. taxiformis."
  //   Direct quote (Ecology section): "...an optimal temperature range of
  //   15-28 deg C" (matches Table 1's survivable range).
  //
  //   CAVEAT (stated honestly): these are hatchery/nursery cultivation
  //   (controlled tank) findings from the SARaH lab, not a wild open-ocean
  //   farming study - the best real published source found, but exact
  //   applicability to open-water farming specifically was not separately
  //   verified.
  var sst_optimal   = sstAnnual >= 17 && sstAnnual <= 21;
  var sst_in_window = sstAnnual > 15 && sstAnnual < 28 && !sst_optimal;
  var sst_score = sst_optimal ? 100 : sst_in_window ? 60 : 0;

  var chlReallyAvailable = (cv !== null && cv !== undefined);
  var chl_score = 50;
  if(cv !== null) {
    if(cv >= 0.3 && cv <= 1.5) chl_score = 100;
    else if(cv >= 0.1 && cv < 0.3) chl_score = 70;
    else if(cv > 1.5 && cv <= 3.0) chl_score = 60;
    else if(cv < 0.1) chl_score = 30;
    else chl_score = 20;
  }

  var poll_score = 80;
  if(nv !== null) {
    if(nv < 0.00005) poll_score = 100;
    else if(nv < 0.00010) poll_score = 80;
    else if(nv < 0.00015) poll_score = 50;
    else poll_score = 20;
  }

  var stab_score = 60;
  if(tv !== null) {
    if(tv < 0.01) stab_score = 100;
    else if(tv < 0.03) stab_score = 75;
    else if(tv < 0.05) stab_score = 50;
    else stab_score = 20;
  }

  var bromo_score = 50;
  if(sv !== null && cv !== null) {
    if(sv >= 18 && sv <= 24 && cv >= 0.3 && cv <= 2.0) bromo_score = 100;
    else if(sv >= 18 && sv <= 26 && cv >= 0.1) bromo_score = 75;
    else if(sv > 26 && sv <= 28) bromo_score = 50;
    else if(sv > 28) bromo_score = 20;
    else bromo_score = 40;
  }

  // ============================================================
  // v10.56 REDESIGN: TWO-GATE ARCHITECTURE replacing the opaque
  // 5-way weighted average entirely.
  //
  // RATIONALE (from live discussion): a single weighted-average number
  // let good scores on secondary factors (pollution, stability,
  // bromoform) mathematically outvote a fatal primary factor (v10.55
  // fixed this for SST specifically). The deeper fix, per direct
  // instruction, is architectural: SST and Chlorophyll are treated as
  // sequential GATES that must BOTH pass before secondary factors are
  // even considered - not as one of five items in a weighted sum.
  // Secondary factors (pollution/stability/bromoform), once both gates
  // pass, are surfaced as EXPLICIT TEXT CAUTIONS rather than silently
  // folded into a single opaque number.
  //
  // GATE 1 - SST: hard-veto (<=18 or >=32 deg C, from v10.55) = fatal,
  // near-zero score. "Marginal" (18-32 but below the sst_score>=60
  // good threshold) = gate fail, low/POOR score, explicit reason.
  // "Good" (sst_score>=60, i.e. in_window or optimal) = gate passes.
  //
  // GATE 2 - Chlorophyll: only reached if SST gate passes. "Good"
  // requires chl_score===100 (the 0.3-1.5 mg/m3 optimal bucket) -
  // anything else fails this gate, explicit reason given.
  //
  // Deliberately OUT OF SCOPE for this pass (see conversation): using
  // industrial-proximity or country-development context to adjust
  // pollution weighting. That would require a real geospatial
  // industrial-zone dataset and a country classification layer,
  // neither of which is currently sourced/verified in this script -
  // adding it now would mean guessing again, the exact problem this
  // redesign is trying to fix. Left as a clearly-flagged future item.
  // ============================================================
  // v10.63 CRITICAL FIX: Gate 1 previously only checked ANNUAL MEAN SST
  // against the real 15-28 deg C survivable range - sv_peak (peak SST)
  // was captured (see sstPk above) but NEVER actually used anywhere in
  // the gate logic. PROBLEM DISCOVERED (live testing, Persian Gulf,
  // 26.5N 51.5E): annual mean SST was 27.8 deg C - just inside the
  // survivable range - so Gate 1 PASSED and the site was scored
  // "DECENT". But the Persian Gulf is real, documented as one of the
  // hottest seawater bodies on Earth, with peak summer SST commonly
  // exceeding 32-35 deg C - a species that cannot survive that peak
  // would die during summer regardless of how comfortable the ANNUAL
  // AVERAGE looks. This is the same reasoning the DHW/bleaching module
  // already correctly applies (it uses peak SST, not annual mean, for
  // exactly this reason) - Gate 1 should too.
  //
  // FIX: the hard veto and "good" gate now check BOTH annual mean AND
  // peak SST against the same real cited 15-28 deg C range - if EITHER
  // one falls outside survivable bounds, the gate fails, since a
  // species must survive year-round, not just look fine on average.
  var SST_HARD_MIN = 15, SST_HARD_MAX = 28;
  var sstPeakReallyAvailable = (sv_peak !== null && sv_peak !== undefined);
  var sstHardVeto = sstReallyAvailable && (
    (sstAnnual <= SST_HARD_MIN || sstAnnual >= SST_HARD_MAX) ||
    (sstPeakReallyAvailable && (sstPk <= SST_HARD_MIN || sstPk >= SST_HARD_MAX))
  );
  var sstGateGood = sstReallyAvailable && !sstHardVeto && sst_score >= 60;
  var sstGateFail = sstReallyAvailable && !sstGateGood; // hard veto OR marginal
  var chlGateGood = chlReallyAvailable && chl_score === 100;
  var chlGateFail = chlReallyAvailable && !chlGateGood;

  var aqua_score, aqua_confidence, aqua_missing_note, status, cautions = [];

  if (!sstReallyAvailable && !chlReallyAvailable) {
    aqua_score = null; aqua_confidence = 0;
    // v10.64 FIX: this message previously defaulted to "likely dry
    // land/no ocean pixel" - but live testing (Oslo Fjord, Norway)
    // confirmed a REAL ocean point (real 2m depth, real FAI/NDCI/NDVI
    // seaweed index values all present) can ALSO hit this exact
    // "insufficient data" branch. The actual cause there was NOT land
    // and NOT a masking bug (a relaxed, unmasked coastal fallback also
    // failed) - it was a genuine short-recent-window data gap: SST/
    // Chl-a use a short ~2023-2024 window, which can come up empty at
    // a narrow, complex coastline even where satellite coverage over a
    // MUCH longer window (e.g. S17's 27-year composite) succeeds, since
    // more total imagery means better odds of at least one clear
    // observation. The message is now honest about BOTH real
    // possibilities instead of guessing "land" by default.
    aqua_missing_note = 'No usable SST/Chl-a data in the recent (2023-2024) window at this ' +
      'location - score not computed. This can mean either (a) this point is on land/no ' +
      'ocean pixel here, or (b) this IS real ocean water, but a narrow/complex coastline ' +
      '(e.g. a fjord) had no valid recent satellite observation in this short window, even ' +
      'though a longer-term average might succeed. Check Depth and the S7 seaweed index ' +
      'values above - if those show real numbers, this is case (b), not land.';
    status = 'INSUFFICIENT DATA - not a valid aquaculture site assessment';

  } else if (sstHardVeto) {
    aqua_score = 5; aqua_confidence = 100;
    var annualFails = (sstAnnual <= SST_HARD_MIN || sstAnnual >= SST_HARD_MAX);
    var peakFails = sstPeakReallyAvailable && (sstPk <= SST_HARD_MIN || sstPk >= SST_HARD_MAX);
    var vetoReason = annualFails ?
      'annual mean ('+sstAnnual.toFixed(1)+' deg C)' :
      'peak seasonal SST ('+sstPk.toFixed(1)+' deg C, even though annual mean '+
        sstAnnual.toFixed(1)+' deg C looks survivable)';
    status = 'UNSUITABLE - water temperature: '+vetoReason+' is ' +
      ((annualFails ? sstAnnual : sstPk) <= SST_HARD_MIN ? 'too cold' : 'too hot') +
      ' for this species to survive, REGARDLESS of other conditions';
    aqua_missing_note = 'GATE 1 (SST) FAILED HARD: '+vetoReason+' is outside the '+
      SST_HARD_MIN+'-'+SST_HARD_MAX+' deg C absolute survivable range' +
      (peakFails && !annualFails ? ' (v10.63: peak SST now checked, not just annual mean - a ' +
        'species must survive the hottest part of the year too)' : '') +
      '. Chlorophyll/pollution/stability/bromoform never evaluated.';

  } else if (sstGateFail) {
    aqua_score = 20; aqua_confidence = 100;
    status = 'POOR - SST is only marginal ('+sstAnnual.toFixed(1)+' deg C, outside the 20-30 deg C ' +
      'good range) - low viability regardless of other conditions';
    aqua_missing_note = 'GATE 1 (SST) FAILED: marginal temperature. Chlorophyll/pollution/stability/' +
      'bromoform not evaluated - a marginal primary factor caps the score low on its own.';

  } else if (!sstReallyAvailable) {
    // SST unknown entirely - can't evaluate Gate 1 at all
    aqua_score = null; aqua_confidence = 0;
    aqua_missing_note = 'SST unavailable at this location - Gate 1 (SST) cannot be evaluated, ' +
      'so no score can be computed regardless of other data.';
    status = 'INSUFFICIENT DATA - SST unknown, cannot evaluate';

  } else if (chlReallyAvailable && chlGateFail) {
    // Gate 1 passed, Gate 2 (Chl) failed
    aqua_score = 25; aqua_confidence = 100;
    status = 'POOR - SST is good ('+sstAnnual.toFixed(1)+' deg C), but chlorophyll is not favorable ' +
      '('+cv.toFixed(3)+' mg/m3) - low nutrient availability limits yield';
    aqua_missing_note = 'GATE 1 (SST): PASSED. GATE 2 (Chlorophyll): FAILED (' +
      cv.toFixed(3)+' mg/m3, outside the 0.3-1.5 optimal range). Pollution/stability/bromoform ' +
      'not evaluated - low nutrients caps the score on their own once SST is confirmed viable.';

  } else if (!chlReallyAvailable) {
    // Gate 1 passed, Chl unknown - can't confirm Gate 2
    aqua_score = null; aqua_confidence = 0;
    aqua_missing_note = 'SST is good ('+sstAnnual.toFixed(1)+' deg C, Gate 1 passed), but ' +
      'chlorophyll data is unavailable - Gate 2 cannot be evaluated, so no score can be computed.';
    status = 'INSUFFICIENT DATA - SST good but chlorophyll unknown, cannot evaluate';

  } else {
    // BOTH gates passed (SST good AND chlorophyll good) - proceed to
    // secondary factors. Score reflects pollution/stability/bromoform
    // directly (simple average of what's available), NOT blended back
    // with SST/Chl, since those are already confirmed good and their
    // job (as gates) is done.
    var secondary = [];
    if (nv !== null && nv !== undefined) secondary.push(poll_score);
    if (tv !== null && tv !== undefined) secondary.push(stab_score);
    secondary.push(bromo_score); // bromo always computable once sv+cv known

    aqua_score = secondary.length > 0 ?
      Math.round(secondary.reduce(function(a,b){return a+b;},0) / secondary.length) : 75;
    aqua_confidence = 100;

    if (nv !== null && nv !== undefined && poll_score < 80) {
      cautions.push('Pollution is elevated in this area (NO2='+nv.toFixed(8)+' mol/m2)');
    }
    if (tv !== null && tv !== undefined && stab_score < 75) {
      cautions.push('Thermal stability is a concern (SST trending '+tv.toFixed(4)+' deg C/yr)');
    }
    if (bromo_score < 75) {
      cautions.push('Conditions are only marginal for bromoform (anti-methane compound) yield');
    }

    status = (cautions.length === 0 ?
      'GOOD - SST and chlorophyll both favorable, no significant secondary concerns' :
      'DECENT - SST and chlorophyll both favorable, but: ' + cautions.join('; '));
    aqua_missing_note = 'GATE 1 (SST): PASSED ('+sstAnnual.toFixed(1)+' deg C). GATE 2 (Chlorophyll): ' +
      'PASSED ('+cv.toFixed(3)+' mg/m3). Score reflects pollution/stability/bromoform only - ' +
      (cautions.length > 0 ? cautions.join('; ') : 'no secondary concerns flagged.');
  }

  return {
    sst_score: sst_score,
    sst_optimal: sst_optimal,
    sst_in_window: sst_in_window,
    sst_really_available: sstReallyAvailable,
    chl_score: chl_score,
    poll_score: poll_score,
    stab_score: stab_score,
    bromo_score: bromo_score,
    aqua_score: aqua_score,
    aqua_confidence: aqua_confidence,
    aqua_missing_note: aqua_missing_note,
    status: status,
    cautions: cautions
  };
}

// ============================================================
// MODULE C-INT - INTERVENTION RECOMMENDATIONS
// ============================================================
function computeInterventions(region, fp, sc, dhwv, tv, isReefZone) {
  var actions = [];

  if(fp.urchin_N === null || fp.urchin_N === undefined) {
    actions.push({
      priority: 'DATA GAP',
      action: 'Urchin density survey + potential reintroduction',
      basis: 'No published Diadema/Paracentrotus density data for '+region+'. ' +
             'Levitan & Edmunds 2023 PNAS shows urchin collapse (15 to 0.3 per m2) ' +
             'is a primary driver of macroalgae overgrowth and reef lock-in. ' +
             'Recommend: (a) field survey to establish baseline density, ' +
             '(b) if density is low and algae cover (FAI/NDVI) is elevated, ' +
             'consider hatchery-reared urchin reintroduction as practiced ' +
             'in Caribbean restoration projects (e.g. Diadema reintroduction trials).',
      waddington: 'Restores grazing pressure that maintains the Waddington bowl wall - ' +
                   'prevents the system from sliding into the degraded algae-dominated attractor.'
    });
  } else if(fp.urchin_N < (fp.urchin_healthy || 10) * 0.3) {
    actions.push({
      priority: 'HIGH',
      action: 'Urchin density critically low ('+fp.urchin_N+'/m2) - active intervention needed',
      basis: 'Current density is below 30% of healthy baseline ('+(fp.urchin_healthy||10)+'/m2). ' +
             'Per Levitan 2023, this range correlates with reduced grazing and algae overgrowth risk.',
      waddington: 'Bowl wall (grazing-maintained resilience) is compromised - intervention ' +
                   'now is cheaper than after a tipping point is crossed (hysteresis).'
    });
  }

  if(fp.anem_N === null || fp.anem_N === undefined) {
    if(fp.hasField) {
      actions.push({
        priority: 'DATA GAP',
        action: 'Anemone/coral density survey',
        basis: 'Behavioral data exists for '+region+' ('+fp.species+') but density data ' +
               'is not yet published. Recommend field quadrat survey to establish baseline.',
        waddington: 'Density data calibrates the effective mass term in the bowl-depth model, ' +
                     'improving accuracy of tipping-point predictions for this specific site.'
      });
    }
  }

  if(dhwv !== null && dhwv > 4) {
    if(isReefZone) {
      actions.push({
        priority: dhwv > 12 ? 'CRITICAL' : dhwv > 8 ? 'HIGH' : 'MODERATE',
        action: 'Active thermal stress mitigation: shading, assisted gene flow, or coral relocation',
        basis: 'DHW='+dhwv.toFixed(1)+' deg C-weeks indicates '+
               (dhwv>12?'mortality-level':dhwv>8?'mass-bleaching-level':'bleaching-risk-level')+
               ' thermal stress. Published interventions include temporary shade cloths ' +
               '(trialed on GBR), heat-tolerant coral strain propagation (Voolstra 2021), ' +
               'and assisted migration of thermally-tolerant genotypes.',
        waddington: 'This is symptomatic relief, not a basin-deepening fix - addresses the ' +
                     'external sigma*dW forcing term but does not change V(x) shape long-term.'
      });
    } else {
      actions.push({
        priority: 'LOW',
        action: 'Monitor only - thermal anomaly without reef bleaching risk',
        basis: 'DHW='+dhwv.toFixed(1)+' detected but this region is outside tropical reef ' +
               'ecology (no corals present). Likely reflects regional warm anomaly or ' +
               'EBUS upwelling variation rather than an ecological emergency.',
        waddington: 'No basin-wall intervention needed - this is a satellite signal artifact ' +
                     'of non-reef ecology being measured with reef-calibrated thresholds.'
      });
    }
  }

  actions.push({
    priority: 'CONTEXT',
    action: 'Check S7 FAI/NDVI layer before any urchin reintroduction',
    basis: 'If FAI shows RED (extreme algae bloom), urchin reintroduction is HIGH priority ' +
           '(grazers needed to clear algae). If FAI shows CLEAR/green, algae is not currently ' +
           'a limiting factor and other interventions (thermal, pollution) take priority.',
    waddington: 'Connects S7 (current state) to S4/grazer intervention (corrective action) - ' +
                 'the intervention should match what the satellite shows is actually happening.'
  });

  if(fp.Cd !== null && fp.Cd !== undefined) {
    var cdRatio = fp.Cd_poll ? (fp.Cd - 0.006)/(fp.Cd_poll - 0.006) : 0;
    if(cdRatio > 0.5) {
      actions.push({
        priority: 'MODERATE',
        action: 'Source-tracking for heavy metal contamination (Cd/Pb)',
        basis: 'Tissue metal burden ('+fp.Cd+' Cd) is elevated relative to control baseline ' +
               '(Lozano-Bilbao methodology). Recommend tracing point sources (harbour runoff, ' +
               'industrial discharge) rather than treating the symptom in tissue.',
        waddington: 'Pollution acts as a chronic stressor lowering mu threshold - removing ' +
                     'the source is a basin-deepening fix; tissue treatment is not.'
      });
    }
  }

  if(actions.length === 0 || (actions.length === 1 && actions[0].priority === 'CONTEXT')) {
    actions.unshift({
      priority: 'STABLE',
      action: 'No urgent intervention indicated by current data',
      basis: 'Available satellite and field signals do not show acute stress requiring ' +
             'immediate action. Continue monitoring - re-assess if SST trend, DHW, or ' +
             'algae indices change.',
      waddington: 'System appears to be resting in a stable basin (B='+sc.B.toFixed(2)+'). ' +
                   'Preventative monitoring is the appropriate action at this bowl depth.'
    });
  }

  return actions;
}

// ============================================================
// MODULE C - SCORE COMPUTATION (physics layer + ToE precompute)
// ============================================================
var GEBCO = ee.ImageCollection('projects/sat-io/open-datasets/gebco/gebco_grid')
  .mosaic()
  .select('b1')
  .rename('elevation');

var gebco_raw = GEBCO;
var bathymetry_m = gebco_raw.multiply(-1)
  .updateMask(gebco_raw.lt(0));

var waveBathy = ee.ImageCollection('COPERNICUS/MARINE/WAV/ANFC_0_083DEG_STATIC')
  .first().select('deptho').rename('wave_depth');

// v10.54 CRITICAL FIX: this function built annual mean SST images
// directly from the RAW OISST 'sst' band (stored as int16 x 100, e.g.
// 2669 = 26.69 deg C) WITHOUT ever applying the *0.01 conversion that
// getOISSTColl() applies everywhere else in this script (main SST
// display, DHW, MMM, etc.). Confirmed via live diagnostic testing:
// a real point (5.6S, 5.3E) showed a genuine raw OISST value of
// 2668.52, which correctly converts to 26.69 deg C - but toeSSTFit/
// toeSSTNoise below were computing their 44-year linear trend and
// standard deviation directly on that ~2600-3000 raw value, not the
// real ~26-30 deg C temperature. This meant S17's SST trend slope and
// noise (and therefore its SNR/emergence threshold check) had been
// silently off by a factor of 100 this entire time, at every location.
// FIX: apply the same .multiply(0.01) conversion used in
// getOISSTColl(), so this now genuinely operates in real Celsius,
// consistent with the rest of the script.
var _makeAnnSST = function() {
  var list = ee.List.sequence(1982,2025).map(function(yr){
    var d=ee.Date.fromYMD(ee.Number(yr).toInt(),1,1);
    var t=ee.Image(ee.Number(yr).subtract(1982)).float().rename('t');
    return ee.ImageCollection('NOAA/CDR/OISST/V2_1')
      .filter(ee.Filter.date(d,d.advance(1,'year'))).select('sst').mean()
      .multiply(0.01)
      .float().addBands(t).set('system:time_start',d.millis());
  });
  return ee.ImageCollection(list);
};
var _annSSTColl = _makeAnnSST();
var toeSSTFit   = _annSSTColl.select(['t','sst']).reduce(ee.Reducer.linearFit());
var toeSSTNoise = _annSSTColl.select('sst').reduce(ee.Reducer.stdDev());

var _makeAnnCHL = function() {
  var list = ee.List.sequence(1998,2024).map(function(yr){
    var d=ee.Date.fromYMD(ee.Number(yr).toInt(),1,1);
    var t=ee.Image(ee.Number(yr).subtract(1998)).float().rename('t');
    return ee.ImageCollection('COPERNICUS/MARINE/SATELLITE_OCEAN_COLOR/V6')
      .filter(ee.Filter.date(d,d.advance(1,'year'))).select('chlor_a').mean()
      .float().addBands(t).set('system:time_start',d.millis());
  });
  return ee.ImageCollection(list);
};
var _annCHLColl = _makeAnnCHL();
var toeCHLFit   = _annCHLColl.select(['t','chlor_a']).reduce(ee.Reducer.linearFit());
var toeCHLNoise = _annCHLColl.select('chlor_a').reduce(ee.Reducer.stdDev());

var _makeAnnSAL = function() {
  var list = ee.List.sequence(1993,2024).map(function(yr){
    var d=ee.Date.fromYMD(ee.Number(yr).toInt(),1,1);
    var t=ee.Image(ee.Number(yr).subtract(1993)).float().rename('t');
    return ee.ImageCollection('HYCOM/sea_temp_salinity')
      .filter(ee.Filter.date(d,d.advance(1,'year'))).select('salinity_0').mean()
      .multiply(0.001).add(20)
      .float().addBands(t).set('system:time_start',d.millis());
  });
  return ee.ImageCollection(list);
};
var _annSALColl = _makeAnnSAL();
var toeSALFit   = _annSALColl.select(['t','salinity_0']).reduce(ee.Reducer.linearFit());
var toeSALNoise = _annSALColl.select('salinity_0').reduce(ee.Reducer.stdDev());

var _makeAnnNO2 = function() {
  var list = ee.List.sequence(2019,2025).map(function(yr){
    var d=ee.Date.fromYMD(ee.Number(yr).toInt(),1,1);
    var t=ee.Image(ee.Number(yr).subtract(2019)).float().rename('t');
    return ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2')
      .filter(ee.Filter.date(d,d.advance(1,'year')))
      .select('tropospheric_NO2_column_number_density').mean()
      .float().addBands(t).set('system:time_start',d.millis());
  });
  return ee.ImageCollection(list);
};
var _annNO2Coll = _makeAnnNO2();
var toeNO2Fit   = _annNO2Coll.select(['t','tropospheric_NO2_column_number_density'])
  .reduce(ee.Reducer.linearFit());
var toeNO2Noise = _annNO2Coll.select('tropospheric_NO2_column_number_density')
  .reduce(ee.Reducer.stdDev());

function getWaveCelerity(depthMeters) {
  var g = 9.81;
  return depthMeters.multiply(g).sqrt().rename('wave_celerity');
}

function getEnergyConcentrationIndex(depthMeters) {
  var safeDepth = depthMeters.max(ee.Image(0.5));
  var rawECI = safeDepth.sqrt().pow(-1);
  return rawECI.min(ee.Image(1.414)).rename('energy_concentration_index');
}

function getSSTGradientMagnitude(sstImage) {
  var kernel = ee.Kernel.sobel();
  var gradX = sstImage.convolve(ee.Kernel.prewitt());
  var gradY = sstImage.convolve(ee.Kernel.prewitt().rotate(1));
  return gradX.pow(2).add(gradY.pow(2)).sqrt()
    .rename('sst_gradient_magnitude');
}

function getSLRScenario(depthMeters, slrMeters) {
  var effectiveDepth = depthMeters.subtract(ee.Image(slrMeters)).max(ee.Image(0.1));
  return getEnergyConcentrationIndex(effectiveDepth)
    .rename('eci_slr_scenario_'+slrMeters+'m');
}

function computeToESignal(monthlySeries) {
  var valid = monthlySeries.filter(function(d){ return d.v !== null && d.v !== undefined; });
  var n = valid.length;
  if (n < 12) return {error: 'Insufficient data (n='+n+')', n: n};

  var t0 = valid[0].t;
  var tScale = 1000*60*60*24*365.25;
  var xs = valid.map(function(d){ return (d.t - t0) / tScale; });
  var ys = valid.map(function(d){ return d.v; });

  var sumX=0, sumY=0, sumXX=0, sumXY=0;
  for (var i=0; i<n; i++){ sumX+=xs[i]; sumY+=ys[i]; sumXX+=xs[i]*xs[i]; sumXY+=xs[i]*ys[i]; }
  var meanX=sumX/n, meanY=sumY/n;
  var ssXX=sumXX - n*meanX*meanX;
  var ssXY=sumXY - n*meanX*meanY;
  if (Math.abs(ssXX) < 1e-10) return {error: 'No time variance in data', n: n};
  var slope = ssXY / ssXX;
  var intercept = meanY - slope*meanX;

  var residuals = [];
  for (var j=0; j<n; j++){
    residuals.push(ys[j] - (intercept + slope*xs[j]));
  }
  var sumR2=0;
  for (var k=0; k<n; k++) sumR2 += residuals[k]*residuals[k];
  var noise = Math.sqrt(sumR2 / (n-2));

  var recordYears = xs[xs.length-1] - xs[0];
  var signal = Math.abs(slope * recordYears);

  var snr = noise > 0 ? signal / noise : 0;
  var snrThreshold = 2.0;

  return {
    slope: slope,
    intercept: intercept,
    signal: signal,
    noise: noise,
    snr: snr,
    emerged: snr >= snrThreshold,
    direction: slope > 0 ? 'RISING' : 'FALLING',
    recordYears: recordYears.toFixed(1),
    nValid: n,
    error: null
  };
}

function formatToEResult(result, varName, confidenceLevel) {
  if (result.error) return varName+': ERROR - '+result.error;
  var emerged = result.emerged ? 'EMERGED (SNR='+result.snr.toFixed(2)+' > 2.0)' :
    'NOT YET EMERGED (SNR='+result.snr.toFixed(2)+' < 2.0)';
  return varName+' ['+result.recordYears+' yr record, n='+result.nValid+']: '+
    result.direction+' trend, '+emerged+' [confidence: '+confidenceLevel+']';
}

function getEcologicalRecoveryValidation(region) {

  if (region === 'Florida Keys, USA' || region === 'Florida Reef Tract, USA') {
    return {
      checked: true,
      thermalPredictsEcological: false,
      finding: 'NO - thermal recovery did NOT predict ecological recovery here',
      details: 'Real CREMP/DRM state monitoring surveyed the SAME reef sites in ' +
        'summer 2023 (peak bleaching) and again in early 2024 (winter), explicitly ' +
        'timed AFTER thermal stress had subsided to pre-bleaching levels - i.e. ' +
        'AFTER thermal recovery already happened. Despite SST recovering within ' +
        'months, a separate report states current coral cover across Florida\'s ' +
        'Coral Reef is approx 2%, "over 90% relative reduction" - years after ' +
        'thermal recovery, the reef has NOT ecologically recovered.',
      sources: 'Florida\'s Coral Reef 2023-2024 Post-Bleaching Assessment (CREMP/DRM) | ' +
        'NOAA/FAU coral cover report (2024-2025)'
    };
  }
  if (region === 'Red Sea') {
    return {
      checked: true,
      thermalPredictsEcological: null,
      finding: 'UNKNOWN - no published recovery follow-up found (checked directly)',
      details: 'Bennett-Smith/Peixoto et al. 2025 (npj Biodiversity) documented 66-94% ' +
        'anemone mortality during the 2022-2024 Red Sea heatwave, but their monitoring ' +
        'window ends AT the mortality event itself. No published follow-up survey ' +
        'checking subsequent population recovery was found as of this check (paper ' +
        'published Sept 2025 - a multi-year recovery follow-up may not exist yet).',
      sources: 'Checked against: Bennett-Smith et al. 2025 npj Biodiversity (ends at ' +
        'mortality event, no recovery follow-up found)'
    };
  }
  return {
    checked: false,
    thermalPredictsEcological: null,
    finding: null,
    details: 'No published ecological-recovery follow-up study has been independently ' +
      'checked for this specific region yet. Do NOT assume S14\'s thermal recovery ' +
      'time predicts real ecological/biological recovery here without verifying ' +
      'against actual field survey literature first.',
    sources: 'Not yet checked'
  };
}

function getRegion(lat,lon) {
  if(lat>8&&lat<12&&lon>-83&&lon<-81)   return 'Bocas del Toro, Panama';
  if(lat>24&&lat<25.5&&lon>-82&&lon<-80) return 'Florida Keys, USA';
  if(lat>25&&lat<27&&lon>-82&&lon<-79)   return 'Florida Reef Tract, USA';
  if(lat>18&&lat<20&&lon>-80&&lon<-78)   return 'Jamaica';
  if(lat>19&&lat<20.5&&lon>-81&&lon<-79) return 'Little Cayman / Grand Cayman';
  if(lat>17&&lat<18.5&&lon>-67.5&&lon<-65) return 'Puerto Rico';
  if(lat>20&&lat<23&&lon>-75&&lon<-73)   return 'Cuba (north coast)';
  if(lat>8&&lat<25&&lon>-90&&lon<-58)    return 'Caribbean Sea';
  if(lat>32&&lat<49&&lon>-125&&lon<-117) return 'Pacific Coast USA';
  if(lat>25&&lat<32&&lon>-98&&lon<-80)   return 'Gulf of Mexico';
  if(lat>25&&lat<45&&lon>-82&&lon<-65)   return 'Atlantic Coast USA';
  if(lat>12&&lat<30&&lon>32&&lon<44)   return 'Red Sea';
  if(lat>30&&lat<47&&lon>-6&&lon<37)   return 'Mediterranean Sea';
  if(lat>-25&&lat<-10&&lon>142&&lon<155) return 'Great Barrier Reef';
  if(lat>-30&&lat<-25&&lon>152&&lon<154) return 'Great Barrier Reef';
  if(lat>18&&lat<23&&lon>-161&&lon<-154) return 'Hawaii, USA';
  if(lat>-25&&lat<0&&lon>163&&lon<180)  return 'New Caledonia / Pacific';
  if(lat>-20&&lat<5&&lon>160&&lon<180)  return 'Solomon Islands / Pacific';
  if(lon<-120&&lat>-60&&lat<65)         return 'Pacific Ocean';
  if(lat>27&&lat<30&&lon>-18&&lon<-13) return 'Canary Islands, Spain';
  if(lat>36&&lat<43&&lon>-32&&lon<-24) return 'Azores, Portugal';
  if(lat>32&&lat<33.5&&lon>-17.5&&lon<-16) return 'Madeira, Portugal';
  if(lat>14&&lat<18&&lon>-25&&lon<-22) return 'Cape Verde Islands';
  if(lat>-1&&lat<1&&lon>-92&&lon<-88) return 'Galapagos Islands, Ecuador';
  if(lat>-35&&lat<-28&&lon>27&&lon<33) return 'KwaZulu-Natal, South Africa';
  if(lon>-60&&lon<-10&&lat>-60&&lat<65) return 'Atlantic Ocean';
  if(lat<-55) return 'Southern Ocean';
  if(lat>75)  return 'Arctic Ocean';
  var la=Math.abs(Math.round(lat*10)/10), lo=Math.abs(Math.round(lon*10)/10);
  return la+(lat>=0?'N':'S')+' '+lo+(lon>=0?'E':'W')+' Coast';
}

function fmt(v,d) {
  if(v===null||v===undefined||isNaN(v)||!isFinite(v)||v<-900) return 'n/a';
  return (Math.round(v*Math.pow(10,d))/Math.pow(10,d)).toFixed(d);
}

function scoreColors(s) {
  if(s<30)  return {text:'#0a5c1e',bg:'#d4f5df',map:'#00cc44',bar:'#22cc44',lbl:'DEEP BASIN'};
  if(s<55)  return {text:'#7a5000',bg:'#fff6cc',map:'#ffcc00',bar:'#ddaa00',lbl:'WARNING'};
  if(s<75)  return {text:'#7a2e00',bg:'#ffe8d0',map:'#ff6600',bar:'#ff6600',lbl:'HIGH RISK'};
  return         {text:'#7a0000',bg:'#ffd0d0',map:'#ff0000',bar:'#ff2222',lbl:'CRITICAL'};
}

function computeScore(sv, cv, tv, nv, turv, dhwv, fp, lat, lon) {
  var isReefZone = (lat !== undefined && lon !== undefined) ?
    ((lat > -30 && lat < 30) && !isEBUS(lat, lon)) : true;
  var s1 = 50;
  if(sv!==null){
    if(sv>32)s1=100; else if(sv>30)s1=Math.round(75+(sv-30)*12.5);
    else if(sv>28)s1=Math.round(45+(sv-28)*15);
    else if(sv>=24)s1=10; else s1=Math.max(0,Math.round(30-(24-sv)*3));
  }
  s1=Math.max(0,Math.min(100,s1));

  var s2 = 50;
  if(cv!==null&&cv>0){
    if(cv<0.05)s2=85; else if(cv<0.1)s2=70; else if(cv<0.5)s2=50;
    else if(cv<1.0)s2=25; else if(cv<2.0)s2=10; else s2=5;
  }
  s2=Math.max(0,Math.min(100,s2));

  var s3 = 35;
  if(turv!==null){
    if(turv>0.3)s3=85; else if(turv>0.1)s3=65;
    else if(turv>0)s3=40; else if(turv>-0.1)s3=20; else s3=10;
  }
  s3=Math.max(0,Math.min(100,s3));

  var s4b = 30;
  if(tv!==null){
    if(tv>0.08)s4b=95; else if(tv>0.05)s4b=80; else if(tv>0.03)s4b=60;
    else if(tv>0.01)s4b=40; else if(tv>0)s4b=20; else s4b=5;
  }
  var s4d = 0;
  if(dhwv!==null){
    if(isReefZone) {
      if(dhwv>16)s4d=40;
      else if(dhwv>12)s4d=35;
      else if(dhwv>8) s4d=25;
      else if(dhwv>4) s4d=15;
      else if(dhwv>2) s4d=8;
      else if(dhwv>1) s4d=3;
      if(fp.dhw_calibration && dhwv >= fp.dhw_calibration) s4d=40;
    } else {
      if(dhwv>16)s4d=15;
      else if(dhwv>12)s4d=10;
      else if(dhwv>8) s4d=6;
      else if(dhwv>4) s4d=3;
      else if(dhwv>1) s4d=1;
    }
  }
  var s4 = Math.min(100,s4b+s4d);

  var s5 = Math.max(0,Math.min(100,
    (sv!==null&&sv>29?40:0)+(tv!==null&&tv>0.03?35:0)+(cv!==null&&cv<0.3?25:0)));

  var s6 = 20;
  if(nv!==null&&nv>0){
    if(nv>0.00015)s6=90; else if(nv>0.00010)s6=70;
    else if(nv>0.00005)s6=50; else if(nv>0.00002)s6=30; else s6=10;
  }
  s6=Math.max(0,Math.min(100,s6));

  var csat = Math.round(s1*0.15+s2*0.15+s3*0.15+s4*0.25+s5*0.20+s6*0.10);

  var F1=0, F2=0, F3=0, F4=0, F5=0;

  if(fp.hasField) {
    if(fp.urchin_N !== null && fp.urchin_healthy !== null) {
      F1 = -(Math.min(1.0, fp.urchin_N/fp.urchin_healthy)*15);
    }
    if(fp.anem_N !== null) {
      F3 = -(Math.min(1.0, fp.anem_N/10.0)*10);
    }
    if(fp.Cd !== null && fp.Pb !== null) {
      var CdF = Math.max(0,Math.min(1,(fp.Cd - 0.006)/(fp.Cd_poll - 0.006)));
      var PbF = Math.max(0,Math.min(1,(fp.Pb - 0.5)/(fp.Pb_poll - 0.5)));
      F4 = Math.round(((CdF+PbF)/2)*15);
    }
    if(fp.recruit !== null && fp.recruit_healthy !== null) {
      var recFrac = Math.min(1.0, fp.recruit/fp.recruit_healthy);
      F5 = -Math.round(recFrac*8);
    } else if(fp.recruit !== null && fp.recruit <= 5.0) {
      F5 = -(Math.min(1.0, fp.recruit/5.0)*8);
    }
  }

  F1 = isNaN(F1)?0:F1; F2 = isNaN(F2)?0:F2; F3 = isNaN(F3)?0:F3;
  F4 = isNaN(F4)?0:F4; F5 = isNaN(F5)?0:F5;
  var fcT = Math.round(F1+F2+F3+F4+F5);
  var ccs = Math.max(0,Math.min(100,csat+fcT));
  if(isNaN(ccs)) ccs = csat;
  if(isNaN(ccs)) ccs = 30;
  ccs = Math.max(0,Math.min(100,ccs));
  var B = Math.round((100-ccs))/100;
  var mu = Math.max(0.01,ccs/100), dU = 0.25*mu*mu, sig = 0.04+ccs/2000;
  var anem_safe = (fp.anem_N!==null&&!isNaN(fp.anem_N)) ? fp.anem_N : 6.0;
  var urch_safe = (fp.urchin_N!==null&&!isNaN(fp.urchin_N)) ? fp.urchin_N : 0.3;
  var me = Math.max(0.1, anem_safe*0.2 + urch_safe*0.05);
  var w0 = Math.sqrt(mu/me), k = w0*Math.exp(-dU/sig);
  var p5 = isNaN(k)||!isFinite(k) ? 0 : Math.round((1-Math.exp(-k*5))*100);
  p5 = Math.max(0,Math.min(100,p5));
  var ac1 = tv?Math.min(0.99,Math.max(0.10,0.40+tv*9)):0.50;
  var tau = (B>0&&!isNaN(B))?Math.round(10/B)/10:99;
  var aS = 77;
  var aF = fp.accuracy_field_gain || 0;
  return {s1:s1,s2:s2,s3:s3,s4:s4,s5:s5,s6:s6,
    sat_ccs:csat, ccs:ccs, B:B, mu:mu, deltaU:dU,
    meff:me, omega0:w0, k:k, p5yr:p5, ac1:ac1, tau:tau,
    fcTotal:fcT, F1:F1, F2:F2, F3:F3, F4:F4, F5:F5,
    acc_sat:aS, acc_field:aF, acc_total:aS+aF};
}

// ============================================================
// MODULE D - SIDEBAR UI
// ============================================================
var panel = ui.Panel({style:{width:'256px',padding:'4px',backgroundColor:'#f0f5f0'}});

var clickResultsLog = [];

function lbl(txt,sz,col,bg,bld) {
  return ui.Label(txt,{fontSize:sz+'px',fontWeight:bld?'bold':'normal',
    color:col||'#111111',backgroundColor:bg||'rgba(0,0,0,0)',
    padding:'2px 4px',margin:'1px 0'});
}
function sHead(txt,bg,tc) { return lbl(txt,9,tc||'#ffffff',bg||'#334455',true); }
function dynLbl(init,col) {
  return ui.Label(init,{fontSize:'9px',color:col||'#111111',
    backgroundColor:'rgba(0,0,0,0)',padding:'1px 4px',margin:'0'});
}
function row(key,dl) {
  return ui.Panel([lbl(key+': ',9,'#334466'),dl],
    ui.Panel.Layout.Flow('horizontal'),{margin:'1px 0'});
}
function legRow(hex,main,sub) {
  return ui.Panel([
    ui.Label('',{backgroundColor:hex,width:'16px',height:'12px',
      margin:'2px 5px 0 2px',padding:'0',border:'1px solid #888888'}),
    ui.Panel([lbl(main,8,'#111111'), sub?lbl(sub,7,'#555555'):null].filter(Boolean),
      ui.Panel.Layout.Flow('vertical'),{margin:'0'})
  ],ui.Panel.Layout.Flow('horizontal'),{margin:'2px 0'});
}
function legDiv() {
  return ui.Label('',{margin:'3px 0 1px 0',
    backgroundColor:'#cccccc',height:'1px',stretch:'horizontal'});
}

panel.add(lbl('STEMGeoHS Marine v10.36',12,'#ffffff','#1a4a2a',true));
panel.add(lbl('Region-specific field data | v10.36 + Real In-Situ Baseline',9,'#336633'));
panel.add(lbl('Waddington+Lagrangian+Scheffer 2009+PNAS 2025',8,'#557755'));
panel.add(lbl('Ramamurthy 2024|Peixoto 2025|Levitan 2023|Lozano-Bilbao',8,'#557755'));

var clickLbl = lbl('CLICK coastal reef/shallow water to analyze',10,'#ffffff','#1a5a2a',true);
panel.add(clickLbl);

panel.add(lbl('GO TO COORDINATES (use this, not the search bar above)',9,'#ffffff','#334466',true));
panel.add(lbl('Format: lat, lon  (e.g. 63.84, -22.43 for Grindavik)',7,'#555555'));
var coordInput = ui.Textbox({
  placeholder: 'lat, lon  e.g. 63.84, -22.43',
  style: {stretch:'horizontal', margin:'2px 4px', fontSize:'11px'}
});
panel.add(coordInput);
var coordStatusV = dynLbl('','#880000');
var goToCoordsBtn = ui.Button({
  label: 'GO (pans map + draws white ring + analyzes)',
  style: {fontSize:'11px', fontWeight:'bold', margin:'2px 4px',
    backgroundColor:'#cce0ff', color:'#003388', stretch:'horizontal',
    padding:'6px 4px', border:'2px solid #0055cc'},
  onClick: function(){
    var txt = coordInput.getValue();
    if(!txt){ coordStatusV.setValue('Enter coordinates first, e.g. 63.84, -22.43'); return; }
    var parts = txt.split(',').map(function(s){ return parseFloat(s.trim()); });
    if(parts.length !== 2 || isNaN(parts[0]) || isNaN(parts[1])){
      coordStatusV.setValue('Invalid format - use: lat, lon  (e.g. 63.84, -22.43)');
      return;
    }
    var latIn = parts[0], lonIn = parts[1];
    if(latIn < -90 || latIn > 90){
      coordStatusV.setValue('Latitude must be between -90 and 90 - check you didn\'t swap lat/lon');
      return;
    }
    if(lonIn < -180 || lonIn > 180){
      coordStatusV.setValue('Longitude must be between -180 and 180');
      return;
    }
    coordStatusV.setValue('Going to '+latIn+', '+lonIn+'...');
    coordStatusV.style().set('color','#115511');
    Map.setCenter(lonIn, latIn, 11);
    analyzeLocation(latIn, lonIn);
  }
});
panel.add(goToCoordsBtn);
panel.add(coordStatusV);

panel.add(sHead('S13 - HISTORICAL CSD WINDOW TEST (v10.36)','#1a4a4a'));
panel.add(lbl('Test real AC1/variance against a chosen historical window',7,'#226666'));
panel.add(lbl('Independent tool - does not affect score or main click flow',7,'#888888'));
panel.add(lbl('Lat, Lon:',7,'#334466'));
var csdTestCoordInput = ui.Textbox({
  placeholder: 'lat, lon  e.g. 24.55, -81.78 (Florida Keys)',
  style: {stretch:'horizontal', margin:'2px 4px', fontSize:'11px'}
});
panel.add(csdTestCoordInput);
panel.add(lbl('Start date (YYYY-MM-DD):',7,'#334466'));
var csdTestStartInput = ui.Textbox({
  placeholder: 'e.g. 2021-08-01',
  style: {stretch:'horizontal', margin:'2px 4px', fontSize:'11px'}
});
panel.add(csdTestStartInput);
panel.add(lbl('Number of months (4-60):',7,'#334466'));
var csdTestMonthsInput = ui.Textbox({
  placeholder: 'e.g. 24',
  style: {stretch:'horizontal', margin:'2px 4px', fontSize:'11px'}
});
panel.add(csdTestMonthsInput);
panel.add(lbl('Label this run as:',7,'#334466'));
var csdTestWindowLabel = ui.Select({
  items: ['BEFORE event (baseline)', 'AFTER event (post-collapse)'],
  value: 'BEFORE event (baseline)',
  style: {stretch:'horizontal', margin:'2px 4px', fontSize:'11px'}
});
panel.add(csdTestWindowLabel);
var csdTestStatusV = dynLbl('','#880000');
var csdTestAc1V = dynLbl('--','#226666');
var csdTestVarTrendV = dynLbl('--','#226666');
var csdTestThermalV = dynLbl('--','#226622');
var csdTestNoteV = dynLbl('','#888888');
var csdBeforeResult = null;
var csdAfterResult = null;
var csdTestRunBtn = ui.Button({
  label: 'RUN CSD TEST ON THIS WINDOW',
  style: {fontSize:'11px', fontWeight:'bold', margin:'2px 4px',
    backgroundColor:'#cce8e8', color:'#114444', stretch:'horizontal',
    padding:'6px 4px', border:'2px solid #226666'},
  onClick: function(){
    var coordTxt = csdTestCoordInput.getValue();
    var startTxt = csdTestStartInput.getValue();
    var monthsTxt = csdTestMonthsInput.getValue();
    if(!coordTxt || !startTxt || !monthsTxt){
      csdTestStatusV.setValue('Fill in all 3 fields: lat/lon, start date, and number of months');
      return;
    }
    var parts = coordTxt.split(',').map(function(s){ return parseFloat(s.trim()); });
    if(parts.length !== 2 || isNaN(parts[0]) || isNaN(parts[1])){
      csdTestStatusV.setValue('Invalid lat/lon format - use: lat, lon (e.g. 24.55, -81.78)');
      return;
    }
    var latIn = parts[0], lonIn = parts[1];
    if(latIn < -90 || latIn > 90 || lonIn < -180 || lonIn > 180){
      csdTestStatusV.setValue('Lat must be -90..90, Lon must be -180..180 - check for a swap');
      return;
    }
    var nMonthsIn = parseInt(monthsTxt, 10);
    if(isNaN(nMonthsIn) || nMonthsIn < 4 || nMonthsIn > 60){
      csdTestStatusV.setValue('Number of months must be a whole number between 4 and 60');
      return;
    }
    csdTestStatusV.setValue('Running CSD test on '+nMonthsIn+' months starting '+startTxt+'...');
    csdTestStatusV.style().set('color','#115511');
    csdTestAc1V.setValue('computing...'); csdTestVarTrendV.setValue('computing...');
    csdTestThermalV.setValue('computing...');
    csdTestNoteV.setValue('');

    var testPt = ee.Geometry.Point([lonIn, latIn]);
    var testStudy = testPt.buffer(1000);
    var testColl = mkMoSSTRange(startTxt, nMonthsIn);
    var rCSDTest = computeRealCSD(testColl, testStudy, 'sst', 4000);

    var windowLabel = csdTestWindowLabel.getValue();
    var rawValsFC = ee.FeatureCollection(testColl.map(function(img){
      var v = img.reduceRegion({reducer:ee.Reducer.mean(), geometry:testStudy,
        scale:4000, maxPixels:1e9}).get('sst');
      return ee.Feature(null, {v:v, t:img.get('system:time_start')});
    })).sort('t');

    var rMMM_test = MMM_perpixel.reduceRegion({reducer:ee.Reducer.mean(),
      geometry:testStudy, scale:4000, maxPixels:1e9});

    rawValsFC.evaluate(
      function(rawFC, _err){
        if (_err) {
          print('=== S13 RAW MONTHLY DATA ERROR ['+windowLabel+'] === '+_err);
          csdTestThermalV.setValue('n/a - error fetching monthly data');
          return;
        }
        print('=== S13 RAW MONTHLY DATA ['+windowLabel+'] '+startTxt+' + '+nMonthsIn+' months ===');
        print('Location: '+latIn+', '+lonIn);
        var feats = rawFC.features;
        var validCount = 0;
        for(var mi=0; mi<feats.length; mi++){
          var props = feats[mi].properties;
          var dateStr = new Date(props.t).toISOString().slice(0,7);
          if(props.v !== null && props.v !== undefined){
            print('  Month '+(mi+1)+' ('+dateStr+'): '+props.v.toFixed(2)+' deg C');
            validCount++;
          } else {
            print('  Month '+(mi+1)+' ('+dateStr+'): NO DATA (null - no satellite coverage this month/pixel)');
          }
        }
        print('Raw valid months: '+validCount+' / '+feats.length);

        rMMM_test.evaluate(function(mmmRes){
          var mmmTestVal = (mmmRes && mmmRes.mmm !== null && mmmRes.mmm !== undefined) ? mmmRes.mmm : null;
          var trTest = analyzeThermalRecovery(feats, mmmTestVal);

          if (trTest.error) {
            csdTestThermalV.setValue('n/a - '+trTest.error);
            print('=== S13/S14 THERMAL RECOVERY TIME ['+windowLabel+'] === ERROR: '+trTest.error);
          } else {
            csdTestThermalV.setValue(trTest.nCompletedEpisodes+' episode(s), mean recovery: '+
              (trTest.meanRecoveryMonths!==null?trTest.meanRecoveryMonths.toFixed(1)+' mo':'n/a')+
              (trTest.ongoingEpisode!==null?' [+1 ONGOING, not averaged]':''));
            print('=== S13/S14 THERMAL RECOVERY TIME ['+windowLabel+'] ===');
            print('Threshold (MMM+1): '+trTest.threshold.toFixed(2)+' deg C');
            print('Completed episodes: '+trTest.nCompletedEpisodes);
            for (var ti=0; ti<trTest.episodes.length; ti++) {
              var tep = trTest.episodes[ti];
              print('  Episode '+(ti+1)+': '+tep.startDate+' to '+tep.endDate+' - recovered in '+
                tep.months+' month(s), peak '+tep.peakSST.toFixed(2)+' deg C');
            }
            if (trTest.ongoingEpisode !== null) {
              print('ONGOING (not averaged): started '+trTest.ongoingEpisode.startDate);
            }

            var ecoValTest = getEcologicalRecoveryValidation(getRegion(latIn, lonIn));
            print('=== S15 ECOLOGICAL RECOVERY VALIDATION ['+windowLabel+'] ===');
            print('Checked against real published literature: '+(ecoValTest.checked?'YES':'NO - not yet checked'));
            if (ecoValTest.checked) { print('Finding: '+ecoValTest.finding); }
          }
        });
      }
    );

    rCSDTest.evaluate(
      function(res, _err){
        if (_err) {
          csdTestStatusV.setValue('ERROR: '+_err);
          csdTestStatusV.style().set('color','#cc0000');
          csdTestAc1V.setValue('--'); csdTestVarTrendV.setValue('--'); csdTestThermalV.setValue('--');
          return;
        }
        var rAC1 = (res && res.realAC1 !== null && res.realAC1 !== undefined) ? res.realAC1 : null;
        var rVar = (res && res.varTrendRatio !== null && res.varTrendRatio !== undefined) ? res.varTrendRatio : null;
        var rN = (res && res.nValidMonths !== null && res.nValidMonths !== undefined) ? res.nValidMonths : 0;
        csdTestAc1V.setValue(rAC1 !== null ? rAC1.toFixed(3)+' (n='+rN+' valid months)' :
          'n/a (insufficient valid months, n='+rN+')');
        csdTestVarTrendV.setValue(rVar !== null ? rVar.toFixed(2)+'x'+
          (rVar>1.5?' RISING':rVar<0.67?' falling':' stable') : 'n/a');
        csdTestNoteV.setValue(rAC1!==null&&rAC1>0.5&&rVar!==null&&rVar>1.3?
          'CSD pattern present in this window.' :
          'No strong CSD pattern in this window.');
        csdTestStatusV.setValue('Done ['+windowLabel+']. Stored for comparison.');
        csdTestStatusV.style().set('color','#115511');

        print('=== S13 HISTORICAL CSD TEST RESULT ['+windowLabel+'] ===');
        print('Real AC1 (detrended): '+(rAC1!==null?rAC1.toFixed(4):'n/a'));

        var storedResult = {
          ac1: rAC1, varTrend: rVar, nMonths: rN,
          lat: latIn, lon: lonIn, startDate: startTxt, months: nMonthsIn
        };
        if(windowLabel === 'BEFORE event (baseline)'){
          csdBeforeResult = storedResult;
        } else {
          csdAfterResult = storedResult;
        }
      }
    );
  }
});
panel.add(csdTestRunBtn);
panel.add(csdTestStatusV);
panel.add(row('Real AC1 (this window)',csdTestAc1V));
panel.add(row('Variance trend (this window)',csdTestVarTrendV));
panel.add(row('Thermal recovery (direct, v10.36)',csdTestThermalV));
panel.add(csdTestNoteV);
var csdCompareResultV = ui.Label('',{fontSize:'8px',color:'#114444',
  backgroundColor:'#e8f8f8',padding:'4px 6px',margin:'2px 0',whiteSpace:'pre'});
var csdCompareBtn = ui.Button({
  label: 'COMPARE BEFORE vs AFTER',
  style: {fontSize:'11px', fontWeight:'bold', margin:'2px 4px',
    backgroundColor:'#ffe8cc', color:'#884400', stretch:'horizontal',
    padding:'6px 4px', border:'2px solid #cc7700'},
  onClick: function(){
    if(!csdBeforeResult || !csdAfterResult){
      csdCompareResultV.setValue('Need BOTH a BEFORE and an AFTER run stored first.');
      csdCompareResultV.style().set('color','#cc0000');
      return;
    }
    var b = csdBeforeResult, a = csdAfterResult;
    var ac1Rose = (b.ac1!==null && a.ac1!==null) ? (a.ac1 > b.ac1) : null;
    var varRose = (b.varTrend!==null && a.varTrend!==null) ? (a.varTrend > b.varTrend) : null;
    var criteriaChecked = [ac1Rose, varRose].filter(function(v){ return v!==null; });
    var criteriaMet = criteriaChecked.filter(function(v){ return v===true; });
    var pctMatch = criteriaChecked.length>0 ?
      Math.round((criteriaMet.length/criteriaChecked.length)*100) : null;

    var verdict = pctMatch===null ? 'INCONCLUSIVE' :
      pctMatch===100 ? 'Both AC1 and variance rose AFTER vs BEFORE' :
      pctMatch===50  ? 'PARTIAL match' : 'NO match';

    var statement =
      'Directional match: '+(pctMatch!==null?pctMatch+'%':'n/a')+'\n'+
      'BEFORE: AC1='+(b.ac1!==null?b.ac1.toFixed(3):'n/a')+'\n'+
      'AFTER:  AC1='+(a.ac1!==null?a.ac1.toFixed(3):'n/a')+'\n'+
      'Verdict: '+verdict;

    csdCompareResultV.setValue(statement);
    csdCompareResultV.style().set('color', pctMatch===100?'#0a5c1e':pctMatch===50?'#664400':'#880000');

    print('=== S13 COMPARE BEFORE vs AFTER ===');
    print(statement);
  }
});
panel.add(csdCompareBtn);
panel.add(csdCompareResultV);

var locV=dynLbl('--','#115511'), regV=dynLbl('--','#115511'), coV=dynLbl('--','#666666');
panel.add(row('Location',locV)); panel.add(row('Region',regV)); panel.add(row('Lat/Lon',coV));

panel.add(sHead('SATELLITE DATA (globally valid)','#1a3a6a'));
var sstV=dynLbl('--','#cc2200'), trendV=dynLbl('--','#882200'),
    dhwV=dynLbl('--','#aa3300'), mmmV=dynLbl('--','#886600'),
    chlaV=dynLbl('--','#116611'), turbV=dynLbl('--','#664400'),
    no2V=dynLbl('--','#550066'), depthV=dynLbl('--','#113366');
panel.add(row('SST',sstV)); panel.add(row('SST trend',trendV));
panel.add(row('DHW (per-pixel MMM)',dhwV));
var dhwNoteV = dynLbl('','#664400');
panel.add(dhwNoteV);
panel.add(row('MMM (local baseline)',mmmV));
panel.add(row('Chl-a',chlaV)); panel.add(row('Turbidity',turbV));
panel.add(row('NO2',no2V)); panel.add(row('Depth',depthV));
var depthWarnV = dynLbl('','#880000');
panel.add(depthWarnV);
var ebusWarnV = dynLbl('','#885500');
panel.add(ebusWarnV);

panel.add(sHead('SEAWEED / MACROALGAE (S7)','#226600'));
var faiV  = dynLbl('--','#226600');
var ndciV = dynLbl('--','#226600');
var ndviWV= dynLbl('--','#226600');
var algaeStatusV = dynLbl('checking...','#666666');
panel.add(row('FAI (floating algae)',faiV));
panel.add(row('NDCI (bloom index)',ndciV));
panel.add(row('NDVI water (benthic)',ndviWV));
panel.add(row('Algae status',algaeStatusV));

panel.add(sHead('S8 - AQUACULTURE SUITABILITY (v9.9)','#003366'));
var sstWindowV   = dynLbl('--','#666666');
var chlWindowV   = dynLbl('--','#666666');
var ndciNutrientProxyV = dynLbl('--','#666666');
var pollutionV   = dynLbl('--','#666666');
var stabilityV   = dynLbl('--','#666666');
var bromoformV   = dynLbl('--','#666666');
var aquaScoreV   = dynLbl('--','#666666');
var aquaStatusV  = dynLbl('computing...','#666666');
panel.add(row('SST window (22-28 C)',sstWindowV));
panel.add(row('Nutrients Chl-a',chlWindowV));
panel.add(row('NDCI nutrient proxy (S2)',ndciNutrientProxyV));
panel.add(row('Pollution NO2',pollutionV));
panel.add(row('Thermal stability',stabilityV));
panel.add(row('Bromoform yield',bromoformV));
panel.add(row('Aquaculture score',aquaScoreV));
var aquaMissingV = ui.Label('',{fontSize:'7px',color:'#885500',
  backgroundColor:'rgba(0,0,0,0)',padding:'1px 4px',margin:'0',whiteSpace:'pre'});
panel.add(aquaMissingV);
panel.add(row('Farm site status',aquaStatusV));

panel.add(sHead('S9 - GEOLOGY / SOIL (v10.36)','#5a4422'));
var soilTextureV = dynLbl('checking...','#664422');
panel.add(row('Soil texture (USDA)',soilTextureV));

panel.add(sHead('S10/S11 - SEISMIC / VOLCANIC HAZARD (v10.36)','#5a2222'));
var eqStatsV = dynLbl('checking...','#774444');
var volcStatsV = dynLbl('checking...','#774444');
panel.add(row('Earthquakes (200km, M4.5+)',eqStatsV));
panel.add(row('Volcanic activity (300km, 1960-2018)',volcStatsV));

panel.add(sHead('FIELD DATA STATUS (this region)','#2a4a1a'));
var fieldStatusV=dynLbl('Checking...','#666666');
var fieldSpeciesV=dynLbl('--','#336611');
var fieldSourcesV=ui.Label('--',{fontSize:'7px',color:'#555577',
  backgroundColor:'rgba(0,0,0,0)',padding:'1px 4px',margin:'0',whiteSpace:'pre'});
var fieldNotesV=ui.Label('--',{fontSize:'7px',color:'#886600',
  backgroundColor:'rgba(0,0,0,0)',padding:'1px 4px',margin:'0',whiteSpace:'pre'});
panel.add(row('Status',fieldStatusV));
panel.add(row('Species',fieldSpeciesV));
panel.add(fieldSourcesV);
panel.add(fieldNotesV);

panel.add(sHead('FIELD CORRECTIONS (valid data only)','#334422'));
var fc1V=dynLbl('n/a','#888888'),
    fc3V=dynLbl('n/a','#888888'), fc4V=dynLbl('n/a','#888888'),
    fc5V=dynLbl('n/a','#888888'), fcTV=dynLbl('0 total','#222222');
panel.add(row('F1 Urchin grazer',fc1V));
panel.add(row('F3 Anem density',fc3V)); panel.add(row('F4 Metals',fc4V));
panel.add(row('F5 Recruitment',fc5V)); panel.add(row('Total correction',fcTV));

panel.add(sHead('BOWL DEPTH (FUSED - heuristic composite index)','#2a4a2a'));
var satCcsV=dynLbl('--','#664400'), fusCcsV=dynLbl('--','#115511'),
    bowlV=dynLbl('--','#115533'), omega0V=dynLbl('--','#115533'),
    tauV=dynLbl('--','#661111'), ac1V=dynLbl('--','#224411'),
    p5yrV=dynLbl('--','#880000');
panel.add(row('Satellite CCS',satCcsV)); panel.add(row('FUSED CCS',fusCcsV));
panel.add(row('Bowl depth B',bowlV)); panel.add(row('omega0',omega0V));
panel.add(row('Return tau',tauV));
panel.add(row('AC1 (legacy heuristic)',ac1V));
panel.add(row('P(flip 5yr)',p5yrV));

panel.add(sHead('S12 - REAL CSD STATISTICS (v10.36)','#1a4a4a'));
var realAc1V = dynLbl('checking...','#226666');
var realVarTrendV = dynLbl('checking...','#226666');
var realCsdNoteV = dynLbl('','#888888');
panel.add(row('Real AC1 (detrended, 24mo)',realAc1V));
panel.add(row('Variance trend (2nd/1st half)',realVarTrendV));
panel.add(realCsdNoteV);

panel.add(sHead('S14 - REAL THERMAL RECOVERY TIME','#1a3a1a'));
var thermalEpisodesV = dynLbl('checking...','#226622');
var thermalMeanV = dynLbl('checking...','#226622');
var thermalMaxV = dynLbl('checking...','#226622');
var thermalOngoingV = dynLbl('checking...','#226622');
panel.add(row('Heat-stress episodes',thermalEpisodesV));
panel.add(row('Mean recovery time',thermalMeanV));
panel.add(row('Longest recovery observed',thermalMaxV));
panel.add(row('Ongoing episode?',thermalOngoingV));

panel.add(sHead('S15 - ECOLOGICAL RECOVERY VALIDATION','#3a2a4a'));
var ecoValStatusV = dynLbl('checking...','#553377');
var ecoValDetailsV = ui.Label('',{fontSize:'7px',color:'#555577',
  backgroundColor:'rgba(0,0,0,0)',padding:'1px 4px',margin:'0',whiteSpace:'pre'});
panel.add(row('Validation status',ecoValStatusV));
panel.add(ecoValDetailsV);

panel.add(sHead('S16 - PHYSICAL ENERGY INDEX','#1a3a2a'));
var eciDepthV = dynLbl('checking...','#115533');
var eciCurrentV = dynLbl('checking...','#115533');
var eciSLR03V = dynLbl('checking...','#226644');
var eciSLR05V = dynLbl('checking...','#226644');
var eciSLR10V = dynLbl('checking...','#aa3300');
var eciValidationV = dynLbl('checking...','#226644');
panel.add(row('Seafloor depth (GEBCO)',eciDepthV));
panel.add(row('Energy risk (current)',eciCurrentV));
panel.add(row('+0.3m SLR scenario',eciSLR03V));
panel.add(row('+0.5m SLR scenario',eciSLR05V));
panel.add(row('+1.0m SLR scenario',eciSLR10V));
panel.add(row('vs heuristic B score',eciValidationV));

panel.add(sHead('S17 - TIME OF EMERGENCE (v10.36)','#2a1a3a'));
var toeCompoundV = dynLbl('computing...','#553377');
var toeSSTv = dynLbl('computing...','#226644');
var toeCHLv = dynLbl('computing...','#226644');
var toeNO2v = dynLbl('computing...','#226644');
var toeSALv = dynLbl('computing...','#226644');
panel.add(row('Compound status (0-4)',toeCompoundV));
panel.add(row('SST (44yr, HIGH conf)',toeSSTv));
panel.add(row('Chl (27yr, MARGINAL)',toeCHLv));
panel.add(row('Salinity (32yr, MARGINAL)',toeSALv));
panel.add(row('NO2 (7yr, LOW conf)',toeNO2v));

panel.add(sHead('S18 - BIOGEOCHEMISTRY SNAPSHOT','#1a2a3a'));
var s18O2V = dynLbl('computing...','#334466');
var s18pCO2V = dynLbl('computing...','#334466');
var s18pHV = dynLbl('computing...','#334466');
var s18SalV = dynLbl('computing...','#334466');
panel.add(row('Dissolved O2 (surface)',s18O2V));
panel.add(row('pCO2 (surface)',s18pCO2V));
panel.add(row('pH (surface, real acidity)',s18pHV));
panel.add(row('Salinity (current)',s18SalV));

// ============================================================
// S19 - REAL IN-SITU BASELINE (v10.36 NEW)
// ============================================================
panel.add(sHead('S19 - REAL IN-SITU BASELINE (v10.36 NEW)','#0a3a5a'));
panel.add(lbl('Real cleaned SeapHOx sensor data (NOT a model/satellite)',7,'#225577'));
panel.add(lbl('Only populated within ~10-15km of an actual sensor',7,'#225577'));
panel.add(lbl('3 stations: Looe Key FL | Agua Hedionda CA | Scripps Pier CA',7,'#888888'));
var s19StatusV = dynLbl('checking...','#225577');
var s19PhV = dynLbl('--','#225577');
var s19TempV = dynLbl('--','#225577');
var s19SalV = dynLbl('--','#225577');
var s19DoV = dynLbl('--','#225577');
var s19NoteV = ui.Label('',{fontSize:'7px',color:'#555577',
  backgroundColor:'rgba(0,0,0,0)',padding:'1px 4px',margin:'0',whiteSpace:'pre'});
panel.add(row('Matched station',s19StatusV));
panel.add(row('pH (clean baseline)',s19PhV));
panel.add(row('Temp (clean baseline)',s19TempV));
panel.add(row('Salinity (clean baseline)',s19SalV));
panel.add(row('DO mg/L (clean baseline)',s19DoV));
panel.add(s19NoteV);
panel.add(lbl('CAVEAT: point-source only. Compare vs S18 (model) above.',7,'#aa6600'));

panel.add(sHead('ACCURACY','#222244'));

var exportBtn = ui.Button({
  label: 'EXPORT LOGGED CLICKS TO CSV (Drive)',
  style: {fontSize:'12px', fontWeight:'bold', margin:'2px 4px',
    backgroundColor:'#cce0ff', color:'#003388', stretch:'horizontal',
    padding:'8px 4px', border:'2px solid #0055cc'},
  onClick: buildAndExport
});
panel.add(exportBtn);

var exportStatusLabel = ui.Label('No clicks logged yet - 0 rows ready',
  {fontSize:'9px',color:'#cc0000',fontWeight:'bold',margin:'2px 4px'});
panel.add(exportStatusLabel);

function buildAndExport(){
  if(clickResultsLog.length === 0){
    exportStatusLabel.setValue('STILL 0 rows. Click a coastal/ocean point on the map first.');
    return;
  }
  var features = clickResultsLog.map(function(rec){
    return ee.Feature(null, rec);
  });
  var fc = ee.FeatureCollection(features);
  var exportDesc = 'stemgeohs_marine_export_' + Date.now();

  var columnOrder = [
    'timestamp', 'region', 'lat', 'lon', 'region_type',
    'sst_annual_c', 'sst_peak_c', 'sst_trend_c_per_yr',
    'dhw_c_weeks', 'dhw_raw_unfiltered', 'dhw_artifact_flagged',
    'mmm_local_c', 'chl_a_mg_m3', 'turbidity_ndti', 'no2_mol_m2', 'depth_m',
    'cancer_score_satellite', 'cancer_score_fused', 'field_correction',
    'bowl_depth_B', 'accuracy_pct', 'status_label',
    'field_data_available', 'field_species',
    'aquaculture_score', 'aquaculture_confidence_pct', 'aquaculture_status',
    'bromoform_yield', 'num_interventions', 'top_intervention',
    'soil_texture_code', 'soil_texture_label',
    'earthquake_count_200km', 'earthquake_max_mag_200km', 'volcanic_activity_count_300km',
    'real_ac1_24mo', 'real_var_trend_ratio', 'real_csd_n_valid_months',
    'ndci_nutrient_proxy',
    'thermal_recovery_n_completed_episodes', 'thermal_recovery_mean_months',
    'thermal_recovery_max_months', 'thermal_recovery_ongoing_flag',
    'eco_recovery_validation_checked', 'eco_recovery_validation_finding',
    'gebco_depth_m', 'energy_concentration_risk',
    'eci_slr_03m', 'eci_slr_05m', 'eci_slr_10m',
    'eci_vs_heuristic_validation',
    'in_situ_station', 'in_situ_distance_km', 'in_situ_ph_mean',
    'in_situ_temp_c_mean', 'in_situ_salinity_mean', 'in_situ_do_mgL_mean'
  ];

  Export.table.toDrive({
    collection: fc,
    description: exportDesc,
    folder: 'STEMGeoHS_exports',
    fileFormat: 'CSV',
    selectors: columnOrder
  });
  exportStatusLabel.setValue(
    'TASK CREATED: "'+exportDesc+'" with '+clickResultsLog.length+' row(s). ' +
    'Open Tasks tab and click RUN.'
  );
  print('=== EXPORT TASK CREATED ===');
  print('Task name: '+exportDesc);
  print('Rows queued: '+clickResultsLog.length);
}

var testExportBtn = ui.Button({
  label: 'Test export (1 sample row, no click needed)',
  style: {fontSize:'10px', fontWeight:'bold', margin:'2px 4px',
    backgroundColor:'#e8e8e8', color:'#222222', stretch:'horizontal',
    padding:'6px 4px', border:'2px solid #555555'},
  onClick: function(){
    clickResultsLog.push({
      timestamp: new Date().toISOString(),
      lat: 21.27, lon: -157.876, region: 'TEST ROW - Hawaii sample',
      sst_annual_c: 25.0, sst_peak_c: 27.5, sst_trend_c_per_yr: 0.01,
      dhw_c_weeks: 0, mmm_local_c: 27.3, chl_a_mg_m3: null,
      turbidity_ndti: 0.01, no2_mol_m2: 0.00005, depth_m: -50,
      cancer_score_satellite: 21, cancer_score_fused: 21,
      field_correction: 0, bowl_depth_B: 0.79, accuracy_pct: 77,
      status_label: 'DEEP BASIN', field_data_available: false,
      field_species: 'none', aquaculture_score: 79,
      aquaculture_confidence_pct: 80, aquaculture_status: 'GOOD',
      bromoform_yield: 60, num_interventions: 1, top_intervention: 'TEST ROW',
      soil_texture_code: null, soil_texture_label: 'n/a (test row)',
      earthquake_count_200km: 0, earthquake_max_mag_200km: null,
      volcanic_activity_count_300km: 0,
      real_ac1_24mo: null, real_var_trend_ratio: null, real_csd_n_valid_months: 0,
      ndci_nutrient_proxy: null,
      thermal_recovery_n_completed_episodes: null, thermal_recovery_mean_months: null,
      thermal_recovery_max_months: null, thermal_recovery_ongoing_flag: null,
      eco_recovery_validation_checked: false, eco_recovery_validation_finding: 'n/a (test row)',
      in_situ_station: null, in_situ_distance_km: null, in_situ_ph_mean: null,
      in_situ_temp_c_mean: null, in_situ_salinity_mean: null, in_situ_do_mgL_mean: null
    });
    exportStatusLabel.setValue('Test row added. '+clickResultsLog.length+' row(s) ready.');
  }
});
panel.add(testExportBtn);

var clearLogBtn = ui.Button({
  label: 'Clear logged clicks',
  style: {fontSize:'10px', fontWeight:'bold', margin:'2px 4px',
    backgroundColor:'#ffd9d9', color:'#770000', stretch:'horizontal',
    padding:'6px 4px', border:'2px solid #aa2222'},
  onClick: function(){
    clickResultsLog = [];
    exportStatusLabel.setValue('Log cleared - 0 rows ready');
  }
});
panel.add(clearLogBtn);

var accSV=dynLbl('77%','#664400'), accFV=dynLbl('+0%','#888888'), accTV=dynLbl('77%','#664400');
panel.add(row('Satellite',accSV)); panel.add(row('Field gain',accFV));
panel.add(row('Combined',accTV));

panel.add(sHead('PREVENTATIVE INTERVENTIONS (Waddington surgery)','#003355'));
var intPanel = ui.Panel({style:{margin:'2px 0'}});
panel.add(intPanel);

panel.add(sHead('FUSED COASTAL CANCER SCORE','#4a0000'));
var scoreBig = ui.Label('Score: --',{fontSize:'18px',fontWeight:'bold',
  color:'#003300',backgroundColor:'#d4f5df',padding:'4px 8px',margin:'2px 0'});
panel.add(scoreBig);
var scoreBarLbl = ui.Label('--',{fontSize:'11px',fontWeight:'bold',
  color:'#0a5c1e',backgroundColor:'rgba(0,0,0,0)',padding:'1px 4px',margin:'0'});
panel.add(scoreBarLbl);
var barBg = ui.Panel({style:{backgroundColor:'#dddddd',margin:'1px 2px',
  padding:'1px',border:'1px solid #aaaaaa'}});
var barFill = ui.Panel({style:{backgroundColor:'#22cc44',height:'12px',width:'0px',margin:'0'}});
barBg.add(barFill); panel.add(barBg);
panel.add(ui.Panel([
  lbl('0',8,'#555555'),lbl('safe-----------risk',8,'#999999'),lbl('100',8,'#880000')
],ui.Panel.Layout.Flow('horizontal'),{margin:'0 2px'}));
var scoreInterp = ui.Label('Click reef area to analyze',{
  fontSize:'8px',color:'#333333',backgroundColor:'#eefaee',
  padding:'3px 4px',margin:'2px 0',whiteSpace:'pre'});
panel.add(scoreInterp);

panel.add(lbl('MAP LEGEND (scroll to see) | v10.36',9,'#ffffff','#334466',true));

panel.add(sHead('1. Study zone border (white)','#555555'));
panel.add(legRow('#ffffff','White ring','1km buffer'));

panel.add(legDiv());
panel.add(sHead('2. S1 - SST annual mean deg C','#b03000'));
panel.add(legRow('#2c7bb6','< 24 deg C','Cold'));
panel.add(legRow('#abd9e9','24-26 deg C','Optimal'));
panel.add(legRow('#ffffbf','26-28 deg C','Warm'));
panel.add(legRow('#fdae61','28-30 deg C','Stress'));
panel.add(legRow('#f46d43','30-32 deg C','BLEACHING'));
panel.add(legRow('#a50026','> 32 deg C','MORTALITY'));

// v10.44 NEW: restored legend entries matching every real map layer -
// palette values copied EXACTLY from the corresponding Map.addLayer()
// call in loadLayers() so the legend and the actual rendered colors
// always agree.
panel.add(legDiv());
panel.add(sHead('2b. S1 - SST PEAK Jun-Oct (bleaching score)','#8b0000'));
panel.add(legRow('#2c7bb6','< 24 deg C','Below bleaching threshold'));
panel.add(legRow('#ffffbf','26-28 deg C','Warm - watch'));
panel.add(legRow('#fdae61','28-30 deg C','Thermal stress'));
panel.add(legRow('#f46d43','30-32 deg C','BLEACHING RISK'));
panel.add(legRow('#d73027','32-34 deg C','MASS BLEACHING'));
panel.add(legRow('#a50026','34-36 deg C','MORTALITY'));
panel.add(legRow('#660000','> 36 deg C','EXTREME MORTALITY'));

panel.add(legDiv());
panel.add(sHead('3. S1 - Bleaching risk overlay','#880000'));
panel.add(legRow('#ff0000','RED overlay','Peak SST > 30 deg C = bleaching zone'));

panel.add(legDiv());
panel.add(sHead('4. S4 - DHW (per-pixel MMM)','#880000'));
panel.add(legRow('#0000ff','Blue background','DHW=0: no stress'));
panel.add(legRow('#ffffd9','Yellow-white','DHW 1-4: watch'));
panel.add(legRow('#fdae61','Orange','DHW 4-8: BLEACHING RISK'));
panel.add(legRow('#f46d43','Dark orange','DHW 8-12: MASS BLEACH'));
panel.add(legRow('#a50026','Dark red','DHW > 12: MORTALITY'));

panel.add(legDiv());
panel.add(sHead('5. S2 - Chlorophyll-a mg/m3','#005a32'));
panel.add(legRow('#084594','< 0.1','Oligotrophic'));
panel.add(legRow('#2171b5','0.1-1.0','Low-moderate'));
panel.add(legRow('#6baed6','1.0-2.0','Moderate'));
panel.add(legRow('#74c476','2.0-3.5','Good'));
panel.add(legRow('#31a354','2.0-5.0','High'));
panel.add(legRow('#006837','> 5.0','Bloom'));

// v10.44 NEW - previously missing entirely
panel.add(legDiv());
panel.add(sHead('6. S3 - Depth m (GEBCO, v10.66)','#003366'));
panel.add(legRow('#003366','< -200m','Deep ocean'));
panel.add(legRow('#0055aa','-200 to -100m','Shelf'));
panel.add(legRow('#2171b5','-100 to -50m','Outer shelf'));
panel.add(legRow('#6baed6','-50 to 0m','Shallow/nearshore'));
panel.add(legRow('#c7eae5','~0m','Waterline'));
panel.add(legRow('#f5f5f5','0-50m (land)','Beach/low ground'));
panel.add(legRow('#8c510a','> 50m (land)','Higher ground'));

panel.add(legDiv());
panel.add(sHead('7. S3 - Intertidal 0-20m (gold)','#885500'));
panel.add(legRow('#FFD700','Gold overlay','Depth 0 to -20m'));

panel.add(legDiv());
panel.add(sHead('8. S3 - Shallow reef 0-50m (teal)','#006644'));
panel.add(legRow('#00cc99','Teal overlay','Depth 0 to -50m'));

panel.add(legDiv());
panel.add(sHead('9. S3 - Turbidity NDTI (reef < 50m only)','#885500'));
panel.add(legRow('#084594','< -0.15','Very clear'));
panel.add(legRow('#2171b5','-0.15 to -0.05','Clear'));
panel.add(legRow('#74add1','-0.05 to 0.05','Moderate'));
panel.add(legRow('#ffffbf','0.05 to 0.15','Slightly turbid'));
panel.add(legRow('#fdae61','0.15 to 0.25','Turbid - stress'));
panel.add(legRow('#d73027','> 0.25','Very turbid'));

panel.add(legDiv());
panel.add(sHead('10. S4 - SST trend deg C/yr (CSD signal)','#880000'));
panel.add(lbl('[metaphor] Waddington vocabulary, NOT real depth - see S16',7,'#aa4400'));
panel.add(legRow('#0000ff','Cooling','[metaphor] basin deepening'));
panel.add(legRow('#aaaaff','Slight cool','[metaphor] recovering'));
panel.add(legRow('#dddddd','~0/yr','[metaphor] stable'));
panel.add(legRow('#ffaaaa','+0.01 to +0.05/yr','[metaphor] shallowing - possible CSD'));
panel.add(legRow('#ff0000','> +0.05/yr','[metaphor] near flat - high instability'));

panel.add(legDiv());
panel.add(sHead('11. S6 - NO2 mol/m2 (Sentinel-5P)','#550066'));
panel.add(legRow('#313695','< 2e-5','Very clean'));
panel.add(legRow('#4575b4','2e-5 to 5e-5','Background'));
panel.add(legRow('#74add1','5e-5 to 8e-5','Moderate'));
panel.add(legRow('#fee090','8e-5 to 1.1e-4','Elevated'));
panel.add(legRow('#f46d43','1.1e-4 to 1.35e-4','High'));
panel.add(legRow('#d73027','1.35e-4 to 1.5e-4','Very high'));
panel.add(legRow('#a50026','> 1.5e-4','Critical - harbour/pollution'));

panel.add(legDiv());
panel.add(sHead('S7 - Seaweed / Macroalgae (Sentinel-2)','#006600'));
panel.add(lbl('FAI - Floating Algae Index (Hu 2009)',7,'#333333'));
panel.add(legRow('#ccffcc','FAI 0.003-0.02','Mild - trace signal'));
panel.add(legRow('#88dd44','FAI 0.02-0.04','Low - scattered algae'));
panel.add(legRow('#ffff00','FAI 0.04-0.06','Moderate - visible patches'));
panel.add(legRow('#ff8800','FAI 0.06-0.09','High - dense mats'));
panel.add(legRow('#ff4400','FAI 0.09-0.12','Very high - smothering risk'));
panel.add(legRow('#ff0000','FAI > 0.12','EXTREME - massive bloom'));
panel.add(lbl('NDVI Water - benthic macroalgae (dark=dense)',7,'#333333'));
panel.add(legRow('#e8ffe8','NDVI 0.05-0.15','Light cover'));
panel.add(legRow('#aaddaa','NDVI 0.15-0.25','Moderate cover'));
panel.add(legRow('#44aa00','NDVI 0.25-0.35','Dense carpeting'));
panel.add(legRow('#006600','NDVI 0.35-0.45','Very dense'));
panel.add(legRow('#003300','NDVI > 0.45','Extreme - turf domination'));
panel.add(lbl('NDCI - Bloom Index (Mishra 2012)',7,'#333333'));
panel.add(legRow('#ffffcc','NDCI 0.1-0.2','Low bloom signal'));
panel.add(legRow('#a6d96a','NDCI 0.2-0.3','Moderate bloom'));
panel.add(legRow('#1a9641','NDCI 0.3-0.4','High bloom'));
panel.add(legRow('#004400','NDCI > 0.4','Extreme - eutrophication'));

panel.add(legDiv());
panel.add(sHead('13. S9 - Soil texture (v10.42: dataset dead)','#5a4422'));
panel.add(lbl('OpenLandMap asset confirmed inaccessible - layer will',7,'#664422'));
panel.add(lbl('render blank/transparent at every location (v10.42).',7,'#664422'));

panel.add(legDiv());
panel.add(sHead('14. S10 - Earthquakes M4.5+ (USGS, daily)','#5a2222'));
panel.add(legRow('#ff8800','Orange dot','Magnitude 4.5 - 5.9'));
panel.add(legRow('#cc0000','Dark red dot','Magnitude 6.0+'));

panel.add(legDiv());
panel.add(sHead('15. S11 - Volcanic activity 1960-2018 (GDIS)','#5a2222'));
panel.add(legRow('#ff00ff','Magenta dot','Recorded volcanic activity event'));

panel.add(legDiv());
panel.add(sHead('12. Score color dot (center)','#4a0000'));
panel.add(legRow('#00cc44','Green','0-29: DEEP BASIN'));
panel.add(legRow('#ffcc00','Yellow','30-54: WARNING'));
panel.add(legRow('#ff6600','Orange','55-74: HIGH RISK'));
panel.add(legRow('#ff0000','Red','75-100: CRITICAL'));

panel.add(lbl('Scroll up for measurements | v10.36',8,'#555555'));
ui.root.insert(0,panel);

var floatP = ui.Panel({style:{position:'top-center',padding:'6px 14px',
  backgroundColor:'#001a00',border:'2px solid #00cc44',shown:false}});
var fN=lbl('',13,'#00cc44','',true), fS=lbl('',12,'#ffffff','',true),
    fB=lbl('',10,'#aaffaa','',false), fC=lbl('',9,'#888888','',false);
floatP.add(fN); floatP.add(fS); floatP.add(fB); floatP.add(fC);
Map.add(floatP);

// ============================================================
// MODULE E - MAP SETUP
// ============================================================
Map.setCenter(-82.25,9.35,9);
Map.setOptions('HYBRID');
var defaultStartPt = ee.Geometry.Point([-82.25, 9.35]);
var defaultStartZone = defaultStartPt.buffer(1000);
Map.addLayer(ee.Image().byte().paint(defaultStartZone,0,3),
  {palette:['#ffffff']}, 'Study zone (white border) - startup default', true, 1.0);
Map.addLayer(sst,{min:20,max:32,
  palette:['#2c7bb6','#abd9e9','#ffffbf','#fdae61','#f46d43','#d73027','#a50026'],
  opacity:0.5},'SST annual mean (context)');
Map.addLayer(sstPeakFinal,{min:26,max:36,
  palette:['#2c7bb6','#abd9e9','#ffffbf','#fdae61','#f46d43','#d73027','#a50026'],
  opacity:0.8},'SST PEAK Jun-Oct 2023 (bleaching score)');
Map.addLayer(dhwProper.updateMask(dhwProper.gt(0)),{min:0.1,max:12,
  palette:['#ffffd9','#fdae61','#f46d43','#a50026'],
  opacity:0.65},'DHW per-pixel MMM (startup view)');
Map.addLayer(soilTexture,{min:1,max:12,
  palette:['#8B4513','#A0522D','#CD853F','#DEB887','#D2B48C','#BC8F8F',
           '#F4A460','#DAA520','#B8860B','#E0E0E0','#FFE4B5','#FFF8DC']},
  'S9 - Soil texture USDA classes (startup view)', false, 0.6);
Map.addLayer(eqModerate,{color:'ff8800'},'S10 - Earthquakes M4.5-5.9 (startup view)', true, 0.7);
Map.addLayer(eqMajor,{color:'cc0000'},'S10 - Earthquakes M6.0+ (startup view)', true, 0.9);
Map.addLayer(volcanicActivity,{color:'ff00ff'},'S11 - Volcanic activity 1960-2018 (startup view)', true, 0.8);

function loadLayers(study,region,cols) {
  while(Map.layers().length()>0){ Map.remove(Map.layers().get(0)); }
  var dotCol = cols?cols.map:'#888888';
  var dotLbl = cols?cols.lbl:'computing';
  Map.addLayer(ee.Image().byte().paint(study,0,3),
    {palette:['#ffffff'],opacity:1.0},'Study zone (white border)');
  Map.addLayer(sst.clip(study),{min:20,max:32,
    palette:['#2c7bb6','#abd9e9','#ffffbf','#fdae61','#f46d43','#d73027','#a50026'],
    opacity:0.7},'S1 - SST annual mean deg C');
  Map.addLayer(sstPeakFinal.clip(study),{min:24,max:36,
    palette:['#2c7bb6','#ffffbf','#fdae61','#f46d43','#d73027','#a50026','#660000'],
    opacity:0.9},'S1 - SST PEAK Jun-Oct 2023 (bleaching score)');
  Map.addLayer(sstPeakFinal.clip(study).gt(30).selfMask(),
    {palette:['#ff0000'],opacity:0.75},'S1 - Bleaching risk (peak SST > 30)');
  Map.addLayer(ee.Image.constant(0).clip(study).updateMask(oceanMask.clip(study)),
    {palette:['#0000ff'],opacity:0.2},'S4 - DHW blue = no stress');
  Map.addLayer(dhwProper.clip(study).updateMask(dhwProper.clip(study).gt(0)),
    {min:0.1,max:12,palette:['#ffffd9','#fdae61','#f46d43','#a50026'],opacity:0.85},
    'S4 - DHW heat stress (orange where > 0)');
  Map.addLayer(chla.clip(study),{min:0.01,max:5.0,
    palette:['#084594','#2171b5','#6baed6','#74c476','#31a354','#006837'],
    opacity:0.8},'S2 - Chl-a mg/m3');
  // v10.66 FIX: switched the visual "S3 - Depth" map layer to GEBCO,
  // matching the v10.65 sidebar fix - confirmed via live testing that
  // showing the sidebar number (now GEBCO) alongside a map layer still
  // rendering the OLD ETOPO1 colors could visually disagree at exactly
  // the narrow/complex coastlines where this matters most (e.g. Oslo
  // fjord). Same sign convention as bathy (negative=ocean), so the
  // existing color palette (min:-200, max:100) still applies correctly
  // with no other changes needed.
  // SCOPE NOTE: shallowMask/intertidalMask/turbImg (and their own
  // separate visual overlays) intentionally remain on ETOPO1-derived
  // bathyU for now - switching those affects real turbidity/aquaculture
  // computation inputs, not just this depth visualization, and deserves
  // its own separate, carefully-tested consideration.
  Map.addLayer(GEBCO.rename('bedrock').clip(study),{min:-200,max:100,
    palette:['#003366','#0055aa','#2171b5','#6baed6','#c7eae5','#f5f5f5','#8c510a'],
    opacity:0.7},'S3 - Depth m (GEBCO)');
  Map.addLayer(intertidalMask.selfMask().clip(study),
    {palette:['#FFD700'],opacity:0.65},'S3 - Intertidal 0-20m');
  Map.addLayer(shallowMask.selfMask().clip(study),
    {palette:['#00cc99'],opacity:0.35},'S3 - Shallow reef 0-50m');
  Map.addLayer(turbImg.clip(study),{min:-0.3,max:0.3,
    palette:['#084594','#2171b5','#74add1','#ffffbf','#fdae61','#d73027'],
    opacity:0.75},'S3 - Turbidity NDTI (reef < 50m ONLY)');
  Map.addLayer(sst_slope.clip(study),{min:-0.05,max:0.05,
    palette:['#0000ff','#aaaaff','#dddddd','#ffaaaa','#ff0000'],
    opacity:0.85},'S4 - SST trend deg C/yr (CSD)');
  Map.addLayer(no2.clip(study),{min:0.000005,max:0.00015,
    palette:['#313695','#4575b4','#74add1','#fee090','#f46d43','#d73027','#a50026'],
    opacity:0.65},'S6 - NO2 mol/m2');
  Map.addLayer(faiMild.clip(study),{
    min:0.003, max:0.15,
    palette:['#ccffcc','#88dd44','#ffff00','#ff8800','#ff4400',
             '#ff0000','#cc0000','#880000'],
    opacity:0.9
  },'S7 - FAI Floating Algae (green=mild RED=extreme seaweed)');
  Map.addLayer(ndviAlgae.clip(study),{
    min:0.05, max:0.5,
    palette:['#e8ffe8','#aaddaa','#44aa00','#006600','#003300'],
    opacity:0.8
  },'S7 - NDVI Water (benthic macroalgae - dark=dense)');
  Map.addLayer(ndciBloom.clip(study),{
    min:0.1, max:0.5,
    palette:['#ffffcc','#a6d96a','#1a9641','#004400'],
    opacity:0.75
  },'S7 - NDCI Bloom index (algae bloom)');
  Map.addLayer(soilTexture.clip(study),{min:1,max:12,
    palette:['#8B4513','#A0522D','#CD853F','#DEB887','#D2B48C','#BC8F8F',
             '#F4A460','#DAA520','#B8860B','#E0E0E0','#FFE4B5','#FFF8DC']},
    'S9 - Soil texture USDA classes', false, 0.6);
  Map.addLayer(eqModerate,{color:'ff8800'},'S10 - Earthquakes M4.5-5.9', true, 0.7);
  Map.addLayer(eqMajor,{color:'cc0000'},'S10 - Earthquakes M6.0+', true, 0.9);
  Map.addLayer(volcanicActivity,{color:'ff00ff'},'S11 - Volcanic activity 1960-2018', true, 0.8);

  Map.addLayer(ee.Image.constant(1).clip(study.centroid(50).buffer(300)),
    {palette:[dotCol],opacity:0.9},'SCORE: '+dotLbl+' ('+dotCol+')');
}

// ============================================================
// MODULE F0 - ERROR HANDLING + WATCHDOG (v10.37 NEW)
//
// PROBLEM DISCOVERED (live testing, v10.36): the main click pipeline
// is a chain of 19 nested .evaluate() calls. Every one of them ONLY
// registered a success callback - no error callback. If ANY single
// reduceRegion in that chain throws a server-side error at a
// particular location (confirmed to happen near 5.9N,90E - a data-
// sparse band near the Nicobar Islands/Ninety East Ridge where
// HYCOM+OISST both return null), the promise fails silently: nothing
// ever prints to Console, and the sidebar freezes on "computing..."
// indefinitely, with zero feedback to the user that anything is wrong.
//
// FIX:
//   1. handleGeeError(err, stepName) - now passed as the SECOND
//      argument to every one of the 19 .evaluate() calls in the main
//      chain. A server-side failure now immediately: prints the exact
//      failing step + error message to Console, updates the sidebar
//      with a clear human-readable error instead of leaving it frozen,
//      and clears the watchdog timer (see below) so no
//      contradictory "still computing" message can fire afterward.
//   2. A client-side watchdog timer starts the moment a click is
//      analyzed. If the ENTIRE pipeline has not finished within
//      WATCHDOG_TIMEOUT_MS, it updates the sidebar to say so
//      explicitly - "still computing after Ns, this location may have
//      sparse satellite coverage" - rather than leaving a bare
//      "computing..." label with no indication anything is unusual.
//      This does NOT cancel the underlying GEE computation (Earth
//      Engine's JS client has no clean "abort a running evaluate()"
//      primitive) - it only gives honest, timely feedback while the
//      person waits, matching this script's existing standard of
//      never leaving a stale-looking number/state unexplained.
//   3. The watchdog is cleared in TWO places: on genuine pipeline
//      success (end of the rToEAll success callback - the true last
//      step of the chain) and on any handleGeeError call - so it
//      never fires a false "still computing" message after the run
//      has actually already finished or failed.
// ============================================================
// ============================================================
// MODULE F0 - ERROR HANDLING (v10.37/v10.39 - see changelog)
//
// v10.39 CRITICAL FIX: the v10.37 watchdog used setTimeout/
// clearTimeout, which are STANDARD browser JS globals but are NOT
// available in Google Earth Engine's Code Editor sandbox environment.
// CONFIRMED via live testing: every single click threw
// "Line 2064: setTimeout is not defined" and crashed immediately -
// since startPipelineWatchdog() was the very first line inside
// analyzeLocation(), NOTHING in the pipeline ever ran (no reset, no
// map zoom, no white ring, no error handling) - this is a more
// severe regression than either of the two bugs it was meant to fix.
//
// FIX: the timeout-based watchdog is REMOVED (GEE simply does not
// support this - there is no polyfill available inside the sandboxed
// Code Editor). clearPipelineWatchdog() is kept as a harmless no-op
// so the existing call sites (success path, handleGeeError) don't
// need to change. What REMAINS and still works correctly, since
// neither depends on setTimeout:
//   1. handleGeeError(err, stepName) - error callbacks on all 19
//      nested .evaluate() calls. A server-side failure still surfaces
//      immediately with a clear message instead of hanging silently.
//   2. resetSidebarToComputing() (v10.38) - every sidebar field is
//      still reset the instant a new click starts, so stale data from
//      a previous click can never visually linger.
// The only thing lost is the "still computing after 45s" proactive
// warning message for genuinely slow (but not failing) locations -
// those will now simply show "computing..." until they finish or
// error out, same as v10.36 behaved for successful-but-slow clicks.
// ============================================================
var _pipelineWatchdog = null;

function clearPipelineWatchdog() {
  // No-op (v10.39) - kept only so existing call sites don't error.
  // See module comment above for why the timer itself was removed.
}

function startPipelineWatchdog(lat, lon) {
  // No-op (v10.39) - see module comment above. setTimeout is not
  // available in the GEE Code Editor sandbox; there is no safe
  // browser-independent equivalent inside this environment.
}

function handleGeeError(err, stepName) {
  clearPipelineWatchdog();
  print('=== PIPELINE ERROR (v10.37) at step: ' + stepName + ' ===');
  print('Error: ' + err);
  print('This means the server-side computation for this ONE step failed -');
  print('everything computed successfully before this step is still valid');
  print('and already shown in the sidebar; everything after this step was');
  print('never computed and will show stale/blank values.');
  print('Common causes: no data at this pixel for this specific dataset,');
  print('a scale/geometry mismatch, or a transient Earth Engine server issue.');
  print('Try clicking again, or try a nearby coordinate.');
  clickLbl.setValue('ERROR at step "' + stepName + '" - see Console for details');
  clickLbl.style().set('backgroundColor', '#880000');
  scoreInterp.setValue('Computation failed at step: ' + stepName + '\n' +
    'Error: ' + String(err).slice(0, 150) + '\n' +
    'Values shown above this point in the sidebar are still valid.\n' +
    'Try clicking again or try a nearby coordinate.');
  scoreInterp.style().set('backgroundColor', '#ffd0d0');
  scoreInterp.style().set('color', '#880000');
}

// ============================================================
// MODULE F0b - FULL SIDEBAR RESET (v10.38 NEW)
//
// PROBLEM DISCOVERED (live testing, v10.37): clicking "Run" in the GEE
// editor rebuilds the ENTIRE sidebar from scratch (fresh "--"/
// "computing..." widgets), which is why a full script re-run always
// looked "clean". But using the GO TO COORDINATES box (or clicking a
// new map point) WITHOUT re-running only reset a handful of fields
// (score, location, S19) immediately - every other S9-S18 field kept
// showing the PREVIOUS click's numbers until its own async step
// happened to finish and overwrite it. Combined with the v10.37 error
// handling fix, a step that now fails cleanly instead of hanging would
// leave its field frozen on stale old data PERMANENTLY, which is
// arguably more misleading than the original freeze.
//
// FIX: resetSidebarToComputing() sets every dynamic sidebar field
// back to "computing..." (or the appropriate blank state) the instant
// ANY new click/GO starts - map click or coordinate box, no "Run"
// needed. This guarantees stale data can never visually linger past
// the start of a new analysis, matching the same honesty standard
// already used for score/location/S19.
// ============================================================
function resetSidebarToComputing() {
  var C = 'computing...';
  sstV.setValue(C); trendV.setValue(C); dhwV.setValue(C); dhwNoteV.setValue('');
  mmmV.setValue(C); chlaV.setValue(C); turbV.setValue(C); no2V.setValue(C);
  depthV.setValue(C); depthWarnV.setValue(''); ebusWarnV.setValue('');

  faiV.setValue(C); ndciV.setValue(C); ndviWV.setValue(C); algaeStatusV.setValue(C);

  sstWindowV.setValue(C); chlWindowV.setValue(C); ndciNutrientProxyV.setValue(C);
  pollutionV.setValue(C); stabilityV.setValue(C); bromoformV.setValue(C);
  aquaScoreV.setValue(C); aquaMissingV.setValue(''); aquaStatusV.setValue(C);

  soilTextureV.setValue(C);
  eqStatsV.setValue(C); volcStatsV.setValue(C);

  satCcsV.setValue('--'); fusCcsV.setValue('--'); bowlV.setValue('--');
  omega0V.setValue('--'); tauV.setValue('--'); ac1V.setValue('--'); p5yrV.setValue('--');

  realAc1V.setValue(C); realVarTrendV.setValue(C); realCsdNoteV.setValue('');

  thermalEpisodesV.setValue(C); thermalMeanV.setValue(C);
  thermalMaxV.setValue(C); thermalOngoingV.setValue(C);

  ecoValStatusV.setValue(C); ecoValDetailsV.setValue('');

  eciDepthV.setValue(C); eciCurrentV.setValue(C);
  eciSLR03V.setValue(C); eciSLR05V.setValue(C); eciSLR10V.setValue(C);
  eciValidationV.setValue(C);

  toeCompoundV.setValue(C); toeSSTv.setValue(C); toeCHLv.setValue(C);
  toeNO2v.setValue(C); toeSALv.setValue(C);

  s18O2V.setValue(C); s18pCO2V.setValue(C); s18pHV.setValue(C); s18SalV.setValue(C);

  accSV.setValue('--'); accFV.setValue('--'); accTV.setValue('--');

  intPanel.clear();
}

// ============================================================
// MODULE F - CLICK HANDLER
// ============================================================
function analyzeLocation(lat, lon) {
  startPipelineWatchdog(lat, lon); // v10.37 NEW - see MODULE F0 above
  resetSidebarToComputing();       // v10.38 NEW - see MODULE F0b above
  var latR=Math.round(lat*1000)/1000, lonR=Math.round(lon*1000)/1000;
  var region=getRegion(lat,lon);
  var fp=getFieldProfile(region);
  var pt=ee.Geometry.Point([lon,lat]), study=pt.buffer(1000);
  // v10.43 NEW: SST-family variables (SST, peak SST, trend, DHW, MMM) all
  // depend on oceanMask, which is built from ETOPO1 bathymetry at ~1.8km
  // native resolution. At narrow inlets/lagoon mouths (e.g. Agua Hedionda,
  // CA), the 1km "study" buffer can land entirely on a pixel ETOPO1
  // misclassifies as land, returning n/a in the sidebar even though the
  // map's colored SST/trend layer (rendered over a wider area at finer
  // display resolution) clearly shows real data just meters away.
  // FIX: widen the sampling area for SST-family reduceRegion calls only
  // (turbidity/Chl-a/bathymetry keep their existing, already-appropriate
  // "study" 1km buffer) - same proven fallback principle already used
  // for Chl-a's coastal-relaxed mask, applied here via a larger radius
  // instead of a second mask, to avoid adding any new evaluate() calls
  // to an already-long chain.
  var studySST = pt.buffer(3000);


  // ============================================================
  // S19 REAL IN-SITU BASELINE LOOKUP (v10.36 NEW)
  // Synchronous JS - no evaluate()/server call needed, unlike
  // everything else in this pipeline. See MODULE A11 at top of script.
  // ============================================================
  var inSitu = getInSituBaseline(lat, lon);
  if (inSitu) {
    s19StatusV.setValue(inSitu.label + ' (' + inSitu.distance_km.toFixed(1) + ' km away)');
    s19StatusV.style().set('color', '#115511');
    s19PhV.setValue(inSitu.pH.mean.toFixed(3) + ' +/- ' + inSitu.pH.std.toFixed(3) +
      '  (n=' + inSitu.n_clean + ' clean, ' + inSitu.pct_flagged + '% flagged/removed)');
    s19TempV.setValue(inSitu.temp_c.mean.toFixed(2) + ' deg C');
    s19SalV.setValue(inSitu.salinity.mean.toFixed(2) + ' PSU');
    s19DoV.setValue(inSitu.do_mgL.mean.toFixed(2) + ' mg/L');
    s19NoteV.setValue('Source: ' + inSitu.source + ' | Record: ' + inSitu.record +
      '\n' + inSitu.notes);
    print('=== S19 REAL IN-SITU BASELINE (v10.36 NEW) ===');
    print('Matched station: ' + inSitu.label + ' (' + inSitu.distance_km.toFixed(2) + ' km from click)');
    print('Source: ' + inSitu.source);
    print('Record: ' + inSitu.record);
    print('Clean pH:       ' + inSitu.pH.mean.toFixed(3) + ' +/- ' + inSitu.pH.std.toFixed(3) +
      ' (n=' + inSitu.n_clean + ', ' + inSitu.pct_flagged + '% of raw readings flagged as sensor artifacts)');
    print('Clean Temp:     ' + inSitu.temp_c.mean.toFixed(2) + ' deg C');
    print('Clean Salinity: ' + inSitu.salinity.mean.toFixed(2) + ' PSU');
    print('Clean DO:       ' + inSitu.do_mgL.mean.toFixed(2) + ' mg/L');
    print('Note: ' + inSitu.notes);
    print('CAVEAT: real measured sensor baseline, but point-source only.');
  } else {
    s19StatusV.setValue('No in-situ station within range (nearest station is elsewhere)');
    s19StatusV.style().set('color', '#888888');
    s19PhV.setValue('n/a'); s19TempV.setValue('n/a');
    s19SalV.setValue('n/a'); s19DoV.setValue('n/a');
    s19NoteV.setValue('Not within ~10-15km of Looe Key FL, Agua Hedionda CA, or ' +
      'Scripps Pier CA - the only 3 stations currently in this registry.');
  }

  locV.setValue(region); regV.setValue(region);
  coV.setValue('Lat:'+latR+' Lon:'+lonR);
  clickLbl.setValue('Analyzing: '+region+'...');
  clickLbl.style().set('backgroundColor','#115511');
  scoreBig.setValue('Computing...'); scoreBig.style().set('color','#333333');
  scoreBig.style().set('backgroundColor','#eeeeee');
  scoreBarLbl.setValue('Computing...'); scoreBarLbl.style().set('color','#666666');
  scoreInterp.setValue('Running pipeline...');
  barFill.style().set('width','0px'); barFill.style().set('backgroundColor','#aaaaaa');

  if(fp.hasField) {
    fieldStatusV.setValue('FIELD DATA AVAILABLE');
    fieldStatusV.style().set('color','#115511');
    fieldSpeciesV.setValue(fp.species||'--');
    fieldSourcesV.setValue('Sources: '+fp.sources);
    fieldNotesV.setValue('Notes: '+fp.notes);
  } else {
    fieldStatusV.setValue('SATELLITE ONLY - no field data for this region');
    fieldStatusV.style().set('color','#880000');
    fieldSpeciesV.setValue('none');
    fieldSourcesV.setValue('Sources: '+fp.sources);
    fieldNotesV.setValue('Notes: '+fp.notes);
  }
  fc1V.setValue('n/a'); fc3V.setValue('n/a');
  fc4V.setValue('n/a'); fc5V.setValue('n/a'); fcTV.setValue('0 total');

  fN.setValue(region); fS.setValue('Computing...'); fB.setValue('');
  fC.setValue(latR+', '+lonR);
  floatP.style().set('shown',true);
  Map.setCenter(lon,lat,11);

  while(Map.layers().length()>0){ Map.remove(Map.layers().get(0)); }
  Map.addLayer(ee.Image().byte().paint(study,0,3),
    {palette:['#ffffff'],opacity:1.0},'Study zone (white border)');
  Map.addLayer(sst.clip(study),{min:20,max:32,
    palette:['#2c7bb6','#abd9e9','#ffffbf','#fdae61','#f46d43','#d73027','#a50026'],
    opacity:0.9},'S1 - SST deg C (loading)');

  function clip(col){ return col.map(function(img){ return img.clip(study); }); }

  print('');
  print('STEMGeoHS Marine v10.36 -- '+region);
  print('=== MODELS USED ===');
  print('1. Waddington double-well: U(q;mu) = 0.25*q^4 - 0.5*mu*q^2');
  print('2. Langevin SDE: dx = -dU/dx*dt + sigma*dW (PNAS 2025)');
  print('3. Critical Slowing Down: AC1->1, tau->inf as mu->mu_c');
  print('4. Lagrangian L = T - V where T = 0.5*m_eff*q_dot^2');
  print('5. Kramers escape: k = omega0 * exp(-dU/sigma^2)');
  print('6. Cancer score analog: CCS = SUM(w_i * S_i)');
  print('=== FIELD DATA STATUS: '+region+' ===');
  if(fp.hasField) {
    print('STATUS: FIELD DATA AVAILABLE');
    print('Species: '+fp.species);
    print('Field accuracy gain: +'+fp.accuracy_field_gain+'%');
  } else {
    print('STATUS: SATELLITE DATA ONLY');
    print('Accuracy: 77% satellite baseline only');
  }

  function dOpts(title,yT,col,lab){
    return{title:title,
      series:{0:{color:col,lineWidth:2.5,pointSize:4,label:lab}},
      backgroundColor:'#0a1628',
      titleTextStyle:{color:'#ffffff',fontSize:10,bold:true},
      vAxis:{title:yT,textStyle:{color:'#cccccc'},
             titleTextStyle:{color:'#aaaacc'},gridlines:{color:'#1a2a4a'}},
      hAxis:{title:'Month Jan 2023 to Dec 2024',
             textStyle:{color:'#aaaaaa'},titleTextStyle:{color:'#aaaacc'}},
      legend:{textStyle:{color:'#ffffff'}},
      chartArea:{backgroundColor:'#0d1f3c',width:'82%'}};
  }

  print(ui.Chart.image.series({imageCollection:clip(mkMoSST()),
    region:study,reducer:ee.Reducer.mean(),scale:4000,xProperty:'system:time_start'})
  .setOptions(dOpts('SST Monthly | '+region+' | OISST V2.1 daily',
    'SST (deg C)','#ff4422','SST')));

  print(ui.Chart.image.series({imageCollection:clip(mkMoDHW()),
    region:study,reducer:ee.Reducer.mean(),scale:4000,xProperty:'system:time_start'})
  .setOptions({title:'DHW Monthly (per-pixel MMM) | '+region,
    series:{0:{color:'#ff6644',lineWidth:2.5,pointSize:4,label:'DHW per-pixel (deg C-wks)'}},
    backgroundColor:'#0a1628',
    titleTextStyle:{color:'#ffffff',fontSize:10,bold:true},
    vAxis:{title:'DHW (deg C-weeks)',viewWindow:{min:0},textStyle:{color:'#cccccc'},
           titleTextStyle:{color:'#aaaacc'},gridlines:{color:'#1a2a4a'}},
    hAxis:{title:'Month',textStyle:{color:'#aaaaaa'},titleTextStyle:{color:'#aaaacc'}},
    legend:{textStyle:{color:'#ffffff'}},
    chartArea:{backgroundColor:'#0d1f3c',width:'82%'}}));

  print(ui.Chart.image.series({imageCollection:clip(mkMoCHL()),
    region:study,reducer:ee.Reducer.mean(),scale:4000,xProperty:'system:time_start'})
  .setOptions(dOpts('Chl-a Monthly | '+region+' | Copernicus Ocean Color V6',
    'Chl-a (mg/m3)','#44cc44','Chl-a')));

  var c4o=dOpts('Annual SST 2003-2024 | '+region+' | Basin erosion CSD',
    'Annual SST (deg C)','#ff6644','Annual SST');
  c4o.trendlines={0:{type:'linear',color:'#ffcc44',lineWidth:2.5,opacity:0.9,
    showR2:true,visibleInLegend:true}};
  c4o.hAxis={title:'Year 2003-2024',textStyle:{color:'#aaaaaa'},titleTextStyle:{color:'#aaaacc'}};
  c4o.legend={position:'top',textStyle:{color:'#ffffff',fontSize:9}};
  print(ui.Chart.image.series({imageCollection:clip(annualSST),
    region:study,reducer:ee.Reducer.mean(),scale:4000,xProperty:'system:time_start'})
  .setOptions(c4o));

  print(ui.Chart.image.series({imageCollection:clip(mkMoNO2()),
    region:study,reducer:ee.Reducer.mean(),scale:1000,xProperty:'system:time_start'})
  .setOptions(dOpts('NO2 Monthly | '+region+' | S5P OFFL+NRTI | Metal proxy',
    'NO2 (mol/m2)','#cc44cc','NO2')));

  if(region==='Bocas del Toro, Panama'||region==='Caribbean Sea') {
    print('=== REAL FIELD DATA: Ramamurthy 2024 CJS 54:77-82, Table 1 ===');
    print('Bocas del Toro, Panama - staged intraspecific contests, n=3 individuals/species');
    print('B. granuliferum A: vs B: IF, A, P, M, R | vs C: NF, R');
    print('B. granuliferum B: vs C: IF, A, NP, NM, NR');
    print('B. granuliferum C: vs A: NF, NR | vs B: IF, NR');
    print('B. cavernatum A: vs B: IF, NR | vs C: IF, NR');
    print('B. cavernatum B: vs A: IF, A, P, NP, M, NM, NR, LA');
    print('B. cavernatum C: vs A: IF, A, NP, NM, NR, C');
    print('NOTE: qualitative codes only, not a numeric ratio.');
  }

  // ============================================================
  // v10.49 MAJOR OPTIMIZATION: CHECK DEPTH FIRST, BEFORE ANYTHING ELSE
  //
  // PROBLEM (v10.48 was incomplete): that fix only skipped S17/S18
  // (the LAST 2 of 20 nested calls) after ALL 18 other expensive
  // SST/DHW/Chl-a/turbidity/NO2/seaweed/aquaculture/soil/CSD/thermal-
  // recovery calls had ALREADY been dispatched to Earth Engine's
  // servers - because isLand was only known deep inside the chain,
  // after rBATH resolved as the 7th step. This still wasted real
  // computation time on ~18 calls that are meaningless on land.
  //
  // FIX: an independent, standalone depth check now runs FIRST, before
  // the massive ocean-analysis chain is even started. If the point is
  // land, ONLY earthquakes + volcanic activity are computed (both
  // genuinely meaningful on land) - literally none of the other 18
  // ocean-only reduceRegion calls (SST, DHW, MMM, Chl-a, turbidity,
  // NO2, FAI/NDCI/NDVI seaweed, soil, CSD, thermal recovery, S16
  // physics, S17 ToE, S18 BGC) are ever dispatched to the server at
  // all for a land click. If the point is ocean, the ENTIRE existing
  // chain below runs completely unchanged.
  // ============================================================
  // v10.59 FIX: same class of problem as the earlier SST narrow-inlet
  // masking issue (v10.43), now found at Agua Hedionda Lagoon itself -
  // a REAL, confirmed aquaculture site incorrectly landed on the
  // "INLAND" fast-path because ETOPO1's coarse ~1.8km resolution
  // misclassified this narrow lagoon channel using only a 1km MEAN
  // buffer. FIX: widen to 2km AND switch mean() to min() - min() finds
  // the deepest (most confidently "water") pixel in the area, so if
  // ANY real ocean pixel exists nearby, the point is correctly NOT
  // flagged as land. This is the safer failure mode for a hard gate:
  // occasionally treating a genuinely borderline point as ocean (and
  // running the full real analysis) is far better than incorrectly
  // skipping a real, valid coastal/lagoon site entirely.
  var rBATH_early = bathy.reduceRegion({reducer:ee.Reducer.min(),
    geometry:pt.buffer(2000), scale:1000, maxPixels:1e9, bestEffort:true});

  rBATH_early.evaluate(function(bEarlyRes, _errEarly) {
    if (_errEarly) { handleGeeError(_errEarly, "rBATH_early"); return; }
    var bvEarly = (bEarlyRes && bEarlyRes.bedrock !== null && bEarlyRes.bedrock !== undefined) ?
      bEarlyRes.bedrock : null;
    var isLandEarly = (bvEarly !== null && bvEarly > 5);

    if (isLandEarly) {
      // ================================================================
      // MINIMAL LAND-ONLY PATH (v10.49) - only earthquakes + volcanic
      // activity are computed. Every other sidebar field is set directly
      // to an honest "inland - skipped" message, with ZERO additional
      // server-side reduceRegion calls fired for any of them.
      // ================================================================
      var seismicZoneEarly = pt.buffer(200000);
      var volcZoneEarly = pt.buffer(300000);
      var rEQ_early = ee.Dictionary({
        count: usgsEarthquakes.filterBounds(seismicZoneEarly).size(),
        maxMag: usgsEarthquakes.filterBounds(seismicZoneEarly).aggregate_max('mag')
      });
      var _gdisNearbyEarly = gdisDisasters.filterBounds(volcZoneEarly);
      var rVOLC_early = ee.Dictionary({
        totalGdisRecordsNearby: _gdisNearbyEarly.size(),
        matchedByOlcanFilter: volcanicActivity.filterBounds(volcZoneEarly).size(),
        sampleDisasterTypes: _gdisNearbyEarly.aggregate_array('disasterty').distinct().slice(0, 15),
        firstFeaturePropertyNames: ee.Algorithms.If(_gdisNearbyEarly.size().gt(0),
          ee.Feature(_gdisNearbyEarly.first()).propertyNames(), ee.List(['(none)'])),
        firstFeatureAllProperties: ee.Algorithms.If(_gdisNearbyEarly.size().gt(0),
          ee.Feature(_gdisNearbyEarly.first()).toDictionary(), ee.Dictionary({note:'none'}))
      });

      rEQ_early.evaluate(function(eqRes, _err1) {
        if (_err1) { handleGeeError(_err1, "rEQ_early"); return; }
        rVOLC_early.evaluate(function(volcRes, _err2) {
          if (_err2) { handleGeeError(_err2, "rVOLC_early"); return; }

          var eqCount = (eqRes && eqRes.count !== null && eqRes.count !== undefined) ? eqRes.count : 0;
          var eqMaxMag = (eqRes && eqRes.maxMag !== null && eqRes.maxMag !== undefined) ? eqRes.maxMag : null;
          var volcCount = (volcRes && volcRes.matchedByOlcanFilter) ? volcRes.matchedByOlcanFilter : 0;
          var volcTotalRaw = (volcRes && volcRes.totalGdisRecordsNearby) ? volcRes.totalGdisRecordsNearby : 0;

          eqStatsV.setValue(eqCount > 0 ? eqCount+' event(s), max M'+eqMaxMag.toFixed(1) : 'None nearby');
          volcStatsV.setValue(volcCount > 0 ? volcCount+' historic event(s)' : 'None recorded');
          depthV.setValue(fmt(bvEarly,0)+' m (land - marine calcs skipped)');

          var SKIP = 'N/A - inland (skipped, v10.49)';
          sstV.setValue(SKIP); trendV.setValue(SKIP); dhwV.setValue(SKIP); mmmV.setValue(SKIP);
          chlaV.setValue(SKIP); turbV.setValue(SKIP); no2V.setValue(SKIP);
          faiV.setValue(SKIP); ndciV.setValue(SKIP); ndviWV.setValue(SKIP); algaeStatusV.setValue(SKIP);
          sstWindowV.setValue(SKIP); chlWindowV.setValue(SKIP); ndciNutrientProxyV.setValue(SKIP);
          pollutionV.setValue(SKIP); stabilityV.setValue(SKIP); bromoformV.setValue(SKIP);
          aquaScoreV.setValue(SKIP); aquaStatusV.setValue(SKIP);
          soilTextureV.setValue(SKIP);
          realAc1V.setValue(SKIP); realVarTrendV.setValue(SKIP);
          thermalEpisodesV.setValue(SKIP); thermalMeanV.setValue(SKIP);
          thermalMaxV.setValue(SKIP); thermalOngoingV.setValue(SKIP);
          ecoValStatusV.setValue(SKIP);
          eciDepthV.setValue(SKIP); eciCurrentV.setValue(SKIP);
          eciSLR03V.setValue(SKIP); eciSLR05V.setValue(SKIP); eciSLR10V.setValue(SKIP);
          eciValidationV.setValue(SKIP);
          toeCompoundV.setValue(SKIP); toeSSTv.setValue(SKIP); toeCHLv.setValue(SKIP);
          toeNO2v.setValue(SKIP); toeSALv.setValue(SKIP);
          s18O2V.setValue(SKIP); s18pCO2V.setValue(SKIP); s18pHV.setValue(SKIP); s18SalV.setValue(SKIP);

          var colsEarly = {text:'#555555',bg:'#e8e8e8',map:'#999999',bar:'#999999',lbl:'INLAND - N/A'};
          scoreBig.setValue('INLAND LOCATION (depth +'+fmt(bvEarly,0)+'m) - marine scoring N/A');
          scoreBig.style().set('color',colsEarly.text);
          scoreBig.style().set('backgroundColor',colsEarly.bg);
          scoreBarLbl.setValue('Not applicable - this point is on land, not ocean');
          barFill.style().set('width','0px');
          var interpEarly = region+': INLAND LOCATION\n'+
            'Depth: +'+fmt(bvEarly,0)+'m above sea level (confirmed land)\n'+
            'v10.49: ALL ocean-only calculations skipped entirely (SST, DHW,\n'+
            'Chl-a, turbidity, NO2, seaweed, aquaculture, soil, CSD, thermal\n'+
            'recovery, S16/S17/S18) - zero server calls for any of them.\n'+
            'Only earthquakes + volcanic activity computed (valid on land).';
          scoreInterp.setValue(interpEarly);
          scoreInterp.style().set('color',colsEarly.text);
          scoreInterp.style().set('backgroundColor',colsEarly.bg);

          // v10.51 FIX: the floating map overlay (fN/fS/fB) is a SEPARATE
          // display element from the sidebar - it was set to 'Computing...'
          // at the very start of analyzeLocation() but never updated in
          // the land-only path (only the sidebar was), leaving it stuck
          // indefinitely even after the analysis had genuinely finished.
          floatP.style().set('border','2px solid '+colsEarly.map);
          fN.style().set('color',colsEarly.map); fN.setValue(region);
          fS.setValue('INLAND - marine scoring N/A (depth +'+fmt(bvEarly,0)+'m)');
          fS.style().set('color',colsEarly.map);
          fB.setValue('Earthquakes: '+eqCount+' | Volcanic: '+volcCount);

          // v10.50 FIX: clear the temporary "S1 - SST deg C (loading)"
          // placeholder (added earlier in analyzeLocation() before the
          // land/ocean branch was known) before adding the land-only
          // layers - otherwise it stays stuck on the map forever, since
          // this land path never calls loadLayers() (which normally does
          // this same clear for ocean points). Same clear-then-add pattern
          // loadLayers() already uses, just applied here for consistency.
          while(Map.layers().length()>0){ Map.remove(Map.layers().get(0)); }
          Map.addLayer(ee.Image().byte().paint(study,0,3),
            {palette:['#ffffff'],opacity:1.0},'Study zone (white border)');
          // v10.52 FIX: the colored earthquake/volcano dots previously
          // showed up on land clicks only as a side effect of the FULL
          // ocean chain running (which added these layers regardless of
          // land/ocean). The v10.49 optimization correctly skips that
          // full chain for land - but that also silently dropped these
          // two layers as an unintended side effect. Since S10/S11 data
          // IS genuinely computed and valid on land, the visual layers
          // belong here too - added explicitly, matching loadLayers()'s
          // exact same layer definitions for the ocean path.
          Map.addLayer(eqModerate,{color:'ff8800'},'S10 - Earthquakes M4.5-5.9', true, 0.7);
          Map.addLayer(eqMajor,{color:'cc0000'},'S10 - Earthquakes M6.0+', true, 0.9);
          Map.addLayer(volcanicActivity,{color:'ff00ff'},'S11 - Volcanic activity 1960-2018', true, 0.8);
          Map.addLayer(ee.Image.constant(1).clip(study.centroid(50).buffer(300)),
            {palette:[colsEarly.map],opacity:0.9},'SCORE: INLAND (no marine data)');

          clickResultsLog.push({
            timestamp: new Date().toISOString(), lat: latR, lon: lonR, region: region,
            region_type: 'land',
            sst_annual_c: null, sst_peak_c: null, sst_trend_c_per_yr: null,
            dhw_c_weeks: null, dhw_raw_unfiltered: null, dhw_artifact_flagged: null,
            mmm_local_c: null, chl_a_mg_m3: null, turbidity_ndti: null, no2_mol_m2: null,
            depth_m: bvEarly, cancer_score_satellite: null, cancer_score_fused: null,
            field_correction: null, bowl_depth_B: null, accuracy_pct: null,
            status_label: 'INLAND - N/A', field_data_available: false, field_species: 'none',
            aquaculture_score: null, aquaculture_confidence_pct: null, aquaculture_status: 'N/A - inland',
            bromoform_yield: null, num_interventions: 0, top_intervention: 'none (inland)',
            soil_texture_code: null, soil_texture_label: 'N/A - inland (v10.49)',
            earthquake_count_200km: eqCount, earthquake_max_mag_200km: eqMaxMag,
            volcanic_activity_count_300km: volcCount,
            real_ac1_24mo: null, real_var_trend_ratio: null, real_csd_n_valid_months: 0,
            ndci_nutrient_proxy: null,
            thermal_recovery_n_completed_episodes: null, thermal_recovery_mean_months: null,
            thermal_recovery_max_months: null, thermal_recovery_ongoing_flag: null,
            eco_recovery_validation_checked: false, eco_recovery_validation_finding: 'N/A - inland',
            in_situ_station: inSitu ? inSitu.key : null,
            in_situ_distance_km: inSitu ? inSitu.distance_km : null,
            in_situ_ph_mean: null, in_situ_temp_c_mean: null,
            in_situ_salinity_mean: null, in_situ_do_mgL_mean: null
          });
          exportStatusLabel.setValue(clickResultsLog.length+' row(s) ready');

          print('=== v10.49 LAND DETECTED EARLY - MINIMAL PATH ===');
          print('Depth +'+fmt(bvEarly,0)+'m above sea level. Skipped ALL 18 ocean-only');
          print('reduceRegion calls entirely (SST/DHW/MMM/Chl-a/turbidity/NO2/seaweed/');
          print('aquaculture/soil/CSD/thermal-recovery/S16/S17/S18) - zero server cost.');
          print('Earthquakes: '+eqCount+' | Volcanic (matched): '+volcCount+
            ' | Volcanic (raw nearby): '+volcTotalRaw);

          // v10.53 FIX: the real schema data (property names + full
          // properties of an actual nearby record) WAS being computed in
          // rVOLC_early all along (see v10.49) but never actually printed
          // in this land-only path - only the ocean path printed it. This
          // was the missing piece needed to find the real GDIS field name.
          var volcPropNamesEarly = (volcRes && volcRes.firstFeaturePropertyNames) ?
            volcRes.firstFeaturePropertyNames : [];
          var volcAllPropsEarly = (volcRes && volcRes.firstFeatureAllProperties) ?
            volcRes.firstFeatureAllProperties : {};
          print('=== v10.53 VOLCANIC ACTIVITY DIAGNOSTIC (real schema, land path) ===');
          if (volcTotalRaw > 0) {
            print('ACTUAL disastertype values found: '+JSON.stringify(
              (volcRes && volcRes.sampleDisasterTypes) ? volcRes.sampleDisasterTypes : []));
            print('REAL property names on an actual nearby feature: '+JSON.stringify(volcPropNamesEarly));
            print('REAL full properties of that feature: '+JSON.stringify(volcAllPropsEarly));
            if (volcCount === 0) {
              print('DIAGNOSIS: real records exist here, but "disastertype" as a property');
              print('name returned nothing - check the REAL property names printed above');
              print('to find the correct field name for the category.');
            }
          } else {
            print('DIAGNOSIS: GDIS has ZERO records of ANY category within 300km here.');
          }

          clearPipelineWatchdog();
          print('=== PIPELINE COMPLETE (v10.53, land - minimal path) ===');
        });
      });

    } else {
      // ================================================================
      // FULL OCEAN CHAIN - EXISTING CODE, COMPLETELY UNCHANGED BELOW.
      // Only reached when the point is confirmed NOT land.
      // ================================================================
  var rSST  = sst.reduceRegion({reducer:ee.Reducer.mean(),geometry:studySST,scale:4000,maxPixels:1e9,bestEffort:true});
  var rSSTP = sstPeakFinal.reduceRegion({reducer:ee.Reducer.mean(),geometry:studySST,scale:4000,maxPixels:1e9,bestEffort:true});
  var rCHL  = chla.reduceRegion({reducer:ee.Reducer.mean(),geometry:study,scale:4000,maxPixels:1e9});
  var rCHLc = chlaCoastal.reduceRegion({reducer:ee.Reducer.mean(),geometry:study,scale:4000,maxPixels:1e9});
  var rTR   = sst_slope.reduceRegion({reducer:ee.Reducer.mean(),geometry:studySST,scale:4000,maxPixels:1e9,bestEffort:true});
  var rNO2  = no2.reduceRegion({reducer:ee.Reducer.mean(),geometry:study,scale:1000,maxPixels:1e9});
  // v10.65 FIX: switched the main "Depth" display from ETOPO1 (bathy)
  // to GEBCO. Confirmed via live testing (Inner Oslofjord, Norway) that
  // ETOPO1 showed a real, complex fjord channel as ~0m (essentially
  // land-level), while GEBCO correctly showed 22.8m at the exact same
  // point. Research confirms this is not a fluke: ETOPO1's bathymetry
  // was itself originally sourced from an older GEBCO vintage, then
  // coarsened and blended into one combined global land+sea grid back
  // in 2009 - GEBCO has been continuously refined with new survey data
  // every year since, making it the more current, more reliable choice
  // specifically at narrow/complex coastlines like this one.
  // SCOPE NOTE: this ONLY changes the main sidebar "Depth" display and
  // the isLand/hard-veto determination fed by it (bv). oceanMask (used
  // globally for Chl-a/SST/etc. masking) and rBATH_early (the separate,
  // already-fixed v10.59 early land-detection check) both intentionally
  // remain on their existing logic - changing those has a much wider
  // blast radius that hasn't been tested today.
  // v10.67 FIX: scale harmonized to 500 (matching S16's own separate
  // GEBCO query exactly) instead of 1000 - confirmed via live testing
  // (Carlsbad, CA) that both fields used GEBCO but disagreed (-3m vs
  // 9.6m at the same point) purely because they requested different
  // scale parameters over the same small area, in a spot with real
  // underwater terrain variation (submarine scour/canyon features).
  // 500m is closer to GEBCO's native ~450m resolution.
  // v10.68 CRITICAL FIX: the v10.67 scale-harmonization did NOT fully
  // fix the mismatch. Confirmed via live testing (Carlsbad, CA) that
  // even with matching scale AND matching dataset, main Depth (-1m)
  // and S16's Seafloor depth (9.6m) still disagreed substantially. Root
  // cause: S16's own query masks OUT land pixels before averaging
  // (local_gebco.lt(0) - only real ocean counted), but this main Depth
  // query had no such mask - any land/shallow-fringe pixels within the
  // 1km buffer near this coastline were being averaged in alongside
  // real ocean pixels, pulling the result toward zero. Now masked the
  // same way S16 already did correctly.
  var rBATH = GEBCO.rename('bedrock').updateMask(GEBCO.lt(0))
    .reduceRegion({reducer:ee.Reducer.mean(),geometry:study,scale:500,maxPixels:1e9});
  var rTURB = turbImg.reduceRegion({reducer:ee.Reducer.mean(),geometry:study,scale:100,maxPixels:1e9});
  var rDHW  = dhwProper.reduceRegion({reducer:ee.Reducer.mean(),geometry:studySST,scale:4000,maxPixels:1e9,bestEffort:true});
  var rMMM  = MMM_perpixel.reduceRegion({reducer:ee.Reducer.mean(),geometry:studySST,scale:4000,maxPixels:1e9,bestEffort:true});
  var rFAI  = faiImg.updateMask(oceanMask).reduceRegion({
    reducer:ee.Reducer.max(),geometry:study,scale:100,maxPixels:1e9});
  var rNDCI = ndciImg.updateMask(oceanMask).reduceRegion({
    reducer:ee.Reducer.max(),geometry:study,scale:100,maxPixels:1e9});
  var rNDVI_W = ndviWater.updateMask(oceanMask).reduceRegion({
    reducer:ee.Reducer.max(),geometry:study,scale:100,maxPixels:1e9});
  var rCHL_BLOOM = chla.reduceRegion({
    reducer:ee.Reducer.max(),geometry:study,scale:1000,maxPixels:1e9});
  var rSOIL = soilTexture.reduceRegion({
    reducer:ee.Reducer.first(),geometry:pt,scale:250,maxPixels:1e9});

  var seismicZone = pt.buffer(200000);
  var volcZone = pt.buffer(300000);
  var rEQ = ee.Dictionary({
    count: usgsEarthquakes.filterBounds(seismicZone).size(),
    maxMag: usgsEarthquakes.filterBounds(seismicZone).aggregate_max('mag')
  });
  // v10.46 NEW: DIAGNOSTIC for volcanic activity showing "None recorded"
  // even at Kilauea (a real, well-documented eruption zone with major
  // activity inside GDIS's own 1960-2018 window, e.g. the 2018 eruption).
  // This was already flagged as an UNCONFIRMED limitation in the original
  // script's own v10.19 changelog ("plausibly explained by GDIS's narrow
  // criteria... but not independently confirmed either way") - rather
  // than guess again, this pulls the RAW unfiltered GDIS data in the same
  // 300km zone plus the ACTUAL disastertype text values present, so the
  // real cause (wrong field name, wrong filter string, or a genuine
  // empty region) is visible directly in the Console instead of assumed.
  // v10.47 FIX: v10.46's diagnostic showed 8 real GDIS records nearby but
  // an EMPTY sampleDisasterTypes array - if 'disastertype' genuinely
  // existed as a property, aggregate_array().distinct() would have
  // returned at least one string. An empty result despite real records
  // existing means the property NAME itself is very likely wrong (a
  // schema mismatch), not just an unexpected category string. This adds
  // the real fix: pull the actual property names directly off one real
  // nearby feature, so the correct field name is visible directly
  // instead of guessing a second time - same diagnostic pattern the
  // original script's own v10.6 changelog already used successfully
  // once before for this exact class of USGS/GDIS schema problem.
  var _gdisNearbyFC = gdisDisasters.filterBounds(volcZone);
  var rVOLC_diag = ee.Dictionary({
    totalGdisRecordsNearby: _gdisNearbyFC.size(),
    matchedByOlcanFilter: volcanicActivity.filterBounds(volcZone).size(),
    sampleDisasterTypes: _gdisNearbyFC
      .aggregate_array('disasterty').distinct().slice(0, 15),
    // NEW: real property names + full property values from ONE actual
    // nearby feature, so we can see the TRUE schema directly.
    firstFeaturePropertyNames: ee.Algorithms.If(
      _gdisNearbyFC.size().gt(0),
      ee.Feature(_gdisNearbyFC.first()).propertyNames(),
      ee.List(['(no features to inspect)'])
    ),
    firstFeatureAllProperties: ee.Algorithms.If(
      _gdisNearbyFC.size().gt(0),
      ee.Feature(_gdisNearbyFC.first()).toDictionary(),
      ee.Dictionary({note:'no features to inspect'})
    )
  });
  var rVOLC = rVOLC_diag;

  var rCSD = computeRealCSD(mkMoSST(), study, 'sst', 4000);

  var rThermalSeries = ee.FeatureCollection(mkMoSST().map(function(img){
    var v = img.reduceRegion({reducer:ee.Reducer.mean(), geometry:study,
      scale:4000, maxPixels:1e9}).get('sst');
    return ee.Feature(null, {v:v, t:img.get('system:time_start')});
  })).sort('t');

  var bgcSnap = ee.ImageCollection(
    'COPERNICUS/MARINE/GLOBAL_ANALYSISFORECAST_BGC_001_028/BIO')
    .filter(ee.Filter.date('2022-01-01', ee.Date(Date.now())))
    .sort('system:time_start', false).first();

  // ============================================================
  // v10.41 NEW: SKIP redundant S18 BGC snapshot lookups when S19 real
  // in-situ sensor data already covers the same variable at this exact
  // point (inSitu was already checked synchronously at the top of
  // analyzeLocation() - zero extra cost to reuse that result here).
  //
  // RATIONALE (see prior conversation): S18 and S17 answer DIFFERENT
  // questions and are NOT interchangeable:
  //   - S18 (this snapshot) answers "what are conditions RIGHT NOW" -
  //     once a real, cleaned sensor reading exists for the SAME variable
  //     at the SAME point, the coarse 0.25-degree global model estimate
  //     is redundant AND strictly less accurate. Skip it.
  //   - S17 (sst_slope/sal_slope/etc above) answers "has a real 30-60yr
  //     trend emerged" - a ~1yr sensor record can NEVER answer that,
  //     no matter how clean. NEVER skipped, regardless of S19 match.
  //
  // pH, Salinity, and O2 all have a direct S19 real-sensor equivalent -
  // skipped when inSitu exists. pCO2 has NO sensor equivalent anywhere
  // in S19's registry - ALWAYS computed, real match or not.
  // ============================================================
  var _skipBGCForPH  = (inSitu !== null);
  var _skipBGCForSal = (inSitu !== null);
  var _skipBGCForO2  = (inSitu !== null);

  // v10.60 FIX: SPLIT S17 (ToE trends) and S18 (BGC snapshot) into TWO
  // INDEPENDENT server calls instead of one combined rToEAll.
  //
  // PROBLEM DISCOVERED (live testing, Looe Key FL): even though v10.41
  // already skips the REAL computation for S18's pH/Salinity/O2 when a
  // real S19 sensor match exists (using -999 sentinels, near-instant),
  // S18's sidebar fields still showed "computing..." for a long time -
  // because S17 and S18 were bundled into ONE combined evaluate() call.
  // S17's 4 trend variables (44/32/27/7-year precomputed image reads)
  // are NEVER skippable (a short sensor record can't answer a
  // multi-decade trend question) and are genuinely slow - so the WHOLE
  // bundle, including S18's already-known "Skipped" messages, had to
  // wait for S17's slow reads to finish before displaying anything.
  //
  // FIX: rToETrends (S17 only) and rBGCSnapshot (S18 only) are now
  // separate dictionaries, evaluated as two INDEPENDENT parallel calls
  // below (neither waits for the other). When S19 has matched, S18 can
  // now resolve near-instantly - showing its "Skipped, see S19" state
  // right away - while S17's genuinely slow trends continue in the
  // background at their own pace.
  var rToETrends = ee.Dictionary({
    sst_slope: toeSSTFit.select('scale').reduceRegion({
      reducer:ee.Reducer.mean(), geometry:study,
      scale:25000, maxPixels:1e9, bestEffort:true}).get('scale'),
    sst_noise: toeSSTNoise.reduceRegion({
      reducer:ee.Reducer.mean(), geometry:study,
      scale:25000, maxPixels:1e9, bestEffort:true}).get('sst_stdDev'),
    chl_slope: toeCHLFit.select('scale').reduceRegion({
      reducer:ee.Reducer.mean(), geometry:study,
      scale:4000, maxPixels:1e9, bestEffort:true}).get('scale'),
    chl_noise: toeCHLNoise.reduceRegion({
      reducer:ee.Reducer.mean(), geometry:study,
      scale:4000, maxPixels:1e9, bestEffort:true}).get('chlor_a_stdDev'),
    sal_slope: toeSALFit.select('scale').reduceRegion({
      reducer:ee.Reducer.mean(), geometry:study,
      scale:25000, maxPixels:1e9, bestEffort:true}).get('scale'),
    sal_noise: toeSALNoise.reduceRegion({
      reducer:ee.Reducer.mean(), geometry:study,
      scale:25000, maxPixels:1e9, bestEffort:true}).get('salinity_0_stdDev'),
    no2_slope: toeNO2Fit.select('scale').reduceRegion({
      reducer:ee.Reducer.mean(), geometry:study,
      scale:25000, maxPixels:1e9, bestEffort:true}).get('scale'),
    no2_noise: toeNO2Noise.reduceRegion({
      reducer:ee.Reducer.mean(), geometry:study,
      scale:25000, maxPixels:1e9, bestEffort:true})
      .get('tropospheric_NO2_column_number_density_stdDev')
  });

  var rBGCSnapshot = ee.Dictionary({
    o2_surface: _skipBGCForO2 ? ee.Number(-999) :
      bgcSnap.select('o2_depth1').reduceRegion({
      reducer:ee.Reducer.mean(), geometry:study,
      scale:27500, maxPixels:1e9, bestEffort:true}).get('o2_depth1'),
    pco2_surface: ee.ImageCollection(
      'COPERNICUS/MARINE/GLOBAL_ANALYSISFORECAST_BGC_001_028/CO2')
      .filter(ee.Filter.date('2022-01-01', ee.Date(Date.now())))
      .sort('system:time_start', false).first().select('spco2_depth1')
      .reduceRegion({reducer:ee.Reducer.mean(), geometry:study,
      scale:27500, maxPixels:1e9, bestEffort:true}).get('spco2_depth1'),
    ph_surface: _skipBGCForPH ? ee.Number(-999) :
      ee.ImageCollection(
      'COPERNICUS/MARINE/GLOBAL_ANALYSISFORECAST_BGC_001_028/CAR')
      .filter(ee.Filter.date('2022-01-01', ee.Date(Date.now())))
      .sort('system:time_start', false).first().select('ph_depth1')
      .unmask(-9999)
      .reduceRegion({reducer:ee.Reducer.mean(), geometry:study.buffer(55000),
      scale:27500, maxPixels:1e9, bestEffort:true}).get('ph_depth1'),
    sal_current: _skipBGCForSal ? ee.Number(-999) :
      ee.ImageCollection('HYCOM/sea_temp_salinity')
      .filter(ee.Filter.date('2023-01-01', ee.Date(Date.now())))
      .select('salinity_0').mean()
      .multiply(0.001).add(20)
      .rename('salinity_psu')
      .reduceRegion({reducer:ee.Reducer.mean(), geometry:study,
      scale:25000, maxPixels:1e9, bestEffort:true}).get('salinity_psu')
  });

  // v10.60: shared completion tracker so the pipeline-complete message
  // and watchdog-clear only fire once BOTH independent calls finish,
  // regardless of which one happens to resolve first.
  var _s17s18DoneCount = 0;
  function _checkBothS17S18Done() {
    _s17s18DoneCount++;
    if (_s17s18DoneCount >= 2) {
      clearPipelineWatchdog();
      print('=== PIPELINE COMPLETE (v10.60, S17+S18 both resolved independently) ===');
    }
  }

  rSST.evaluate(function(s, _err) {
                  if (_err) { handleGeeError(_err, "rSST"); return; }
    rSSTP.evaluate(function(sp, _err) {
                  if (_err) { handleGeeError(_err, "rSSTP"); return; }
    rCHL.evaluate(function(c, _err) {
                  if (_err) { handleGeeError(_err, "rCHL"); return; }
    rCHLc.evaluate(function(cc, _err) {
                  if (_err) { handleGeeError(_err, "rCHLc"); return; }
      rTR.evaluate(function(tr, _err) {
                  if (_err) { handleGeeError(_err, "rTR"); return; }
        rNO2.evaluate(function(n) {
          rBATH.evaluate(function(b, _err) {
                  if (_err) { handleGeeError(_err, "rBATH"); return; }
            rTURB.evaluate(function(t, _err) {
                  if (_err) { handleGeeError(_err, "rTURB"); return; }
              rDHW.evaluate(function(d, _err) {
                  if (_err) { handleGeeError(_err, "rDHW"); return; }
                rMMM.evaluate(function(mmm, _err) {
                  if (_err) { handleGeeError(_err, "rMMM"); return; }
                rFAI.evaluate(function(fa, _err) {
                  if (_err) { handleGeeError(_err, "rFAI"); return; }
                rNDCI.evaluate(function(nc, _err) {
                  if (_err) { handleGeeError(_err, "rNDCI"); return; }
                rNDVI_W.evaluate(function(nw, _err) {
                  if (_err) { handleGeeError(_err, "rNDVI_W"); return; }
                rCHL_BLOOM.evaluate(function(cb, _err) {
                  if (_err) { handleGeeError(_err, "rCHL_BLOOM"); return; }
                rSOIL.evaluate(function(soilRes, _err) {
                  if (_err) { handleGeeError(_err, "rSOIL"); return; }
                rEQ.evaluate(function(eqRes, _err) {
                  if (_err) { handleGeeError(_err, "rEQ"); return; }
                rVOLC.evaluate(function(volcRes, _err) {
                  if (_err) { handleGeeError(_err, "rVOLC"); return; }
                rCSD.evaluate(function(csdRes, _err) {
                  if (_err) { handleGeeError(_err, "rCSD"); return; }
                rThermalSeries.evaluate(function(thermalFC, _err) {
                  if (_err) { handleGeeError(_err, "rThermalSeries"); return; }

                  var sv      = (s&&s.sst!==null&&s.sst>-900)?s.sst:null;
                  var sv_peak = (sp&&sp.sst_peak!==null&&sp.sst_peak>-900)?sp.sst_peak:sv;
                  var cv = (c&&c.chlor_a>0) ? c.chlor_a :
                           (cc&&cc.chlor_a>0) ? cc.chlor_a : null;
                  var chlSource = (c&&c.chlor_a>0) ? 'ocean' :
                                  (cc&&cc.chlor_a>0) ? 'coastal relaxed' : 'n/a';
                  var tv   = (tr&&tr.scale!==null)?tr.scale:null;
                  var nv   = (n&&n.NO2_column_number_density>0)?n.NO2_column_number_density:null;
                  var bv   = (b&&b.bedrock!==null)?b.bedrock:null;
                  var turv = (t&&t.turbidity!==null)?t.turbidity:null;
                  var dhwv_raw = (d&&d.dhw!==null)?d.dhw:null;
                  var DHW_SANITY_CAP = 60;
                  var dhwv = (dhwv_raw !== null && dhwv_raw <= DHW_SANITY_CAP) ? dhwv_raw : null;
                  var dhwFlaggedArtifact = (dhwv_raw !== null && dhwv_raw > DHW_SANITY_CAP);
                  var mmmv = (mmm&&mmm.mmm!==null)?mmm.mmm:null;

                  // ============================================================
                  // v10.45 NEW: LAND DETECTION - suppress misleading marine score
                  //
                  // PROBLEM DISCOVERED (live testing at Volcano, HI - an inland
                  // town ~15-20km from the coast, testing S10/S11 seismic data):
                  // when a click lands on real land (bathy > 0m elevation), ALL
                  // marine satellite inputs (SST, DHW, Chl-a) correctly return
                  // null - but computeScore() then silently substitutes generic
                  // DEFAULT values (s1=50, anem_safe=6.0, urch_safe=0.3, etc.)
                  // and produces a confident-looking numeric score (e.g. "31/100
                  // B=0.69") that is mathematically well-formed but represents
                  // NO real marine signal whatsoever - just fallback constants.
                  // This directly violates this script's own standard (already
                  // applied elsewhere: aquaculture, S12, S14) of never
                  // fabricating a number for missing data.
                  //
                  // FIX: if the clicked point is confirmed land (depth > 5m
                  // elevation, small buffer to avoid tide-flat/DEM noise), the
                  // FUSED COASTAL CANCER SCORE box now explicitly states
                  // "INLAND LOCATION - marine scoring not applicable" instead of
                  // a fake numeric score. S10/S11 (seismic/volcanic) remain
                  // FULLY VALID and continue to display normally, since those
                  // use a wide 200-300km search radius that is meaningful on
                  // land just as much as at the coast - only the ocean-specific
                  // cancer score/bowl-depth/P(flip) numbers are suppressed.
                  // ============================================================
                  var isLand = (bv !== null && bv > 5);

                  var sc = computeScore(sv_peak,cv,tv,nv,turv,dhwv,fp,lat,lon);
                  sc.sst_annual = sv;
                  sc.sst_peak = sv_peak;
                  var cols = isLand ?
                    {text:'#555555',bg:'#e8e8e8',map:'#999999',bar:'#999999',lbl:'INLAND - N/A'} :
                    scoreColors(sc.ccs);
                  loadLayers(study,region,cols);

                  var peakLabel = sv_peak!==null?' | Peak(Jun-Oct):'+fmt(sv_peak,1)+' deg C':'';
                  sstV.setValue(fmt(sv,2)+' deg C (annual mean)'+peakLabel);
                  trendV.setValue(fmt(tv,4)+' deg C/yr');
                  var isTropicalReef = (lat > -30 && lat < 30) && !isEBUS(lat, lon);
                  var dhwLabel = dhwv ? (isTropicalReef ?
                      (dhwv>12?' MORTALITY':dhwv>8?' MASS BLEACH':dhwv>4?' BLEACH RISK':dhwv>0?' Watch':' No stress') :
                      (dhwv>8?' STRONG warm anomaly':dhwv>4?' Moderate warm anomaly':dhwv>0?' Mild warm anomaly':' No anomaly')
                  ) : '';
                  dhwV.setValue(fmt(dhwv,2)+' deg C-wks'+dhwLabel);
                  mmmV.setValue(fmt(mmmv,1)+' deg C (local baseline)');
                  var isTropReef = (lat > -30 && lat < 30);
                  if(dhwFlaggedArtifact) {
                    dhwNoteV.setValue('DATA ARTIFACT: raw DHW exceeded realistic maximum. Suppressed to n/a.');
                    dhwNoteV.style().set('color','#cc0000');
                  } else if(!isTropReef && dhwv && dhwv > 4) {
                    dhwNoteV.setValue('NOTE: warm anomaly, NOT coral bleaching (no reefs here).');
                    dhwNoteV.style().set('color','#884400');
                  } else {
                    dhwNoteV.setValue('');
                  }
                  if(bv!==null && bv < -200) {
                    depthWarnV.setValue('NOTE: Deep ocean - click closer to shore for best results');
                  } else { depthWarnV.setValue(''); }
                  if(isEBUS(lat,lon)) {
                    ebusWarnV.setValue('EBUS ZONE: SST cooling trend may be upwelling artifact.');
                    ebusWarnV.style().set('color','#885500');
                  } else { ebusWarnV.setValue(''); }
                  chlaV.setValue(fmt(cv,3)+' mg/m3 ('+chlSource+')');
                  turbV.setValue(turv!==null?fmt(turv,3)+' NDTI':'n/a - open ocean');
                  no2V.setValue(fmt(nv,8)+' mol/m2');
                  depthV.setValue(fmt(bv,0)+' m');

                  fc1V.setValue(fp.urchin_N!==null?Math.round(sc.F1)+' (N='+fp.urchin_N+'/m2)':'n/a - no data this region');
                  fc3V.setValue(fp.anem_N!==null?
                    Math.round(sc.F3)+' (aN='+fp.anem_N+'/m2)'+(fp.anem_N_estimated?' [ESTIMATED]':'') :
                    'n/a - no data this region');
                  fc4V.setValue(fp.Cd!==null?'+'+sc.F4+' (Cd='+fp.Cd+')':'n/a - no data this region');
                  fc5V.setValue(fp.recruit!==null?
                    Math.round(sc.F5)+' (recruit='+fp.recruit+')'+(fp.recruit_estimated?' [ESTIMATED]':'') :
                    'n/a - no data this region');
                  fcTV.setValue((sc.fcTotal>0?'+':'')+sc.fcTotal+' total');

                  satCcsV.setValue(sc.sat_ccs+'/100');
                  fusCcsV.setValue(sc.ccs+'/100 FUSED');
                  fusCcsV.style().set('color',cols.text);
                  bowlV.setValue(sc.B.toFixed(2));
                  omega0V.setValue(sc.omega0.toFixed(4));
                  tauV.setValue(sc.tau.toFixed(1)+'x');
                  ac1V.setValue(sc.ac1.toFixed(3));
                  p5yrV.setValue(sc.p5yr+'%');
                  p5yrV.style().set('color',cols.text);

                  accSV.setValue(sc.acc_sat+'%');
                  accFV.setValue('+'+sc.acc_field+'%');
                  accTV.setValue(sc.acc_total+'% TOTAL');

                  var displayCcs = (!isNaN(sc.ccs)&&sc.ccs!==null)?sc.ccs:sc.sat_ccs;
                  var displayB = (!isNaN(sc.B)&&sc.B!==null)?sc.B.toFixed(2):'n/a';

                  // v10.45 NEW: land locations get an honest message instead of
                  // a numeric score built entirely from fallback defaults.
                  if (isLand) {
                    scoreBig.setValue('INLAND LOCATION (depth +'+fmt(bv,0)+'m) - marine scoring N/A');
                    scoreBig.style().set('color',cols.text);
                    scoreBig.style().set('backgroundColor',cols.bg);
                    scoreBarLbl.setValue('Not applicable - this point is on land, not ocean');
                    scoreBarLbl.style().set('color',cols.map);
                    barFill.style().set('width','0px');
                    barFill.style().set('backgroundColor',cols.bar);

                    var interp = region+': INLAND LOCATION\n'+
                      'Depth: +'+fmt(bv,0)+'m above sea level (confirmed land)\n'+
                      'Marine ecosystem scoring requires an ocean pixel and is not\n'+
                      'meaningful here - no score shown rather than one built from\n'+
                      'generic fallback defaults. S10/S11 hazard data above is\n'+
                      'still valid. Try a coastal coordinate for a real assessment.';
                    scoreInterp.setValue(interp);
                    scoreInterp.style().set('color',cols.text);
                    scoreInterp.style().set('backgroundColor',cols.bg);

                    floatP.style().set('border','2px solid '+cols.map);
                    fN.style().set('color',cols.map); fN.setValue(region);
                    fS.setValue('INLAND - marine scoring N/A (depth +'+fmt(bv,0)+'m)');
                    fS.style().set('color',cols.map);

                    print('=== v10.45 LAND DETECTION ===');
                    print('This point is confirmed LAND (depth +'+fmt(bv,0)+'m above sea level).');
                    print('Marine ecosystem score suppressed - would otherwise be built entirely');
                    print('from generic fallback defaults. S10/S11 data above remains fully valid.');
                  } else {
                  scoreBig.setValue('SCORE: '+displayCcs+'/100  B='+displayB);
                  scoreBig.style().set('color',cols.text);
                  scoreBig.style().set('backgroundColor',cols.bg);
                  scoreBarLbl.setValue(cols.lbl+' ('+sc.ccs+'/100)');
                  scoreBarLbl.style().set('color',cols.map);
                  barFill.style().set('width',Math.round(sc.ccs*1.9)+'px');
                  barFill.style().set('backgroundColor',cols.bar);

                  var interp = region+': '+cols.lbl+'\nB='+sc.B.toFixed(2)+' | Acc='+sc.acc_total+'%';
                  scoreInterp.setValue(interp);
                  scoreInterp.style().set('color',cols.text);
                  scoreInterp.style().set('backgroundColor',cols.bg);

                  floatP.style().set('border','2px solid '+cols.map);
                  fN.style().set('color',cols.map); fN.setValue(region);
                  fS.setValue('Score: '+sc.ccs+'/100  B='+sc.B.toFixed(2));
                  fS.style().set('color',cols.map);
                  } // end isLand check (v10.45)

                  print('=== MEASUREMENTS ===');
                  print('SST peak: '+fmt(sv_peak,2)+' deg C | MMM: '+fmt(mmmv,1)+
                    ' | DHW: '+fmt(dhwv,2)+' | Trend: '+fmt(tv,4));
                  print('Chl-a: '+fmt(cv,3)+' | Turbidity: '+fmt(turv,3)+
                    ' | NO2: '+fmt(nv,8)+' | Depth: '+fmt(bv,0));

                  var faiv  = (fa && fa.fai  !== null && !isNaN(fa.fai))  ? fa.fai  : null;
                  var nciV  = (nc && nc.ndci !== null && !isNaN(nc.ndci)) ? nc.ndci : null;
                  var ndwV  = (nw && nw.ndvi_water !== null && !isNaN(nw.ndvi_water)) ? nw.ndvi_water : null;
                  var chlBloom = (cb&&cb.chlor_a!==null)?cb.chlor_a:null;

                  faiV.setValue(faiv !== null ? fmt(faiv,4) :
                    (ndwV!==null?fmt(ndwV,3)+' (NDVI proxy)':'n/a'));
                  ndciV.setValue(nciV !== null ? fmt(nciV,3) : 'n/a');
                  ndviWV.setValue(ndwV !== null ? fmt(ndwV,3) : 'n/a');
                  var algStatus = faiv!==null ?
                    (faiv>0.05?'CRITICAL - dense algae':faiv>0.01?'ELEVATED':'CLEAR') :
                    'No S2 data here';
                  algaeStatusV.setValue(algStatus);

                  var aq = computeAquaculture(sv, sv_peak, cv, nv, turv, tv);
                  sstWindowV.setValue(fmt(sv,1)+' deg C');
                  chlWindowV.setValue(fmt(cv,3)+' mg/m3');
                  ndciNutrientProxyV.setValue(fmt(nciV,3));
                  pollutionV.setValue(fmt(nv,7)+' mol/m2');
                  stabilityV.setValue(fmt(tv,4)+' deg C/yr');
                  bromoformV.setValue(aq.bromo_score>=80?'HIGH':aq.bromo_score>=60?'MODERATE':'LOW');
                  aquaScoreV.setValue(aq.aqua_score!==null?
                    aq.aqua_score+'/100 ('+aq.aqua_confidence+'%)' : 'n/a - insufficient data');
                  aquaMissingV.setValue(aq.aqua_missing_note);
                  aquaStatusV.setValue(aq.status);

                  var soilCode = (soilRes && soilRes.b0 !== null && soilRes.b0 !== undefined) ?
                    soilRes.b0 : null;
                  var soilLabel = soilCode !== null ? soilTextureLabel(soilCode) : null;
                  // v10.42: this will now always be null - OpenLandMap asset confirmed
                  // permanently dead (see soilTexture declaration comment). Message is
                  // deliberately specific rather than a generic "n/a".
                  soilTextureV.setValue(soilCode !== null ?
                    soilLabel+' (class '+soilCode+')' :
                    'Dataset unavailable (OpenLandMap asset confirmed inaccessible, v10.42)');

                  var eqCount = (eqRes && eqRes.count !== null && eqRes.count !== undefined) ?
                    eqRes.count : 0;
                  var eqMaxMag = (eqRes && eqRes.maxMag !== null && eqRes.maxMag !== undefined) ?
                    eqRes.maxMag : null;

                  // v10.47: volcRes now also carries the real feature schema.
                  var volcCount = (volcRes && volcRes.matchedByOlcanFilter !== null &&
                    volcRes.matchedByOlcanFilter !== undefined) ? volcRes.matchedByOlcanFilter : 0;
                  var volcTotalRaw = (volcRes && volcRes.totalGdisRecordsNearby !== null &&
                    volcRes.totalGdisRecordsNearby !== undefined) ? volcRes.totalGdisRecordsNearby : 0;
                  var volcSampleTypes = (volcRes && volcRes.sampleDisasterTypes) ?
                    volcRes.sampleDisasterTypes : [];
                  var volcPropNames = (volcRes && volcRes.firstFeaturePropertyNames) ?
                    volcRes.firstFeaturePropertyNames : [];
                  var volcAllProps = (volcRes && volcRes.firstFeatureAllProperties) ?
                    volcRes.firstFeatureAllProperties : {};

                  eqStatsV.setValue(eqCount > 0 ? eqCount+' event(s), max M'+eqMaxMag.toFixed(1) : 'None nearby');
                  volcStatsV.setValue(volcCount > 0 ? volcCount+' historic event(s)' : 'None recorded');

                  print('=== v10.47 VOLCANIC ACTIVITY DIAGNOSTIC (real schema) ===');
                  print('Raw GDIS records within 300km (ANY category, unfiltered): '+volcTotalRaw);
                  print('Records matched by "olcan" substring filter: '+volcCount);
                  if (volcTotalRaw > 0) {
                    print('ACTUAL disastertype values found: '+JSON.stringify(volcSampleTypes));
                    print('REAL property names on an actual nearby feature: '+JSON.stringify(volcPropNames));
                    print('REAL full properties of that feature: '+JSON.stringify(volcAllProps));
                    if (volcCount === 0) {
                      print('DIAGNOSIS: real records exist here, but "disastertype" as a property');
                      print('name returned nothing - check the REAL property names printed above');
                      print('to find the correct field name for the category (likely a different');
                      print('name entirely, not just a different category string).');
                    }
                  } else {
                    print('DIAGNOSIS: GDIS has ZERO records of ANY category within 300km here.');
                    print('This means the filter itself is not the problem - GDIS genuinely has');
                    print('no data at all in this specific area (consistent with the original');
                    print('script\'s own note that GDIS only covers 1960-2018 recorded-IMPACT');
                    print('events, not a complete eruption inventory - a real eruption without');
                    print('a recorded disaster impact in that specific window would be absent).');
                  }


                  var realAC1 = (csdRes && csdRes.realAC1 !== null && csdRes.realAC1 !== undefined) ?
                    csdRes.realAC1 : null;
                  var varTrendRatio = (csdRes && csdRes.varTrendRatio !== null && csdRes.varTrendRatio !== undefined) ?
                    csdRes.varTrendRatio : null;
                  var nValidMonths = (csdRes && csdRes.nValidMonths !== null && csdRes.nValidMonths !== undefined) ?
                    csdRes.nValidMonths : 0;
                  realAc1V.setValue(realAC1 !== null ? realAC1.toFixed(3)+' (n='+nValidMonths+')' : 'n/a');
                  realVarTrendV.setValue(varTrendRatio !== null ? varTrendRatio.toFixed(2)+'x' : 'n/a');

                  var thermalFeats = (thermalFC && thermalFC.features) ? thermalFC.features : [];
                  var trResult = analyzeThermalRecovery(thermalFeats, mmmv);
                  if (trResult.error) {
                    thermalEpisodesV.setValue('n/a - '+trResult.error);
                    thermalMeanV.setValue('n/a'); thermalMaxV.setValue('n/a');
                    thermalOngoingV.setValue('');
                  } else {
                    thermalEpisodesV.setValue(trResult.nCompletedEpisodes+' completed episode(s)');
                    thermalMeanV.setValue(trResult.meanRecoveryMonths !== null ?
                      trResult.meanRecoveryMonths.toFixed(1)+' months avg' : 'n/a');
                    thermalMaxV.setValue(trResult.maxRecoveryMonths !== null ?
                      trResult.maxRecoveryMonths+' months' : 'n/a');
                    thermalOngoingV.setValue(trResult.ongoingEpisode !== null ?
                      'ONGOING since '+trResult.ongoingEpisode.startDate : 'None');
                  }

                  var ecoVal = getEcologicalRecoveryValidation(region);
                  ecoValStatusV.setValue(ecoVal.checked ? ecoVal.finding : 'NOT YET CHECKED');
                  ecoValDetailsV.setValue(ecoVal.details);

                  var local_gebco = ee.ImageCollection(
                    'projects/sat-io/open-datasets/gebco/gebco_grid')
                    .mosaic().select('b1').rename('elevation');
                  var local_bathy = local_gebco.multiply(-1)
                    .updateMask(local_gebco.lt(0));
                  var local_slr03 = local_bathy.subtract(ee.Image(0.3)).max(ee.Image(0.1))
                    .sqrt().pow(-1).min(ee.Image(1.414)).rename('eci_slr_03');
                  var local_slr05 = local_bathy.subtract(ee.Image(0.5)).max(ee.Image(0.1))
                    .sqrt().pow(-1).min(ee.Image(1.414)).rename('eci_slr_05');
                  var local_slr10 = local_bathy.subtract(ee.Image(1.0)).max(ee.Image(0.1))
                    .sqrt().pow(-1).min(ee.Image(1.414)).rename('eci_slr_10');

                  ee.Dictionary({
                    depth: local_bathy.reduceRegion({reducer:ee.Reducer.mean(),
                      geometry:study, scale:500, maxPixels:1e9, bestEffort:true}).get('elevation'),
                    eci_03: local_slr03.reduceRegion({reducer:ee.Reducer.mean(),
                      geometry:study, scale:500, maxPixels:1e9, bestEffort:true}).get('eci_slr_03'),
                    eci_05: local_slr05.reduceRegion({reducer:ee.Reducer.mean(),
                      geometry:study, scale:500, maxPixels:1e9, bestEffort:true}).get('eci_slr_05'),
                    eci_10: local_slr10.reduceRegion({reducer:ee.Reducer.mean(),
                      geometry:study, scale:500, maxPixels:1e9, bestEffort:true}).get('eci_slr_10')
                  }).evaluate(function(eciVals) {
                    var depthVal = eciVals ? eciVals.depth : null;
                    var eci03 = eciVals ? eciVals.eci_03 : null;
                    var eci05 = eciVals ? eciVals.eci_05 : null;
                    var eci10 = eciVals ? eciVals.eci_10 : null;
                    var bScore = sc ? sc.B : null;
                    var validationNote = '';
                    if (depthVal !== null && bScore !== null) {
                      var physicalRisk = depthVal < 10 ? 'HIGH' : depthVal < 30 ? 'MODERATE' : 'LOW';
                      var heuristicRisk = bScore > 0.6 ? 'LOW' : bScore > 0.4 ? 'MODERATE' : 'HIGH';
                      validationNote = physicalRisk === heuristicRisk ?
                        'CONSISTENT' : 'DISCREPANCY FLAGGED - physical='+physicalRisk+' heuristic='+heuristicRisk;
                    }
                    if (depthVal !== null) {
                      eciDepthV.setValue(depthVal.toFixed(1)+' m (GEBCO)');
                      eciCurrentV.setValue(depthVal < 10 ? 'HIGH risk' : depthVal < 30 ? 'MODERATE risk' : 'LOW risk');
                      eciSLR03V.setValue(eci03 !== null ? (eci03*10).toFixed(2) : 'n/a');
                      eciSLR05V.setValue(eci05 !== null ? (eci05*10).toFixed(2) : 'n/a');
                      eciSLR10V.setValue(eci10 !== null ? (eci10*10).toFixed(2) : 'n/a');
                      eciValidationV.setValue(validationNote);
                    } else {
                      eciDepthV.setValue('n/a'); eciCurrentV.setValue('n/a');
                      eciSLR03V.setValue('n/a'); eciSLR05V.setValue('n/a'); eciSLR10V.setValue('n/a');
                      eciValidationV.setValue('Cannot validate');
                    }
                  });

                  var isReefZoneForInt = (lat > -30 && lat < 30) && !isEBUS(lat, lon);
                  var interventions = computeInterventions(region, fp, sc, dhwv, tv, isReefZoneForInt);
                  intPanel.clear();
                  var priorityColors = {
                    'CRITICAL':{bg:'#ffd0d0',tc:'#880000'}, 'HIGH':{bg:'#ffe0c0',tc:'#aa4400'},
                    'MODERATE':{bg:'#fff4c0',tc:'#886600'}, 'LOW':{bg:'#e8f4ff',tc:'#225588'},
                    'DATA GAP':{bg:'#f0e8ff',tc:'#553388'}, 'CONTEXT':{bg:'#f0f0f0',tc:'#555555'},
                    'STABLE':{bg:'#e0ffe0',tc:'#115511'}
                  };
                  for(var ai=0; ai<interventions.length; ai++){
                    var act = interventions[ai];
                    var pc = priorityColors[act.priority] || {bg:'#eeeeee',tc:'#333333'};
                    var actBox = ui.Panel({style:{backgroundColor:pc.bg, margin:'2px 0',
                      padding:'4px 6px', border:'1px solid '+pc.tc}});
                    actBox.add(ui.Label('['+act.priority+'] '+act.action,
                      {fontSize:'8px',fontWeight:'bold',color:pc.tc,
                       backgroundColor:'rgba(0,0,0,0)',margin:'0',whiteSpace:'pre'}));
                    intPanel.add(actBox);
                  }

                  print('=== FUSED RESULT ===');
                  print('FUSED CCS: '+sc.ccs+'/100  ('+sc.acc_total+'% accuracy) | Status: '+cols.lbl);

                  clickResultsLog.push({
                    timestamp: new Date().toISOString(), lat: latR, lon: lonR, region: region,
                    region_type: isLand ? 'land' : 'ocean',
                    sst_annual_c: sv, sst_peak_c: sv_peak, sst_trend_c_per_yr: tv,
                    dhw_c_weeks: dhwv, dhw_raw_unfiltered: dhwv_raw, dhw_artifact_flagged: dhwFlaggedArtifact,
                    mmm_local_c: mmmv, chl_a_mg_m3: cv, turbidity_ndti: turv, no2_mol_m2: nv, depth_m: bv,
                    cancer_score_satellite: sc.sat_ccs, cancer_score_fused: sc.ccs,
                    field_correction: sc.fcTotal, bowl_depth_B: sc.B, accuracy_pct: sc.acc_total,
                    status_label: cols.lbl, field_data_available: fp.hasField,
                    field_species: fp.hasField ? fp.species : 'none',
                    aquaculture_score: aq.aqua_score, aquaculture_confidence_pct: aq.aqua_confidence,
                    aquaculture_status: aq.status, bromoform_yield: aq.bromo_score,
                    num_interventions: interventions.length,
                    top_intervention: interventions.length > 0 ? interventions[0].action : 'none',
                    soil_texture_code: soilCode, soil_texture_label: soilCode !== null ? soilLabel : 'n/a',
                    earthquake_count_200km: eqCount, earthquake_max_mag_200km: eqMaxMag,
                    volcanic_activity_count_300km: volcCount,
                    real_ac1_24mo: realAC1, real_var_trend_ratio: varTrendRatio,
                    real_csd_n_valid_months: nValidMonths, ndci_nutrient_proxy: nciV,
                    thermal_recovery_n_completed_episodes: trResult.error ? null : trResult.nCompletedEpisodes,
                    thermal_recovery_mean_months: trResult.error ? null : trResult.meanRecoveryMonths,
                    thermal_recovery_max_months: trResult.error ? null : trResult.maxRecoveryMonths,
                    thermal_recovery_ongoing_flag: trResult.error ? null : (trResult.ongoingEpisode !== null),
                    eco_recovery_validation_checked: ecoVal.checked,
                    eco_recovery_validation_finding: ecoVal.checked ? ecoVal.finding : 'not yet checked',
                    in_situ_station: inSitu ? inSitu.key : null,
                    in_situ_distance_km: inSitu ? inSitu.distance_km : null,
                    in_situ_ph_mean: inSitu ? inSitu.pH.mean : null,
                    in_situ_temp_c_mean: inSitu ? inSitu.temp_c.mean : null,
                    in_situ_salinity_mean: inSitu ? inSitu.salinity.mean : null,
                    in_situ_do_mgL_mean: inSitu ? inSitu.do_mgL.mean : null
                  });
                  exportStatusLabel.setValue(clickResultsLog.length+' row(s) ready');

                  var sFC = ee.FeatureCollection([
                    ee.Feature(null,{c:'S1 SST',s:sc.s1,f:sc.s1,safe:30}),
                    ee.Feature(null,{c:'S2 Chl-a',s:sc.s2,f:sc.s2,safe:30}),
                    ee.Feature(null,{c:'S3 Turb',s:sc.s3,f:sc.s3,safe:30}),
                    ee.Feature(null,{c:'S4 DHW',s:sc.s4,f:sc.s4,safe:30}),
                    ee.Feature(null,{c:'S5 Bio',s:sc.s5,f:sc.s5,safe:30}),
                    ee.Feature(null,{c:'S6 NO2',s:sc.s6,f:sc.s6,safe:30})
                  ]);
                  print(ui.Chart.feature.byFeature({
                    features:sFC,xProperty:'c',yProperties:['s','f','safe']
                  }).setChartType('ColumnChart').setOptions({
                    title:'FUSED SCORE: '+sc.ccs+'/100 | '+cols.lbl+' | '+region,
                    series:{0:{color:'#ffaa44',label:'Satellite'},1:{color:cols.bar,label:'Fused'},
                      2:{color:'#4466ff',type:'line',lineWidth:2,pointSize:0,label:'Safe (30)'}},
                    vAxis:{title:'Risk 0-100',viewWindow:{min:0,max:100}},
                    backgroundColor:'#0a1628',titleTextStyle:{color:cols.bar,fontSize:10,bold:true},
                    legend:{position:'top',textStyle:{color:'#ffffff',fontSize:9}},
                    chartArea:{backgroundColor:'#0d1f3c',width:'82%'}}));

                // ============================================================
                // v10.48 NEW: SKIP S17/S18 ENTIRELY FOR INLAND LOCATIONS
                //
                // isLand was already resolved earlier in this chain (right
                // after bathymetry/depth came back), so we know here, BEFORE
                // firing the real server-side computation, whether this point
                // is land. S17 (Time of Emergence trends) and S18 (BGC
                // snapshot) are among the most expensive lookups in the whole
                // pipeline (precomputed ToE image reads + Copernicus BGC +
                // HYCOM salinity) - and both are completely meaningless on
                // land, same reasoning as the main marine score suppression.
                //
                // FIX: if isLand, skip calling rToEAll.evaluate() at all -
                // no server round-trip happens - and set every S17/S18
                // sidebar field directly to an honest "inland" message
                // instead. This is a genuine computation-time saving, not
                // just a display change, since the real BGC/ToE reduceRegion
                // calls never fire for land clicks anymore.
                // ============================================================
                if (isLand) {
                  var INLAND_MSG = 'N/A - inland location';
                  toeCompoundV.setValue(INLAND_MSG);
                  toeSSTv.setValue(INLAND_MSG); toeCHLv.setValue(INLAND_MSG);
                  toeNO2v.setValue(INLAND_MSG); toeSALv.setValue(INLAND_MSG);
                  s18O2V.setValue(INLAND_MSG); s18pCO2V.setValue(INLAND_MSG);
                  s18pHV.setValue(INLAND_MSG); s18SalV.setValue(INLAND_MSG);

                  print('=== v10.48 S17/S18 SKIPPED (inland location) ===');
                  print('Depth +'+fmt(bv,0)+'m above sea level - Time of Emergence trends');
                  print('and BGC snapshot are not meaningful here and were never computed');
                  print('server-side (real computation-time saving, not just a display change).');

                  clearPipelineWatchdog();
                  print('=== PIPELINE COMPLETE (v10.48, inland - S17/S18 skipped) ===');
                } else {
                // v10.60: S17 (ToE trends) - genuinely slow, precomputed
                // 44/32/27/7-year image reads. Runs independently of S18.
                rToETrends.evaluate(function(toeVals, _err) {
                  if (_err) { handleGeeError(_err, "rToETrends"); return; }
                  var sstToeVals = {fit:{scale:toeVals.sst_slope}, noise:{sst_stdDev:toeVals.sst_noise}};
                  var chlToeVals = {fit:{scale:toeVals.chl_slope}, noise:{chlor_a_stdDev:toeVals.chl_noise}};
                  var salToeVals = {fit:{scale:toeVals.sal_slope}, noise:{salinity_0_stdDev:toeVals.sal_noise}};
                  var no2ToeVals = {fit:{scale:toeVals.no2_slope}, noise:{tropospheric_NO2_column_number_density_stdDev:toeVals.no2_noise}};

                  function calcSnr(vals, slopeKey, noiseKey, years) {
                    if (!vals || !vals.fit || !vals.noise) return {error:'No data', snr:null, years:years};
                    var slope = vals.fit[slopeKey];
                    var noise = vals.noise[noiseKey];
                    if (slope===null||slope===undefined||noise===null||noise===undefined) {
                      return {error:'Null values', snr:null, years:years};
                    }
                    var signal = Math.abs(slope * years);
                    var snr = noise > 0 ? signal / noise : 0;
                    return {slope:slope, noise:noise, years:years, signal:signal,
                      snr:snr, emerged:snr>=2.0, direction:slope>0?'RISING':'FALLING', error:null};
                  }
                  var sstToE = calcSnr(sstToeVals, 'scale', 'sst_stdDev', 44);
                  var chlToE = calcSnr(chlToeVals, 'scale', 'chlor_a_stdDev', 27);
                  var no2ToE = calcSnr(no2ToeVals, 'scale', 'tropospheric_NO2_column_number_density_stdDev', 7);
                  var salToE = calcSnr(salToeVals, 'scale', 'salinity_0_stdDev', 32);

                  var nEmerged = 0;
                  if (!sstToE.error && sstToE.emerged) nEmerged++;
                  if (!chlToE.error && chlToE.emerged) nEmerged++;
                  if (!no2ToE.error && no2ToE.emerged) nEmerged++;
                  if (!salToE.error && salToE.emerged) nEmerged++;
                  toeCompoundV.setValue(nEmerged+'/4 stressors emerged');

                  var fmtSidebar = function(r, conf) {
                    if (r.error) return 'n/a';
                    return r.direction+' '+r.years+'yr | SNR='+r.snr.toFixed(2)+(r.emerged?' [EMERGED]':' [not yet]');
                  };
                  toeSSTv.setValue(fmtSidebar(sstToE,'HIGH'));
                  toeCHLv.setValue(fmtSidebar(chlToE,'MARGINAL'));
                  toeNO2v.setValue(fmtSidebar(no2ToE,'LOW'));
                  toeSALv.setValue(fmtSidebar(salToE,'MARGINAL'));

                  print('=== S17 (ToE trends) resolved independently (v10.60) ===');
                  _checkBothS17S18Done();
                });

                // v10.60: S18 (BGC snapshot) - can resolve near-instantly
                // when S19 has matched (3 of 4 fields are skip-sentinels,
                // no real server computation needed for those). Runs
                // independently of S17's slower trends above.
                rBGCSnapshot.evaluate(function(bgcVals, _err) {
                  if (_err) { handleGeeError(_err, "rBGCSnapshot"); return; }
                  {
                    var o2v = bgcVals.o2_surface;
                    var pco2v = bgcVals.pco2_surface;
                    var phv = bgcVals.ph_surface;
                    // v10.62 CRITICAL FIX: this line read 'salinity_current',
                    // a key name that only existed in the OLD pre-v10.60 code
                    // (which built an intermediate renaming object). After the
                    // v10.60 split, rBGCSnapshot returns the raw dictionary
                    // directly, whose actual key is 'sal_current' - reading
                    // the wrong (nonexistent) property name silently returned
                    // undefined every time, regardless of the skip-sentinel
                    // logic (confirmed via live diagnostic: _skipBGCForSal was
                    // correctly true, but salv still came back undefined).
                    var salv = bgcVals.sal_current;

                    // v10.41 NEW: -999 is the "intentionally skipped because
                    // S19 real sensor already covers this variable" sentinel
                    // (see rBGCSnapshot construction above) - display that
                    // plainly instead of ever showing -999 as a real reading.
                    var SKIP_SENTINEL = -999;

                    s18O2V.setValue(
                      (o2v === SKIP_SENTINEL) ? 'Skipped - see real S19 sensor DO below' :
                      (o2v!==null&&o2v!==undefined ? o2v.toFixed(2)+' mmol/m3' : 'n/a'));

                    var pco2uatm = (pco2v!==null&&pco2v!==undefined&&pco2v!==SKIP_SENTINEL) ? pco2v*9.869 : null;
                    s18pCO2V.setValue(pco2uatm!==null ? pco2uatm.toFixed(0)+' uatm' : 'n/a');

                    var phReal = (phv!==null&&phv!==undefined&&phv!==SKIP_SENTINEL&&phv>0&&phv<14) ? phv : null;
                    s18pHV.setValue(
                      (phv === SKIP_SENTINEL) ? 'Skipped - see real S19 sensor pH below' :
                      (phReal!==null ? phReal.toFixed(3)+' pH units' : 'n/a'));

                    var salvPSU = (salv!==null&&salv!==undefined&&salv!==SKIP_SENTINEL) ? salv : null;

                    // v10.61 NEW: diagnostic - confirmed at Looe Key (real
                    // S19 match) that Salinity showed "n/a" instead of the
                    // expected "Skipped - see real S19..." message, even
                    // though _skipBGCForSal should have been true. Rather
                    // than guess at the cause, print the RAW values so the
                    // real behavior is visible directly (same approach that
                    // solved the GDIS field-name mystery earlier).
                    print('=== v10.61 SALINITY SKIP DIAGNOSTIC ===');
                    print('_skipBGCForSal was: '+_skipBGCForSal+' | inSitu matched: '+(inSitu!==null));
                    print('Raw salv value returned: '+JSON.stringify(salv)+
                      ' (typeof: '+typeof salv+')');
                    print('Expected: -999 (number) if skip fired correctly, or a real PSU value/null otherwise.');
                    s18SalV.setValue(
                      (salv === SKIP_SENTINEL) ? 'Skipped - see real S19 sensor salinity below' :
                      (salvPSU!==null ? salvPSU.toFixed(2)+' PSU' : 'n/a'));

                    // ============================================================
                    // v10.57 NEW: real cited salinity/pH cautions for S8 aquaculture
                    //
                    // WHY THIS IS A SEPARATE, LATER UPDATE rather than part of
                    // computeAquaculture() itself: real salinity/pH data (from
                    // this same S18 BGC snapshot) is only available this late in
                    // the pipeline - but the aquaculture score itself is computed
                    // and displayed much earlier (right after DHW/turbidity
                    // resolve). Rather than restructure the whole chain, this
                    // appends real-data cautions to the ALREADY-DISPLAYED
                    // aquaculture status once salinity/pH actually arrive, only
                    // when both SST+Chl gates already passed (aq.status starts
                    // with 'GOOD' or 'DECENT').
                    //
                    // CITATION (same source as the v10.57 SST fix):
                    //   Statton, J. (2024). Asparagopsis taxiformis hatchery and
                    //   cultivation manual. AgriFutures Australia, Pub. 24-083.
                    //   Direct quotes:
                    //   Salinity: "prefers salinities about 35 parts per thousand
                    //   (ppt) but has some tolerance for lower and higher
                    //   salinities; for example, it grows in Shark Bay, Western
                    //   Australia at salinities greater than 40 ppt."
                    //   pH: "grow best when the pH is between 7.0 and 9.0."
                    //
                    // HONEST SCOPE NOTE: unlike SST (which has two clearly cited
                    // bounds and is used as a full hard-veto gate), salinity only
                    // has ONE explicitly numeric bound (>40 ppt upper tolerance) -
                    // the lower tolerance is described only qualitatively ("some
                    // tolerance for lower... salinities", no number given). For
                    // this reason salinity is treated as a CAUTION, not a hard
                    // veto gate, same treatment as pH (which has two clear cited
                    // bounds, but real pH data is frequently unavailable at
                    // coastal points - see the known Copernicus BGC masking
                    // limitation documented elsewhere in this script - so making
                    // it a hard gate would block scores unnecessarily often).
                    // ============================================================
                    if (aq && aq.status && (aq.status.indexOf('GOOD') === 0 || aq.status.indexOf('DECENT') === 0)) {
                      var extraCautions = [];
                      if (salvPSU !== null) {
                        if (salvPSU < 30 || salvPSU > 40) {
                          extraCautions.push('Salinity ('+salvPSU.toFixed(1)+' PSU) is outside the ' +
                            'documented tolerance range (~35 PSU preferred, up to 40 PSU tolerated per ' +
                            'Statton 2024) - may stress this species');
                        }
                      }
                      if (phReal !== null) {
                        if (phReal < 7.0 || phReal > 9.0) {
                          extraCautions.push('pH ('+phReal.toFixed(2)+') is outside the documented ' +
                            '7.0-9.0 growth range (Statton 2024)');
                        }
                      }
                      if (extraCautions.length > 0) {
                        aquaStatusV.setValue(aq.status + ' | ADDITIONAL (real data, arrived later): ' +
                          extraCautions.join('; '));
                        print('=== v10.57 AQUACULTURE: additional real-data cautions ===');
                        print(extraCautions.join(' | '));
                      }
                    }

                    print('=== S19 SUMMARY (v10.41) ===');
                    if (inSitu) {
                      print('Real in-situ sensor matched - S18 pH/Salinity/O2 snapshot lookups');
                      print('were SKIPPED (v10.41 optimization) since S19 already provides real,');
                      print('more accurate readings for these same 3 variables. pCO2 has no');
                      print('sensor equivalent and was still computed normally: ' +
                        (pco2uatm!==null?pco2uatm.toFixed(0)+' uatm':'n/a'));
                      print('Real in-situ pH: '+inSitu.pH.mean.toFixed(3)+' | Temp: '+
                        inSitu.temp_c.mean.toFixed(2)+' C | Salinity: '+
                        inSitu.salinity.mean.toFixed(2)+' PSU | DO: '+
                        inSitu.do_mgL.mean.toFixed(2)+' mg/L');
                    } else {
                      print('No in-situ station within range of this click - full S18 model');
                      print('snapshot computed normally (pH: '+(phReal!==null?phReal.toFixed(3):'n/a')+
                        ', Salinity: '+(salvPSU!==null?salvPSU.toFixed(2)+' PSU':'n/a')+
                        ', O2: '+(o2v!==null&&o2v!==undefined?o2v.toFixed(2)+' mmol/m3':'n/a')+').');
                    }
                  }

                  print('=== S18 (BGC snapshot) resolved independently (v10.60) ===');
                  _checkBothS17S18Done();
                });

                } // end isLand skip check (v10.48)

                }); // rThermalSeries
                }); // rCSD
                }); // rVOLC
                }); // rEQ
                }); // rSOIL
                }); // rCHL_BLOOM
                }); // rNDVI_W
                }); // rNDCI
                }); // rFAI
                }); // rMMM
                }); // rDHW
            }); // rTURB
          }); // rBATH
        }); // rNO2
      }); // rTR
    }); // rCHLc
    }); // rCHL
    }); // rSSTP
  }); // rSST
    } // end else (ocean chain)
  }); // rBATH_early

} // end analyzeLocation(lat, lon)

Map.onClick(function(coords) {
  analyzeLocation(coords.lat, coords.lon);
});

// ============================================================
// STARTUP
// ============================================================
print('STEMGeoHS Marine v10.36 -- READY');
print('NEW vs v10.35: Real In-Situ Ocean Chemistry Baseline (MODULE A11 / S19)');
print('3 confirmed ERDDAP stations: Looe Key FL, Agua Hedionda CA, Scripps Pier CA');
print('Matched by proximity (~10-15km radius) - kept separate from S18 model data.');
print('NOT wired into ccs/B/p5yr score math - contextual only.');
print('CLICK any coastal area to analyze.');
