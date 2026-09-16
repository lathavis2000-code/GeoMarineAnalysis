// CI smoke tests for the STEMGeoHS Marine Earth Engine script.
//
// SCOPE, stated honestly: this file cannot run Earth Engine. Everything that
// depends on a real reduceRegions/evaluate result is UNTESTED here and always
// will be. What it does cover is the part of the tool that is plain JavaScript
// and has historically broken: the file parsing at all, loading without a
// ReferenceError, staying inside the ES5 subset the Code Editor supports, the
// statistics, the score guards, and whether click handlers are actually wired.
//
// A v10.164 release shipped a ReferenceError that broke FIND SWEET SPOT
// outright, and v10.149/v10.158 shipped ES6 methods that only failed in the
// browser. Both classes are caught below. Run: node tests/ci.js
'use strict';

var fs = require('fs');
var path = require('path');
var vm = require('vm');
var es5 = require('./es5.js');
var stub = require('./stub.js');

var ROOT = path.join(__dirname, '..');
var failures = [];
var checks = 0;

function ok(name, cond, detail) {
  checks++;
  if (cond) { console.log('  PASS  ' + name); }
  else { console.log('  FAIL  ' + name + (detail ? '  -> ' + detail : '')); failures.push(name); }
}
function section(t) { console.log('\n== ' + t + ' =='); }

// WHICH FILE IS THE TOOL. The repo root holds several GeoMarineAnalysisV10_*.js
// - V10_150, _152 and _153 are historical snapshots kept alongside the live
// V10_161 - so a glob picks the wrong one or refuses to choose. The README is
// the repo's own declaration of which file to paste into the Code Editor, and
// it has been re-pointed on every rename (V10_155 -> _156 -> _161), so it is
// the authority here. Resolving through it also means a README that falls out
// of step with the code fails CI instead of quietly misleading a new user.
function findToolFile() {
  var readme = path.join(ROOT, 'README.md');
  if (!fs.existsSync(readme)) {
    console.error('FATAL: README.md is missing, so the current tool file cannot be identified.');
    process.exit(2);
  }
  var named = {}, m, re = /GeoMarineAnalysisV10_\d+\.js/g;
  var txt = fs.readFileSync(readme, 'utf8');
  while ((m = re.exec(txt)) !== null) named[m[0]] = 1;
  var names = Object.keys(named);
  if (names.length !== 1) {
    console.error('FATAL: README.md must name exactly one GeoMarineAnalysisV10_*.js; it names ' +
      names.length + ': ' + JSON.stringify(names));
    process.exit(2);
  }
  var p = path.join(ROOT, names[0]);
  if (!fs.existsSync(p)) {
    console.error('FATAL: README.md points at ' + names[0] + ', which does not exist.');
    process.exit(2);
  }
  // Exactly one file may declare TOOL_VERSION. The header calls out three
  // consecutive releases that shipped a stale version marker because the number
  // was typed in more than one place; two live files would be that failure mode
  // at the file level.
  var declaring = fs.readdirSync(ROOT).filter(function (f) {
    return /^GeoMarineAnalysisV10_\d+\.js$/.test(f) &&
      /^var TOOL_VERSION\s*=/m.test(fs.readFileSync(path.join(ROOT, f), 'utf8'));
  });
  if (declaring.length !== 1 || declaring[0] !== names[0]) {
    console.error('FATAL: TOOL_VERSION must be declared by exactly one file, and it must be the one ' +
      'README.md names (' + names[0] + '). Declaring files: ' + JSON.stringify(declaring));
    process.exit(2);
  }
  return p;
}

var TOOL = findToolFile();
var SRC = fs.readFileSync(TOOL, 'utf8');
console.log('tool file: ' + path.basename(TOOL) + '  (' + SRC.length + ' bytes)');

section('1. Parse');
var parsed = true;
try { new vm.Script(SRC, { filename: TOOL }); }
catch (e) { parsed = false; ok('parses as JavaScript', false, e.message); }
if (parsed) ok('parses as JavaScript', true);

section('2. ES5 subset (Earth Engine Code Editor has no ES6)');
var es6 = es5.scan(SRC);
ok('no ES6 constructs in executable code', es6.length === 0,
  es6.slice(0, 8).map(function (f) { return 'line ' + f.line + ': ' + f.rule; }).join('; '));

section('3. Loads under stubbed ee/ui/Map');
var h = stub.install(global);
var loaded = true;
try { vm.runInThisContext(SRC, { filename: TOOL }); }
catch (e) { loaded = false; ok('whole file executes without throwing', false, e.message); }
if (loaded) ok('whole file executes without throwing', true);
if (!loaded) { report(); return; }

section('4. Version identity (declared once, printed once, and they agree)');
var decls = SRC.split('\n').filter(function (l) { return /^var TOOL_VERSION\s*=/.test(l); });
ok('exactly one TOOL_VERSION declaration', decls.length === 1, decls.length + ' found');
ok('startup banner names the running version',
  h.prints.length > 0 && h.prints[0].indexOf(TOOL_VERSION) !== -1,
  'banner was: ' + (h.prints[0] || '(none)'));
ok('a changelog entry exists for the running version',
  h.prints.some(function (p) { return p.indexOf(TOOL_VERSION.replace(/^v/, 'v')) === 0; }),
  'no startup entry begins with ' + TOOL_VERSION);

section('5. Seasonal Mann-Kendall');
function synth(nYears, driftPerYear, seed) {
  var rows = [], s = seed || 1;
  for (var y = 2010; y < 2010 + nYears; y++) {
    for (var m = 1; m <= 12; m++) {
      s = (s * 1103515245 + 12345) % 2147483648;
      rows.push({ m: y + '-' + (m < 10 ? '0' : '') + m, v: (y - 2010) * driftPerYear + Math.sin(m) * 5 + (s / 2147483648 - 0.5) });
    }
  }
  return rows;
}
var mkNull = seasonalMannKendall(synth(12, 0, 7));
var mkUp = seasonalMannKendall(synth(12, 0.5, 7));
var mkDown = seasonalMannKendall(synth(12, -0.5, 7));
ok('null series is not significant', mkNull.p > 0.05, 'p=' + mkNull.p);
ok('rising series is significant and positive', mkUp.p < 0.05 && mkUp.S > 0, 'p=' + mkUp.p + ' S=' + mkUp.S);
ok('falling series is significant and negative', mkDown.p < 0.05 && mkDown.S < 0, 'p=' + mkDown.p + ' S=' + mkDown.S);
ok('Sen slope recovers the injected drift', Math.abs(mkUp.slope - 0.5) < 0.05, 'slope=' + mkUp.slope);
ok('a seasonal cycle alone produces no trend', Math.abs(mkNull.tau) < 0.15, 'tau=' + mkNull.tau);
var tied = [];
for (var ty = 2010; ty < 2022; ty++) for (var tm = 1; tm <= 12; tm++) tied.push({ m: ty + '-' + (tm < 10 ? '0' : '') + tm, v: Math.sin(tm) });
var mkTied = seasonalMannKendall(tied);
ok('fully-tied series reports zero variance, not "insufficient seasons"',
  mkTied.error && mkTied.error.indexOf('zero variance') === 0, mkTied.error);
ok('single point reports insufficient seasons',
  seasonalMannKendall([{ m: '2019-01', v: 1 }]).error.indexOf('insufficient seasons') === 0);
ok('empty input does not throw', seasonalMannKendall([]).error !== null);
ok('NaN values are dropped, not propagated',
  !isNaN(seasonalMannKendall(synth(12, 0.5, 7).concat([{ m: '2030-01', v: NaN }])).tau));

section('6. Acoustic sites and climatology gates');
var nSites = Object.keys(ACOUSTIC_SITES).length;
ok('at least one acoustic site is registered', nSites >= 1, nSites + ' sites');
Object.keys(ACOUSTIC_SITES).forEach(function (k) {
  var s = ACOUSTIC_SITES[k];
  ok('site ' + k + ' has series/varKind/varUnit',
    Array.isArray(s.series) && s.series.length > 0 && !!s.varKind && !!s.varUnit);
  ok('site ' + k + ' months are unique and sorted', (function () {
    var seen = {}, prev = '';
    for (var i = 0; i < s.series.length; i++) {
      var m = s.series[i].m;
      if (seen[m] || m <= prev) return false;
      seen[m] = 1; prev = m;
    }
    return true;
  })());
});
// Regression on the v10.172 error. The point is NOT that every site passes -
// FK02 has 20 valid months and is correctly refused - but that the gate's
// decision matches the gate's own rule: nTotal against CLIM_MIN_TOTAL_SAMPLES,
// and calendar months holding AT LEAST ONE sample against
// CLIM_MIN_DISTINCT_MONTHS. v10.172 asserted a refusal by counting months with
// >=CLIM_MIN_SAMPLES_PER_MONTH instead, which is an internal quality counter
// and gates nothing. Asserting "all sites pass" would have re-hidden that, and
// would break the moment a legitimately-thin site is added.
Object.keys(ACOUSTIC_SITES).forEach(function (k) {
  var series = ACOUSTIC_SITES[k].series;
  var tv = series.map(function (r) { return { t: Date.parse(r.m + '-15T00:00:00Z'), v: r.v }; });
  var c = computeUsableClimatology(tv);
  var cal = {}, distinctAny = 0;
  series.forEach(function (r) { var m = parseInt(r.m.substring(5, 7), 10); cal[m] = (cal[m] || 0) + 1; });
  for (var m2 = 1; m2 <= 12; m2++) if ((cal[m2] || 0) >= 1) distinctAny++;
  var shouldPass = series.length >= CLIM_MIN_TOTAL_SAMPLES && distinctAny >= CLIM_MIN_DISTINCT_MONTHS;
  ok('climatology gate on ' + k + ' matches its own rule (expect ' + (shouldPass ? 'accept' : 'refuse') + ')',
    c.ok === shouldPass,
    'ok=' + c.ok + ' nTotal=' + c.nTotalSamples + ' nDistinct=' + c.nDistinctMonths +
    ' (>=1-sample months=' + distinctAny + ')');
  // And the gate must never be decided by the quality counter.
  var wellSampled = 0;
  for (var m3 = 1; m3 <= 12; m3++) if ((cal[m3] || 0) >= CLIM_MIN_SAMPLES_PER_MONTH) wellSampled++;
  if (shouldPass && wellSampled < CLIM_MIN_DISTINCT_MONTHS) {
    ok('gate on ' + k + ' is not decided by CLIM_MIN_SAMPLES_PER_MONTH', c.ok === true,
      'wellSampled=' + wellSampled + ' would fail a >=3-years reading, but the real gate accepts');
  }
});

section('7. Acoustic correction guards (A1)');
var fk = computeAcousticCorrection(24.43313, -81.93068);
ok('FK01 matches and is eligible', fk.matched === true && fk.eligible === true);
ok('FK01 contributes nothing on a null result', fk.A1 === 0 && fk.applied === false, 'A1=' + fk.A1);
ok('FK01 earns no accuracy for a null result', fk.accEarned === 0);
var sb = computeAcousticCorrection(42.470793, -70.24294);
ok('SB02 matches on distance', sb.matched === true);
ok('SB02 is REFUSED on variable kind', sb.eligible === false && sb.A1 === 0);
var away = computeAcousticCorrection(9.175, -81.981);
ok('a site with no coverage is neither matched nor eligible',
  away.matched === false && away.eligible === false && away.A1 === 0);

section('8. A1 applied path (no real datum reaches it, so it is synthesised)');
function fakeSite(drift, lat, lon) {
  var rows = [];
  for (var y = 2015; y < 2023; y++) for (var m = 1; m <= 12; m++)
    rows.push({ m: y + '-' + (m < 10 ? '0' : '') + m, v: 2.0 + (y - 2015) * drift + Math.sin(m) * 0.3, nh: 700, dep: '01' });
  return { label: 'SYNTH', lat: lat, lon: lon, radius_km: 20, source: 'synthetic', record: 'synthetic',
    band: 'test', varLabel: 'diel ratio (synthetic)', varKind: 'diel_ratio', varUnit: 'dB', series: rows };
}
ACOUSTIC_SITES.__test_falling = fakeSite(-0.15, 1.0, 1.0);
ACOUSTIC_SITES.__test_rising = fakeSite(0.15, 2.0, 2.0);
ACOUSTIC_SITES.__test_extreme = fakeSite(-5.0, 4.0, 4.0);
ACOUSTIC_SITES.__test_thin = { label: 'SYNTH thin', lat: 3, lon: 3, radius_km: 20, source: 's', record: 'r',
  band: 'b', varLabel: 'diel ratio (thin)', varKind: 'diel_ratio', varUnit: 'dB',
  series: [{ m: '2019-01', v: 2.0 }, { m: '2020-01', v: 1.0 }, { m: '2019-02', v: 2.0 }, { m: '2020-02', v: 1.0 }] };
var fall = computeAcousticCorrection(1.0, 1.0);
var rise = computeAcousticCorrection(2.0, 2.0);
var extreme = computeAcousticCorrection(4.0, 4.0);
var thin = computeAcousticCorrection(3.0, 3.0);
ok('falling biophony RAISES stress (positive A1)', fall.applied === true && fall.A1 > 0, 'A1=' + fall.A1);
ok('rising biophony LOWERS stress (negative A1)', rise.applied === true && rise.A1 < 0, 'A1=' + rise.A1);
ok('an applied term earns its accuracy credit', fall.accEarned === ACOUSTIC_ACC_GAIN);
ok('|A1| never exceeds the cap', Math.abs(extreme.A1) <= ACOUSTIC_MAX_POINTS, 'A1=' + extreme.A1);
ok('a series below the design floor does not vote', thin.applied === false && thin.A1 === 0, thin.reason);

section('9. Score inertness (A1 must move nothing where it is not measured)');
var fp0 = { hasField: false, urchin_N: null, anem_N: null, Cd: null, Pb: null, recruit: null };
[['Bocas', 9.175, -81.981], ['Hawaii', 19.711, -156.053], ['SB02', 42.470793, -70.24294], ['FK01', 24.43313, -81.93068]]
  .forEach(function (c) {
    var sc = computeScore(29.1, 0.6, 0.04, 8e-5, 0.12, 6, fp0, c[1], c[2]);
    var expected = Math.max(0, Math.min(100, sc.sat_ccs + sc.fcTotal));
    ok(c[0] + ': fused score unchanged by the acoustic term',
      sc.acousticA1 === 0 && sc.ccs === expected, 'A1=' + sc.acousticA1 + ' ccs=' + sc.ccs + ' expected=' + expected);
  });
var scApplied = computeScore(29.1, 0.6, 0.04, 8e-5, 0.12, 6, fp0, 1.0, 1.0);
ok('an applied A1 does move the fused score',
  scApplied.ccs === Math.max(0, Math.min(100, scApplied.sat_ccs + scApplied.fcTotal + scApplied.acousticA1)) &&
  scApplied.acousticA1 > 0, 'ccs=' + scApplied.ccs + ' A1=' + scApplied.acousticA1);
delete ACOUSTIC_SITES.__test_falling; delete ACOUSTIC_SITES.__test_rising;
delete ACOUSTIC_SITES.__test_extreme; delete ACOUSTIC_SITES.__test_thin;

section('10. Insufficient-data contract keeps one shape');
var okScore = computeScore(29.1, 0.6, 0.04, 8e-5, 0.12, 6, fp0, 24.43313, -81.93068);
var noScore = computeScore(null, null, null, null, null, null, fp0, 24.43313, -81.93068);
ok('no satellite input yields a null score, not a number', noScore.ccs === null && noScore.insufficientData === true);
['acousticA1', 'acousticApplied', 'acousticEligible', 'acousticMatched', 'acousticNote', 'acc_acoustic']
  .forEach(function (k) {
    ok('both returns carry ' + k, (k in okScore) && (k in noScore));
  });

section('11. S20 replay window (v10.176) and its baseline-overlap disclosure');
ok('replay floor reaches the acoustic record', S20_REPLAY_MIN_YEAR <= 2018,
  'floor=' + S20_REPLAY_MIN_YEAR + ', earliest acoustic month is 2018-11');
ok('MMM baseline constants drive the collection, not a literal',
  typeof MMM_BASELINE_START_YEAR === 'number' && typeof MMM_BASELINE_END_YEAR === 'number' &&
  MMM_BASELINE_END_YEAR > MMM_BASELINE_START_YEAR,
  MMM_BASELINE_START_YEAR + '-' + MMM_BASELINE_END_YEAR);
// The disclosure exists because the widened window now overlaps the baseline.
// If a future change moves the baseline entirely before the replay floor, the
// overlap is gone and the note is dead code - this asserts the relationship
// rather than the sentence, so it fails loudly either way.
var overlaps = S20_REPLAY_MIN_YEAR <= MMM_BASELINE_END_YEAR;
ok('replay window overlapping the MMM baseline is disclosed in the panel text', !overlaps ||
  SRC.indexOf('BASELINE OVERLAP') !== -1,
  'window ' + S20_REPLAY_MIN_YEAR + '-' + S20_REPLAY_MAX_YEAR + ' vs baseline ' +
  MMM_BASELINE_START_YEAR + '-' + MMM_BASELINE_END_YEAR);
// Drive the real handler across the boundaries. ee is re-pointed so evaluate()
// yields a DHW and the SUCCESS branch - the only one carrying the note - runs.
(function () {
  var savedEe = global.ee;
  function ch() {
    var f = function () { return ch(); };
    return new Proxy(f, {
      get: function (t, k) {
        if (k === 'evaluate') return function (cb) { cb({ dhw: 9.5, latestDate: null }, null); };
        if (k === 'then' || k === 'toJSON' || k === Symbol.toPrimitive || k === Symbol.iterator) return undefined;
        return ch();
      },
      apply: function () { return ch(); }
    });
  }
  global.ee = ch();
  global.lastClickLat = 24.43313; global.lastClickLon = -81.93068;
  function run(m) {
    speciesDateInput.setValue(m); speciesRiskStatusV.setValue('');
    try { speciesRiskBtn._cb(); } catch (e) { return 'THREW:' + e.message; }
    return String(speciesRiskStatusV.getValue());
  }
  var below = run((S20_REPLAY_MIN_YEAR - 1) + '-12');
  var floor = run(S20_REPLAY_MIN_YEAR + '-01');
  var above = run((S20_REPLAY_MAX_YEAR + 1) + '-01');
  var ceil = run(S20_REPLAY_MAX_YEAR + '-12');
  var inBase = run(MMM_BASELINE_END_YEAR + '-08');
  var outBase = run((MMM_BASELINE_END_YEAR + 1) + '-08');
  ok('a month below the floor is refused', below.indexOf('must be') !== -1, below.slice(0, 70));
  ok('the floor month itself is accepted', floor.indexOf('must be') === -1, floor.slice(0, 70));
  ok('a month above the ceiling is refused', above.indexOf('must be') !== -1, above.slice(0, 70));
  ok('the ceiling month itself is accepted', ceil.indexOf('must be') === -1, ceil.slice(0, 70));
  ok('a replay inside the baseline carries the overlap disclosure',
    inBase.indexOf('BASELINE OVERLAP') !== -1, inBase.slice(0, 90));
  ok('a replay after the baseline carries no disclosure',
    outBase.indexOf('BASELINE OVERLAP') === -1, outBase.slice(0, 90));
  ok('the disclosure names the self-suppressing direction',
    inBase.indexOf('SELF-SUPPRESSING') !== -1 && inBase.indexOf('UNDERSTATED') !== -1);
  ok('no disclosure line is long enough to clip in the 256px panel',
    inBase.split('\n').every(function (l) { return l.length <= 56; }),
    'longest=' + Math.max.apply(null, inBase.split('\n').map(function (l) { return l.length; })));
  global.ee = savedEe;
})();

section('12. Click path is wired (not just the functions)');
['s21VerdictV', 's21bVerdictV'].forEach(function (name) {
  global[name].setValue('__SENTINEL__');
});
try { analyzeLocation(24.43313, -81.93068); } catch (e) { /* stubs die later in the pipeline; S21 runs first */ }
ok('a map click reaches S21', s21VerdictV.getValue() !== '__SENTINEL__');
ok('a map click reaches S21b', s21bVerdictV.getValue() !== '__SENTINEL__');

section('13. S21 panel text belongs to the site on screen');
// The defect this section exists for has now shipped FOUR times (v10.173
// BIOPHONY TREND, v10.177 caveat, v10.178 STATE VARIABLE and TEST). Every
// instance passed node --check, the ES5 gate and every assertion above,
// because none of them is a syntax or type error: they are true sentences
// bound to a shared widget and rendered at a site they do not describe.
// The invariant that catches the whole class: after clicking site X, no S21
// label may contain a number or a phrase that belongs only to site Y.
(function () {
  var S21_LABELS = ['s21VerdictV', 's21VarV', 's21WhyV', 's21CaveatV',
                    's21RangeV', 's21SeasonV', 's21NV'];
  // Values that appear in exactly one site's record, and nowhere else.
  var FINGERPRINTS = [
    {site: 'FK01', at: [24.43313, -81.93068], marks: ['0.2951', '1.40 dB', '28 pooled']},
    {site: 'FK02', at: [24.4888, -81.66632],  marks: ['0.0297', '0.93 dB', '20 pooled']},
    {site: 'SB02', at: [42.470793, -70.24294], marks: ['0.3366', '13.27 dB', '44 pooled']}
  ];

  function textAt(lat, lon) {
    try { updateS21Acoustic(lat, lon); } catch (e) { return '__THREW__ ' + e.message; }
    return S21_LABELS.map(function (n) {
      return global[n] ? String(global[n].getValue()) : '';
    }).join('\n');
  }

  FINGERPRINTS.forEach(function (mine) {
    var txt = textAt(mine.at[0], mine.at[1]);
    ok(mine.site + ': its own figures are on screen',
      mine.marks.some(function (m) { return txt.indexOf(m) !== -1; }),
      'none of ' + mine.marks.join('/') + ' found');
    FINGERPRINTS.forEach(function (other) {
      if (other.site === mine.site) return;
      var leaked = other.marks.filter(function (m) { return txt.indexOf(m) !== -1; });
      ok(mine.site + ': no ' + other.site + ' figure leaks into its panel',
        leaked.length === 0, 'leaked: ' + leaked.join(', '));
    });
  });

  // The ratio claim is a SAFETY claim - it says the series is immune to
  // per-deployment gain drift. It is only true of a ratio. Asserting it at a
  // band-level site tells the reader SB02 is protected while S21b, in the same
  // panel, reports its drift as unresolved.
  var sbTxt = textAt(42.470793, -70.24294);
  ok('SB02 is never told it cancels a calibration offset',
    sbTxt.indexOf('cancels any constant per-deployment') === -1);
  ok('SB02 is told the opposite, explicitly',
    sbTxt.indexOf('NOT a ratio') !== -1 && sbTxt.indexOf('cancels NO') !== -1);
  var fkTxt = textAt(24.43313, -81.93068);
  ok('FK01 keeps the ratio argument, which is true there',
    fkTxt.indexOf('cancels any constant per-deployment') !== -1);

  // Reef language at a 42 N bank was one of the four false statements in v10.177.
  ok('no reef language reaches SB02',
    sbTxt.toLowerCase().indexOf('reef') === -1);

  // A click with no hydrophone in range must leave nothing behind from the last.
  updateS21Acoustic(0, 0);
  var empty = S21_LABELS.map(function (n) {
    return global[n] ? String(global[n].getValue()) : '';
  }).join('\n');
  var stale = [];
  FINGERPRINTS.forEach(function (f) {
    f.marks.forEach(function (m) { if (empty.indexOf(m) !== -1) stale.push(f.site + ':' + m); });
  });
  ok('a no-coverage click clears every site-specific figure',
    stale.length === 0, 'stale: ' + stale.join(', '));
})();

function report() {
  console.log('\n' + '-'.repeat(60));
  if (failures.length === 0) {
    console.log(checks + ' checks passed.');
    process.exit(0);
  }
  console.log(failures.length + ' of ' + checks + ' checks FAILED:');
  failures.forEach(function (f) { console.log('  - ' + f); });
  process.exit(1);
}
report();
