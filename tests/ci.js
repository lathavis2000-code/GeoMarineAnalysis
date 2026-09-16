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

  // The fingerprint checks above catch a figure from the WRONG site. They do
  // not catch a sentence that contradicts its own live number - which is how
  // the fifth instance shipped: the note read "Quality varies: 12/12 calendar
  // months have 3+ years, so the rest lean on the fitted harmonic" at SB02,
  // where there is no rest, while S21b two rows below said "0 fully imputed".
  // Half the sentence was computed and half was not.
  //
  // Assert the BRANCH MARKER, not a phrase. The first version of this test
  // grepped for 'leans on the fitted harmonic' and failed on the correct code,
  // because that string also occurs inside the negation "so no month leans on
  // the fitted harmonic". A substring cannot tell a claim from its denial -
  // which is the same mistake, in a test, as the defect it guards against.
  [{site: 'FK01', at: [24.43313, -81.93068],  clears: true,  thin: 4, absent: 1},
   // FK02 has 20 valid months against CLIM_MIN_TOTAL_SAMPLES=26, so it never
   // reaches the quality sentence at all. That is correct, not a gap.
   {site: 'FK02', at: [24.4888, -81.66632],   clears: false, thin: 7, absent: 1},
   {site: 'SB02', at: [42.470793, -70.24294], clears: true,  thin: 0, absent: 0}
  ].forEach(function (s) {
    updateS21Acoustic(s.at[0], s.at[1]);
    var note = String(s21NoteV.getValue());
    var varies  = note.indexOf('Quality VARIES') !== -1;
    var uniform = note.indexOf('Quality is UNIFORM') !== -1;
    var blocked = note.indexOf('does NOT clear the climatology gates') !== -1;

    ok(s.site + ': the note takes the branch its gate arithmetic requires',
      s.clears ? (!blocked && (varies || uniform)) : (blocked && !varies && !uniform),
      'clears=' + s.clears + ' blocked=' + blocked + ' varies=' + varies + ' uniform=' + uniform);

    if (s.clears) {
      var shortfall = (s.thin + s.absent) > 0;
      ok(s.site + ': VARIES only when some month is short of its own years',
        varies === shortfall, 'thin=' + s.thin + ' absent=' + s.absent + ' varies=' + varies);
      ok(s.site + ': UNIFORM only when every month carries its own years',
        uniform === !shortfall, 'thin=' + s.thin + ' absent=' + s.absent + ' uniform=' + uniform);
      // The exact wording that shipped: "the rest" with nothing left over.
      ok(s.site + ': never says "the rest" of a complete climatology',
        !(!shortfall && note.indexOf('the rest') !== -1));
      // A count of 1 must not take a plural verb.
      ok(s.site + ': agrees in number with the count it prints',
        note.indexOf(' 1\n  have ') === -1 && note.indexOf(' 1\n  lean ') === -1,
        note.replace(/\n/g, ' | ').slice(0, 160));

      // S21b describes FK01's climatology too, and says it "leans on the
      // harmonic fit for 5 of its 12 months and imputes one outright" - all
      // five non-well-sampled months, one of them wholly. v10.179's note
      // counted only the four PARTIAL ones, so the tool printed 4 and 5 for
      // the same set on the same screen. Neither was false; they defined
      // "lean" differently and said nothing about it. The headline count must
      // match the convention the other panel already uses.
      if (shortfall) {
        var head = note.match(/the other (\d+)\s*\n\s*leans? on the fitted harmonic/);
        ok(s.site + ': the leaner count matches S21b (every month short of 3 years)',
          !!head && parseInt(head[1], 10) === (s.thin + s.absent),
          'headline=' + (head ? head[1] : 'not found') +
          ' expected=' + (s.thin + s.absent));
        // And the split must add back up to the headline.
        // [\s\S] not . - the split wraps a line ("their own thin\n  samples"),
        // and JS dot does not cross a newline without the s flag. The first
        // version of this check used .*? and failed on correct output.
        var split = note.match(/- (\d+) blended with [\s\S]*?, and (\d+) with none at all/);
        if (s.thin > 0 && s.absent > 0) {
          ok(s.site + ': the split sums to the headline',
            !!split && (parseInt(split[1], 10) + parseInt(split[2], 10)) === (s.thin + s.absent),
            split ? split[1] + '+' + split[2] : 'split not found');
        }
      }
    }
  });

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

section('14. Panel version stamps cannot go stale');
// CHANGELOG.md states the rule and why it exists: "The running version is
// declared once, in code, as TOOL_VERSION ... three consecutive rounds of this
// file shipped with a stale marker because the number was typed in more than
// one place." Section 4 already checks the declaration. This checks the other
// half: a version typed into a PANEL HEADER.
//
// v10.178 stamped the S21 header with the then-current version. One release
// later the header said v10.178 on screen while the tool ran v10.179. A marker
// that must be hand-edited on every release will go stale, so the convention is
// that a stamp names the version a panel was INTRODUCED in - a fact about the
// past, which cannot drift.
(function () {
  var heads = SRC.match(/sHead\('[^']*'/g) || [];

  var running = heads.filter(function (h) { return h.indexOf(TOOL_VERSION) !== -1; });
  ok('no panel header carries the RUNNING version',
    running.length === 0, running.join(' | '));

  // Compare componentwise, NOT as floats. The first version of this check used
  // parseFloat and reported v10.67 as NEWER than v10.180, because 10.67 > 10.18
  // as a decimal. A dotted version is not a number - which is the same species
  // of error as everything else this file guards: a representation quietly
  // standing in for the thing it represents.
  function vparts(v) {
    return v.replace(/^v/, '').split('.').map(function (n) { return parseInt(n, 10); });
  }
  function older(a, b) {
    var x = vparts(a), y = vparts(b);
    for (var i = 0; i < Math.max(x.length, y.length); i++) {
      if ((x[i] || 0) !== (y[i] || 0)) return (x[i] || 0) < (y[i] || 0);
    }
    return false;
  }
  ok('componentwise version compare beats parseFloat',
    older('v10.67', 'v10.180') && older('v10.172', 'v10.180') &&
    !older('v10.180', 'v10.67') && !older('v10.180', 'v10.180'));

  var stamped = heads.filter(function (h) { return /\(v10\.\d+\)/.test(h); });
  var newer = stamped.filter(function (h) {
    var m = h.match(/\(v(10\.\d+)\)/);
    return !(m && older('v' + m[1], TOOL_VERSION));
  });
  ok('every panel version stamp predates the running version',
    newer.length === 0, newer.join(' | '));

  ok('the S21 stamp names the release that introduced it',
    SRC.indexOf("sHead('S21 - REAL PASSIVE-ACOUSTIC BIOPHONY (v10.172)'") !== -1);
})();

section('15. Study scripts LOAD, not just parse');
// node --check proves a file is syntactically valid. It says nothing about an
// identifier that is never defined, or one that is read before it is assigned.
// sb02_sar_export.js shipped with both: buildPersistence() read a collection
// named SCENES that exists nowhere in the file, and the block calling it ran at
// module level ~115 lines before s1Joined was assigned. Both passed node --check
// and the ES5 gate, and the script died in the Code Editor on line 778 before
// CHECK 1 could print. Executing the file under the stub catches that class.
//
// Each script runs in its OWN vm context so it cannot clobber the globals the
// sections above rely on.
(function () {
  var dir = path.join(ROOT, 'studies');
  if (!fs.existsSync(dir)) { ok('no studies/ directory to load', true); return; }

  var scripts = [];
  (function walk(d) {
    fs.readdirSync(d).forEach(function (f) {
      var full = path.join(d, f);
      if (fs.statSync(full).isDirectory()) return walk(full);
      if (/\.js$/.test(f)) scripts.push(full);
    });
  })(dir);

  ok('studies/ contains at least one script to load', scripts.length > 0);

  scripts.forEach(function (file) {
    var rel = path.relative(ROOT, file);
    var sandbox = { console: { log: function () {}, error: function () {} } };
    vm.createContext(sandbox);
    stub.install(sandbox);
    var threw = null;
    try {
      vm.runInContext(fs.readFileSync(file, 'utf8'), sandbox, { filename: rel });
    } catch (e) {
      threw = e.constructor.name + ': ' + e.message;
    }
    ok(rel + ' executes without throwing', threw === null, threw);

    if (threw === null && /sb02_sar_export/.test(rel)) {
      var exps = sandbox.Export.exports;
      // STAGE A writes one persistence asset per relative orbit.
      var assets = exps.filter(function (e) { return e.dest === 'toAsset'; });
      ok(rel + ': STAGE A exports one asset per relative orbit',
        assets.length === 3, 'got ' + assets.length);

      // PR #19: the write path appends _<PARAM_SET_ID> and the read path did
      // not, so asset mode asked for images STAGE A never wrote. Both build the
      // id from one expression now - assert the written id still carries it.
      ok(rel + ': every persistence assetId carries the parameter set id',
        assets.every(function (e) { return /_p\d+$/.test(String(e.opts.assetId)); }),
        assets.map(function (e) { return String(e.opts.assetId); }).join(' | '));

      // And the asset id must end with the same suffix as its description, so
      // the two cannot drift apart the way they did before.
      ok(rel + ': assetId and description agree',
        assets.every(function (e) {
          return String(e.opts.assetId).indexOf(String(e.opts.description)) !== -1;
        }),
        assets.map(function (e) {
          return String(e.opts.description) + ' -> ' + String(e.opts.assetId);
        }).join(' | '));

      ok(rel + ': the placeholder asset root is gone',
        !assets.some(function (e) { return /CHANGE_ME/.test(String(e.opts.assetId)); }));
    }
  });
})();

section('16. The ES5 gate itself works');
// The gate is the only thing standing between an ES2015 library method and a
// Code Editor runtime failure, and it had three defects at once: no Math.* rule
// (which is how Math.log10 reached a study script), line numbers that drifted
// after the first block comment (236 lines off in the SAR file), and no CLI at
// all - `node tests/es5.js <file>` printed nothing and exited 0 for any input.
// A gate nobody tests is a gate nobody can trust, so test it.
(function () {
  // Each sample is one line so the expected line number is unambiguous.
  var SAMPLES = [
    ['const declaration', 'const a = 1;'],
    ['let declaration', 'let b = 2;'],
    ['arrow function', 'var f = function () { return 0; }; var g = x => x;'],
    ['Object.assign', 'var o = Object.assign({}, {});'],
    ['Array.prototype.includes', 'var c = [1].includes(1);'],
    ['Array.prototype.find/findIndex', 'var d = [1].find(function (x) { return x; });'],
    ['String.prototype.repeat', 'var e = "x".repeat(3);'],
    ['Math.* ES2015', 'var h = Math.log10(10);'],
    ['Math.* ES2015', 'var i2 = Math.trunc(1.5);'],
    ['Map/Set constructor', 'var j = new Map();'],
    ['template literal', 'var k = `x`;']
  ];
  SAMPLES.forEach(function (sample) {
    var rule = sample[0], src = sample[1];
    var hits = es5.scan(src);
    ok('gate flags: ' + rule,
      hits.some(function (h) { return h.rule === rule; }),
      'scan returned ' + JSON.stringify(hits));
  });

  // Math.* specifically, because its absence is what let Math.log10 ship.
  ok('gate has a Math.* rule at all',
    es5.scan('Math.log10(1);').length > 0);
  // ES5 Math must NOT trip it.
  ok('gate does not flag ES5 Math',
    es5.scan('Math.log(1); Math.sqrt(4); Math.max(1,2); Math.LN10;').length === 0,
    JSON.stringify(es5.scan('Math.log(1); Math.sqrt(4); Math.max(1,2); Math.LN10;')));

  // Line numbers must survive a block comment and a multi-line string - the
  // stripper used to blank those newlines, so every later finding was reported
  // at the wrong line.
  var padded = '/*\n' + new Array(40).join('filler\n') + '*/\n' +
               'var s = "a\\\nb";\n' +
               'const late = 1;\n';
  var lateLine = padded.split('\n').indexOf('const late = 1;') + 1;
  var found = es5.scan(padded).filter(function (h) { return h.rule === 'const declaration'; });
  ok('line numbers survive a block comment and a multi-line string',
    found.length === 1 && found[0].line === lateLine,
    'reported ' + (found.length ? found[0].line : 'nothing') + ', expected ' + lateLine);

  // Clean ES5 must stay clean - a gate that cries wolf gets ignored.
  ok('gate passes plain ES5',
    es5.scan('var x = 1;\nfunction f(a) { return a * 2; }\n').length === 0);
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
